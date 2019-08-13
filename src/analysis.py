#!/usr/bin/env python

"""
cll-time_course
"""

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')

import os
import random
import string

import GPy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from scipy.stats import zscore
from sklearn.preprocessing import LabelEncoder
from statsmodels.stats.multitest import multipletests

from peppy import Project
from ngs_toolkit.atacseq import ATACSeqAnalysis


# graphics settings
sns.set_style("white")
plt.rcParams['svg.fonttype'] = 'none'

# random seed
SEED = int("".join(
    LabelEncoder()
    .fit(list(string.ascii_uppercase))
    .transform(list("BOCKLAB")).astype(str)))
random.seed(SEED)
np.random.seed(SEED)


def main():
    # Start project and analysis objects
    analysis = ATACSeqAnalysis(
        from_pep=os.path.join("metadata", "project_config.yaml"))
    cell_types = list(set([
        sample.cell_type for sample in analysis.samples]))

    # Get accessibility matrix with normalization
    # generate consensus peak set from all samples, cell types together
    analysis.get_consensus_sites()

    # annotate consensus peak set
    analysis.get_peak_gene_annotation()
    analysis.get_peak_genomic_location()
    analysis.get_peak_chromatin_state(
        chrom_state_file=get_chrom_state_file())
    analysis.calculate_peak_support(region_type="summits")

    # get coverage values for each region in each sample
    analysis.measure_coverage()
    analysis.coverage = analysis.coverage.loc[
        ~analysis.coverage.index.str.contains("chrX|chrY"), :]

    # data normalization
    analysis.matrix_norm = data_normalization(analysis)

    # annotate matrix with sample metadata and save
    analysis.annotate_samples(
        numerical_attributes=analysis.numerical_attributes)

    # annotate matrix with feature metadata
    analysis.annotate_features()
    analysis.to_pickle()

    #
    # Unsupervised analysis
    # # for all samples jointly
    analysis.unsupervised_analysis()

    # # for each cell type separately
    for cell_type in cell_types:
        analysis.unsupervised_analysis(
            samples=[s for s in analysis.samples if (s.cell_type == cell_type)],
            attributes_to_plot=['patient_id', 'timepoint', 'batch'],
            output_prefix="accessibility_{}_only".format(cell_type))

    #
    # Time Series Analysis
    matrix_file = os.path.abspath(os.path.join(
        "results", analysis.name + ".matrix_norm.csv"))
    gp_output_dir = os.path.join(analysis.results_dir, "gp_fits")
    mohgp_output_dir = os.path.join(analysis.results_dir, "mohgp_fits")
    fits_file = os.path.join(gp_output_dir, prefix + ".GP_fits.all_cell_types.csv")
    prefix = "accessibility.qnorm_pcafix_cuberoot"

    # fit GPs with varying and constant kernel to detect variable regulatory elements
    fit_gaussian_processes(
        analysis.matrix_norm,
        cell_types=cell_types,
        matrix_file=matrix_file,
        prefix=prefix)  # wait for jobs to complete

    analysis.fits = gather_gaussian_processes(
        analysis.matrix_norm,
        matrix_file=matrix_file,
        prefix=prefix)
    analysis.fits.to_csv(fits_file, index=True)

    visualize_gaussian_process_fits(
        analysis, analysis.fits, output_dir=gp_output_dir, prefix=prefix)

    # cluster variable regulatory elements with hierarchical mixtures of GPs (MOHGP)
    analysis.assignments = dict()
    for alpha in [0.05]:  # 0.01
        clust_prefix = prefix + ".p={}".format(alpha)
        fit_MOHGP(
            analysis.matrix_norm,
            matrix_file=matrix_file,
            fits_file=fits_file,
            cell_types=cell_types,
            n_clust=[3, 5, 4, 4, 4, 3],
            prefix=clust_prefix,
            output_dir=mohgp_output_dir,
            alpha=alpha)
        # wait for jobs to finish

        analysis.assignments[alpha] = gather_MOHGP(
            analysis.matrix_norm, analysis.fits,
            cell_types=cell_types,
            n_clust=[3, 5, 4, 4, 4, 3],
            prefix=clust_prefix,
            fits_dir=mohgp_output_dir,
            alpha=alpha,
            posterior_threshold=0.8)
        analysis.assignments[alpha].to_csv(os.path.join(
            mohgp_output_dir, clust_prefix + ".GP_fits.mohgp_clusters.csv"), index=True)
        analysis.assignments[alpha] = pd.read_csv(os.path.join(
            mohgp_output_dir, clust_prefix + ".GP_fits.mohgp_clusters.csv"), index_col=0)

        visualize_clustering_assignments(
            analysis.matrix_norm, analysis.assignments[alpha], prefix=clust_prefix,
            output_dir=mohgp_output_dir)

        # export clusters at gene level
        clusters_to_signatures(analysis.assignments[alpha])

    # Not included in publication:
    # # Linear time
    # gp_linear(analysis)
    # # Randomized samples
    # gp_random(analysis)

    #
    # Plot distribution of clusters per cell type dependent on their dynamic pattern
    alpha = 0.05
    assignments = analysis.assignments[alpha]
    prefix = "accessibility.qnorm_pcafix_cuberoot" + ".p={}".format(alpha)
    l_assignments = cluster_stats(assignments, prefix=prefix)

    #
    # Get enrichments of region clusters
    enrichments_dir = os.path.join(analysis.results_dir, "cluster_enrichments")
    l_assignments['comparison_name'] = l_assignments['cell_type'] + \
        "_" + l_assignments['cluster_name'].astype(str)

    # In addition, get enrichments for all of changing regions for each cell type
    _a = l_assignments.copy()
    _a['comparison_name'] = _a['comparison_name'].str.replace("_.*", "")
    l_assignments = l_assignments.append(_a)

    differential_enrichment(
        analysis,
        l_assignments,
        data_type="ATAC-seq",
        output_dir=enrichments_dir,
        output_prefix=prefix,
        genome="hg19",
        directional=False,
        max_diff=10000,
        sort_var="D",
        as_jobs=True
    )

    # get HOMER consensus motifs
    from ngs_toolkit.general import homer_combine_motifs
    comparison_dirs = [os.path.join("results", "cluster_enrichments", x) + ".all" for x in l_assignments['comparison_name'].unique()]
    output_dir = os.path.join("results", "cluster_enrichments_nostringent")

    for dir_ in comparison_dirs:
        cpus = 8
        genome = "hg19"
        combined_motifs = os.path.join(output_dir, "homerMotifs.filtered.motifs")
        cmd = (
            "findMotifsGenome.pl {bed} {genome}r {dir} -p {cpus} -nomotif -mknown {motif_file}"
            .format(bed=os.path.join(dir_, prefix + "_regions.bed"),
                    genome=genome, cpus=cpus, dir=dir_,
                    motif_file=combined_motifs))
        subprocess.call("sbatch -J homer.{d} -o {dir}.homer.log -p shortq -c 8 --mem 20000"
                        .format(d=os.path.basename(dir_), dir=dir_)
                        .split(" ")
                        + ['--wrap', cmd])

    homer_combine_motifs(
        comparison_dirs, output_dir,
        region_prefix=prefix,
        reduce_threshold=0.6, match_threshold=10, info_value=0.6,
        p_value_threshold=1e-10, fold_enrichment=None,
        cpus=8, run=1, as_jobs=True, genome="hg19",
        motif_database=None, known_vertebrates_TFs_only=False)

    # Not included in publication
    # # Inspect remaining regions
    # collect_differential_enrichment(
    #     l_assignments,
    #     directional=False,
    #     steps=['homer_consensus'],
    #     data_type="ATAC-seq",
    #     output_dir=enrichments_dir,
    #     output_prefix=prefix,
    #     permissive=True)

    # for simple, label in [(False, ""), (True, "-simple")]:
    #     enrichment_table = pd.read_csv(os.path.join(
    #         enrichments_dir, prefix + ".meme_ame.csv"))
    #     if simple:
    #         enrichment_table = enrichment_table[enrichment_table['comparison_name'].str.contains("_")]
    #     plot_differential_enrichment(
    #         enrichment_table,
    #         "motif",
    #         data_type="ATAC-seq",
    #         direction_dependent=False,
    #         output_dir=enrichments_dir,
    #         comp_variable="comparison_name",
    #         output_prefix=prefix + label,
    #         top_n=5)
    #     enrichment_table = pd.read_csv(os.path.join(
    #         enrichments_dir, prefix + ".lola.csv"))
    #     if simple:
    #         enrichment_table = enrichment_table[enrichment_table['comparison_name'].str.contains("_")]
    #     plot_differential_enrichment(
    #         enrichment_table,
    #         "lola",
    #         data_type="ATAC-seq",
    #         direction_dependent=False,
    #         output_dir=enrichments_dir,
    #         comp_variable="comparison_name",
    #         output_prefix=prefix + label,
    #         top_n=5)
    #     enrichment_table = pd.read_csv(os.path.join(
    #         enrichments_dir, prefix + ".homer_consensus.csv"))
    #     if simple:
    #         enrichment_table = enrichment_table[enrichment_table['comparison_name'].str.contains("_")]
    #     plot_differential_enrichment(
    #         enrichment_table,
    #         "homer_consensus",
    #         data_type="ATAC-seq",
    #         direction_dependent=False,
    #         output_dir=enrichments_dir,
    #         comp_variable="comparison_name",
    #         output_prefix=prefix + label,
    #         top_n=1)
    #     enrichment_table = pd.read_csv(os.path.join(
    #         enrichments_dir, prefix + ".enrichr.csv"))
    #     if simple:
    #         enrichment_table = enrichment_table[enrichment_table['comparison_name'].str.contains("_")]
    #     plot_differential_enrichment(
    #         enrichment_table,
    #         "enrichr",
    #         data_type="ATAC-seq",
    #         direction_dependent=False,
    #         output_dir=enrichments_dir,
    #         comp_variable="comparison_name",
    #         output_prefix=prefix + label,
    #         top_n=5)

    # each cell type independently
    for cell_type in cell_types:
        enrichment_table = pd.read_csv(os.path.join(
            enrichments_dir, prefix + ".lola.csv"))
        enrichment_table = enrichment_table.loc[
            enrichment_table['comparison_name'].str.contains("{}_".format(cell_type)), :]
        plot_differential_enrichment(
            enrichment_table,
            "lola",
            data_type="ATAC-seq",
            direction_dependent=False,
            output_dir=enrichments_dir,
            comp_variable="comparison_name",
            output_prefix=prefix + ".{}_only".format(cell_type),
            top_n=25)

        enrichment_table = pd.read_csv(os.path.join(
            enrichments_dir, prefix + ".enrichr.csv"))
        enrichment_table = enrichment_table.loc[
            enrichment_table['comparison_name'].str.contains("{}_".format(cell_type)), :]
        plot_differential_enrichment(
            enrichment_table,
            "enrichr",
            data_type="ATAC-seq",
            direction_dependent=False,
            output_dir=enrichments_dir,
            comp_variable="comparison_name",
            output_prefix=prefix + ".{}_only".format(cell_type),
            top_n=5)

    # Figure 2c, 3c
    prefix = "accessibility.qnorm_pcafix_cuberoot.p=0.05"
    enrichment_table = pd.read_csv(os.path.join(enrichments_dir, prefix + ".lola.csv"))
    plot_lola_enrichments(
        enrichment_table[enrichment_table['comparison_name'].str.contains("_")],
        top_n=5, comp_variable='comparison_name')

    # Figure 2d, 3d-e
    enrichment_table = pd.read_csv(os.path.join(enrichments_dir, prefix + ".enrichr.csv"))
    plot_enrichments(
        enrichment_table[enrichment_table['comparison_name'].str.contains("_")],
        top_n=5, comp_variable='comparison_name')

    # TF analysis
    transcription_factor_accessibility(analysis)

    # CLL indentity analysis
    correlation_to_bcell(analysis)
    differentiation_assessment(analysis)

    # Microenvironment analysis
    cytokine_receptor_repertoire(analysis)
    cytokine_interplay(analysis)

    # Cross-data type comparison
    scrna_comparison(analysis)

    # Off-target analysis
    off_target_signature(analysis)

    # Response prediction at day 0 from ATAC-seq
    pre_treatment_correlation_with_response(analysis)


def get_chrom_state_file():
    from ngs_toolkit.utils import download_gzip_file
    import pandas as pd

    url = (
        "https://egg2.wustl.edu/roadmap/data/byFileType/"
        + "chromhmmSegmentations/ChmmModels/coreMarks/jointModel/"
        + "final/E032_15_coreMarks_dense.bed.gz"
    )
    chrom_state_file = os.path.abspath("E032_15_coreMarks_hg19_dense.bed")
    download_gzip_file(url, chrom_state_file)

    # Test
    assert os.path.exists(chrom_state_file)
    assert os.stat(chrom_state_file).st_size > 0
    return chrom_state_file


def print_name(function):
    """
    Decorator to print the name of the decorated function.
    """
    import inspect

    def wrapper(obj, *args, **kwargs):
        print(function.__name__)
        return function(obj, *args, **kwargs)
    return wrapper


def data_normalization(analysis):
    """
    Perform normalization of a chromatin accessibility coverage matrix
    with quantile normalization followed by PCA-based batch correction and
    cube-root tranformation.
    """
    from ngs_toolkit.utils import normalize_quantiles_r

    counts_qnorm = pd.DataFrame(
        normalize_quantiles_r(analysis.matrix_raw.values),
        index=analysis.matrix_raw.index,
        columns=analysis.matrix_raw.columns
    )
    # normalize batch effect
    cell_types = list(set([sample.cell_type for sample in analysis.samples]))
    counts_qnorm_pcafix = subtract_principal_component_by_attribute(
        counts_qnorm.T, attributes=cell_types[:-1], pcs=[1] * 5).T
    # cube-root transform
    sign = (counts_qnorm_pcafix >= 0).astype(int).replace(0, -1)
    return sign * np.absolute(counts_qnorm_pcafix) ** (1 / 3.)


def subtract_principal_component_by_attribute(df, pcs=[1], attributes=["ATAC-seq"]):
    """
    Given a matrix `df` (n_samples, n_variables), remove `pc{i}` (1-based) from matrix.
    """
    from sklearn.decomposition import PCA

    assert len(pcs) == len(attributes), "Length of PCs and attributes must be identical."

    X2 = pd.DataFrame(index=df.index, columns=df.columns)
    for pc, attr in zip(pcs, attributes):
        pc -= 1
        sel = df.index[df.index.str.contains(attr)]
        print(pc, attr, sel)
        X = df.loc[sel, :]

        # PCA
        pca = PCA()
        X_hat = pca.fit_transform(X)

        # Remove PC
        X2.loc[sel, :] = X - np.outer(X_hat[:, pc], pca.components_[pc, :])
    for sample in df.index:
        if X2.loc[sample, :].isnull().all():
            X2.loc[sample, :] = df.loc[sample, :]
    return X2


def fit_gaussian_processes(
        matrix,
        matrix_file,
        cell_types=None,
        linear=False, randomize=False,
        prefix="accessibiliy",
        output_dir="/scratch/users/arendeiro/cll-time_course/gp_fit_job/",
        chunks=2000, total_job_lim=800, refresh_time=10, in_between_time=0.01,
        partition="longq", cpus=2, mem=4000,
        library="GPy"):
    """
    Estimate temporal variability of regulatory elements by comparing
    the fit of a Gaussian Process (GP) regression model with a
    variable kernel with and another with a static kernel.

    This is done in parallel across jobs in a cluster with the job manager slurm.
    """
    def count_jobs_running(cmd="squeue", sep="\n"):
        """
        Count running jobs on a cluster by invoquing a command that lists the jobs.
        """
        import subprocess
        return subprocess.check_output(cmd).split(sep).__len__()

    def submit_job_if_possible(cmd, total_job_lim=800, refresh_time=10, in_between_time=5):
        """
        If less than `total_job_lim` jobs are running, submit, else, try again in `refresh_time` seconds.
        """
        import time
        import os

        submit = count_jobs_running() < total_job_lim
        while not submit:
            time.sleep(refresh_time)
            submit = count_jobs_running() < total_job_lim
        os.system(cmd)
        time.sleep(in_between_time)

    if cell_types is None:
        cell_types = matrix.columns.get_level_values("cell_type").drop_duplicates()

    # setup output dirs
    fits_dir = os.path.join(output_dir, "fits")
    log_dir = os.path.join(output_dir, "log")
    for _dir in [output_dir, fits_dir, log_dir]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    # get ranges of matrix features
    r = np.arange(0, matrix.shape[0], chunks)

    for cell_type in tqdm(cell_types, desc="cell_type", dynamic_ncols=True):
        for start, end in tqdm(zip(r, r[1:]) + [(r[-1], matrix.shape[0])], desc="chunk", dynamic_ncols=True):
            range_name = "{}-{}".format(start, end)
            name = ".".join([prefix, cell_type, range_name, library])
            log = os.path.join(output_dir, "log", name + ".log")
            job = " ".join([
                "python -u /home/arendeiro/jobs/gp_fit_job.py",
                "--data-range {}".format(range_name),
                "--range-delim -",
                "--matrix-file {}".format(matrix_file),
                "--cell-type {}".format(cell_type),
                "--library {}".format(library),
                "--output-prefix {}".format(name),
                "--output-dir {}".format(fits_dir)]
                + ["--linear" if linear else ""]
                + ["--randomize" if randomize else ""])
            cmd = " ".join([
                "sbatch",
                "-J {}".format(name),
                "-o {}".format(log),
                "-p {}".format(partition),
                "-c {}".format(cpus),
                "--mem {}".format(mem),
                "--wrap",
                "'{}'".format(job)])

            if not os.path.exists(os.path.join(fits_dir, name + ".csv")):
                submit_job_if_possible(cmd, total_job_lim=total_job_lim, refresh_time=refresh_time, in_between_time=in_between_time)


def gather_gaussian_processes(
        matrix,
        matrix_file,
        cell_types=None,
        prefix="accessibility",
        output_dir="/scratch/users/arendeiro/cll-time_course/gp_fit_job/",
        chunks=2000, total_job_lim=800, refresh_time=10, in_between_time=0.01,
        partition="longq", cpus=2, mem=4000,
        library="GPy",
        permissive=False):
    """
    Collect the output of distributed Gaussian Process fitting procedure.
    """
    if cell_types is None:
        cell_types = matrix.columns.get_level_values("cell_type").drop_duplicates()

    # Collect output of parallel jobs
    fits = pd.DataFrame()
    r = np.arange(0, matrix.shape[0], chunks)[1:]
    for cell_type in tqdm(cell_types, desc="cell_type"):
        for start, end in tqdm(zip(r, r[1:]) + [(r[-1], matrix.shape[0])], desc="chunk"):
            range_name = "{}-{}".format(start, end)
            name = ".".join([prefix, cell_type, range_name, library])
            try:
                df = pd.read_csv(os.path.join(output_dir, "fits", name + ".csv"), index_col=0)
            except IOError:
                if permissive:
                    print("GP fit file not found: '{}'".format(name))
                    continue
                else:
                    raise

            df['cell_type'] = cell_type

            fits = fits.append(df)

    # correct p-values
    fits['q_value'] = np.concatenate(fits.groupby("cell_type", sort=False)['p_value'].apply(lambda x: multipletests(x, method="fdr_bh")[1]))

    return fits


def visualize_gaussian_process_fits(analysis, fits, output_dir, prefix="accessibility"):
    """
    Visualize statistics of Gaussian Process fitting in isolation as well as their relationship.
    Plot some examples of temporaly dynamic regions.
    """

    # # Visualize the relationship between the fits and parameters
    # g = sns.PairGrid(fits.drop("cell_type", axis=1).sample(n=2000))
    # g.map(plt.scatter, alpha=0.2, s=2, rasterized=True)
    # g.savefig(os.path.join(output_dir, prefix + ".fits.parameters.pairwise.all_cell_types.svg"), dpi=300, bbox_inches="tight")

    # Plot likelihood relationships
    g = sns.FacetGrid(data=fits, col="cell_type", col_wrap=2)
    g.map(plt.scatter, "White", "RBF", alpha=0.1, s=2, rasterized=True)
    g.savefig(os.path.join(output_dir, prefix + ".fits.RBF_vs_White.cell_types.svg"), dpi=300, bbox_inches="tight")

    g = sns.FacetGrid(data=fits, col="cell_type", col_wrap=2)
    g.map(plt.scatter, "White", "D", alpha=0.1, s=2, rasterized=True)
    g.savefig(os.path.join(output_dir, prefix + ".fits.D_vs_White.cell_types.svg"), dpi=300, bbox_inches="tight")

    n_cell_types = len(fits['cell_type'].drop_duplicates())
    n_row = n_col = int(np.ceil(np.sqrt(n_cell_types)))
    fig, axis = plt.subplots(n_row, n_col, figsize=(n_col * 4, n_row * 4))
    axis = axis.flatten()
    for i, cell_type in enumerate(fits['cell_type'].drop_duplicates()):
        f = fits[fits['cell_type'] == cell_type].head(60000)
        d = axis[i].scatter(f['White'], f["D"], c=f["mean_posterior_std"], cmap="BuGn", edgecolor='grey', alpha=0.5, s=5, rasterized=True)
        plt.colorbar(d, ax=axis[i], label='STD of posterior mean')
        axis[i].set_xlabel("log L(Data|Constant)")
        axis[i].set_ylabel("D statistic\n(2 * [log L(Data|Varying) - log L(Data|Constant)])")
        axis[i].set_title(cell_type)
    fig.savefig(os.path.join(output_dir, prefix + ".fits.D_vs_White.mean_posterior_std.cell_types.svg"), dpi=300, bbox_inches="tight")

    # Let's rank regions
    n_top = 6
    e = fits.sort_values("p_value").head(n_top)  # ~fits['cell_type'].str.contains("NK")
    examples = e.index
    example_ct = e['cell_type']
    example_acc = analysis.matrix_norm.loc[examples, :]
    cell_types = example_acc.columns.get_level_values("cell_type").drop_duplicates()
    example_acc['cell_type'] = example_ct

    # Plot some of the top examples
    n_col = len(cell_types)

    fig, axis = plt.subplots(n_top, n_col, figsize=(n_col * 3, n_top * 3), sharex=True, sharey=False)
    for i, region in enumerate(examples):
        for j, cell_type in enumerate(cell_types):
            samples = example_acc.columns.get_level_values("cell_type") == cell_type
            cur = example_acc.loc[
                region,
                samples
            ].drop_duplicates().squeeze()
            x = np.log2(1 + cur.index.get_level_values("timepoint").str.replace("d", "").astype(int).values)
            axis[i, j].scatter(x, cur.values, alpha=0.8, s=5)

            # Let's fit again the DPs just to demonstrate
            kernel = GPy.kern.RBF(input_dim=1) + GPy.kern.Bias(input_dim=1)

            model = GPy.models.GPRegression(x.reshape(-1, 1), cur.values.reshape(-1, 1), kernel)
            model.optimize()
            model.plot_f(ax=axis[i, j], lower=2.5, upper=97.5, legend=None, plot_density=True, plot_data=False, color="red")

    for i, ax in enumerate(axis[:, 0]):
        ax.set_ylabel(examples[i])
        # ax.set_title(example_acc.iloc[i]['cell_type'].squeeze())
    for i, ax in enumerate(axis[0, :]):
        ax.set_title(cell_types[i])
    fig.savefig(os.path.join(output_dir, prefix + ".top_variable.scatter.all_samples.svg"), dpi=300, bbox_inches="tight")


@print_name
def fit_MOHGP(
        matrix,
        matrix_file,
        fits_file,
        linear=False,
        cell_types=None,
        n_clust=4,
        prefix="accessibility",
        output_dir="/scratch/users/arendeiro/cll-time_course/mohgp_fit_job/",
        alpha=0.05,
        partition="shortq", cpus=24, mem=48000):
    """
    Cluster temporaly variable regulatory elements with
    a Hierarchical Mixture of Gaussian Processes (MOHGP).
    """
    # Fit MOHGPs in parallel jobs per cell type
    library = "GPy"

    if cell_types is None:
        cell_types = matrix.columns.get_level_values("cell_type").drop_duplicates()
    if type(n_clust) is int:
        n_clust = [n_clust] * len(cell_types)

    # setup output dirs
    log_dir = os.path.join(output_dir, "log")
    for _dir in [output_dir, log_dir]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    for cell_type, n in zip(cell_types, n_clust):
        name = ".".join([prefix, cell_type, "GPclust"])
        log = os.path.join(output_dir, "log", name + ".log")
        job = " ".join([
            "python -u ~/jobs/gpclust_job.py",
            "--matrix-file {}".format(matrix_file),
            "--fits-file {}".format(fits_file),
            "--n_clust {}".format(n),
            "--cell-type {}".format(cell_type),
            "--output-prefix {}".format(name),
            "--alpha {}".format(alpha),
            "--matrix-header-range 9",
            "--output-dir {}".format(output_dir)]
            + ["--linear" if linear else ""])
        cmd = " ".join([
            "sbatch",
            "-J {}".format(name),
            "-o {}".format(log),
            "-p {}".format(partition),
            "-c {}".format(cpus),
            "--mem {}".format(mem),
            "--wrap '{}'".format(job)])
        os.system(cmd)


@print_name
def gather_MOHGP(
        matrix, fits, prefix="accessibility",
        fits_dir="/scratch/users/arendeiro/cll-time_course/mohgp_fit_job/",
        cell_types=None, n_clust=4, alpha=0.05, posterior_threshold=0.8):
    """
    Cluster temporaly variable regulatory elements with
    a Hierarchical Mixture of Gaussian Processes (MOHGP).
    """
    if cell_types is None:
        cell_types = matrix.columns.get_level_values("cell_type").drop_duplicates()
    if type(n_clust) is int:
        n_clust = [n_clust] * cell_types.shape[0]

    assignments = pd.DataFrame()
    for cell_type, n in zip(cell_types, n_clust):
        print(cell_type)

        # Get variable regions for cell type
        variable = fits[(fits['p_value'] < alpha) & (fits['cell_type'] == cell_type)].index

        # Read in the posterior probabilities matrix (Phi) of cluster assignments
        name = ".".join([prefix, cell_type, "GPclust"])
        phi = pd.DataFrame(
            np.fromfile(os.path.join(fits_dir, name + "." + cell_type + ".mohgp.posterior_probs_phi.npy"))
            .reshape((len(variable), n)),
            index=variable)

        # Apply threshold probability to filter out some regions
        phi['assigned'] = (phi > posterior_threshold).any(axis=1)

        # Assign to cluster
        phi['cluster_assignment'] = phi.loc[:, phi.dtypes == float].apply(np.argmax, axis=1).to_frame(name="cluster").squeeze()

        ass = fits.loc[phi.index]
        ass = ass.loc[ass['cell_type'] == cell_type]

        assignments = assignments.append(ass.join(phi).reset_index(), ignore_index=True)

    return assignments.set_index("index")


@print_name
def visualize_clustering_assignments(
        matrix, assignments, prefix="accessibility",
        output_dir="results/mohgp_fit_job/"):
    """
    Visualize discovered feature clusters and observe their temporal dynamics across samples.
    """

    all_signal = pd.DataFrame()
    for cell_type in matrix.columns.get_level_values("cell_type").drop_duplicates():
        print(cell_type)

        regions = assignments[assignments['cell_type'] == cell_type].sort_values("cluster_assignment")
        regions_assigned = regions[regions['assigned'] == True]

        # Plot
        matrix2 = matrix.loc[regions.index, matrix.columns.get_level_values(
            "cell_type") == cell_type].astype(float)
        matrix2 = matrix2.T.groupby(['patient_id', "timepoint"]).mean().T
        tp = pd.Series(matrix2.columns.get_level_values("timepoint").str.replace(
            "d", "").astype(int), index=matrix2.columns).sort_values()


        background = matrix.loc[:, matrix.columns.get_level_values(
            "cell_type") == cell_type].astype(float)
        background = background.T.groupby(['patient_id', "timepoint"]).mean().T


        sample_display_names = cell_type + " - " + matrix2.loc[:, tp.index].columns.get_level_values(
            'patient_id') + ", " + matrix2.loc[:, tp.index].columns.get_level_values('timepoint')

        # all variable regions with assignments and threshold mask
        col_colors = [[plt.get_cmap("Paired")(i) for i in regions['cluster_assignment']], [
            plt.get_cmap("binary")(i) for i in regions['assigned'].astype(float)]]
        cm = sns.color_palette('inferno', tp.max() + 1)
        row_colors = [cm[t] for t in tp]

        for col_cluster, label in [(True, "clustered"), (False, "ordered")]:
            g = sns.clustermap(
                matrix2.loc[:, tp.index].T,
                col_colors=col_colors, row_colors=row_colors,
                row_cluster=False, col_cluster=col_cluster, z_score=1,
                xticklabels=False, yticklabels=sample_display_names,
                rasterized=True, figsize=(8, 0.2 * matrix2.shape[1]), metric="correlation", robust=True)
            g.ax_heatmap.set_xticklabels(
                g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(
                g.ax_heatmap.get_yticklabels(), rotation=0)
            g.ax_col_dendrogram.set_rasterized(True)
            g.savefig(os.path.join(
                output_dir, prefix + "." + cell_type +
                ".mohgp.fitted_model.clustermap.cluster_labels.{}.with_posterior_probs.svg".format(label)),
                dpi=300, bbox_inches="tight")

            # only variable and with assignments above threshold
            g = sns.clustermap(
                matrix2.loc[regions_assigned.index, tp.index].T,
                col_colors=[plt.get_cmap("Paired")(i)
                            for i in regions_assigned['cluster_assignment']], row_colors=row_colors,
                row_cluster=False, col_cluster=col_cluster, z_score=1,
                xticklabels=False, yticklabels=sample_display_names,
                rasterized=True, figsize=(8, 0.2 * matrix2.shape[1]), metric="correlation", robust=True)
            g.ax_heatmap.set_xticklabels(
                g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(
                g.ax_heatmap.get_yticklabels(), rotation=0)
            g.ax_col_dendrogram.set_rasterized(True)
            g.savefig(os.path.join(
                output_dir, prefix + "." + cell_type +
                ".mohgp.fitted_model.clustermap.cluster_labels.{}.only_posterior_above_threshold.svg".format(label)),
                dpi=300, bbox_inches="tight")

            # only variable and with assignments above threshold: mean per timepoint
            matrix_mean = matrix2.loc[
                regions_assigned.index, tp.index].T.groupby(level="timepoint").mean()
            background_mean = background.loc[
                :, tp.index].T.groupby(level="timepoint").mean()
            g = sns.clustermap(
                matrix_mean,
                col_colors=[plt.get_cmap("Paired")(i)
                            for i in regions_assigned['cluster_assignment']], row_colors=row_colors,
                row_cluster=False, col_cluster=col_cluster, z_score=1, xticklabels=False, yticklabels=True,
                rasterized=True, figsize=(8, 0.2 * matrix_mean.shape[0]), metric="correlation", robust=True)
            g.ax_heatmap.set_xticklabels(
                g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(
                g.ax_heatmap.get_yticklabels(), rotation=0)
            g.ax_col_dendrogram.set_rasterized(True)
            g.savefig(os.path.join(
                output_dir, prefix + "." + cell_type +
                ".mohgp.fitted_model.mean_acc.clustermap.cluster_labels.{}.only_posterior_above_threshold.svg".format(label)),
                dpi=300, bbox_inches="tight")

            # only variable and with assignments above threshold: mean per timepoint, mean per cluster
            cluster_matrix_mean = matrix_mean.T.join(regions_assigned[['cluster_assignment']]).groupby('cluster_assignment').mean().T
            cluster_background_mean = background_mean.mean(axis=1)
            cluster_matrix_mean_norm = cluster_matrix_mean.copy()
            for i in cluster_matrix_mean.columns:
                cluster_matrix_mean_norm.loc[:, i] = cluster_matrix_mean.loc[:, i] / cluster_background_mean

            for z_score, label2, cbar_label in [(None, "", "Mean accessibility"), (1, "z_score.", "Mean accessibility\n(Z-score)")]:
                g = sns.clustermap(
                    cluster_matrix_mean,
                    row_cluster=False, col_cluster=col_cluster, z_score=z_score, xticklabels=True, yticklabels=True,
                    rasterized=True, square=True, metric="correlation", robust=True, cbar_kws={"label": cbar_label})
                g.ax_heatmap.set_xticklabels(
                    g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
                g.ax_heatmap.set_yticklabels(
                    g.ax_heatmap.get_yticklabels(), rotation=0)
                g.ax_col_dendrogram.set_rasterized(True)
                g.savefig(os.path.join(
                    output_dir, prefix + "." + cell_type +
                    ".mohgp.fitted_model.mean_acc.clustermap.cluster_labels.{}.{}only_posterior_above_threshold.svg".format(label, label2)),
                    dpi=300, bbox_inches="tight")

                g = sns.clustermap(
                    cluster_matrix_mean_norm,
                    row_cluster=False, col_cluster=col_cluster, z_score=z_score, xticklabels=True, yticklabels=True,
                    rasterized=True, square=True, metric="correlation", robust=True, cbar_kws={"label": cbar_label})
                g.ax_heatmap.set_xticklabels(
                    g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
                g.ax_heatmap.set_yticklabels(
                    g.ax_heatmap.get_yticklabels(), rotation=0)
                g.ax_col_dendrogram.set_rasterized(True)
                g.savefig(os.path.join(
                    output_dir, prefix + "." + cell_type +
                    ".mohgp.fitted_model.mean_acc.clustermap.cluster_labels.total_norm{}.{}only_posterior_above_threshold.svg".format(label, label2)),
                    dpi=300, bbox_inches="tight")

        cluster_matrix_mean_norm['cell_type'] = cell_type
        all_signal = all_signal.append(cluster_matrix_mean_norm)

        # Cluster patterns
        # this is a mock of the MOHGP underlying posterior.
        cs = regions_assigned['cluster_assignment'].drop_duplicates(
        ).dropna().shape[0]

        fig, axis = plt.subplots(1, cs, figsize=(
            cs * 4, 1 * 4), sharex=True, sharey=True)
        for i, cluster in enumerate(regions_assigned['cluster_assignment'].drop_duplicates().sort_values()):
            cluster_regions = regions_assigned[regions_assigned['cluster_assignment'] == cluster].index

            X = matrix2.loc[cluster_regions, tp.index].mean(
                axis=0).T.reset_index()
            X['time'] = np.log2(
                1 + X["timepoint"].str.replace("d", "").astype(int).values)
            # d = X.groupby('time').mean().squeeze().to_frame(name="mean")
            # d['upper_q'] = X.groupby('time').quantile(.975)
            # d['lower_q'] = X.groupby('time').quantile(.025)

            kernel = GPy.kern.RBF(1.0, variance=0.5) + \
                GPy.kern.Bias(1.0, variance=0.05)
            m = GPy.models.GPRegression(X=X[['time']], Y=X[[0]], kernel=kernel)
            m.optimize()
            m.plot([0 - 0.5, max(X['time']) + 0.5], ax=axis[i], legend=None)

            # axis[i].set_ylim((1.5, 3.5))
            axis[i].set_title("Cluster {}\n(n = {})".format(
                1 + cluster, cluster_regions.shape[0]))
        axis[0].set_ylabel("Chromatin accessibility")
        for ax in axis:
            ax.set_xlabel("Time (log2)")
        fig.savefig(os.path.join(
            output_dir, prefix + "." + cell_type +
            ".".join([
                ".mohgp", "fitted_model", "mean_acc", "clustermap", "cluster_labels",
                "only_posterior_above_threshold", "variable", "cluster_means.svg"])), dpi=300, bbox_inches="tight")


def clusters_to_signatures(assignments):
    """
    Save gene sets in ATAC-seq dynamic clusters to disk to use with scRNA-seq data.
    """

    for i, cell_type in enumerate(assignments['cell_type'].drop_duplicates()):
        for cluster in range(5):
            ass = assignments.loc[
                (assignments['cell_type'] == cell_type) &
                (assignments['cluster_assignment'] == cluster), :].sort_values("p_value")
            if ass.empty:
                continue
            print(cell_type, cluster)

            # get genes in cluster
            all_g = analysis.coverage_annotated.loc[ass.index, "gene_name"].str.split(",").apply(pd.Series).stack().drop_duplicates()
            top_g = analysis.coverage_annotated.loc[ass.head(200).index, "gene_name"].str.split(",").apply(pd.Series).stack().drop_duplicates()
            all_g_coding = all_g[
                (~all_g.str.startswith("LINC"))
                & (~all_g.str.startswith("LOC"))
                & (~all_g.str.startswith("MIR"))].str.replace("-AS1", "")
            top_g_coding = top_g[
                (~top_g.str.startswith("LINC"))
                & (~top_g.str.startswith("LOC"))
                & (~top_g.str.startswith("MIR"))].str.replace("-AS1", "")

            all_g.to_csv(os.path.join(
                "results", "single_cell_RNA", "ATAC-seq_signatures", "atac-seq.{}.cluster_{}.gene_level.csv"
                .format(cell_type, cluster)),
                index=False)
            top_g.to_csv(os.path.join(
                "results", "single_cell_RNA", "ATAC-seq_signatures", "atac-seq.{}.cluster_{}.gene_level.top200.csv"
                .format(cell_type, cluster)),
                index=False)
            all_g_coding.to_csv(os.path.join(
                "results", "single_cell_RNA", "ATAC-seq_signatures", "atac-seq.{}.cluster_{}.gene_level.only_coding.csv"
                .format(cell_type, cluster)),
                index=False)
            top_g_coding.to_csv(os.path.join(
                "results", "single_cell_RNA", "ATAC-seq_signatures", "atac-seq.{}.cluster_{}.gene_leveltop200.only_coding.csv"
                .format(cell_type, cluster)),
                index=False)


def gp_linear(
        analysis,
        prefix="accessibility.qnorm_pcafix_cuberoot.linear_time",
        gp_output_dir="{results_dir}/gp_fits.linear_time",
        mohgp_output_dir="{results_dir}/mohgp_fits.linear_time",
        output_dir="{results_dir}/gp_fits"):
    """
    GP analysis with linear time.
    """

    if "{results_dir}" in gp_output_dir:
        gp_output_dir = gp_output_dir.format(results_dir=analysis.results_dir)
    if "{results_dir}" in mohgp_output_dir:
        mohgp_output_dir = mohgp_output_dir.format(results_dir=analysis.results_dir)

    # Try with linear time
    fit_gaussian_processes(
        analysis.matrix_norm,
        cell_types=cell_types,
        matrix_file=matrix_file,
        linear_time=True,
        prefix=prefix)  # wait for jobs to complete
    fits = gather_gaussian_processes(
        analysis.matrix_norm,
        matrix_file=matrix_file,
        prefix=prefix)
    fits.to_csv(os.path.join(gp_output_dir, prefix + ".GP_fits.all_cell_types.csv"), index=True)
    fits = pd.read_csv(os.path.join(gp_output_dir, prefix + ".GP_fits.all_cell_types.csv"), index_col=0)

    visualize_gaussian_process_fits(analysis, fits, output_dir=gp_output_dir, prefix=prefix)

    # cluster variable regulatory elements with hierarchical mixtures of GPs (MOHGP)

    for alpha in [0.01, 0.05]:
        prefix = prefix + ".p={}".format(alpha)
        fit_MOHGP(
            analysis.matrix_norm,
            linear=True,
            matrix_file=os.path.abspath(os.path.join(
                "results", analysis.name + ".accessibility.annotated_metadata.csv")),
            fits_file=os.path.join(
                gp_output_dir, prefix + ".GP_fits.all_cell_types.csv"),
            cell_types=cell_types,
            n_clust=[3, 5, 4, 4, 4, 3],
            prefix=prefix, output_dir=mohgp_output_dir, alpha=alpha)
        # wait for jobs to finish
        assignments = gather_MOHGP(
            analysis.matrix_norm, fits,
            cell_types=cell_types,
            n_clust=[3, 4, 4, 4, 4, 3],
            prefix=prefix, fits_dir=mohgp_output_dir, alpha=alpha, posterior_threshold=0.8)
        assignments.to_csv(os.path.join(mohgp_output_dir, prefix + ".GP_fits.linear_time.mohgp_clusters.csv"), index=True)

        visualize_clustering_assignments(
            analysis.matrix_norm, assignments, prefix=prefix,
            output_dir=mohgp_output_dir)


def gp_random(
        analysis,
        prefix="accessibility.qnorm_pcafix_cuberoot.random",
        output_dir="{results_dir}/gp_fits"):
    """
    GP analysis with permuted time.
    """

    if "{results_dir}" in output_dir:
        output_dir = output_dir.format(results_dir=analysis.results_dir)

    all_fits = pd.DataFrame()
    for i in range(30):
        # randomize times to test robustness
        fit_gaussian_processes(
            analysis.matrix_norm,
            cell_types=['CLL'],
            matrix_file=os.path.abspath(os.path.join(
                "results", analysis.name + ".accessibility.annotated_metadata.csv")),
            prefix=prefix + ".random_{}".format(i),
            randomize=True)  # wait for jobs to complete
        fits = gather_gaussian_processes(
            analysis.matrix_norm,
            cell_types=['CLL'],
            matrix_file=os.path.abspath(os.path.join(
                "results", analysis.name + ".accessibility.annotated_metadata.csv")),
            prefix=prefix + "_{}".format(i),
            permissive=True)
        fits['iteration'] = i
        all_fits = all_fits.append(fits.reset_index(), ignore_index=True)

    all_fits.to_csv(os.path.join(output_dir, prefix + ".GP_fits.all_cell_types.csv"), index=False)
    all_fits = pd.read_csv(os.path.join(output_dir, prefix + ".GP_fits.all_cell_types.csv"))

    fig, axis = plt.subplots(1)
    axis.scatter(all_fits["White"], all_fits["D"], alpha=0.5, s=2, rasterized=True)
    fig.savefig(os.path.join(output_dir, prefix + ".fits.D_vs_White.cell_types.svg"), dpi=300, bbox_inches="tight")

    fig, axis = plt.subplots(1, 1, figsize=(1 * 4, 1 * 4))
    d = axis.scatter(all_fits['White'], all_fits["D"], c=all_fits["mean_posterior_std"], cmap="BuGn", edgecolor='grey', alpha=0.5, s=5, rasterized=True)
    plt.colorbar(d, ax=axis, label='STD of posterior mean')
    axis.set_xlabel("log L(Data|Constant)")
    axis.set_ylabel("D statistic\n(2 * [log L(Data|Varying) - log L(Data|Constant)])")
    fig.savefig(os.path.join(output_dir, prefix + ".fits.D_vs_White.mean_posterior_std.cell_types.svg"), dpi=300, bbox_inches="tight")

    # Let's rank regions
    n_top = 6
    e = all_fits.sort_values("p_value").set_index("index").head(n_top)  # ~fits['cell_type'].str.contains("NK")
    examples = e.index
    example_ct = e['cell_type']
    example_acc = analysis.matrix_norm.loc[examples, :]
    cell_types = example_acc.columns.get_level_values("cell_type").drop_duplicates()
    example_acc['cell_type'] = example_ct

    # Plot some of the top examples
    n_col = len(cell_types)

    fig, axis = plt.subplots(n_top, n_col, figsize=(n_col * 3, n_top * 3), sharex=True, sharey=False)
    for i, region in enumerate(examples):
        for j, cell_type in enumerate(cell_types):
            samples = example_acc.columns.get_level_values("cell_type") == cell_type
            cur = example_acc.loc[
                region,
                samples
            ].drop_duplicates().squeeze()
            x = np.log2(1 + cur.index.get_level_values("timepoint").str.replace("d", "").astype(int).values)
            axis[i, j].scatter(x, cur.values, alpha=0.8, s=5)

            # Let's fit again the DPs just to demonstrate
            kernel = GPy.kern.RBF(input_dim=1) + GPy.kern.Bias(input_dim=1)

            model = GPy.models.GPRegression(x.reshape(-1, 1), cur.values.reshape(-1, 1), kernel)
            model.optimize()
            model.plot_f(ax=axis[i, j], lower=2.5, upper=97.5, legend=None, plot_density=True, plot_data=False, color="red")

    for i, ax in enumerate(axis[:, 0]):
        ax.set_ylabel(examples[i])
        # ax.set_title(example_acc.iloc[i]['cell_type'].squeeze())
    for i, ax in enumerate(axis[0, :]):
        ax.set_title(cell_types[i])
    fig.savefig(os.path.join(output_dir, prefix + ".top_variable.scatter.all_samples.svg"), dpi=300, bbox_inches="tight")

    # How recurrent is it that the same randomized reg.element is variable
    sns.distplot(-np.log10(all_fits.groupby("index")['p_value'].min()))


def cluster_stats(assignments, prefix):
    import itertools

    def overlap((a, b), func=max):
        """
        Return overlap of A and B sets as the maximum of either intersection (percentage).
        """
        return (
            func(
                len(set(a).intersection(set(b))),
                len(set(b).intersection(set(a))))
            /
            float(func(len(a), len(b)))
        ) * 100

    cluster_labels = pd.DataFrame({
        "Bcell": {0: "down", 1: "up", 2: "sine"},
        # "Bulk": {0: "down", 1: "up2", 2: "up1", 3: "other"},
        "Bulk": {0: "down", 1: "up", 2: "up", 3: "other"},
        "CD4": {0: "down", 1: "up", 2: "sine", 3: "other"},
        "CD8": {0: "up", 1: "down", 2: "sine", 3: "other"},
        "CLL": {0: "down", 1: "up", 2: "middle", 3: "sine"},
        "Mono": {0: "down", 1: "up", 2: "sine"},
        "NK": {0: "up", 1: "down"},
    })
    cluster_labels = pd.DataFrame(cluster_labels)
    cluster_labels.index.name = "cluster_assignment"
    assignments = pd.merge(
        assignments.reset_index(),
        pd.melt(cluster_labels.reset_index(),
            id_vars=['cluster_assignment'], var_name="cell_type", value_name="cluster_name")).set_index("index")

    cluster_counts = assignments.reset_index().groupby(['cell_type', 'cluster_name'])['index'].count().to_frame(name="count")
    counts_fraction = (cluster_counts / cluster_counts.groupby(level=['cell_type']).sum()) * 100.

    fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4), tight_layout=True)
    sns.factorplot(
        data=counts2.reset_index(), x="cell_type", y="count",
        hue="cluster_name", order=counts2.groupby(level="cell_type").sum().sort_values('count', ascending=False).index,
        kind="bar", ax=axis[0])
    sns.factorplot(
        data=counts_fraction.reset_index(), x="cell_type", y="count",
        hue="cluster_name", order=counts_fraction.reorder_levels([1, 0]).loc['down'].sort_values('count', ascending=False).index,
        kind="bar", ax=axis[1])
    axis[0].set_ylabel("Count")
    axis[1].set_ylabel("% of total")
    sns.despine(fig)
    fig.savefig(os.path.join("results", prefix + ".mohgp.cluster_name.counts.barplot.svg"), bbox_inches="tight")

    # Physical overlap
    a = assignments.sort_values(['cell_type', 'cluster_name'])
    g = a.groupby(['cell_type', 'cluster_name']).groups.values()
    n = map(lambda x: "_".join([str(i) for i in x]), a.groupby(['cell_type', 'cluster_name']).groups.keys())
    l = len(n)
    comb = pd.DataFrame(np.array(map(overlap, itertools.product(g, repeat=2))).reshape((l, l)))
    ns = pd.DataFrame(np.array(map(lambda x: ":".join(x), itertools.product(n, repeat=2))).reshape((l, l)))
    comb.index = [x[0] for x in ns[0].str.split(":")]
    comb.columns = [x[1] for x in ns.loc[0].str.split(":")]
    comb = comb.sort_index(axis=0).sort_index(axis=1)

    fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4), tight_layout=True)
    sns.heatmap(data=comb, cmap="inferno", cbar_kws={"label": "Percentage overlap"}, square=True, ax=axis[0])
    axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=90, fontsize="x-small", ha="left", va="center")
    axis[0].set_yticklabels(axis[0].get_yticklabels(), rotation=0, fontsize="x-small", ha="right", va="center")
    comb2 = comb.copy()
    np.fill_diagonal(comb2.values, np.nan)
    sns.heatmap(data=comb2, cmap="inferno", cbar_kws={"label": "Percentage overlap (no diagonal)"}, square=True, ax=axis[1])
    axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=90, fontsize="x-small", ha="left", va="center")
    axis[1].set_yticklabels(axis[1].get_yticklabels(), rotation=0, fontsize="x-small", ha="right", va="center")
    fig.savefig(os.path.join("results", prefix + ".mohgp.cluster_name.overlap.heatmap.svg"), dpi=300)

    return assignments


def plot_lola_enrichments(
        enrichment_table, top_n=20,
        comp_variable='comparison_name',
        output_dir="{results_dir}/cluster_enrichments"):
    import scipy

    if "{results_dir}" in output_dir:
        output_dir = output_dir.format(results_dir=analysis.results_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Each cell type separately, plot enrichments of the different clusters
    for cell_type in enrichment_table[comp_variable].drop_duplicates().str.replace("_.*", "").unique():
        enr = enrichment_table[enrichment_table[comp_variable].str.contains(cell_type)]

        # get a unique label for each lola region set
        enr["label"] = (
            enr["description"].astype(str) + ", " +
            enr["cellType"].astype(str) + ", " +
            enr["tissue"].astype(str) + ", " +
            enr["antibody"].astype(str) + ", " +
            enr["treatment"].astype(str))
        enr["label"] = (
            enr["label"]
            .str.replace("nan", "").str.replace("None", "")
            .str.replace(", , ", "").str.replace(", $", ""))

        # Replace inf values with
        enr["pValueLog"] = enr["pValueLog"].replace(
            np.inf,
            enr.loc[enr["pValueLog"] != np.inf, "pValueLog"].max()
        )
        enr = enr[~enr['pValueLog'].isnull()]

        # Normalize enrichments per dataset with Z-score prior to comparing various region sets
        for comparison_name in enr['comparison_name'].drop_duplicates():
            mask = enr['comparison_name'] == comparison_name
            enr.loc[mask, "z_p"] = scipy.stats.zscore(enr.loc[mask, "pValueLog"])

        # Plot top_n terms of each comparison in barplots
        top_data = enr.set_index("label").groupby(comp_variable)["pValueLog"].nlargest(top_n).reset_index()

        n = len(enr[comp_variable].drop_duplicates())
        n_side = int(np.ceil(np.sqrt(n)))

        # pivot table
        lola_pivot = pd.pivot_table(enr,
            values="z_p", columns=comp_variable, index="label").fillna(0)
        lola_pivot = lola_pivot.replace(np.inf, lola_pivot[lola_pivot != np.inf].max().max())

        top = enr.set_index('label').groupby(comp_variable)['pValueLog'].nlargest(top_n)
        top_terms = top.index.get_level_values('label').unique()

        # plot clustered heatmap
        shape = lola_pivot.loc[top_terms, :].shape
        g = sns.clustermap(
            lola_pivot.loc[top_terms, :], figsize=(max(6, 0.12 * shape[1]), max(6, 0.12 * shape[0])), square=True,
            cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"}, metric="correlation")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, "cluster_enrichments.{}".format(cell_type) + ".lola.cluster_specific.svg"), bbox_inches="tight", dpi=300)


    # All cell types and clusters together
    all_samples = enrichment_table.index
    no_cll = enrichment_table[~enrichment_table[comp_variable].str.contains("CLL")].index

    for mask, label in [
        (all_samples, "all_samples"),
        (no_cll, "no_cll")
    ]:
        enr = enrichment_table.copy().loc[mask]

        # get a unique label for each lola region set
        enr["label"] = (
            enr["description"].astype(str) + ", " +
            enr["cellType"].astype(str) + ", " +
            enr["tissue"].astype(str) + ", " +
            enr["antibody"].astype(str) + ", " +
            enr["treatment"].astype(str))
        enr["label"] = (
            enr["label"]
            .str.replace("nan", "").str.replace("None", "")
            .str.replace(", , ", "").str.replace(", $", ""))

        # Replace inf values with
        enr["pValueLog"] = enr["pValueLog"].replace(
            np.inf,
            enr.loc[enr["pValueLog"] != np.inf, "pValueLog"].max()
        )
        enr = enr[~enr['pValueLog'].isnull()]

        # Normalize enrichments per dataset with Z-score prior to comparing various region sets
        for comparison_name in enr['comparison_name'].drop_duplicates():
            mask = enr['comparison_name'] == comparison_name
            enr.loc[mask, "z_p"] = scipy.stats.zscore(enr.loc[mask, "pValueLog"])

        # Plot top_n terms of each comparison in barplots
        top_data = enr.set_index("label").groupby(comp_variable)["pValueLog"].nlargest(top_n).reset_index()

        n = len(enr[comp_variable].drop_duplicates())
        n_side = int(np.ceil(np.sqrt(n)))

        # pivot table
        lola_pivot = pd.pivot_table(enr,
            values="z_p", columns=comp_variable, index="label").fillna(0)
        lola_pivot = lola_pivot.replace(np.inf, lola_pivot[lola_pivot != np.inf].max().max())

        top = enr.set_index('label').groupby(comp_variable)['pValueLog'].nlargest(top_n)
        top_terms = top.index.get_level_values('label').unique()

        # plot clustered heatmap
        shape = lola_pivot.loc[top_terms, :].shape
        g = sns.clustermap(
            lola_pivot.loc[top_terms, :], figsize=(max(6, 0.12 * shape[1]), max(6, 0.12 * shape[0])), square=True,
            cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"}, metric="correlation")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, "cluster_enrichments.{}.top{}".format(label, top_n) + ".lola.cluster_specific.svg"), bbox_inches="tight", dpi=300)

        # plot sorted heatmap
        g = sns.clustermap(
            lola_pivot.loc[top_terms, :], figsize=(max(6, 0.12 * shape[1]), max(6, 0.12 * shape[0])), square=True,
            cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"}, metric="correlation", col_cluster=False)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, "cluster_enrichments.{}.top{}.sorted".format(label, top_n) + ".lola.cluster_specific.svg"), bbox_inches="tight", dpi=300)

        # plot correlation
        g = sns.clustermap(
            lola_pivot.loc[top_terms, :].corr(), figsize=(max(6, 0.12 * shape[1]), max(6, 0.12 * shape[1])), square=True,
            cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, "cluster_enrichments.{}.top{}.corr".format(label, top_n) + ".lola.cluster_specific.svg"), bbox_inches="tight", dpi=300)


def plot_enrichments(
        enrichment_table, top_n=20,
        comp_variable='comparison_name',
        output_dir="{results_dir}/cluster_enrichments"):
    import scipy

    if "{results_dir}" in output_dir:
        output_dir = output_dir.format(results_dir=analysis.results_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # enrichment_table["description"] = enrichment_table["description"].str.decode("utf-8")
    enrichment_table["log_p_value"] = (-np.log10(enrichment_table["p_value"])).replace({np.inf: 300})

    # Cell types separately
    for cell_type in enrichment_table[comp_variable].str.split("_").apply(lambda x: x[0]).unique():
        for gene_set_library in enrichment_table["gene_set_library"].unique():
            print(cell_type, gene_set_library)
            if gene_set_library == "Epigenomics_Roadmap_HM_ChIP-seq":
                continue

            enr = enrichment_table[
                (enrichment_table['gene_set_library'] == gene_set_library) &
                (enrichment_table[comp_variable].str.contains(cell_type))]

            # Normalize enrichments per dataset with Z-score prior to comparing various region sets
            for comparison_name in enr[comp_variable].drop_duplicates():
                mask = enr[comp_variable] == comparison_name
                enr.loc[mask, "z_log_p_value"] = scipy.stats.zscore(enr.loc[mask, "log_p_value"])

            # Plot top_n terms of each comparison in barplots
            top_data = (
                enr[enr["gene_set_library"] == gene_set_library]
                .set_index("description")
                .groupby(comp_variable)
                ["log_p_value"]
                .nlargest(top_n)
                .reset_index())

            # pivot table
            enrichr_pivot = pd.pivot_table(
                enr[enr["gene_set_library"] == gene_set_library],
                values="z_log_p_value", columns="description", index=comp_variable).fillna(0)

            top = enr[enr["gene_set_library"] == gene_set_library].set_index('description').groupby(comp_variable)['z_log_p_value'].nlargest(top_n)
            top_terms = top.index.get_level_values('description').unique()

            # plot sorted heatmap
            g = sns.clustermap(enrichr_pivot[list(set(top_terms))].T, figsize=(max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
                cbar_kws={"label": "-log10(p-value) of enrichment\nof differential genes"}, metric="correlation", square=True, col_cluster=False)
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, "cluster_enrichments.enrichr.{}.{}.cluster_specific.sorted.svg".format(cell_type, gene_set_library)), bbox_inches="tight", dpi=300)

    # All cell types and clusters together
    for gene_set_library in enrichment_table["gene_set_library"].unique():
        print(gene_set_library)
        if gene_set_library == "Epigenomics_Roadmap_HM_ChIP-seq":
            continue

        enr = enrichment_table[enrichment_table['gene_set_library'] == gene_set_library]

        all_samples = enrichment_table.index
        no_cll = enrichment_table[~enrichment_table[comp_variable].str.contains("CLL")].index

        for mask, label in [
            (all_samples, "all_samples"),
            (no_cll, "no_cll")
        ]:
            enr2 = enrichment_table.copy().loc[mask]

            # Normalize enrichments per dataset with Z-score prior to comparing various region sets
            for comparison_name in enr2['comparison_name'].drop_duplicates():
                mask = enr2['comparison_name'] == comparison_name
                enr2.loc[mask, "z_log_p_value"] = scipy.stats.zscore(enr2.loc[mask, "log_p_value"])

            # Plot top_n terms of each comparison in barplots
            top_data = (
                enr2[enr2["gene_set_library"] == gene_set_library]
                .set_index("description")
                .groupby(comp_variable)
                ["log_p_value"]
                .nlargest(top_n)
                .reset_index())

            # pivot table
            enrichr_pivot = pd.pivot_table(
                enr2[enr2["gene_set_library"] == gene_set_library],
                values="z_log_p_value", columns="description", index=comp_variable).fillna(0)

            top = enr2[enr2["gene_set_library"] == gene_set_library].set_index('description').groupby(comp_variable)['z_log_p_value'].nlargest(top_n)
            top_terms = top.index.get_level_values('description').unique()

            # plot clustered heatmap
            shape = enrichr_pivot[list(set(top_terms))].shape
            g = sns.clustermap(enrichr_pivot[list(set(top_terms))].T, figsize=(max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
                cbar_kws={"label": "-log10(p-value) of enrichment\nof differential genes"}, metric="correlation", square=True)
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, "cluster_enrichments.all_cell_types_clusters" + ".enrichr.{}.{}.top{}.cluster_specific.svg".format(gene_set_library, label, top_n)), bbox_inches="tight", dpi=300)

            # plot sorted heatmap
            g = sns.clustermap(enrichr_pivot[list(set(top_terms))].T, figsize=(max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
                cbar_kws={"label": "-log10(p-value) of enrichment\nof differential genes"}, metric="correlation", square=True, col_cluster=False)
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, "cluster_enrichments.all_cell_types_clusters" + ".sorted.enrichr.{}.{}.top{}.cluster_specific.svg".format(gene_set_library, label, top_n)), bbox_inches="tight", dpi=300)


    # pivot table
    enr = enrichment_table

    # Normalize enrichments per dataset with Z-score prior to comparing various region sets
    for comparison_name in enr['comparison_name'].drop_duplicates():
        mask = enr['comparison_name'] == comparison_name
        enr.loc[mask, "z_log_p_value"] = scipy.stats.zscore(enr.loc[mask, "log_p_value"])

    # Plot top_n terms of each comparison in barplots
    top_data = (
        enr
        .set_index("description")
        .groupby(comp_variable)
        ["log_p_value"]
        .nlargest(top_n)
        .reset_index())

    # pivot table
    enrichr_pivot = pd.pivot_table(
        enr,
        values="z_log_p_value", columns="description", index=comp_variable).fillna(0)

    top = enr.set_index('description').groupby(comp_variable)['z_log_p_value'].nlargest(top_n)
    top_terms = top.index.get_level_values('description').unique()

    # plot correlation
    shape = enrichr_pivot.shape
    g = sns.clustermap(
        enrichr_pivot.T.corr(), square=True,
        cbar_kws={"label": "Correlation in enrichment\nof differential regions"})
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.fig.savefig(os.path.join(output_dir, "cluster_enrichments.all_cell_types_clusters.corr" + ".enrichr.svg"), bbox_inches="tight", dpi=300)


def specific_cll_enrichments():
    """
    Let's plot region- and gene-level accessibility of most prominent enriched terms in CLL
    temporaly dynamic regions.
    """
    import scipy

    def stringify(x):
        return (
            x.str.replace('[', '')
            .str.replace(']', '')
            .str.replace(' ', '')
            .str.split(',')
            .apply(pd.Series).stack()
            .drop_duplicates().values)

    def cluster_genes(X, metric="correlation"):
        import scipy
        Z = scipy.cluster.hierarchy.linkage(X, metric=metric)
        d = scipy.cluster.hierarchy.dendrogram(Z, labels=X.index, no_plot=True)
        return d['ivl']

    # Get gene level values
    g_acc = get_gene_level_accessibility(analysis)
    g_acc_red = g_acc.loc[:, g_acc.columns.get_level_values("cell_type") == "CLL"].T.groupby(level=["cell_type", "timepoint"]).mean().T

    # Load enrichments
    comp_variable ='comparison_name'
    enrichment_table = pd.read_csv(os.path.join("results_deconvolve", "GPclust", "gp_fit_job.sorted.GPclust.enrichr.csv"))
    enrichment_table = enrichment_table[(enrichment_table[comp_variable].str.contains("CLL_"))]
    enrichment_table["log_p_value"] = (-np.log10(enrichment_table["p_value"])).replace({np.inf: 300})
    enr = enrichment_table[enrichment_table["gene_set_library"].str.contains("GO_Biological_Process_201")]

    ang_genes = stringify(enr.loc[enr['description'].str.contains("angiogenesis"), "genes"])
    negreg_genes = stringify(enr.loc[enr['description'].str.contains("negative regulation"), "genes"])
    ik_genes = stringify(enr.loc[enr['description'].str.contains("I-kappaB"), "genes"])
    mapk_genes = stringify(enr.loc[enr['description'].str.contains("MAPK"), "genes"])

    # Plot heatmaps
    fig, axis = plt.subplots(1, 4, figsize=(4 * 4, 1 * 6))
    sns.heatmap(
        g_acc_red.loc[ang_genes].apply(scipy.stats.zscore, axis=1)
        .loc[cluster_genes(g_acc_red.loc[ang_genes].apply(scipy.stats.zscore, axis=1).dropna())], ax=axis[0], vmin=-2, vmax=2)
    sns.heatmap(
        g_acc_red.loc[negreg_genes].apply(scipy.stats.zscore, axis=1)
        .loc[cluster_genes(g_acc_red.loc[negreg_genes].apply(scipy.stats.zscore, axis=1).dropna())], ax=axis[1], vmin=-2, vmax=2)
    sns.heatmap(
        g_acc_red.loc[ik_genes].apply(scipy.stats.zscore, axis=1)
        .loc[cluster_genes(g_acc_red.loc[ik_genes].apply(scipy.stats.zscore, axis=1).dropna())], ax=axis[2], vmin=-2, vmax=2)
    sns.heatmap(
        g_acc_red.loc[mapk_genes].apply(scipy.stats.zscore, axis=1)
        .loc[cluster_genes(g_acc_red.loc[mapk_genes].apply(scipy.stats.zscore, axis=1).dropna())], ax=axis[3], vmin=-2, vmax=2)

    for ax in axis:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize="xx-small")
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize="xx-small")
    fig.savefig(os.path.join("results", "specific_enrichments.enrichr.CLL" + ".svg"), bbox_inches="tight", dpi=300)

    import urllib3
    from bs4 import BeautifulSoup

    targets = pd.DataFrame()
    for url in ["http://www.bu.edu/nf-kb/gene-resources/target-genes/",
                "http://www.bu.edu/nf-kb/gene-resources/human-genes/"]:
        http_pool = urllib3.connection_from_url(url)
        r = http_pool.urlopen('GET', url)
        soup = BeautifulSoup(r.data, "lxml")
        tables = soup.findAll("table")
        for table in tables:
            if table.findParent("table") is None:
                name = table.find_previous_sibling()
                df = pd.read_html(str(table), header=0)[0]
                df = df.loc[:, ~df.isnull().all()]
                df['section'] = name.getText()
                df = df.rename(columns={"Acc #": "Access #"})
                targets = targets.append(df)

    targets.to_csv(os.path.join("metadata", "nfkb_genes.csv"), index=False)

    targets = targets[~targets['Human Gene Name'].isnull()]


    fig, axis = plt.subplots(4, 4, figsize=(4 * 4, 4 * 4))
    axis = axis.flatten()
    for i, section in enumerate(targets['section'].drop_duplicates()):
        genes = targets.loc[targets['section'] == section, 'Human Gene Name'].drop_duplicates().dropna()
        x = g_acc_red.loc[genes].dropna()
        if x.shape[0] == 0:
            continue
        sns.heatmap(
            x.apply(scipy.stats.zscore, axis=1)
            .loc[cluster_genes(x.apply(scipy.stats.zscore, axis=1))], ax=axis[i], vmin=-2, vmax=2)

    for ax in axis:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize="xx-small")
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize="xx-small")


def transcription_factor_accessibility(
        analysis,
        output_dir="{results_dir}/tf_accessibility",
        bed_dir="/home/arendeiro/resources/regions/LOLACore/hg19/encode_tfbs/regions/"):
    from glob import glob
    import pybedtools

    if "{results_dir}" in output_dir:
        output_dir = output_dir.format(results_dir=analysis.results_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    tfs = [
        "NFKB",
        "BATF", "PU1", "EBF1", "IKZF1", "PAX5", "ATF2", "ATF3", "POU2F2", "NFATC1", "FOXM1", "NFIC", "MTA3", "MEF2A", "MEF2C", "TCF3",
        "RUNX3", "BCL11A", "BCL3", "JUN", "FOS", "IRF1", "IRF3", "IRF4",
        "POL2", "CTCF", "SP1", "SP2", "ELK1", "ELK4", "STAT1", "STAT3", "STAT5A",
        "GATA1", # "REST",
        "NANOG", "POU5F", # "SOX2"

        # from scRNA-seq
        "GTF2B", "MBD4", "TAF1",
    ]

    all_res = pd.DataFrame()
    for factor_name in tfs:
        print(factor_name)
        # get consensus TF regions from LOLA database
        transcription_factor_regions_sets = glob(
            os.path.join(bed_dir, "*{}*".format(factor_name.capitalize())))[:5]
        bt = pybedtools.BedTool(transcription_factor_regions_sets[0]).sort()
        for fn in transcription_factor_regions_sets[1:]:
            bt = bt.cat(pybedtools.BedTool(fn).sort().merge())
        bt = bt.merge()
        bt.saveas(os.path.join("data", "external", "TF_binding_sites" + factor_name + ".bed"))

        # get regions overlapping with TF sites
        transcription_factor_r = analysis.sites.intersect(bt, wa=True).to_dataframe(names=['chrom', 'start', 'end'])
        transcription_factor_r.index = transcription_factor_r['chrom'] + ":" + transcription_factor_r['start'].astype(str) + "-" + transcription_factor_r['end'].astype(str)
        transcription_factor_a = analysis.matrix_norm.loc[transcription_factor_r.index].dropna()

        # group regions by quantile of accessibility across all experiments
        lower = 0.0
        upper = 1.0
        n_groups = 10
        r = np.arange(lower, upper + (upper / n_groups), upper / n_groups)
        mean = transcription_factor_a.mean(axis=1)

        res = pd.DataFrame()
        for l_quantile, u_quantile in zip(r, r[1:]):
            i = mean[(mean.quantile(l_quantile) > mean) & (mean < mean.quantile(u_quantile))].index

            m = transcription_factor_a.loc[i, :].mean()
            m.index = m.index.get_level_values("sample_name")
            m['upper_quantile'] = u_quantile
            res = res.append(m, ignore_index=True)
        res = res.set_index('upper_quantile')
        i = pd.DataFrame(map(pd.Series, res.columns.str.split("_")))
        i[2] = i[2].str.replace('d', '').astype(int)
        res.columns = pd.MultiIndex.from_arrays(i[[3, 2, 1]].values.T, names=['cell_type', 'timepoint', 'patient_id'])

        res = res.sort_index(axis=1, level=['cell_type', 'timepoint'], sort_remaining=False)

        d = res.dropna().T.groupby(level=['cell_type', 'timepoint']).mean().mean(1).reset_index()
        d["transcription_factor"] = factor_name
        d["binding_sites"] = transcription_factor_a.shape[0]
        d = d.rename(columns={0: "accessibility"})
        all_res = all_res.append(d, ignore_index=True)


        g = sns.clustermap(res.dropna().T, col_cluster=False, z_score=1, rasterized=True, figsize=(res.shape[0] * 0.12, res.shape[1] * 0.12), row_cluster=False)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.savefig(os.path.join(output_dir, "{}_binding.per_quantile.sorted.svg".format(factor_name)), dpi=300, bbox_inches="tight")

        res_mean = res.dropna().T.groupby(level=['cell_type', 'timepoint']).mean()
        g = sns.clustermap(res_mean.dropna(), col_cluster=False, z_score=1, rasterized=False, square=True, row_cluster=False)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.savefig(os.path.join(output_dir, "{}_binding.per_quantile.mean_patients.sorted.svg".format(factor_name)), dpi=300, bbox_inches="tight")


    # Get "background" accessibility per cell type per timepoint
    # group regions by quantile of accessibility across all experiments
    lower = 0.0
    upper = 1.0
    n_groups = 10
    r = np.arange(lower, upper + (upper / n_groups), upper / n_groups)
    mean = analysis.matrix_norm.mean(axis=1)

    res = pd.DataFrame()
    for l_quantile, u_quantile in zip(r, r[1:]):
        i = mean[(mean.quantile(l_quantile) > mean) & (mean < mean.quantile(u_quantile))].index

        m = analysis.matrix_norm.loc[i, :].mean()
        m.index = m.index.get_level_values("sample_name")
        m['upper_quantile'] = u_quantile
        res = res.append(m, ignore_index=True)
    res = res.set_index('upper_quantile')
    i = pd.DataFrame(map(pd.Series, res.columns.str.split("_")))
    i[2] = i[2].str.replace('d', '').astype(int)
    res.columns = pd.MultiIndex.from_arrays(i[[3, 2, 1]].values.T, names=['cell_type', 'timepoint', 'patient_id'])

    res = res.sort_index(axis=1, level=['cell_type', 'timepoint'], sort_remaining=False)
    d = res.dropna().T.groupby(level=['cell_type', 'timepoint']).mean().mean(1).reset_index()
    d["transcription_factor"] = "background"
    d["binding_sites"] = analysis.matrix_norm.shape[0]
    d = d.rename(columns={0: "accessibility"})
    all_res = all_res.append(d, ignore_index=True)

    # Save
    all_res.to_csv(os.path.join(output_dir, "all_factor_binding.csv"), index=False)
    all_res = pd.read_csv(os.path.join(output_dir, "all_factor_binding.csv"))


    # Normalize to background
    for cell_type in all_res['cell_type'].drop_duplicates():
        for tf in all_res['transcription_factor'].drop_duplicates():
            s = all_res.loc[(all_res['cell_type'] == cell_type) & (all_res['transcription_factor'] == tf)].set_index(['cell_type', 'timepoint'])['accessibility']
            b = all_res.loc[(all_res['cell_type'] == cell_type) & (all_res['transcription_factor'] == "background")].set_index(['cell_type', 'timepoint'])['accessibility']
            all_res.loc[(all_res['cell_type'] == cell_type) & (all_res['transcription_factor'] == tf), 'norm_accessibility'] = (s - b).values  # ((s - b) / b.std()).values

    all_res.to_csv(os.path.join(output_dir, "all_factor_binding.normalized.csv"), index=False)
    all_res = pd.read_csv(os.path.join(output_dir, "all_factor_binding.normalized.csv"))

    # Plot
    g = sns.factorplot(data=all_res[all_res['transcription_factor'] == "background"], x="timepoint", y="accessibility", hue="cell_type")
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.cell_type_hue.lineplots.background_only.svg"), dpi=300, bbox_inches="tight")

    g = sns.lmplot(
        data=all_res[all_res['transcription_factor'] != "background"],
        x="timepoint", y="norm_accessibility",
        hue="cell_type", col="transcription_factor",
        ci=75,
        col_wrap=5, sharey=False)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.cell_type_hue.lineplots.free_y.svg"), dpi=300, bbox_inches="tight")

    g = sns.lmplot(
        data=all_res[all_res['transcription_factor'] != "background"],
        x="timepoint", y="norm_accessibility",
        hue="cell_type", col="transcription_factor",
        ci=75, col_wrap=5)
    g.savefig(os.path.join(output_dir, "_all_factor_binding.mean_patients.cell_type_hue.lineplots.svg"), dpi=300, bbox_inches="tight")


    p = pd.pivot_table(all_res[all_res['transcription_factor'] != "background"], index="transcription_factor", columns=["cell_type", "timepoint"], values="norm_accessibility")

    g = sns.clustermap(p, col_cluster=False, rasterized=False, square=True, cbar_kws={"label": "Mean accessibility"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.clustermap.sorted.svg"), dpi=300, bbox_inches="tight")
    g = sns.clustermap(p, col_cluster=False, z_score=0, rasterized=False, square=True, cbar_kws={"label": "Mean accessibility\n(Z-score)"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.clustermap.sorted.z_score.svg"), dpi=300, bbox_inches="tight")

    g = sns.clustermap((p.T - p.mean(1)).T, col_cluster=False, rasterized=False, square=True, cbar_kws={"label": "Accessibility\n(signal minus mean)"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.clustermap.sorted.minus_mean.svg"), dpi=300, bbox_inches="tight")

    g = sns.clustermap((p.T / p.mean(1)).T, col_cluster=False, rasterized=False, square=True, cbar_kws={"label": "Accessibility\n(signal over mean)"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.clustermap.sorted.over_mean.svg"), dpi=300, bbox_inches="tight")

    p_z = p.copy()
    for cell_type in p.columns.get_level_values("cell_type").unique():
        # p_z.loc[:, p_z.columns.get_level_values("cell_type") == cell_type] = scipy.stats.zscore(
        #    p_z.loc[:, p_z.columns.get_level_values("cell_type") == cell_type], axis=1)
        p_z.loc[:, p_z.columns.get_level_values("cell_type") == cell_type] = p_z.loc[:, p_z.columns.get_level_values("cell_type") == cell_type].apply(scipy.stats.zscore, axis=1)

    g = sns.clustermap(p_z, col_cluster=False, rasterized=False, square=True, cbar_kws={"label": "Mean accessibility\n(Z-score)"}, robust=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.clustermap.sorted.z_score_per_cell_type.svg"), dpi=300, bbox_inches="tight")

    g = sns.clustermap(p_z - p_z.mean(0), col_cluster=False, rasterized=False, square=True, cbar_kws={"label": "Mean accessibility\n(Z-score)"}, robust=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.clustermap.sorted.z_score_per_cell_type.minus_mean.svg"), dpi=300, bbox_inches="tight")

    for cell_type in p.columns.get_level_values("cell_type").unique():
        pp = p.loc[:, p.columns.get_level_values("cell_type") == cell_type]
        g = sns.clustermap(pp, col_cluster=False, rasterized=False, square=True, cbar_kws={"label": "Mean accessibility"})
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.clustermap.{}_only.sorted.svg".format(cell_type)), dpi=300, bbox_inches="tight")
        g = sns.clustermap(pp, col_cluster=False, z_score=0, rasterized=False, square=True, cbar_kws={"label": "Mean accessibility\n(Z-score)"}, vmin=-5, vmax=5)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.clustermap.{}_only.sorted.z_score.svg".format(cell_type)), dpi=300, bbox_inches="tight")

        pp_z = p_z.loc[:, p_z.columns.get_level_values("cell_type") == cell_type]
        g2 = sns.clustermap(pp_z - pp_z.mean(0), col_cluster=False, rasterized=False, square=True, cbar_kws={"label": "Mean accessibility\n(Z-score)"}, robust=True, metric="correlation")
        g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
        g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
        g2.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.clustermap.{}_only.sorted.z_score.minus_mean.svg".format(cell_type)), dpi=300, bbox_inches="tight")

    # plot mean TF accesibility per cell type
    tf_mean = p.T.groupby(['cell_type']).mean().T
    g = sns.clustermap(tf_mean.iloc[g2.dendrogram_row.reordered_ind, :], col_cluster=False, row_cluster=False, rasterized=False, square=True, cbar_kws={"label": "Mean accessibility\n(Z-score)"}, robust=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.clustermap.tf_global_accessibility_per_cell_type.svg"), dpi=300, bbox_inches="tight")
    g = sns.clustermap(tf_mean.iloc[g2.dendrogram_row.reordered_ind, :],
        col_cluster=False, row_cluster=False, rasterized=False, square=True, cbar_kws={"label": "Mean accessibility\n(Z-score)"}, robust=True, z_score=1, cmap="PuOr_r")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.clustermap.tf_global_accessibility_per_cell_type.zscore.svg"), dpi=300, bbox_inches="tight")

    # plot binding sites per TF
    tf_n = all_res[['transcription_factor', "binding_sites"]].drop_duplicates().set_index('transcription_factor').drop('background', axis=0)
    tf_n = tf_n.loc[tf_mean.iloc[g2.dendrogram_row.reordered_ind, :].index, :]

    g = sns.clustermap(tf_n,
        col_cluster=False, row_cluster=False, rasterized=False, square=True, cbar_kws={"label": "TFBS number per TF"}, cmap="Greens")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.clustermap.tfbs_per_tf.svg"), dpi=300, bbox_inches="tight")


    # investigate expression of TFs
    # try to get it again based on cross-patient differences
    expr = pd.read_csv(os.path.join(
        "results", "single_cell_RNA", "13_4_Overtime_nUMI_Cutoff", "SigGenes_overTime.tsv"), sep="\t")
    # expr2 = expr.loc[expr['qvalue'] < 0.05]

    m_fc = expr.loc[(expr['gene'].isin(tfs)) & (expr['cellType'] == "CLL"), :].groupby('gene')['logFC'].apply(lambda x: max(abs(x)))
    s = (expr.loc[(expr['gene'].isin(tfs)) & (expr['cellType'] == "CLL"), :].groupby('gene')['logFC'].mean() > 0).astype(int).replace(0, -1)
    fc = (m_fc * s).sort_values()
    p = expr.loc[(expr['gene'].isin(tfs)) & (expr['cellType'] == "CLL"), :].groupby('gene')['pvalue'].min().sort_values()

    fig, axis = plt.subplots(1, 1, figsize=(3, 3))
    axis.scatter(fc.loc[p.index], -np.log10(p), alpha=0.5)
    for tf in p.index:
        axis.text(fc.loc[tf], -np.log10(p.loc[tf]), tf, color="black")
    axis.set_xlim(-fc.abs().max() - 1, fc.abs().max() + 1)
    axis.set_xlabel("log2(fold-change)")
    axis.set_ylabel("-log10(min(p-value))")
    fig.savefig(os.path.join(output_dir, "scRNA-seq_expression.tfs.min_patients.signed_fold_change.svg"), dpi=300, bbox_inches="tight")



    fig, axis = plt.subplots(1)
    sns.heatmap(expr.loc[["RELB", "REL", "NFKB1", "NFKB2"]].T.groupby(['cell_type', 'timepoint']).mean(), ax=axis)
    axis.set_yticklabels(axis.get_yticklabels(), rotation=0)
    axis.set_xticklabels(axis.get_xticklabels(), rotation=90)



def correlation_to_bcell(analysis):
    cll = analysis.matrix_norm.loc[:, (analysis.matrix_norm.columns.get_level_values("cell_type") == "CLL")]
    bcell = analysis.matrix_norm.loc[:, (analysis.matrix_norm.columns.get_level_values("cell_type") == "Bcell")]

    m_cll = cll.mean(axis=1)
    m_cll.name = "CLL"
    m_bcell = bcell.mean(axis=1)
    m_bcell.name = "Bcell"
    cll.columns = cll.columns.get_level_values("sample_name")
    cll = cll.join(m_cll)
    cll = cll.join(m_bcell)

    c = cll.corr()
    c_m = pd.melt(c.reset_index(), id_vars=['index'])
    res = pd.DataFrame()
    for sample in c_m['index']:
        if sample in ['CLL', "Bcell"]: continue

        a1 = c_m.loc[(c_m['index'] == sample) & (c_m['variable'] == "Bcell"), "value"].squeeze()
        a2 = c_m.loc[(c_m['index'] == sample) & (c_m['variable'] == "CLL"), "value"].squeeze()

        res = res.append(pd.Series([a1, a2] + sample.split("_")[1:-1], index=['bcell', 'cll', 'patient_id', 'timepoint']), ignore_index=True)

    res = res.drop_duplicates()

    res['timepoint'] = res['timepoint'].str.replace("d", "").astype(int)
    res['ratio'] = res['bcell'] / res['cll']

    g = sns.factorplot(data=res, x='timepoint', y='ratio', col="patient_id")
    g.savefig(os.path.join(output_dir, "correlation_to_bcell.ratio_Bcell_CLL.lineplots.svg".format(factor_name)), dpi=300, bbox_inches="tight")

    diff = res.groupby(['patient_id']).apply(lambda x: np.log2(x.loc[x['timepoint'].argmax(), 'ratio'].squeeze() / x.loc[x['timepoint'] == 0, 'ratio'].squeeze()))
    g = sns.barplot(data=diff.reset_index(), x='patient_id', y=0)
    g.set_ylim(diff.min(), diff.max())

    fc = np.log2(bcell.mean(axis=1) / cll.mean(axis=1)).sort_values()
    g = sns.clustermap(cll.loc[fc.abs().sort_values().tail(2000).index, :], yticklabels=False, metric="correlation", z_score=1)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)

    # from scipy.optimize import linprog

    # res = np.empty((0, 2))
    # for i in tqdm(cll.index):
    #     b = bcell.loc[i, :].mean()
    #     c = cll.loc[i, :].mean()
    #     c = [b, c]
    #     A = [[1, 1]]
    #     b = [1]
    #     x0_bounds = (0, None)
    #     x1_bounds = (0, None)

    #     res = np.vstack([res, linprog(c, A_ub=A, b_ub=b, bounds=(x0_bounds, x1_bounds)).x])


def differentiation_assessment(analysis):
    # Load quantification matrix from major analysis
    accessibility = analysis.matrix_norm
    accessibility.columns = accessibility.columns.get_level_values("sample_name")

    for sample_set in ["fzhao"]:
        # Make new config file pointing to respective annotation sheet
        new_annot = "annotation.{}_samples.csv".format(sample_set)
        new_config = os.path.join("metadata", "project_config.{}_samples.yaml".format(sample_set))
        c = open(os.path.join("metadata", "project_config.yaml"), 'r').read().replace("annotation.csv", new_annot)
        with open(new_config, 'w') as handle:
            handle.write(c)

        # Start new analysis
        prj = Project(new_config)
        for sample in prj.samples:
            sample.filtered = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
        tmp_analysis = ATACSeqAnalysis(name="cll-time_course.{}_samples".format(sample_set), prj=prj, samples=prj.samples)
        # Get consensus peak set from major tmp_analysis
        tmp_analysis.set_consensus_sites(os.path.join("results", "cll-time_course_peak_set.bed"))
        # Get coverage values for each peak in each sample
        tmp_analysis.measure_coverage(tmp_analysis.samples)
        tmp_analysis.coverage.to_csv(os.path.join(
            "results", "analysis.cll-time_course.{}_samples".format(sample_set) + ".coverage.csv"))
        tmp_analysis.coverage = tmp_analysis.coverage.loc[:, ~tmp_analysis.coverage.columns.str.contains("D199")]
        # Normalize cell types jointly (quantile normalization)
        tmp_analysis.normalize_coverage_quantiles(samples=[s for s in tmp_analysis.samples if "D199" not in s.name])
        tmp_analysis.coverage_qnorm.to_csv(os.path.join(
            "results", "analysis.cll-time_course.{}_samples".format(sample_set) + ".coverage_qnorm.csv"))
        tmp_analysis.to_pickle()

        # Join matrix with CLL samples from major tmp_analysis
        q = accessibility.join(tmp_analysis.coverage_qnorm).drop(['chrom', 'start', 'end'], axis=1)

        # Compute correlation
        for region_set in ['all', 'differential', 'cll_differential']:
            if region_set == "all":
                c = q.corr()
            elif region_set == "differential":
                diff = analysis.assignments.index.unique()
                c = q.loc[diff, :].corr()
            elif region_set == "cll_differential":
                diff = analysis.assignments.loc[analysis.assignments['cell_type'] == "CLL"].index.unique()
                c = q.loc[diff, :].corr()

            g = sns.clustermap(tmp_analysis.coverage_qnorm.drop(['chrom', 'start', 'end'], axis=1).corr(), xticklabels=False, figsize=(8, 8), rasterized=True)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.savefig(os.path.join("results", "differentiation_assessment.{}_only.map.{}.svg".format(sample_set, region_set)), dpi=200)

            g = sns.clustermap(c, xticklabels=False, figsize=(8, 8), rasterized=True)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.savefig(os.path.join("results", "differentiation_assessment.{}_together.map.{}.svg".format(sample_set, region_set)), dpi=200)

            mask = c.index.str.contains("d_CLL") | (c.index.str.contains("CB") & c.index.str.contains("NaiveBcell|chedBcell|Tcell"))
            g = sns.clustermap(c.loc[mask, mask], xticklabels=False, figsize=(8, 8), rasterized=True)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.savefig(os.path.join("results", "differentiation_assessment.{}_together.relevant.map.{}.svg".format(sample_set, region_set)), dpi=200)

            c2 = pd.melt(c.loc[mask, mask].reset_index(), "index")
            annot = pd.DataFrame(map(pd.Series, c2['index'].str.split("_"))).iloc[:, [1, 2, 3]]
            annot[2] = annot[2].str.replace("d|ND", "").astype(int)
            c2.index = pd.MultiIndex.from_arrays(annot.T.values, names=["patient_id", "timepoint", "cell_type"])

            c3 = c2[c2['variable'].str.contains("ND") & ~c2.index.get_level_values("cell_type").str.contains("ND") & ~c2['index'].str.contains("ND")]
            c3['healty_subset'] = [x[3] for x in c3['variable'].str.split("_")]
            g = sns.factorplot(data=c3.reset_index(), x="timepoint", y="value", col="patient_id", hue="healty_subset", sharey=False)
            g.savefig(os.path.join("results", "differentiation_assessment.{}_correlation_to_normal.time_dependent.{}.svg".format(sample_set, region_set)), dpi=300)

            g = sns.factorplot(data=c3.reset_index(), x="timepoint", y="value", col="patient_id", hue="healty_subset", sharey=False, sharex=False)
            g.savefig(os.path.join("results", "differentiation_assessment.{}_correlation_to_normal.time_dependent.{}.svg".format(sample_set, region_set)), dpi=300)

            data = pd.pivot_table(data=c3.reset_index(), index=['patient_id', 'timepoint', 'cell_type'], columns='healty_subset', values="value")
            data_z = pd.DataFrame(zscore(data, axis=0), index=data.index, columns=data.columns)
            data_z = (data - data.mean(0)) / data.std(0)

            g = sns.factorplot(
                data=pd.melt(data_z.reset_index(), id_vars=['patient_id', 'timepoint', 'cell_type'], var_name="healty_subset"),
                x="timepoint", y="value", col="patient_id", hue="healty_subset", sharey=False, sharex=False)
            g.savefig(os.path.join("results", "differentiation_assessment.{}_correlation_to_normal.time_dependent.zscore.{}.svg".format(sample_set, region_set)), dpi=300)

            # Get ratios along time for each patient
            data['switched_to_naive'] = np.log2(data['SwitchedBcell'] / data['NaiveBcell'])
            data['switched_to_unswitched'] = np.log2(data['SwitchedBcell'] / data['UnswitchedBcell'])
            data['unswitched_to_naive'] = np.log2(data['UnswitchedBcell'] / data['NaiveBcell'])
            data['Bcell_mean'] = data[['SwitchedBcell', 'UnswitchedBcell', 'NaiveBcell']].mean(axis=1)
            data['b_to_t'] = np.log2(data['Bcell_mean'] / data['Tcell'])
            data['cd8_to_cd4'] = np.log2(data['CD8Tcell'] / data['CD4Tcell'])

            g = sns.factorplot(
                data=pd.melt(data.loc[:, data.columns.str.contains("_to_")].reset_index(), id_vars=['patient_id', 'timepoint', 'cell_type'], var_name="healty_subset"),
                x="timepoint", y="value", col="patient_id", hue="healty_subset", sharey=True, sharex=False)
            for ax in g.axes.flatten():
                ax.axhline(0, linestyle="--", color="black")
            g.savefig(os.path.join("results", "differentiation_assessment.{}_correlation_to_normal.time_dependent.ratios.fold_change.{}.svg".format(sample_set, region_set)), dpi=300)

            g = sns.factorplot(
                data=pd.melt(data.loc[:, data.columns.str.contains("_to_")].reset_index(), id_vars=['patient_id', 'timepoint', 'cell_type'], var_name="healty_subset"),
                x="timepoint", y="value", hue="healty_subset", sharey=True, sharex=False)
            for ax in g.axes.flatten():
                ax.axhline(0, linestyle="--", color="black")
            g.savefig(os.path.join("results", "differentiation_assessment.{}_correlation_to_normal.time_dependent.ratios.fold_change.mean_patients.{}.svg".format(sample_set, region_set)), dpi=300)

            from ngs_toolkit.graphics import radar_plot

            data2 = data.reset_index()

            f = radar_plot(
                data2[data2['patient_id'] != "KI"],
                subplot_var="patient_id", group_var="timepoint",
                radial_vars=["NaiveBcell", "SwitchedBcell", "UnswitchedBcell"],
                cmap="inferno", scale_to_max=False)
            f.savefig(os.path.join("results", "differentiation_assessment.{}_correlation_to_normal.time_dependent.radar_plot.{}.svg".format(sample_set, region_set)), dpi=300)

            f = radar_plot(
                data2[data2['patient_id'] != "KI"],
                subplot_var="patient_id", group_var="timepoint",
                radial_vars=["NaiveBcell", "SwitchedBcell", "UnswitchedBcell"],
                cmap="inferno", scale_to_max=True)
            f.savefig(os.path.join("results", "differentiation_assessment.{}_correlation_to_normal.time_dependent.radar_plot.scale_to_max.{}.svg".format(sample_set, region_set)), dpi=300)


def get_gene_level_accessibility(analysis):
    """
    Get gene-level measurements of chromatin accessibility.
    """
    assert hasattr(analysis, "gene_annotation")
    acc = analysis.matrix_norm.copy()

    g = analysis.gene_annotation['gene_name'].str.split(",").apply(pd.Series).stack()
    g.index = g.index.droplevel(1)
    g.name = "gene_name"
    acc2 = analysis.matrix_norm.join(g).drop("gene_name", axis=1)
    acc2.index = analysis.matrix_norm.join(g).reset_index().set_index(['index', 'gene_name']).index
    acc2.columns = analysis.matrix_norm.columns
    return acc2.groupby(level=['gene_name']).mean()


def get_cll_gene_expression():
    df = pd.read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81274/suppl/GSE81274_CLL_expression_matrix.csv.gz")
    df = df.drop("ensembl_transcript_id", axis=1).groupby("ensembl_gene_id").max()

    return np.log2((df / df.sum(axis=0)) * 1e6)


def cytokine_receptor_repertoire(
        analysis,
        output_dir="{results_dir}/cytokine_receptor"):
    """
    """
    from bioservices.kegg import KEGG
    import requests
    from ngs_toolkit.graphics import barmap

    assert hasattr(analysis, "gene_annotation")

    if "{results_dir}" in output_dir:
        output_dir = output_dir.format(results_dir=analysis.results_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get ligand-receptor info
    lr = pd.read_excel("https://images.nature.com/original/nature-assets/ncomms/2015/150722/ncomms8866/extref/ncomms8866-s3.xlsx", 1)
    lr = lr.loc[
        ~lr['Pair.Evidence'].str.contains("EXCLUDED"),
        ['Pair.Name', 'Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol']]
    lr.columns = ['pair', 'ligand', 'receptor']
    lr_genes = lr['ligand'].unique().tolist() + lr['receptor'].unique().tolist()

    # Get cell adhesion molecules (CAMs)
    # Get HGNC id to gene symbol mapping
    url_query = "".join([
        """http://grch37.ensembl.org/biomart/martservice?query=""",
        """<?xml version="1.0" encoding="UTF-8"?>""",
        """<!DOCTYPE Query>""",
        """<Query  virtualSchemaName = "default" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >""",
        """<Dataset name = "hsapiens_gene_ensembl" interface = "default" >""",
        """<Attribute name = "ensembl_gene_id" />""",
        """<Attribute name = "external_gene_name" />""",
        """<Attribute name = "hgnc_id" />""",
        """<Attribute name = "hgnc_symbol" />""",
        """</Dataset>""",
        """</Query>"""])
    req = requests.get(url_query, stream=True)
    mapping = pd.DataFrame(
        (x.strip().split(",") for x in list(req.iter_lines())),
        columns=["ensembl_gene_id", "gene_name", "hgnc_id", "hgnc_symbol"]).replace("", np.nan)

    # Get HGNC ids from genes in CAM pathway
    k = KEGG()
    k.organism = "hsa"
    k.get("hsa04514 ")
    genes = list()
    for kegg_gene in k.parse(k.get("hsa04514"))['GENE']:
        g = k.parse(k.get("hsa:" + kegg_gene))['DBLINKS']
        if 'HGNC' in g:
            genes.append(g['HGNC'])

    cam_genes = mapping.loc[mapping['hgnc_id'].isin(genes), 'hgnc_symbol'].tolist()

    # Bring gene annotation as part of the accessibility index
    acc = analysis.matrix_norm.copy()
    acc.index = analysis.matrix_norm.join(analysis.gene_annotation).reset_index().set_index(['index', 'gene_name']).index

    # order by cell type, timepoint, patient
    acc = acc.T.reset_index().sort_values(by=['cell_type', 'timepoint', 'patient_id'], axis=0).set_index(acc.columns.names).T

    # Get all reg. elements which are associated to each gene
    full_acc = acc.loc[acc.index.get_level_values("gene_name").isin(lr_genes), :]
    red_acc = full_acc.groupby(level="gene_name").mean()

    full_acc_time = full_acc.T.groupby(level=["cell_type", "timepoint"]).mean().T
    red_acc_time = red_acc.T.groupby(level=["cell_type", "timepoint"]).mean().T

    for variable, label1 in [("red", "gene_level")]:  # ("full", "region_level"),
        for ext, label2 in [("", ""), ("_time", ".timepoint_mean")]:
            for z, label3 in [(None, ""), (1, ".zscore")]:
                m = eval(variable + "_" + "acc" + ext)

                w, h = m.shape[0] * 0.01, m.shape[1] * 0.12

                g = sns.clustermap(
                    m.T, col_cluster=True, row_cluster=True, z_score=z,
                    figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
                g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.sig_only{}.clustermap{}.svg".format(label1, label2, label3)), dpi=300, bbox_inches="tight")

                g = sns.clustermap(
                    m.T,
                    col_cluster=True, row_cluster=False, z_score=z,
                    figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
                g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.sig_only{}.clustermap{}.sorted.svg".format(label1, label2, label3)), dpi=300, bbox_inches="tight")

    # Differential regulatory elements with changing ligand or receptor
    # get variable regions for all cell types
    variable_regions = analysis.assignments.index
    full_acc_sig = full_acc.loc[full_acc.index.get_level_values("index").isin(variable_regions.tolist())]
    red_acc_sig = full_acc_sig.groupby(level="gene_name").mean()
    full_acc_time_sig = full_acc_sig.T.groupby(level=["cell_type", "timepoint"]).mean().T
    red_acc_time_sig = red_acc_sig.T.groupby(level=["cell_type", "timepoint"]).mean().T

    for variable, label1 in [("full", "region_level"), ("red", "gene_level")]:
        for ext, label2 in [("", ""), ("_time", ".timepoint_mean")]:
            for z, label3 in [(None, ""), (1, ".zscore")]:
                print(label1, label2, label3)
                m = eval(variable + "_" + "acc" + ext + "_sig")

                w, h = m.shape[0] * 0.01, m.shape[1] * 0.12

                # g = sns.clustermap(
                #     m.T, col_cluster=True, row_cluster=True, z_score=z,
                #     figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
                # g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
                # g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.sig_only{}.clustermap{}.svg".format(label1, label2, label3)), dpi=300, bbox_inches="tight")

                # g = sns.clustermap(
                #     m.T, col_cluster=True, row_cluster=False, z_score=z,
                #     figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
                # g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
                # g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.sig_only{}.clustermap{}.sorted.svg".format(label1, label2, label3)), dpi=300, bbox_inches="tight")

                # Without CLL
                m2 = m.loc[:, m.columns.get_level_values("cell_type") != "CLL"]
                g = sns.clustermap(
                    m2.T, col_cluster=True, row_cluster=True, z_score=z,
                    figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
                g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.sig_only{}.clustermap{}.noCLL.svg".format(label1, label2, label3)), dpi=300, bbox_inches="tight")

                g = sns.clustermap(
                    m2.T, col_cluster=True, row_cluster=False, z_score=z,
                    figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
                g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.sig_only{}.clustermap{}.sorted.noCLL.svg".format(label1, label2, label3)), dpi=300, bbox_inches="tight")


    # for each cell type
    cmap_mean = plt.get_cmap("PuOr_r")
    for cell_type in [u'Bcell', u'CD4', u'CD8', u'Mono', u'NK']:
        variable_regions = analysis.assignments[(analysis.assignments['cell_type'] == cell_type)].index.drop_duplicates().tolist()
        for variable, label1 in [("red", "gene_level")]:  # ("full", "region_level")
            for ext, label2 in [("", ""), ("_time", ".timepoint_mean")]:
                for row_cluster, label3 in [(False, "sorted")]:  # (True, "clustered")
                    for z, label4 in [(1, ".zscore")]:  # (None, ""),
                        print(cell_type, label1, label2, label3, label4)
                        m = eval(variable + "_" + "acc" + ext + "_sig")

                        if variable == "full":
                            m = m.loc[variable_regions, m.columns.get_level_values("cell_type") == cell_type]
                        else:
                            variable_genes = full_acc.index[full_acc.index.get_level_values("index").isin(variable_regions)].get_level_values("gene_name").drop_duplicates().tolist()
                            m = m.loc[m.index.get_level_values("gene_name").isin(variable_genes), m.columns.get_level_values("cell_type") == cell_type].groupby("gene_name").mean()

                        w, h = m.shape[0] * 0.1, m.shape[1] * 1.2

                        if variable == "red":
                            square = True
                            xticklabels = True
                        else:
                            square = False
                            xticklabels = False

                        mean = m.mean(1)
                        norm = matplotlib.colors.Normalize(vmin=-abs(mean).max(), vmax=abs(mean).max())

                        g = sns.clustermap(
                            m.T, col_cluster=True, row_cluster=row_cluster, z_score=z, square=square,
                            figsize=(h, w), robust=True, xticklabels=xticklabels, rasterized=True, col_colors=cmap_mean(norm(mean)))
                        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize=h / 3.)
                        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="center", fontsize=w / 3.)
                        g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.{}.sig_only{}.clustermap.{}{}.svg".format(
                            cell_type, label1, label2, label3, label4)), dpi=300, bbox_inches="tight")

                        if variable == "red":
                            if (label4 == ".zscore") and (label3 == "sorted"):
                                fig = barmap(
                                    g.data2d, z_score=None, square=square, ylims=(-1, 1), figsize=(g.fig.get_figwidth(), g.fig.get_figheight()))
                                for ax in fig.axes:
                                    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="left", fontsize=h / 6.)
                                    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="center", fontsize=w / 3.)
                                fig.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.{}.sig_only{}.barmap.{}{}.svg".format(cell_type, label1, label2, label3, label4)), dpi=300, bbox_inches="tight")

                            # only receptors, cytokines
                            m2 = m.loc[
                                (m.index.str.startswith("CD") | m.index.str.startswith("CCL") | m.index.str.startswith("CCR") | m.index.str.startswith("CXC")
                                    | m.index.str.startswith("IL") | m.index.str.startswith("ITG") | m.index.str.startswith("IFN")
                                    | m.index.str.startswith("VCAM") | m.index.str.startswith("NCAM") | m.index.str.startswith("ICAN")) &
                                (~m.index.str.startswith("CDH")), :]

                            w, h = m2.shape[0] * 0.1, m2.shape[1] * 1.2

                            mean = m2.mean(1)
                            norm = matplotlib.colors.Normalize(vmin=-abs(mean).max(), vmax=abs(mean).max())

                            g2 = sns.clustermap(
                                m2.T, col_cluster=True, row_cluster=row_cluster, z_score=z, square=square,
                                figsize=(h, w), robust=True, xticklabels=xticklabels, rasterized=True,
                                col_colors=cmap_mean(norm(mean)))
                            g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize=h / 3.)
                            g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90, ha="center", fontsize=w)
                            g2.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.{}.sig_only{}.clustermap.{}{}.only_CD.svg".format(
                                cell_type, label1, label2, label3, label4)), dpi=300, bbox_inches="tight")

                            if (label4 == ".zscore") and (label3 == "sorted"):
                                fig = barmap(
                                    g2.data2d, z_score=None, square=square, ylims=(-1, 1), figsize=(g.fig.get_figwidth(), g.fig.get_figheight()))
                                for ax in fig.axes:
                                    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="left", fontsize=h / 6.)
                                    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="center", fontsize=w)
                                fig.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.{}.sig_only{}.barmap.{}{}.only_CD.svg".format(cell_type, label1, label2, label3, label4)), dpi=300, bbox_inches="tight")

                            # sorted alphabetically
                            g2 = sns.clustermap(
                                m2.T, col_cluster=False, row_cluster=row_cluster, z_score=z, square=square,
                                figsize=(h, w), robust=True, xticklabels=xticklabels, rasterized=True,
                                col_colors=cmap_mean(norm(mean)))
                            g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize=h / 3.)
                            g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90, ha="center", fontsize=w)
                            g2.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.{}.sig_only{}.clustermap.{}{}.only_CD.sorted_name.svg".format(
                                cell_type, label1, label2, label3, label4)), dpi=300, bbox_inches="tight")

                            if (label4 == ".zscore") and (label3 == "sorted"):
                                fig = barmap(
                                    g2.data2d, z_score=None, square=square, ylims=(-1, 1), figsize=(g.fig.get_figwidth(), g.fig.get_figheight()))
                                for ax in fig.axes:
                                    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="left", fontsize=h / 6.)
                                    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="center", fontsize=w)
                                fig.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.{}.sig_only{}.barmap.{}{}.only_CD.sorted_name.svg".format(cell_type, label1, label2, label3, label4)), dpi=300, bbox_inches="tight")


                            # only cam pathway
                            m3 = m.loc[(m.index.isin(cam_genes)), :]

                            w, h = m3.shape[0] * 0.1, m3.shape[1] * 1.2

                            mean = m3.mean(1)
                            norm = matplotlib.colors.Normalize(vmin=-abs(mean).max(), vmax=abs(mean).max())

                            g2 = sns.clustermap(
                                m3.T, col_cluster=True, row_cluster=row_cluster, z_score=z, square=square,
                                figsize=(h, w), robust=True, xticklabels=xticklabels, rasterized=True,
                                col_colors=cmap_mean(norm(mean)))
                            g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize=h / 3.)
                            g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90, ha="center")
                            g2.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.{}.sig_only{}.clustermap.{}{}.only_CAM.svg".format(
                                cell_type, label1, label2, label3, label4)), dpi=300, bbox_inches="tight")

                            if (label4 == ".zscore") and (label3 == "sorted"):
                                fig = barmap(
                                    g2.data2d, z_score=None, square=square, ylims=(-1, 1), figsize=(g.fig.get_figwidth(), g.fig.get_figheight()))
                                for ax in fig.axes:
                                    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="left", fontsize=h / 6.)
                                    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="center", fontsize=w)
                                fig.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.{}.sig_only{}.barmap.{}{}.only_CAM.svg".format(cell_type, label1, label2, label3, label4)), dpi=300, bbox_inches="tight")




        red_acc_time_sig.to_csv(os.path.join(output_dir, "ligand-receptor_repertoire.{}.gene_level.sig_only.timepoint_mean.clustermap.csv".format(cell_type)), index=True)

    cll_expr = get_cll_gene_expression()


    # Only genes belonging to CAM pathway
    red_acc_time_sig_specific = red_acc_time_sig.loc[
        (red_acc_time_sig.index.isin(cam_genes)) |
        (red_acc_time_sig.index.str.startswith("VEG"))]
    w, h = red_acc_time_sig_specific.shape[0] * 0.12, red_acc_time_sig_specific.shape[1] * 0.12
    g = sns.clustermap(
        red_acc_time_sig_specific.T, col_cluster=True, row_cluster=False, z_score=1,
        figsize=(w, h), robust=True, xticklabels=True, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="left", fontsize="xx-small")
    g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.all_cell_types.gene_level.sig_only.cam_pathway.timepoint_mean.clustermap.svg"), dpi=300, bbox_inches="tight")

    red_acc_time_sig_specific = red_acc_time_sig.loc[
        (red_acc_time_sig.index.str.startswith("TGF")) |
        (red_acc_time_sig.index.str.startswith("WNT")) |
        (red_acc_time_sig.index.str.startswith("NOTCH")) |
        (red_acc_time_sig.index.str.startswith("TNF")) |
        (red_acc_time_sig.index.str.startswith("IL")) |
        (red_acc_time_sig.index.str.startswith("TLR"))]
    w, h = red_acc_time_sig_specific.shape[0] * 0.12, red_acc_time_sig_specific.shape[1] * 0.12
    g = sns.clustermap(
        red_acc_time_sig_specific.T, col_cluster=True, row_cluster=False, z_score=1,
        figsize=(w, h), robust=True, xticklabels=True, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="left", fontsize="xx-small")
    g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.all_cell_types.gene_level.sig_only.pathways.timepoint_mean.clustermap.svg"), dpi=300, bbox_inches="tight")

    red_acc_time_sig_specific = red_acc_time_sig.loc[
        (red_acc_time_sig.index.str.startswith("SLC")) |
        (red_acc_time_sig.index.str.startswith("CD")) |
        (red_acc_time_sig.index.str.startswith("IL")) |
        (red_acc_time_sig.index.str.startswith("CXC")) |
        (red_acc_time_sig.index.str.startswith("CC")) |
        (red_acc_time_sig.index.str.startswith("CCL"))]
    w, h = red_acc_time_sig_specific.shape[0] * 0.12, red_acc_time_sig_specific.shape[1] * 0.12
    g = sns.clustermap(
        red_acc_time_sig_specific.T, col_cluster=True, row_cluster=False, z_score=1,
        figsize=(w, h), robust=True, xticklabels=True, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="left", fontsize="xx-small")
    g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.all_cell_types.gene_level.sig_only.cytokine-receptors.timepoint_mean.clustermap.svg"), dpi=300, bbox_inches="tight")



    red_acc_time_sig_specific = red_acc_time_sig.loc[
        (red_acc_time_sig.index.str.startswith("TGF")) |
        (red_acc_time_sig.index.str.startswith("WNT")) |
        (red_acc_time_sig.index.str.startswith("NOTCH")) |
        (red_acc_time_sig.index.str.startswith("TNF")) |
        (red_acc_time_sig.index.str.startswith("IL")) |
        (red_acc_time_sig.index.str.startswith("TLR")) |
        (red_acc_time_sig.index.str.startswith("SLC")) |
        (red_acc_time_sig.index.str.startswith("CD")) |
        (red_acc_time_sig.index.str.startswith("IL")) |
        (red_acc_time_sig.index.str.startswith("CXC")) |
        (red_acc_time_sig.index.str.startswith("CC")) |
        (red_acc_time_sig.index.str.startswith("CCL"))]

    m = red_acc_time_sig_specific.loc[:, red_acc_time_sig_specific.columns.get_level_values("cell_type") == "CLL"]
    g = sns.clustermap(
        red_acc_time_sig_specific.T, col_cluster=True, row_cluster=False, z_score=1,
        robust=True, xticklabels=True, rasterized=True)

    for cell_type in acc.columns.levels[3]:
        m = red_acc_time_sig_specific.loc[:, red_acc_time_sig_specific.columns.get_level_values("cell_type") == cell_type]
        # w, h = m.shape[0] * 0.3, m.shape[1] * 0.12
        fig = barmap(
            m.T.iloc[:, g.dendrogram_col.reordered_ind], z_score=0,
            figsize=(4, 4), ylims=(-3, 3))
        fig.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.gene_level.sig_only.all_molecules.timepoint_mean.barmap.allclust.svg".format(cell_type)), dpi=300, bbox_inches="tight")



    variable = analysis.assignments.index

    full_acc_sig = full_acc.loc[
        full_acc.index.get_level_values("index").isin(variable.tolist()),
        :]
    red_acc_sig = full_acc_sig.groupby(level="gene_name").mean()
    full_acc_time_sig = full_acc_sig.T.groupby(level=["cell_type", "timepoint"]).mean().T
    red_acc_time_sig = red_acc_sig.T.groupby(level=["cell_type", "timepoint"]).mean().T

    w, h = red_acc_sig.shape[0] * 0.12, red_acc_sig.shape[1] * 0.12
    g = sns.clustermap(
        red_acc_sig.T, col_cluster=True, row_cluster=True, z_score=1,
        figsize=(w, h), square=True, robust=True, xticklabels=True, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="left", fontsize="xx-small")
    g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.all_cell_types.gene_level.sig_only.clustermap.svg"), dpi=300, bbox_inches="tight")

    # same as above, by same order but average per timepoint
    g = sns.clustermap(
        red_acc_time_sig.T.iloc[:, g.dendrogram_col.reordered_ind], col_cluster=False, row_cluster=False, z_score=1,
        figsize=(w, h), square=True, robust=True, xticklabels=True, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="left", fontsize="xx-small")
    g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.all_cell_types.gene_level.sig_only.timepoint_mean.clustermap.sorted.svg".format(cell_type)), dpi=300, bbox_inches="tight")


def cytokine_interplay(analysis):
    import scipy

    # Get ligand-receptor info
    ligand_receptor = pd.read_excel("https://images.nature.com/original/nature-assets/ncomms/2015/150722/ncomms8866/extref/ncomms8866-s3.xlsx", 1)
    ligand_receptor = ligand_receptor.loc[
        ~ligand_receptor['Pair.Evidence'].str.contains("EXCLUDED"),
        ['Pair.Name', 'Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol']]
    ligand_receptor.columns = ['pair', 'ligand', 'receptor']
    ligand_receptor_genes = ligand_receptor['ligand'].unique().tolist() + ligand_receptor['receptor'].unique().tolist()

    # create receptor->ligand mapping
    receptor_ligand_dict = ligand_receptor[['receptor', 'ligand']].set_index('receptor').to_dict()['ligand']
    ligand_receptor_dict = ligand_receptor[['receptor', 'ligand']].set_index('ligand').to_dict()['receptor']

    # expression
    ligand_receptor_expr = pd.read_excel("https://images.nature.com/original/nature-assets/ncomms/2015/150722/ncomms8866/extref/ncomms8866-s5.xlsx", 1, index_col=[0, 1, 2])
    ligand_receptor_expr = np.log2(1 + ligand_receptor_expr.drop("F5.PrimaryCells.Expression_Max", axis=1).dropna())
    ligand_receptor_expr = ligand_receptor_expr.loc[~(ligand_receptor_expr == 0).all(axis=1), :]

    # Get differential regulatory elements with changing ligand or receptor
    # get variable regions for all cell types
    variable = analysis.assignements.index

    # Get variable genes from those regions
    acc = analysis.matrix_norm.copy()
    g = analysis.gene_annotation['gene_name'].str.split(",").apply(pd.Series).stack()
    g.index = g.index.droplevel(1)
    g.name = "gene_name"
    variable_genes = list(set(g.loc[g.index.isin(variable)]))

    # Get gene-level accessibility averaged per cell type per timepoint
    acc = analysis.matrix_norm.loc[variable]
    acc.index = analysis.matrix_norm.join(analysis.gene_annotation).reset_index().set_index(['index', 'gene_name']).index

    acc_mean_gene = get_gene_level_accessibility(analysis)
    acc_mean_gene = acc_mean_gene.T.groupby(['cell_type', 'timepoint']).mean().T

    expression_threshold = 0.5
    fig, axis = plt.subplots(2, 2, figsize=(4 * 2, 6 * 2))
    # Get all ligands whose receptor is expressed in CD8 or CD4 cells
    for i, cell_type in enumerate(["CD8+ T cells", "CD4+ T cells"]):
        cell_type_acc = cell_type.split("+")[0]

        expressed = ligand_receptor_expr.loc[ligand_receptor_expr[cell_type] > expression_threshold, cell_type]
        expressed_receptors = expressed[
            expressed.index.get_level_values("ApprovedSymbol")
            .isin(ligand_receptor['receptor'].drop_duplicates())].index.get_level_values("ApprovedSymbol")

        expressed_receptors = expressed_receptors[expressed_receptors.isin(variable_genes)]
        ligands = [receptor_ligand_dict[r] for r in expressed_receptors]

        # Get accessibility of ligands in CLL
        cll_ligand_acc = acc_mean_gene.loc[
            ligands,
            (acc_mean_gene.columns.get_level_values("cell_type") == "CLL") &
            (acc_mean_gene.columns.get_level_values("timepoint") <= 150)
        ]

        # Get accessibility of receptors in other cell type
        other_receptor_acc = acc_mean_gene.loc[
            expressed_receptors,
            (acc_mean_gene.columns.get_level_values("cell_type") == cell_type_acc) &
            (acc_mean_gene.columns.get_level_values("timepoint") <= 150)
        ]

        other_receptor_acc = other_receptor_acc.loc[~cll_ligand_acc.isnull().all(1).values].apply(scipy.stats.zscore, axis=1)
        cll_ligand_acc = cll_ligand_acc.dropna().apply(scipy.stats.zscore, axis=1)

        # Cluster ligands in CLL
        dist_mat = scipy.cluster.hierarchy.distance.squareform(scipy.cluster.hierarchy.distance.pdist(cll_ligand_acc, metric="correlation"))
        linkage = scipy.cluster.hierarchy.complete(np.triu(dist_mat))
        leave_order = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True)['leaves']
        cll_ligand_acc = cll_ligand_acc.iloc[leave_order]
        other_receptor_acc = other_receptor_acc.iloc[leave_order]

        # Plot the accessibility of ligands of the expressed receptors in CLL
        sns.heatmap(
            cll_ligand_acc,
            robust=True, ax=axis[i, 0])
        axis[i, 0].set_xlabel("Time", ha="center", va="center")
        axis[i, 0].set_xticklabels(axis[i, 0].get_xticklabels(), rotation=90, ha="right", fontsize="x-small")
        axis[i, 0].set_ylabel("CLL", ha="center", va="center")
        axis[i, 0].set_yticklabels(axis[i, 0].get_yticklabels(), rotation=0, ha="right", fontsize="xx-small")

        sns.heatmap(
            other_receptor_acc,
            robust=True, ax=axis[i, 1])
        axis[i, 1].set_xlabel("Time", ha="center", va="center")
        axis[i, 1].set_xticklabels(axis[i, 1].get_xticklabels(), rotation=90, ha="right", fontsize="x-small")
        axis[i, 1].set_ylabel(cell_type, ha="center", va="center")
        axis[i, 1].set_yticklabels(axis[i, 1].get_yticklabels(), rotation=0, ha="right", fontsize="xx-small")

    axis[0, 0].set_title("Ligands")
    axis[0, 1].set_title("Receptors")


def scrna_comparison(analysis):
    from ngs_toolkit.utils import normalize_quantiles_r
    output_dir = os.path.join("results", "atac_scrna_comparison")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    prefix = "accessibility.qnorm_pcafix_cuberoot.p=0.05"


    # Gene expression level
    # check signatures in scRNA-seq data
    # get scRNA-seq data grouped by cell type, patient, timepoint
    df = pd.read_csv("results/single_cell_RNA/MeanMatrix.csv", index_col=0)
    df.columns = pd.MultiIndex.from_arrays(pd.Series(df.columns).str.split("_").apply(pd.Series).values.T, names=['patient_id', 'timepoint', 'cell_type'])

    df2 = df.T.groupby(['cell_type', 'timepoint']).mean().T
    df2 = df2.loc[:,
        (~df2.columns.get_level_values("cell_type").isin(["NurseLikeCell", "NA"])) &
        (df2.columns.get_level_values("timepoint").isin(["d0", "d30"]))]
    df2 = df2.loc[~(df2.sum(axis=1) == 0), :]

    # try quantile normalization
    to_norm = df.loc[:, ~df.columns.get_level_values("cell_type").isin(["NurseLikeCell", "NA"])]
    tpm_qnorm = pd.DataFrame(
        normalize_quantiles_r(to_norm.values),
        index=to_norm.index,
        columns=to_norm.columns
    )
    analysis.scrna_cell_mean = to_norm
    analysis.scrna_cell_mean_qnorm = tpm_qnorm
    from ngs_toolkit.general import unsupervised_analysis


    # Pairwise correlations
    g = sns.clustermap(
        to_norm.corr(), xticklabels=False, annot=True,  # yticklabels=sample_display_names,
        cmap="Spectral_r", figsize=(0.2 * to_norm.shape[1], 0.2 * to_norm.shape[1]), cbar_kws={"label": "Pearson correlation"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize='xx-small')
    g.ax_heatmap.set_xlabel(None, visible=False)
    g.ax_heatmap.set_ylabel(None, visible=False)

    # Pairwise correlations
    g = sns.clustermap(
        tpm_qnorm.corr(), xticklabels=False, annot=True,  # yticklabels=sample_display_names,
        cmap="Spectral_r", figsize=(0.2 * tpm_qnorm.shape[1], 0.2 * tpm_qnorm.shape[1]), cbar_kws={"label": "Pearson correlation"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize='xx-small')
    g.ax_heatmap.set_xlabel(None, visible=False)
    g.ax_heatmap.set_ylabel(None, visible=False)

    g = sns.clustermap(tpm_qnorm.loc[['IKZF1', 'NFATC1', 'PAX5']])
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)

    p = tpm_qnorm.loc[:, tpm_qnorm.columns.get_level_values("cell_type") == "CLL"]
    sns.swarmplot(
        data=pd.melt(p.loc[['IKZF1', 'NFATC1', 'PAX5']].T.reset_index(),
                     id_vars=['patient_id', 'timepoint', 'cell_type']),
        x='variable', y='value', hue="timepoint")

    #

    mean_signatures = pd.DataFrame()
    for i, cell_type in enumerate(analysis.assignments['cell_type'].drop_duplicates()):
        for cluster in range(5):
            ass = analysis.assignments.loc[
                (analysis.assignments['cell_type'] == cell_type) &
                (analysis.assignments['cluster_assignment'] == cluster), :].sort_values("p_value")
            if ass.empty:
                continue

            # get genes in cluster
            top_g = analysis.coverage_annotated.loc[ass.head(200).index, "gene_name"].str.split(",").apply(pd.Series).stack().drop_duplicates()
            top_g_coding = top_g[(~top_g.str.startswith("LINC")) & (~top_g.str.startswith("LOC")) & (~top_g.str.startswith("MIR"))].str.replace("-AS1", "")

            sig = df2.loc[top_g_coding].mean()
            sig = sig.reset_index()
            sig['ATAC-seq_cell_type'] = cell_type
            sig['cluster'] = cluster
            mean_signatures = mean_signatures.append(sig, ignore_index=True)

    mean_signatures.to_csv(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.mean_expression.csv"), index=False)
    mean_signatures = pd.read_csv(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.mean_expression.csv"))

    # Plot signature change with time within cell types for each ATAC-seq cluster
    piv = pd.pivot_table(
        data=mean_signatures,
        index=['cell_type', 'timepoint'],
        columns=['ATAC-seq_cell_type', 'cluster'])

    # remove effect of different expression mean per timepoint/cell type
    piv2 = (piv.T / df2.mean()).T
    # remove effect of different signatures having more/less expressed genes
    piv3 = (piv2 / piv2.mean(0)) ** 2

    piv3.to_csv(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.mean_expression.normalized.csv"), index=True)

    fig, axis = plt.subplots(1, 5, figsize=(6, 3), sharey=False, gridspec_kw={'width_ratios': [5, 4, 4, 3, 2]})
    for i, cell_type in enumerate(piv3.index.levels[0]):
        p = piv3.loc[piv3.index.get_level_values("cell_type") == cell_type, piv3.columns.get_level_values("ATAC-seq_cell_type") == cell_type]
        p.columns = p.columns.droplevel(0)
        p = pd.melt(p.reset_index(), id_vars=['cell_type', 'timepoint'])
        sns.barplot(data=p, x='cluster', y='value', hue='timepoint', ax=axis[i])
        axis[i].set_ylim((0.5, p['value'].max()))
        axis[i].set_title(cell_type)
    fig.savefig(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.mean_expression.normalized.barplot.svg"), bbox_inches="tight")

    # get measurements per patient
    df2 = df.loc[
        :,
        (~df.columns.get_level_values("cell_type").isin(["NurseLikeCell", "NA"])) &
        (df.columns.get_level_values("timepoint").isin(["d0", "d30"]))]
    df2 = df2.loc[~(df2.sum(axis=1) == 0), :]


    # try to get it again based on cross-patient differences
    scrna_diff = pd.read_csv(os.path.join(
        "results", "single_cell_RNA", "13_4_Overtime_nUMI_Cutoff", "SigGenes_overTime.tsv"), sep="\t")
    scrna_diff2 = scrna_diff.loc[scrna_diff['qvalue'] < 0.05]
    scrna_diff2.groupby('cellType')['gene'].nunique()
    scrna_diff2['log_pvalue'] = -np.log10(scrna_diff2['pvalue'])
    scrna_diff2 = scrna_diff2.rename(
        columns={"logFC": "log2FoldChange", "cellType": "comparison_name", "gene": "gene_name"})



    # observe some ATAC-seq signatures at the RNA level
    fig, axis = plt.subplots(1, 5, figsize=(3 * 5, 3 * 1))
    axis = axis.flatten()
    for i, cell_type in enumerate(["CD4", "CD8", "CLL", "Mono", "NK"]):
        for cluster in range(5):
            ass = assignments.loc[
                (assignments['cell_type'] == cell_type) &
                (assignments['cluster_assignment'] == cluster), :].sort_values("p_value")
            if ass.empty:
                continue

            print(cell_type, cluster)

            # get genes in cluster
            top_g = analysis.coverage_annotated.loc[ass.head(200).index, "gene_name"].str.split(",").apply(pd.Series).stack().drop_duplicates()
            top_g_coding = top_g[(~top_g.str.startswith("LINC")) & (~top_g.str.startswith("LOC")) & (~top_g.str.startswith("MIR"))].str.replace("-AS1", "").drop_duplicates()

            # # get respective cell type values
            # p = df2.loc[top_g_coding, df2.columns.get_level_values("cell_type") == cell_type].dropna().T
            # # keep only patients with two timepoints
            # p = p[p.groupby(level=['patient_id']).count() > 1].dropna()
            # # remove genes with no variation
            # p = p.loc[:, p.sum(0) != 0]
            # p = p.loc[:, p.std(0) > 0]

            # g = sns.clustermap(p, xticklabels=True, metric="correlation", rasterized=True, row_cluster=False, figsize=(p.shape[1] * 0.05, p.shape[0] * 0.12))
            # g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            # g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=3)
            # g.savefig(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.{}.cluster_{}.clustermap.svg".format(cell_type, cluster)), bbox_inches="tight", dpi=300)

            # g = sns.clustermap(p, z_score=1, xticklabels=True, metric="correlation", rasterized=True, row_cluster=False, figsize=(p.shape[1] * 0.05, p.shape[0] * 0.12), robust=True)
            # g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            # g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=3)
            # g.savefig(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.{}.cluster_{}.clustermap.zscore.svg".format(cell_type, cluster)), bbox_inches="tight", dpi=300)


            # intersect with differential expressed
            diff_genes = scrna_diff2.loc[(scrna_diff2['comparison_name'] == cell_type) & (scrna_diff2['qvalue'] < 0.05), 'gene_name']
            top_g_coding = top_g_coding[top_g_coding.isin(diff_genes)]
            # get respective cell type values
            p = df2.loc[top_g_coding, df2.columns.get_level_values("cell_type") == cell_type].dropna().T
            # keep only patients with two timepoints
            p = p[p.groupby(level=['patient_id']).count() > 1].dropna()
            # remove genes with no variation
            p = p.loc[:, p.sum(0) != 0]
            p = p.loc[:, p.std(0) > 0]

            # sort by timepoint
            p = p.sort_index(level="timepoint", axis=0)

            if p.empty or p.shape[1] < 2:
                continue

            g = sns.clustermap(p, xticklabels=True, metric="correlation", rasterized=False, row_cluster=False, figsize=(p.shape[1] * 0.05, p.shape[0] * 0.12))
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=3)
            g.savefig(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.{}.cluster_{}.diff_rna.clustermap.svg".format(cell_type, cluster)), bbox_inches="tight", dpi=300)

            g = sns.clustermap(p, z_score=1, xticklabels=True, metric="correlation", rasterized=False, row_cluster=False, figsize=(p.shape[1] * 0.05, p.shape[0] * 0.12), robust=True)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=3)
            g.savefig(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.{}.cluster_{}.diff_rna.clustermap.zscore.svg".format(cell_type, cluster)), bbox_inches="tight", dpi=300)

            fc = np.log2(p.groupby('timepoint').mean().iloc[1] / p.groupby('timepoint').mean().iloc[0])

            axis[i].scatter(fc.rank() / fc.shape[0], fc, label=cluster, alpha=0.8)
            # sns.distplot(fc, ax=axis[i])
            axis[i].set_title(cell_type)
            axis[i].legend()
    fig.savefig(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.diff_rna.rank_fold-change.scatter.svg"), bbox_inches="tight", dpi=300)




    # Get averaged gene expression of each ATAC-seq cluster for each single cell
    from ngs_toolkit.scrna import seurat_rdata_to_pandas
    # get Seurat-normalized gene expression and metadata from R
    rdata_file = os.path.join("results", "single_cell_RNA", "10_Seurat_raw", "inclDay30_noIGHLK_negbinom", "inclDay30_noIGHLK.RData")
    expression, metadata = seurat_rdata_to_pandas(rdata_file, "pbmc")
    metadata = metadata.rename(columns={"nGene": "genes_covered", "nUMI": "umis_detected", "CellType": "assigned_cell_type", "sample": "sample_id"})
    metadata["patient_id"] = metadata['sample_id'].str.replace(r"\d", "").str.split("_").apply(lambda x: x[0])
    metadata["timepoint"] = metadata['sample_id'].str.split("_").apply(lambda x: x[-1]).str.replace("d", "").astype(int)

    cell_signatures = pd.DataFrame()
    for i, cell_type in enumerate(["CD4", "CD8", "CLL", "NK", "Mono"]):
        for cluster in range(5):
            ass = analysis.assignments.loc[
                (analysis.assignments['cell_type'] == cell_type) &
                (analysis.assignments['cluster_assignment'] == cluster), :].sort_values("p_value")
            if ass.empty:
                continue

            # get genes in cluster
            top_g = analysis.coverage_annotated.loc[ass.index, "gene_name"].str.split(",").apply(pd.Series).stack().drop_duplicates()
            top_g_coding = top_g[(~top_g.str.startswith("LINC")) & (~top_g.str.startswith("LOC")) & (~top_g.str.startswith("MIR"))].str.replace("-AS1", "")

            # get cells of respective cell type
            for patient in metadata['patient_id'].unique():
                for timepoint in metadata['timepoint'].unique():
                    print(cell_type, cluster, patient, timepoint)
                    cells = metadata.loc[
                        (metadata['assigned_cell_type'] == cell_type) &
                        (metadata['patient_id'] == patient) &
                        (metadata['timepoint'] == timepoint)].index
                    pat_cells = metadata.loc[
                        (metadata['assigned_cell_type'] == cell_type) &
                        (metadata['patient_id'] == patient) &
                        (metadata['timepoint'] == timepoint)].index
                    if len(cells) < 1:
                        continue
                    exp = expression.loc[top_g_coding, cells].dropna()
                    pat_exp = expression.loc[:, pat_cells].dropna()
                    exp = (exp.mean(axis=0) - pat_exp.mean().mean()) / pat_exp.std().std()
                    # exp = exp.mean(axis=0)

                    exp = exp.to_frame(name="signature")
                    exp['cell_type'] = cell_type
                    exp['cluster'] = cluster
                    exp['patient'] = patient
                    exp['timepoint'] = timepoint
                    cell_signatures = cell_signatures.append(exp.reset_index(), ignore_index=True)


    fig, axis = plt.subplots(3, 5, figsize=(4 * 5, 4 * 3), sharey=True)
    axis = iter(axis.flatten())
    for cell_type in ['CD4', 'CD8', 'CLL']:
        for cluster in cell_signatures['cluster'].unique():
            ax = axis.next()
            ax.set_title("{} - cluster {}".format(cell_type, cluster))
            p = cell_signatures[(cell_signatures['cell_type'] == cell_type) & (cell_signatures['cluster'] == cluster)]
            if p.empty: continue
            sns.barplot(data=p, x="patient", y='signature', hue='timepoint', ax=ax)
    fig.savefig(os.path.join(output_dir, "_scRNA-seq_signal_on_ATAC-seq_clusters.diff_atac.expression_single_cells.std_patient.barplot.by_cluster.svg"), bbox_inches="tight", dpi=300)

    fig, axis = plt.subplots(3, 4, figsize=(4 * 4, 4 * 3), sharey=True)
    axis = iter(axis.flatten())
    for cell_type in ['CD4', 'CD8', 'CLL']:
        for patient in cell_signatures['patient'].unique():
            ax = axis.next()
            ax.set_title("{} - patient {}".format(cell_type, patient))
            p = cell_signatures[(cell_signatures['cell_type'] == cell_type) & (cell_signatures['patient'] == patient)]
            sns.barplot(data=p, x="cluster", y='signature', hue='timepoint', ax=ax)
    fig.savefig(os.path.join(output_dir, "_scRNA-seq_signal_on_ATAC-seq_clusters.diff_atac.expression_single_cells.std_patient.barplot.by_patient.svg"), bbox_inches="tight", dpi=300)


    p = cell_signatures[(cell_signatures['cell_type'].isin(["CD4", "CD8", "CLL"]) & (cell_signatures['cluster'].isin([0, 1])))]
    grid = sns.FacetGrid(data=p, row="cell_type", col="cluster", hue='timepoint')
    grid.map(sns.violinplot, "patient", 'signature', alpha=0.2)
    grid.savefig(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.diff_atac.expression_single_cells.barplot.svg"), bbox_inches="tight", dpi=300)


    # # make heatmap of expression in CLL cluster regions
    # e = expression.T.join(metadata)
    # ee = e.groupby(['assigned_cell_type', 'patient_id', 'timepoint']).mean().T
    # diff_regions = analysis.assignments.loc[
    #     (analysis.assignments['cell_type'] == "CLL") &
    #     (analysis.assignments['cluster_assignment'].isin([0, 1])), :].index
    # genes = analysis.coverage_annotated.loc[diff_regions, "gene_name"].str.split(",").apply(pd.Series).stack().drop_duplicates()


    # grid = sns.clustermap(
    #     ee.loc[genes, ee.columns.get_level_values("assigned_cell_type") == "CLL"].dropna(),
    #     z_score=1, cmap="RdBu_r", metric="correlation", rasterized=True,
    #     vmax=1,
    #     robust=True,
    #     xticklabels=True, yticklabels=False
    # )
    # # grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize=4)
    # grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
    # grid.savefig(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.diff_atac.expression_mean.clustermap.svg"), bbox_inches="tight", dpi=300)




    # Get enrichments of ATAC-seq region clusters
    atac_enr = pd.read_csv(os.path.join(
        analysis.results_dir, "cluster_enrichments", prefix + ".enrichr.csv"))
    atac_enr['data_type'] = "ATAC-seq"
    atac_enr.loc[~atac_enr['comparison_name'].str.contains("_"), 'comparison_name'] = atac_enr.loc[~atac_enr['comparison_name'].str.contains("_"), 'comparison_name'] + "_all"
    atac_enr['cell_type'] = atac_enr['comparison_name'].str.split("_").apply(lambda x: x[0])
    atac_enr['cluster'] = atac_enr['comparison_name'].str.split("_").apply(lambda x: x[1])


    # Get enrichments of scRNA-seq over time
    scrna_enr = pd.read_csv(os.path.join(
        "results", "single_cell_RNA", "13_4_Overtime_nUMI_Cutoff", "All_Enrich_.tsv"), sep="\t")
    scrna_enr['data_type'] = "scRNA-seq"
    scrna_enr = scrna_enr.rename(columns={
        "database": "gene_set_library", "category": "description",
        "grp": "comparison_name",
        "pval": "p_value", "zScore": "z_score", 'combinedScore': "combined_score"})
    scrna_enr['cell_type'] = scrna_enr['comparison_name'].str.split("_").apply(lambda x: x[0])
    scrna_enr['patient'] = scrna_enr['comparison_name'].str.split("_").apply(lambda x: x[1])
    scrna_enr['timepoint'] = scrna_enr['comparison_name'].str.split("_").apply(lambda x: x[2])
    scrna_enr['direction'] = scrna_enr['comparison_name'].str.split("_").apply(lambda x: x[3])
    scrna_enr['cell_type'] = scrna_enr['cell_type'].replace("NK", "NKcell")

    # see the enrichemnts
    plot_differential_enrichment(
        scrna_enr,
        "enrichr",
        data_type="RNA-seq",
        direction_dependent=True,
        output_dir=output_dir,
        comp_variable="comparison_name",
        output_prefix="scRNA-seq.differential_enrichment",
        top_n=5)


    # try to get it again based on cross-patient differences
    scrna_diff = pd.read_csv(os.path.join(
        "results", "single_cell_RNA", "13_4_Overtime_nUMI_Cutoff", "SigGenes_overTime.tsv"), sep="\t")
    scrna_diff2 = scrna_diff.loc[scrna_diff['qvalue'] < 0.05]
    scrna_diff2.groupby('cellType')['gene'].nunique()
    scrna_diff2['log_pvalue'] = -np.log10(scrna_diff2['pvalue'])
    scrna_diff2 = scrna_diff2.rename(
        columns={"logFC": "log2FoldChange", "cellType": "comparison_name", "gene": "gene_name"})

    scrna_diff['intercept'] = 1
    analysis.expression = scrna_diff[['gene', 'intercept']].drop_duplicates().set_index("gene")

    for label1, index in [
            # ("", scrna_diff2.index),
            (".noRP", ~scrna_diff2['gene_name'].str.contains("RPL|RP-|RPS|MT-|HLA"))]:
        for label2, max_diff in [
                ("", 10000),
                (".max200pvalue", 200)]:

            prefix = "scrna_diff.cross_patient.enrichments{}{}".format(label1, label2)

            # differential_enrichment(
            #     analysis,
            #     scrna_diff2.loc[index, :].set_index("gene_name"),
            #     data_type="RNA-seq",
            #     output_dir=output_dir,
            #     output_prefix=prefix,
            #     genome="hg19",
            #     directional=True,
            #     max_diff=max_diff,
            #     sort_var="pvalue",
            #     as_jobs=True)
            # # wait for jobs to finish

            collect_differential_enrichment(
                scrna_diff2.loc[index, :].set_index("gene_name"),
                directional=True,
                data_type="RNA-seq",
                output_dir=output_dir,
                output_prefix=prefix,
                permissive=True)

            enrichment_table = pd.read_csv(os.path.join(
                output_dir, prefix + ".enrichr.csv"))

            plot_differential_enrichment(
                enrichment_table,
                "enrichr",
                data_type="RNA-seq",
                direction_dependent=True,
                output_dir=output_dir,
                comp_variable="comparison_name",
                output_prefix=prefix,
                top_n=10)

    # try to match scRNA-seq and ATAC-seq

    atac_enr['p_value'] = -np.log10(atac_enr['p_value'])
    scrna_enr['p_value'] = -np.log10(scrna_enr['p_value'])

    cell_types = (atac_enr['cell_type'].drop_duplicates()[atac_enr['cell_type'].drop_duplicates().isin(scrna_enr['cell_type'].unique().tolist())])
    gene_set_libraries = (atac_enr['gene_set_library'].drop_duplicates()[atac_enr['gene_set_library'].drop_duplicates().isin(scrna_enr['gene_set_library'].unique().tolist())])
    for cell_type in cell_types:
        at = atac_enr[atac_enr['cell_type'] == cell_type]
        sc = scrna_enr[scrna_enr['cell_type'] == cell_type]

        for gene_set_library in gene_set_libraries:
            print((cell_type, gene_set_library))
            a = at.loc[at['gene_set_library'] == gene_set_library, :]
            s = sc.loc[sc['gene_set_library'] == gene_set_library, :]

            g_s = (
                s.groupby(['description', 'timepoint', 'direction'])
                [['combined_score', 'p_value']]
                .apply(max).reset_index())
            g_s['data_type'] = "scRNA-seq"
            g_s['comparison_name'] = "all " + g_s['timepoint'] + " " + g_s['direction']
            g_s.columns = g_s.columns.get_level_values(0)

            q = pd.concat([
                a[['data_type', 'comparison_name', 'description', 'combined_score', 'p_value']],
                s[['data_type', 'comparison_name', 'description', 'combined_score', 'p_value']],
                g_s[['data_type', 'comparison_name', 'description', 'combined_score', 'p_value']]])

            for metric in ['p_value', 'combined_score']:
                p = pd.pivot_table(data=q, index="description", columns=['data_type', 'comparison_name'], values=metric)

                g = sns.clustermap(p.fillna(0).T, z_score=0, xticklabels=True, metric="correlation", rasterized=True)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=3)
                g.savefig(os.path.join(output_dir, "comparison.{}.{}.{}.clustermap.z_score.svg".format(cell_type, gene_set_library, metric)), dpi=300, bbox_inches="tight")

                # Only scRNA-seq grouped
                p = pd.pivot_table(data=q.loc[(q['comparison_name'].str.contains("all")) | (q['data_type'] == "ATAC-seq"), :], index="description", columns=['data_type', 'comparison_name'], values=metric)

                g = sns.clustermap(p.fillna(0).T, z_score=0, xticklabels=True, metric="correlation", rasterized=True)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=3)
                g.savefig(os.path.join(output_dir, "comparison.{}.{}.{}.clustermap.z_score.only_atac_scrnamean.svg".format(cell_type, gene_set_library, metric)), dpi=300, bbox_inches="tight")

                comparisons = q.loc[q['data_type'] == "ATAC-seq", "comparison_name"].drop_duplicates().sort_values()
                groups = q.loc[q['data_type'] == "scRNA-seq", "comparison_name"].drop_duplicates().sort_values()
                groups = groups[groups.str.contains("all")]

                keywords = ["BCR", "AP-1", "NFKB", "IL-", "CXCR", "p38", "TGF", "WNT"]
                n_top = 10

                fig, axis = plt.subplots(len(groups), len(comparisons), figsize=(3 * len(comparisons), 3 * len(groups)))
                for i, comparison_name in enumerate(comparisons):
                    a = q.loc[(q['data_type'] == "ATAC-seq") & (q['comparison_name'] == comparison_name)].rename(columns={metric: "ATAC-seq"})

                    for j, group in enumerate(groups):
                        b = q.loc[(q['data_type'] == "scRNA-seq") & (q['comparison_name'] == group)].rename(columns={metric: "scRNA-seq"})
                        joint = a.set_index("description")[["ATAC-seq"]].join(b.set_index("description")[["scRNA-seq"]]).dropna()
                        joint.index = joint.index.str.replace("_Homo.*", "")
                        axis[j, i].scatter(joint["ATAC-seq"], joint['scRNA-seq'], alpha=0.7)

                        j2 = joint[joint.index.str.contains("|".join(keywords))].dropna()
                        for k, row in j2.iterrows():
                            axis[j, i].text(row["ATAC-seq"], row['scRNA-seq'], row.name, fontsize="small")
                        j3 = (joint['ATAC-seq'] + joint['scRNA-seq']).sort_values().tail(n_top).index
                        for k, row in joint.loc[j3, :].iterrows():
                            axis[j, i].text(row["ATAC-seq"], row['scRNA-seq'], row.name, fontsize="small")

                for i, ax in enumerate(axis[0, :]):
                    ax.set_title(comparisons.iloc[i])
                for i, ax in enumerate(axis[:, 0]):
                    ax.set_ylabel(groups.iloc[i])
                fig.savefig(os.path.join(output_dir, "comparison.{}.{}.{}.scatter.svg".format(cell_type, gene_set_library, metric)), dpi=300, bbox_inches="tight")


        # Highlight 1:
        # TF binding loss in CLL

        # get ATAC-seq LOLA enrichments
        # curate to get per TF value across cell lines etc...
        prefix = "accessibility.qnorm_pcafix_cuberoot.p=0.05"
        lola = pd.read_csv(os.path.join(
            analysis.results_dir, "cluster_enrichments", prefix + ".lola.csv"))
        lola = lola[lola['collection'] == "encode_tfbs"]
        lola['antibody'] = (
            lola['antibody']
            .str.replace("_\(.*", "").str.replace("\(.*", "")
            .str.replace("eGFP-", "").str.replace("HA-", "")

            .str.replace("-C20", "").str.replace("-N19", "")
            .str.replace("_C9B9", "").str.replace("-4H8", "")
            .str.replace("-110", "")

            .str.replace("alpha", "").str.replace("gamma", "")

            .str.replace("-", "")
        ).str.upper()

        # get maximum across ENCODE ChIP-seq sets
        lola_red = lola.groupby(['comparison_name', 'direction', 'antibody'])['pValueLog', 'logOddsRatio'].mean()
        lola_red['data_type'] = "ATAC-seq"
        lola_red = lola_red.rename(columns={'pValueLog': 'log_p_value', 'logOddsRatio': 'combined_score'})

        # get scRNA-seq Enrichr enrichments
        label1 = ".noRP"
        label2 = ""
        prefix = "scrna_diff.cross_patient.enrichments{}{}".format(label1, label2)
        enrichment_table = pd.read_csv(os.path.join(output_dir, prefix + ".enrichr.csv"))

        # from ENCODE_ChIP-seq database
        encode = enrichment_table[enrichment_table['gene_set_library'] == 'ENCODE_TF_ChIP-seq_2015']

        # curate
        encode['tf'] = (
            encode['description']
            .str.replace("_.*", "")
            .str.replace("eGFP-", "")
            .str.replace("HA-", "")
            .str.replace("phospho.*", ""))
        # get maximum across ENCODE ChIP-seq sets
        encode['log_p_value'] = -np.log10(encode['p_value'])
        encode_red = encode.groupby(['comparison_name', 'direction', 'tf'])['log_p_value', 'combined_score'].mean()
        encode_red['data_type'] = "scRNA-seq"

        # join and compare
        joint = pd.concat([lola_red, encode_red])
        for metric in ['log_p_value', 'combined_score']:

            # scatter

            ## not directional
            joint2 = joint.groupby(['data_type', 'comparison_name', 'antibody']).mean()
            for cmap in ['YlOrBr', "summer", "coolwarm"]:
                fig, axis = plt.subplots(1, 1, figsize=(4, 4))
                piv = pd.pivot_table(data=joint2, index="antibody", columns=['data_type', 'comparison_name'], values=metric)
                top = piv[[('scRNA-seq', "CLL"), ('ATAC-seq', "CLL")]].dropna().apply(zscore, axis=0).mean(axis=1)
                axis.scatter(piv[('scRNA-seq', "CLL")], piv[('ATAC-seq', "CLL")], alpha=0.5, color=plt.get_cmap(cmap)(top))
                for i in top.sort_values().tail(20).index:
                    axis.text(piv.loc[i, ('scRNA-seq', "CLL")], piv.loc[i, ('ATAC-seq', "CLL")], i)
                axis.set_xlabel("ATAC-seq ChIP-seq database overlap ({})".format(metric))
                axis.set_ylabel("scRNA-seq gene set overlap ({})".format(metric))
                fig.savefig(os.path.join(output_dir, "lola_vs_encode_enrichments.{}.CLL_joint_non_directional.scatter.{}.svg".format(metric, cmap)), dpi=300, bbox_inches="tight")

            ## directional
            fig, axis = plt.subplots(1, 1, figsize=(4, 4))
            piv = pd.pivot_table(data=joint, index="antibody", columns=['data_type', 'comparison_name', 'direction'], values=metric)
            top = piv[[('ATAC-seq', "CLL_down", "all"), ('scRNA-seq', "CLL", "down")]].dropna().apply(zscore, axis=0).mean(axis=1)
            axis.scatter(piv[('ATAC-seq', "CLL_down", "all")], piv[('scRNA-seq', "CLL", "down")], alpha=0.5, color=plt.get_cmap("YlGnBu")(top))
            for i in top.sort_values().tail(20).index:
                axis.text(piv.loc[i, ('ATAC-seq', "CLL_down", "all")], piv.loc[i, ('scRNA-seq', "CLL", "down")], i)
            fig.savefig(os.path.join(output_dir, "lola_vs_encode_enrichments.{}.CLL_down.scatter.svg".format(metric)), dpi=300, bbox_inches="tight")

            fig, axis = plt.subplots(1, 1, figsize=(4, 4))
            piv = pd.pivot_table(data=joint, index="antibody", columns=['data_type', 'comparison_name', 'direction'], values=metric)
            top = piv[[('ATAC-seq', "CLL_up", "all"), ('scRNA-seq', "CLL", "up")]].dropna().apply(zscore, axis=0).mean(axis=1)
            axis.scatter(piv[('ATAC-seq', "CLL_up", "all")], piv[('scRNA-seq', "CLL", "up")], alpha=0.5, color=plt.get_cmap("YlOrRd")(top))
            for i in top.sort_values().tail(20).index:
                axis.text(piv.loc[i, ('ATAC-seq', "CLL_up", "all")], piv.loc[i, ('scRNA-seq', "CLL", "up")], i)
            fig.savefig(os.path.join(output_dir, "lola_vs_encode_enrichments.{}.CLL_up.scatter.svg".format(metric)), dpi=300, bbox_inches="tight")

            # barplots
            l = lola_red.loc[('CLL_down', 'all')].sort_values('log_p_value', ascending=False).head(20).reset_index()
            e = encode_red.loc[('CLL', 'down')].sort_values('log_p_value', ascending=False).head(20).reset_index()
            fig, axis = plt.subplots(1, 2, figsize=(4, 3))
            sns.barplot(data=l, x="log_p_value", y="antibody", orient="horizontal", color=sns.color_palette("colorblind")[0], ax=axis[0])
            sns.barplot(data=e, x="log_p_value", y="tf", orient="horizontal", color=sns.color_palette("colorblind")[0], ax=axis[1])
            sns.despine(fig)
            fig.savefig(os.path.join(output_dir, "lola_vs_encode_enrichments.{}.CLL_down.barplot.svg".format(metric)), dpi=300, bbox_inches="tight")

            # heatmaps
            piv = pd.pivot_table(
                data=joint.reset_index().dropna(),
                index='antibody', columns=['data_type', 'comparison_name', 'direction'],
                values=metric, fill_value=0).replace(np.inf, 600)

            for label, z_score in [("", None), ("z0", 0)]:
                g = sns.clustermap(piv, z_score=z_score, rasterized=True)
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=3)
                g.savefig(os.path.join(output_dir, "lola_vs_encode_enrichments.{}.heatmap.{}.svg".format(metric, label)), dpi=300, bbox_inches="tight")


        # Report number of enriched TFs in
        reported_tfs = dict()
        metric = "log_p_value"
        piv = pd.pivot_table(data=joint, index="antibody", columns=['data_type', 'comparison_name', 'direction'], values=metric)

        for direction in ['down', 'up']:
            # get down/upregulated in ATAC-seq down/up cluster and scRNA-seq down/upregualted genes
            piv_cll_direction = piv.loc[
                :,
                (
                    (piv.columns.get_level_values("data_type") == "ATAC-seq") &
                    (piv.columns.get_level_values("comparison_name") == "CLL_{}".format(direction)))
                |
                (
                    (piv.columns.get_level_values("data_type") == "scRNA-seq") &
                    (piv.columns.get_level_values("comparison_name") == "CLL") &
                    (piv.columns.get_level_values("direction") == direction))
                ].dropna()

            # scale between 0-1
            standard_score = lambda x: (x - x.min()) / (x.max() - x.min())
            # now report enriched the ones above a line where the sum is higher than 0.70
            reported_tfs[direction] = piv_cll_direction.loc[piv_cll_direction.apply(standard_score, axis=0).sum(axis=1) > 0.7, :].index.tolist()
        reported_tfs['both'] = [x for x in reported_tfs['down'] if x in reported_tfs['up']]


        # follow with ATAC-seq in-depth analysis
        tfs = g2.data2d.index # tf_m.iloc[g2.dendrogram_row.reordered_ind].index
        cell_type = "CLL"
        c = df2.loc[tfs, df2.columns.get_level_values("cell_type") == cell_type]
        c = c.loc[:, c.columns.get_level_values("patient_id") != "FE"]

        # expression of TFs in aggregated scRNA-seq data
        fig, axis = plt.subplots(1, 2, figsize=(2 * (c.shape[1] * 0.05), c.shape[0] * 0.12))
        sns.heatmap(c,
            xticklabels=True, rasterized=True, square=True, ax=axis[0])
        sns.heatmap(
            pd.DataFrame(scipy.stats.zscore(c, 1), index=c.index, columns=c.columns),
            xticklabels=True, yticklabels=False, rasterized=False, square=True, ax=axis[1])
        for ax in axis:
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize="xx-small")
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize="xx-small")
        fig.savefig(os.path.join(output_dir, "scRNA-seq.TF_expression.sorted_like_ATAC-seq.clustermap.svg"), bbox_inches="tight", dpi=300)

        c = c.sort_index(level='timepoint', axis=1)

        fig, axis = plt.subplots(1, 2, figsize=(2 * (c.shape[1] * 0.05), c.shape[0] * 0.12))
        sns.heatmap(c,
            xticklabels=True, rasterized=True, square=True, ax=axis[0])
        sns.heatmap(
            pd.DataFrame(scipy.stats.zscore(c, 1), index=c.index, columns=c.columns),
            xticklabels=True, yticklabels=False, rasterized=False, square=True, ax=axis[1])
        for ax in axis:
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize="xx-small")
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize="xx-small")
        fig.savefig(os.path.join(output_dir, "scRNA-seq.TF_expression.sorted_like_ATAC-seq.clustermap.sorted_timepoint.svg"), bbox_inches="tight", dpi=300)

        fig, axis = plt.subplots(1, 1, figsize=(1 * (c.shape[1] * 0.05), c.shape[0] * 0.12))
        sns.heatmap(c.mean(1).to_frame(),
            xticklabels=True, rasterized=True, square=True, ax=axis)
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0, fontsize="xx-small")
        axis.set_xticklabels(axis.get_xticklabels(), rotation=90, fontsize="xx-small")
        fig.savefig(os.path.join(output_dir, "scRNA-seq.TF_expression.sorted_like_ATAC-seq.clustermap.mean_timepoints_patients.svg"), bbox_inches="tight", dpi=300)



        # Highlight 2:
        # CD8 cells activation vs apoptosis reduction

        # from WikiPathways_2016 database
        # get ATAC enrichments
        prefix = "accessibility.qnorm_pcafix_cuberoot.p=0.05"
        atac = pd.read_csv(os.path.join(analysis.results_dir, "cluster_enrichments", prefix + ".enrichr.csv"))
        atac = atac[atac['gene_set_library'] == "NCI-Nature_2016"]
        atac['log_p_value'] = -np.log10(atac['p_value'])
        atac['description'] = atac['description'].str.replace("_Homo .*", "").str.replace("_Mus .*", "")
        atac['data_type'] = "ATAC-seq"

        # get scRNA-seq Enrichr enrichments
        label1 = ".noRP"
        label2 = ""
        prefix = "scrna_diff.cross_patient.enrichments{}{}".format(label1, label2)
        enrichment_table = pd.read_csv(os.path.join(output_dir, prefix + ".enrichr.csv"))
        rna = enrichment_table[enrichment_table['gene_set_library'] == 'NCI-Nature_2016']
        rna['log_p_value'] = -np.log10(rna['p_value'])
        rna['description'] = rna['description'].str.replace("_Homo .*", "").str.replace("_Mus .*", "")
        rna['data_type'] = "scRNA-seq"

        # barplots
        a = atac.loc[atac['comparison_name'] == 'CD8_up', :].reset_index().sort_values('log_p_value', ascending=False).head(20)
        r = rna.loc[(rna['comparison_name'] == 'CD8') & (rna['direction'] == 'up'), :].reset_index().sort_values('log_p_value', ascending=False).head(20)

        fig, axis = plt.subplots(1, 2, figsize=(4, 3))
        sns.barplot(data=a, x="log_p_value", y="description", orient="horizontal", color=sns.color_palette("colorblind")[0], ax=axis[0])
        sns.barplot(data=r, x="log_p_value", y="description", orient="horizontal", color=sns.color_palette("colorblind")[0], ax=axis[1])
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "atac_vs_rna_enrichr_enrichments.NCI-Nature.CD8_up.barplot.svg"), dpi=300, bbox_inches="tight")

        # scatter
        joint = atac.loc[atac['comparison_name'] == 'CD8_up', :].rename(columns={'log_p_value': 'ATAC-seq'}).groupby('description')[['ATAC-seq']].mean()
        joint = joint.join(rna.loc[(rna['comparison_name'] == 'CD8') & (rna['direction'] == 'up'), :].rename(columns={'log_p_value': 'scRNA-seq'}).groupby('description')[['scRNA-seq']].mean())

        fig, axis = plt.subplots(1, figsize=(3, 3))
        axis.scatter(joint['ATAC-seq'], joint['scRNA-seq'], alpha=0.5, s=10, color=plt.get_cmap("YlOrRd")(joint.dropna().mean(1)))
        for term in joint.dropna().mean(1).sort_values().tail(20).index:
            axis.text(joint['ATAC-seq'].loc[term].mean(), joint['scRNA-seq'].loc[term].mean(), term)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "atac_vs_rna_enrichr_enrichments.NCI-Nature.CD8_up.scatter.svg"), dpi=300, bbox_inches="tight")

        fig, axis = plt.subplots(1, figsize=(1, 2))
        sns.barplot(
            data=(joint.dropna().apply(scipy.stats.zscore, 0).mean(1)).sort_values(ascending=True).tail(10).reset_index(),
            x=0, y="description", orient="horizontal", palette="YlOrRd", ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "atac_vs_rna_enrichr_enrichments.NCI-Nature.CD8_up.joint_zscores.mean.barplot.svg"), dpi=300, bbox_inches="tight")

        fig, axis = plt.subplots(1, figsize=(1, 2))
        sns.barplot(
            data=(joint.dropna().apply(scipy.stats.zscore, 0).mean(1) / joint.dropna().apply(scipy.stats.zscore, 0).std(1)).sort_values(ascending=True).tail(10).reset_index(),
            x=0, y="description", orient="horizontal", palette="YlOrRd", ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "atac_vs_rna_enrichr_enrichments.NCI-Nature.CD8_up.joint_zscores.mean_over_std.barplot.svg"), dpi=300, bbox_inches="tight")

        # plot expression of some genes enriched

        df = pd.read_csv("results/single_cell_RNA/MeanMatrix.csv", index_col=0)
        df.columns = pd.MultiIndex.from_arrays(pd.Series(df.columns).str.split("_").apply(pd.Series).values.T, names=['patient_id', 'timepoint', 'cell_type'])
        # take genes of top 20 pathways together
        diff_path_genes = set(r['genes'].sum().replace('][', ', ').replace('[', '').replace(']', '').split(', '))
        # expr = scrna_diff2.loc[~scrna_diff2['gene_name'].str.contains("RPL|RP-|RPS|MT-|HLA"), :]

        diff_path_genes_expr = df.loc[
            diff_path_genes,
            (df.columns.get_level_values("cell_type") == "CD8") &
            (df.columns.get_level_values("timepoint").isin(["d0", "d30"]))].T.sort_index(level="timepoint", axis=0)
        w, h = (diff_path_genes_expr.shape[1]) * 0.05, (diff_path_genes_expr.shape[0]) * 0.12,
        g = sns.clustermap(diff_path_genes_expr, figsize=(w, h), square=True, row_cluster=False, z_score=None, metric="correlation")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)



        # Highlight 3:
        # Cytokine/receptor repertoire changes across cell types

        df = pd.read_csv("results/single_cell_RNA/MeanMatrix.csv", index_col=0)
        df.columns = pd.MultiIndex.from_arrays(pd.Series(df.columns).str.split("_").apply(pd.Series).values.T, names=['patient_id', 'timepoint', 'cell_type'])
        df2 = df.T.groupby(['cell_type', 'timepoint']).mean().T

        df2 = df2.loc[:,
            (~df2.columns.get_level_values("cell_type").isin(["NurseLikeCell", "NA"])) &
            (df2.columns.get_level_values("timepoint").isin(["d0", "d30"]))]
        df2 = df2.loc[~(df2.sum(axis=1) == 0), :]

        g = sns.clustermap(
            df2.loc[tfs, :].dropna().T,
            metric="correlation", z_score=1, row_cluster=False,
            col_colors=[
                plt.get_cmap("summer")(np.log10(1 + df.loc[tfs, :].dropna().T.min())),
                plt.get_cmap("summer")(np.log10(1 + df.loc[tfs, :].dropna().T.mean())),
                plt.get_cmap("summer")(np.log10(1 + df.loc[tfs, :].dropna().T.max()))])
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)




        # include FACS


        # Highlight 4:
        # Within and between cell type temporal dynamics



def off_target_signature(analysis):
    import itertools
    from scipy.stats import pearsonr
    import scipy

    def overlap((a, b), func=max):
        """
        Return overlap of A and B sets as the {maximum} of either intersection (percentage).
        """
        a = set(a)
        b = set(b)
        return (
            func(
                len(a.intersection(b)),
                len(b.intersection(a)))
            /
            float(func(len(a), len(b)))
        ) * 100

    # Gene expression
    expr = pd.read_csv("results/single_cell_RNA/MeanMatrix.csv", index_col=0)
    c = pd.Series(expr.columns).str.split("_").apply(pd.Series)
    c[1] = c[1].str.replace('d', '').astype(int)
    c = c[[1, 2, 0]]
    expr.columns = pd.MultiIndex.from_arrays(c.values.T, names=['timepoint', 'cell_type', 'patient_id'])
    expr = expr.sort_index(level=['timepoint', 'cell_type'], axis=1)

    # Z-score values within each cell type
    expr_z = expr.copy()
    for cell_type in expr.columns.levels[1]:
        expr_z.loc[:, expr_z.columns.get_level_values("cell_type") == cell_type] = scipy.stats.zscore(
            expr_z.loc[:, expr_z.columns.get_level_values("cell_type") == cell_type], axis=1)

    # Gene expression diff
    expr_diff = pd.read_csv(os.path.join(
        "results", "single_cell_RNA", "13_4_Overtime_nUMI_Cutoff", "SigGenes_overTime.tsv"), sep="\t").rename(
            columns={"cellType": "cell_type", "logFC": "log_fold_change"}).set_index("gene")
    expr_diff['intercept'] = 1
    expr_diff['patient_id'] = expr_diff['patient'].str.split("_").apply(lambda x: x[0])
    expr_diff['timepoint'] = expr_diff['patient'].str.split("_").apply(lambda x: x[1]).str.replace("d", "").astype(int)
    expr_diff['direction'] = (expr_diff['log_fold_change'] > 0).astype(int).replace(0, -1)

    expr_diff2 = expr_diff.loc[expr_diff['qvalue'] < 0.05]


    # plot scatter of fold changes across cell types
    expr_diff_red = expr_diff[expr_diff['timepoint'] == 30].reset_index().groupby(['cell_type', 'gene'])['log_fold_change'].mean().reset_index(level="cell_type")
    combs = list(itertools.combinations(expr_diff_red['cell_type'].unique(), 2))
    n = int(np.ceil(np.sqrt(len(combs))))
    fig, axis = plt.subplots(n, n, figsize=(n * 3, n * 3), sharex=True, sharey=True)
    axis = axis.flatten()
    for i, (c1, c2) in enumerate(combs):
        j = expr_diff_red.loc[expr_diff_red['cell_type'] == c1, "log_fold_change"].to_frame(name=c1)
        j = j.join(expr_diff_red.loc[expr_diff_red['cell_type'] == c2, "log_fold_change"].to_frame(name=c2))
        axis[i].scatter(j[c1], j[c2], alpha=0.3, s=2, rasterized=True)
        axis[i].set_xlabel(c1)
        axis[i].set_ylabel(c2)
        axis[i].axhline(0, linestyle="--", color="black", alpha=0.75, zorder=-20)
        axis[i].axvline(0, linestyle="--", color="black", alpha=0.75, zorder=-20)

        r = pearsonr(j.dropna()[c1], j.dropna()[c2])[0]
        axis[i].text(1, -1, "r = {0:.3f}".format(r))
    fig.savefig(os.path.join(analysis.results_dir, "offtarget_signature.fold_change.scatter.svg"), dpi=300, bbox_inches="tight")

    # Fold change correlations at day 30
    combs = list(itertools.permutations(expr_diff_red['cell_type'].unique(), 2))
    corrs = list()
    for i, (c1, c2) in enumerate(combs):
        j = expr_diff_red.loc[expr_diff_red['cell_type'] == c1, "log_fold_change"].to_frame(name=c1)
        j = j.join(expr_diff_red.loc[expr_diff_red['cell_type'] == c2, "log_fold_change"].to_frame(name=c2))
        r = pearsonr(j.dropna()[c1], j.dropna()[c2])[0]
        corrs.append([c1, c2, r])

    fig, axis = plt.subplots(1, 1, figsize=(4 * 1, 4), tight_layout=True)
    sns.heatmap(pd.pivot_table(data=pd.DataFrame(corrs), index=0, columns=1), cmap="RdBu_r", vmin=-0.5, vmax=0.5, ax=axis, square=True)
    fig.savefig(os.path.join(analysis.results_dir, "offtarget_signature.fold_change.correlation.heatmap.svg"), dpi=300)

    # plot number of shared genes across cell types
    piv = pd.pivot_table(data=expr_diff2, index='gene', columns='cell_type', values='direction')

    a = expr_diff2.sort_values(['cell_type', 'direction'])
    g = a.groupby(['cell_type', 'direction']).groups.values()
    n = map(lambda x: "_".join([str(i) for i in x]), a.groupby(['cell_type', 'direction']).groups.keys())
    l = len(n)
    comb = pd.DataFrame(np.array(map(overlap, itertools.product(g, repeat=2))).reshape((l, l)))
    ns = pd.DataFrame(np.array(map(lambda x: ":".join(x), itertools.product(n, repeat=2))).reshape((l, l)))
    comb.index = [x[0] for x in ns[0].str.split(":")]
    comb.columns = [x[1] for x in ns.loc[0].str.split(":")]
    comb = comb.sort_index(axis=0).sort_index(axis=1)

    fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4), tight_layout=True)
    sns.heatmap(data=comb, cmap="inferno", cbar_kws={"label": "Percentage overlap"}, square=True, ax=axis[0])
    axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=90, fontsize="x-small", ha="center", va="top")
    axis[0].set_yticklabels(axis[0].get_yticklabels(), rotation=0, fontsize="x-small", ha="right", va="center")
    comb2 = comb.copy()
    np.fill_diagonal(comb2.values, np.nan)
    sns.heatmap(data=comb2, cmap="inferno", cbar_kws={"label": "Percentage overlap (no diagonal)"}, square=True, ax=axis[1])
    axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=90, fontsize="x-small", ha="center", va="top")
    axis[1].set_yticklabels(axis[1].get_yticklabels(), rotation=0, fontsize="x-small", ha="right", va="center")
    fig.savefig(os.path.join(analysis.results_dir, "offtarget_signature.fold_change.overlap.heatmap.svg"), dpi=300)


    # get genes which are diff acrross patietns and cell types
    piv = pd.pivot_table(data=expr_diff2, index='gene', columns='cell_type', values='direction')
    common = piv.loc[(piv['CLL'] == 1) & (piv.sum(1) > 3)]
    common = piv[piv.sum(axis=1) > 1]
    common = piv[piv.sum(axis=1) < -1]
    common = piv[piv.sum(axis=1).abs() > 2]

    # piv = pd.pivot_table(data=expr_diff2, index='gene', columns=['cell_type', 'patient_id', 'timepoint'], values='direction')
    # # filter on agreement across patients
    # piv2 = piv.T.groupby(['cell_type', 'timepoint']).sum().abs()
    # piv3 = piv[(piv2 > 1).any()]
    # s = piv3.sum(1).sort_values()
    # common = s[s.abs() > 10]
    # plt.scatter(s.rank(), s, alpha=0.1, s=5)

    # common.to_csv(os.path.join("results", "offtarget_signature.csv"), index=True)
    # common.to_clipboard(index=True)


    g = sns.clustermap(expr.loc[common.index].T, xticklabels=False, z_score=1, rasterized=True, row_cluster=False)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=5)
    # g.savefig(os.path.join(analysis.results_dir, "offtarget_signature.svg"), dpi=300, bbox_inches="tight")

    # get patients which have at least two timepoints per cell type
    f = expr_z.columns.to_frame()
    f['intercept'] = 1
    c = pd.pivot_table(f, index=['cell_type', 'patient_id'], columns='timepoint')
    i = pd.melt(c[c.sum(1) > 1].reset_index(), id_vars=['cell_type', 'patient_id']).dropna()[['timepoint', 'cell_type', 'patient_id']]
    expr_z = expr_z.loc[:, pd.MultiIndex.from_arrays(i.values.T, names=['timepoint', 'cell_type', 'patient_id'])]

    g = sns.clustermap(
        expr_z.loc[common.index, ~expr_z.columns.get_level_values("cell_type").isin(["NA", "NurseLikeCell"])].dropna().T,
        xticklabels=False, z_score=None, rasterized=True, row_cluster=True, metric="correlation")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=5)

    g = sns.clustermap(
        expr_z.loc[common.index, ~expr_z.columns.get_level_values("cell_type").isin(["NA", "NurseLikeCell"])].dropna().T.sort_index(axis=0, level=['timepoint', 'cell_type']),
        xticklabels=False, z_score=None, rasterized=True, row_cluster=False, metric="correlation")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=5)

    g = sns.clustermap(
        expr_z.loc[common.index, ~expr_z.columns.get_level_values("cell_type").isin(["NA", "NurseLikeCell"])].dropna().T.sort_index(axis=0, level=['timepoint', 'cell_type']),
        xticklabels=False, z_score=None, rasterized=True, row_cluster=False, metric="correlation")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=5)


    g = sns.clustermap(
        expr_z.loc[common.index, expr_z.columns.get_level_values("cell_type").isin(["CLL", "CD4", "CD8", "Mono"])].dropna().T.sort_index(axis=0, level=['cell_type', 'timepoint']),
        xticklabels=False, z_score=None, rasterized=True, row_cluster=False, metric="correlation")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=5)


    e = expr_z.loc[common.index, expr_z.columns.get_level_values("cell_type").isin(["CLL", "CD4", "CD8", "Mono"])].dropna().T.sort_index(axis=0, level=['timepoint', 'cell_type'])
    diff = e[e.index.get_level_values('timepoint') >= 120].mean() - e[e.index.get_level_values('timepoint') == 0].mean()
    norm = matplotlib.colors.Normalize(vmin=-diff.abs().max(), vmax=diff.abs().max())
    g = sns.clustermap(
        e, xticklabels=True, z_score=None, rasterized=False, row_cluster=False, metric="correlation", square=True,
        col_colors=[plt.get_cmap("PuOr_r")(norm(diff))])
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=3)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
    g.savefig(os.path.join(analysis.results_dir, "offtarget_signature.svg"), dpi=300, bbox_inches="tight")


    g = sns.clustermap(
        e.groupby(['cell_type', 'timepoint']).mean(), xticklabels=True, z_score=None, rasterized=False, row_cluster=False, metric="correlation", square=True,
        col_colors=[plt.get_cmap("PuOr_r")(norm(diff))])
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=3)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
    g.savefig(os.path.join(analysis.results_dir, "offtarget_signature.mean.svg"), dpi=300, bbox_inches="tight")


    # come up with a ibrutinib signature score
    diff = e[e.index.get_level_values('timepoint') >= 120].mean() - e[e.index.get_level_values('timepoint') == 0].mean()
    up = diff[diff > 0].index
    down = diff[diff <= 0].index

    sig = e.loc[:, up].mean(axis=1) / e.loc[:, down].mean(axis=1)

    g = sns.clustermap(
        e.groupby(['cell_type', 'timepoint']).mean(),
        col_colors=[plt.get_cmap("RdBu_r")(diff)],
        z_score=1, cmap="RdBu_r",
        xticklabels=True, rasterized=False, row_cluster=False, metric="correlation", square=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=3)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
    g.savefig("/home/arendeiro/plot.svg", bbox_inches="tight", dpi=300)


    # Make a few plots illustrating marker expression
    genes = ['CXCR4', "CD44", "FGR", "ZFPL2"]
    cell_types = ['CD4', "CD8", 'CLL', "NK", "Mono"]

    d = pd.melt(e.loc[:, genes].reset_index(),
        var_name="gene", id_vars=["cell_type", "patient_id", "timepoint"])
    d = d[d['patient_id'] == 'PT']

    fig, axis = plt.subplots(1, len(genes), figsize=(4 * len(genes), 4))
    for i, gene in enumerate(genes):
        sns.barplot(data=d[d['gene'] == gene], x="timepoint", y="value", hue="cell_type", ax=axis[i])
        axis[i].set_title(gene)
    fig.savefig("/home/arendeiro/plot2.svg", bbox_inches="tight", dpi=300)



    df = pd.read_csv(os.path.join("results", "offtarget_signature.enrichment.WikiPathways.tsv"), sep='\t')
    df["Term"] = df["Term"].str.replace("_.*", "")
    df2 = df.sort_values('Adjusted P-value').head(20)

    fig, axis = plt.subplots(1, figsize=(5.0 / 2., 8.5 / 2.))
    sns.barplot(x=-np.log10(df2['Adjusted P-value']), y=df2['Term'], ax=axis, orient="horizontal", estimator=max)
    # axis.scatter(x=-np.log10(df2['Adjusted P-value']), y=-np.log10(df2['Adjusted P-value']).rank(), s=df2['odds_ratio'], alpha=0.6)
    fig.savefig(os.path.join("offtarget_signature.enrichment.WikiPathways.barplot2.svg"), dpi=300, bbox_inches="tight")

    fig, axis = plt.subplots(1, figsize=(5.0 / 2., 8.5 / 2.))
    sns.barplot(x=df2['odds_ratio'], y=df2['Term'], ax=axis, orient="horizontal", estimator=max)
    # axis.scatter(x=-np.log10(df2['Adjusted P-value']), y=-np.log10(df2['Adjusted P-value']).rank(), s=df2['odds_ratio'], alpha=0.6)
    fig.savefig(os.path.join("offtarget_signature.enrichment.WikiPathways.barplot2.odds.svg"), dpi=300, bbox_inches="tight")

    # Calculate odds ratio:
    df['d'] = df['Overlap'].str.split("/").apply(lambda x: int(x[0]))
    df['c'] = df['Overlap'].str.split("/").apply(lambda x: int(x[1]))
    df['b'] = common.shape[0] - df['d']
    df['a'] = 20000 - (df['b'] + df['c'] + df['d'])
    df['odds_ratio'] = df.apply(lambda x: fisher_exact([[x['a'], x['b']], [x['c'], x['d']]])[0], axis=1)

    # reduce redundant terms
    df['Term'] = df['Term'].str.lower()
    df2 = df.groupby(['Term'])[['odds_ratio', 'Adjusted P-value']].max()

    fig, axis = plt.subplots(1, figsize=(3, 3))
    axis.scatter(x=df2['odds_ratio'], y=-np.log10(df2['Adjusted P-value']), s=3, alpha=0.6)
    for term in df2.sort_values('Adjusted P-value').head(10).index:
        axis.text(x=df2.loc[term, 'odds_ratio'], y=-np.log10(df2.loc[term, 'Adjusted P-value']), s=term, fontsize=6, ha="right")
    axis.set_xlabel("Odds ratio")
    axis.set_ylabel("-log(q-value)")
    sns.despine(fig)
    fig.savefig(os.path.join("offtarget_signature.enrichment.WikiPathways.scatter.svg"), dpi=300, bbox_inches="tight")






    # Get a cell type mixture for each patient, timepoint
    mix = e.groupby(['patient_id', 'timepoint']).mean().T

    a = (mix.loc[up, mix.columns.get_level_values("timepoint") == 0].mean() /
            mix.loc[down, mix.columns.get_level_values("timepoint") == 0].mean())
    b = (mix.loc[up, mix.columns.get_level_values("timepoint").isin([30, 120])].mean() /
            mix.loc[down, mix.columns.get_level_values("timepoint").isin([30, 120])].mean())
    # select earliest timepoint per patient
    q = b.reset_index().groupby('patient_id').apply(lambda x: x['timepoint'].min())
    b = b.loc[zip(q.index, q)]
    t, p = scipy.stats.ttest_rel(b, a)


    # Get validation cohort
    akh = pd.read_csv(
        os.path.join("..", "cll-ibrutinib", "results", "cll-ibrutinib_AKH-RNA.expression_counts.gene_level.quantile_normalized.log2_tpm.annotated_metadata.csv"),
        header=range(22), index_col=0).sort_index(level=['patient_id', 'timepoint_name'], axis=1)


    a = (akh.loc[up, akh.columns.get_level_values("timepoint_name") == "before_Ibrutinib"].mean() /
            akh.loc[down, akh.columns.get_level_values("timepoint_name") == "before_Ibrutinib"].mean())
    b = (akh.loc[up, akh.columns.get_level_values("timepoint_name") == "after_Ibrutinib"].mean() /
            akh.loc[down, akh.columns.get_level_values("timepoint_name") == "after_Ibrutinib"].mean())
    t, p = scipy.stats.ttest_rel(b, a)

    data = a.reset_index()[['patient_id', 'timepoint_name', 0]].append(b.reset_index()[['patient_id', 'timepoint_name', 0]])
    fig, axis = plt.subplots(1, 1, figsize=(3, 3))
    sns.pointplot(x="timepoint_name", y=0, hue="patient_id", data=data, axis=axis)
    axis.text(0.9, 0.9, "t = {:3f}\np = {:3f}".format(t, p))
    fig.savefig(os.path.join(analysis.results_dir, "offtarget_signature.akh.validation.paired_pointplot.svg"), dpi=300, bbox_inches="tight")

    data = a.reset_index()[['patient_id', 'timepoint_name', 0]].append(b.reset_index()[['patient_id', 'timepoint_name', 0]])
    data.sort_values(0)

    tp = (data.loc[:, 'timepoint_name'] == 'before_Ibrutinib').sum()
    tn = (data.loc[:, 'timepoint_name'] == 'after_Ibrutinib').sum()
    tpr = list()
    fpr = list()
    for threshold in np.linspace(data[0].min(), data[0].max(), 20):
        p = (data.loc[data[0] <= threshold, 'timepoint_name'] == 'before_Ibrutinib').sum()
        f = (data.loc[data[0] <= threshold, 'timepoint_name'] == 'after_Ibrutinib').sum()
        tpr.append(float(p) / tp)
        fpr.append(float(f) / tn)
    from sklearn import metrics
    auc = metrics.auc(fpr, tpr)

    fig, axis = plt.subplots(1, figsize=(3, 3))
    axis.plot(fpr, tpr)
    axis.plot((0, 1), (0, 1), linestyle="--", color="black", alpha=0.5)
    axis.text(0.7, 0.3, "AUC = {:3f}".format(auc))
    axis.set_xlabel("FPR")
    axis.set_ylabel("TPR")
    fig.savefig(os.path.join(analysis.results_dir, "offtarget_signature.akh.validation.roc_auc.svg"), dpi=300, bbox_inches="tight")


    g2 = sns.clustermap(
        akh.sort_index(level='time_since_treatment', axis=1).loc[diff.index].dropna(), xticklabels=True, z_score=0,
        rasterized=True, col_cluster=True, metric="correlation", row_cluster=True)
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.savefig(os.path.join(analysis.results_dir, "offtarget_signature.akh.svg"), dpi=300, bbox_inches="tight")
    g2 = sns.clustermap(
        akh.sort_index(level='time_since_treatment', axis=1).loc[diff[g.dendrogram_col.reordered_ind].index].dropna(), xticklabels=True, z_score=0,
        rasterized=True, col_cluster=False, metric="correlation", row_cluster=False)
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.savefig(os.path.join(analysis.results_dir, "offtarget_signature.akh.sorted_time.svg"), dpi=300, bbox_inches="tight")
    g2 = sns.clustermap(
        akh.sort_index(level='time_since_treatment', axis=1).loc[diff[g.dendrogram_col.reordered_ind].index].dropna(), xticklabels=True, z_score=0,
        rasterized=True, col_cluster=True, metric="correlation", row_cluster=False)
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.savefig(os.path.join(analysis.results_dir, "offtarget_signature.akh.clustered.svg"), dpi=300, bbox_inches="tight")
    g2 = sns.clustermap(
        akh.sort_index(level=['patient_id', 'timepoint_name'], axis=1).loc[diff[g.dendrogram_col.reordered_ind].index].dropna(), xticklabels=True, z_score=0,
        rasterized=True, col_cluster=False, metric="correlation", row_cluster=False)
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.savefig(os.path.join(analysis.results_dir, "offtarget_signature.akh.sorted_patient.svg"), dpi=300, bbox_inches="tight")
    g2 = sns.clustermap(
        akh.sort_index(level='timepoint_name', axis=1).loc[diff[g.dendrogram_col.reordered_ind].index].dropna(), xticklabels=True, z_score=0,
        rasterized=True, col_cluster=False, metric="correlation", row_cluster=False)
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.savefig(os.path.join(analysis.results_dir, "offtarget_signature.akh.sorted_timepoint.svg"), dpi=300, bbox_inches="tight")




    # come up with a ibrutinib signature score
    expr = pd.read_csv("results/single_cell_RNA/MeanMatrix.csv", index_col=0)
    c = pd.Series(expr.columns).str.split("_").apply(pd.Series)
    c[1] = c[1].str.replace('d', '').astype(int)
    c = c[[1, 2, 0]]
    expr.columns = pd.MultiIndex.from_arrays(c.values.T, names=['timepoint', 'cell_type', 'patient_id'])
    expr = expr.sort_index(level=['timepoint', 'cell_type'], axis=1)
    # Z-score values within each cell type
    expr_z = expr.copy()
    for cell_type in expr.columns.levels[1]:
        expr_z.loc[:, expr_z.columns.get_level_values("cell_type") == cell_type] = scipy.stats.zscore(
            expr_z.loc[:, expr_z.columns.get_level_values("cell_type") == cell_type], axis=1)

    e = expr_z.loc[common.index, expr_z.columns.get_level_values("cell_type").isin(["CLL", "CD4", "CD8", "Mono"])].dropna().T.sort_index(axis=0, level=['timepoint', 'cell_type'])
    diff = e[e.index.get_level_values('timepoint') >= 120].mean() - e[e.index.get_level_values('timepoint') == 0].mean()
    up = diff[diff > 0].index
    down = diff[diff <= 0].index
    # sig = e.loc[:, up].mean(axis=1) / e.loc[:, down].mean(axis=1)

    # Get second validation cohort
    landau = pd.read_csv(os.path.join("landau.2017.cll_ibrutinib_RNA-seq.csv"), index_col=0)
    landau.columns = pd.MultiIndex.from_arrays(
        np.array(landau.columns.str.split("_").tolist()).T,
        names=['patient_id', 'timepoint'])

    # # get only patients with timepoint 1 and 3
    # ts = landau.T.reset_index().groupby('patient_id')['timepoint'].sum()
    # landau = landau.loc[:, ts.loc[ts.str.contains("T01") & ts.str.contains("T02") & ts.str.contains("T03")].index]

    res = pd.DataFrame()
    for timepoint in landau.columns.levels[1]:
        r = (landau.loc[up, landau.columns.get_level_values("timepoint") == timepoint].mean() /
             landau.loc[down, landau.columns.get_level_values("timepoint") == timepoint].mean())
        r.name = timepoint
        res = res.append(r.reset_index(level=1, drop=True))
    res = res.T

    fig, axis = plt.subplots(1, 2, figsize=(2 * 3, 3), sharey=True)
    t, p = scipy.stats.ttest_rel(*res[['T03', 'T01']].dropna().T.values.tolist())
    fig.suptitle('Landau et al, cohort' + "\nt = {:.1f}\np = {:.1E}".format(t, p), fontsize=12)
    t, p = scipy.stats.ttest_rel(*res[['T02', 'T01']].dropna().T.values.tolist())

    sns.pointplot(
        x="timepoint", y="value", hue="index",
        data=res[['T01', 'T02']].reset_index().melt(id_vars="index", var_name="timepoint"),
        ax=axis[0], palette="tab20")
    axis[0].text(0.9, 0.9, "t = {:.1f}\np = {:.1E}".format(t, p))
    axis[0].set_title("1 month vs pre-theraphy")

    t, p = scipy.stats.ttest_rel(*res[['T03', 'T02']].dropna().T.values.tolist())
    sns.pointplot(
        x="timepoint", y="value", hue="index",
        data=res[['T02', 'T03']].reset_index().melt(id_vars="index", var_name="timepoint"),
        ax=axis[1], palette="tab20")
    axis[1].text(0.9, 0.9, "t = {:.1f}\np = {:.1E}".format(t, p))
    axis[1].set_title("6 month vs pre-theraphy")
    for ax in axis:
        ax.set_ylabel("Ibrutinib signature score")
    fig.savefig(os.path.join(analysis.results_dir, "offtarget_signature.landau.validation.paired_pointplot.svg"), dpi=300, bbox_inches="tight")

    # ROC curve
    from sklearn import metrics
    data = res.reset_index().melt(id_vars="index", var_name="timepoint")
    data = data.sort_values("timepoint").dropna()

    fig, axis = plt.subplots(1, figsize=(3, 3))
    axis.plot((0, 1), (0, 1), linestyle="--", color="black", alpha=0.5)
    axis.set_xlabel("FPR")
    axis.set_ylabel("TPR")
    for c, label in enumerate(["T03-T01","T02-T01","T03-T02"]):
        print(label)
        t2, t1 = label.split("-")

        tp = (data.loc[:, 'timepoint'] == t1).sum()
        tn = (data.loc[:, 'timepoint'] == t2).sum()
        tpr = list()
        fpr = list()
        # for threshold in np.linspace(data["value"].min(), data["value"].max(), 20):
        for threshold in data["value"].sort_values().unique():
            p = (data.loc[data["value"] <= threshold, 'timepoint'] == t1).sum()
            n = (data.loc[data["value"] <= threshold, 'timepoint'] == t2).sum()
            tpr.append(float(p) / tp)
            fpr.append(float(n) / tn)
        auc = metrics.auc(fpr, tpr)

        axis.plot(fpr, tpr, label="{}; AUC = {:3f}".format(label, auc), color=sns.color_palette()[c])

        # randomize and plot
        r_fpr = list()
        r_tpr = list()
        r_aucs = list()
        for i in range(100):
            data2 = data.loc[data['timepoint'].isin([t1, t2])].copy()
            data2['value'] = data2['value'].sample(frac=1).values
            tp = (data2.loc[:, 'timepoint'] == t1).sum()
            tn = (data2.loc[:, 'timepoint'] == t2).sum()
            tpr = list()
            fpr = list()
            # for threshold in np.linspace(data2["value"].min(), data2["value"].max(), 20):
            for threshold in data2["value"].sort_values().unique():
                p = (data2.loc[data2["value"] <= threshold, 'timepoint'] == t1).sum()
                n = (data2.loc[data2["value"] <= threshold, 'timepoint'] == t2).sum()
                tpr.append(float(p) / tp)
                fpr.append(float(n) / tn)
            auc = metrics.auc(fpr, tpr)
            r_aucs.append(auc)
            # axis.plot(fpr, tpr, color=sns.color_palette()[c], alpha=0.1, linestyle=":")
            r_fpr.append(fpr)
            r_tpr.append(tpr)

        axis.plot(
            np.array(r_fpr).mean(0), np.array(r_tpr).mean(0),
            color=sns.color_palette()[c], linestyle=":",
            label="Random {}; AUC = {:3f}".format(label, np.mean(r_aucs)))
    axis.legend()
    fig.savefig(os.path.join(analysis.results_dir, "offtarget_signature.landau.validation.roc_auc.svg"), dpi=300, bbox_inches="tight")

    g2 = sns.clustermap(
        landau.sort_index(level='timepoint', axis=1).loc[diff.index].dropna(),
        xticklabels=True, z_score=0, center=0, cmap="RdBu_r",
        rasterized=True, col_cluster=True, metric="correlation", row_cluster=True)
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.savefig(os.path.join(analysis.results_dir, "offtarget_signature.landau.svg"), dpi=300, bbox_inches="tight")

    g2 = sns.clustermap(
        landau.sort_index(level='timepoint', axis=1).loc[diff.index].dropna(),
        xticklabels=True, z_score=0, center=0, cmap="RdBu_r",
        rasterized=True, col_cluster=False, metric="correlation", row_cluster=True)
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.savefig(os.path.join(analysis.results_dir, "offtarget_signature.landau.sorted_time.svg"), dpi=300, bbox_inches="tight")

    g2 = sns.clustermap(
        landau.sort_index(level='timepoint', axis=1).loc[diff.index].dropna(),
        xticklabels=True, z_score=0, center=0, cmap="RdBu_r",
        rasterized=True, col_cluster=True, metric="correlation", row_cluster=True)
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.savefig(os.path.join(analysis.results_dir, "offtarget_signature.landau.clustered.svg"), dpi=300, bbox_inches="tight")
    g2 = sns.clustermap(
        landau.sort_index(level=['patient_id', 'timepoint'], axis=1).loc[diff.index].dropna(),
        xticklabels=True, z_score=0, center=0, cmap="RdBu_r",
        rasterized=True, col_cluster=False, metric="correlation", row_cluster=True)
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.savefig(os.path.join(analysis.results_dir, "offtarget_signature.landau.sorted_patient.svg"), dpi=300, bbox_inches="tight")
    g2 = sns.clustermap(
        landau.sort_index(level='timepoint', axis=1).loc[diff.index].dropna(),
        xticklabels=True, z_score=0, center=0, cmap="RdBu_r",
        rasterized=True, col_cluster=False, metric="correlation", row_cluster=True)
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.savefig(os.path.join(analysis.results_dir, "offtarget_signature.landau.sorted_timepoint.svg"), dpi=300, bbox_inches="tight")

    # # Investigate cell cycle changes in single cells
    # # stdscore = lambda x: (x - x.min()) / (x.max() - x.min())

    # # load predictions
    # # cc = pd.read_table(os.path.join("results", "single_cell_RNA", "50_cell_cycle", "all_cells.txt"))
    # cc = pd.read_table(os.path.join("data", "cell_cycle_assignments.all_cells.txt"))
    # cc = cc.join(
    #     cc['sample']
    #     .str.replace("FE1_", "")
    #     .str.replace("KI1_", "")
    #     .str.split("_").apply(pd.Series).rename(columns={0: "patient_id", 1: "timepoint"}))
    # cc['patient_id'] = cc['patient_id'].str.replace(r'\d', "")
    # cc['timepoint'] = cc['timepoint'].str.replace("d", "").astype(int)

    # # plot distribution of scores
    # cc['mean'] = (cc.loc[:, "G2M.Score"] + cc.loc[:, "S.Score"]) / 2
    # cc['fc'] = np.log2(cc.loc[:, "G2M.Score"] / cc.loc[:, "S.Score"])

    # phases = sorted(cc["Phase"].unique())
    # fig, axis = plt.subplots(1, 2, figsize=(8, 4))
    # for i, phase in enumerate(phases):
    #     axis[0].scatter(
    #         cc.loc[cc["Phase"] == phase, "S.Score"],
    #         cc.loc[cc["Phase"] == phase, "G2M.Score"],
    #         c=plt.get_cmap("Paired")(i), label=phase, s=3, alpha=0.1)
    #     axis[1].scatter(
    #         cc.loc[cc["Phase"] == phase, "mean"],
    #         cc.loc[cc["Phase"] == phase, "fc"],
    #         c=plt.get_cmap("Paired")(i), label=phase, s=3, alpha=0.1)

    # axis[0].plot((-0.2, 1), (-0.2, 1), linestyle="--", color="black", alpha=0.5)
    # axis[1].axhline(0, linestyle="--", color="black", alpha=0.5)
    # axis[1].axvline(0, linestyle="--", color="black", alpha=0.5)
    # axis[0].set_xlabel("S score")
    # axis[0].set_ylabel("G2M score")
    # axis[1].set_xlabel("Mean")
    # axis[1].set_ylabel("log2 fold change (G2M / S)")
    # axis[0].set_xlim(-0.2, 1)
    # axis[0].set_ylim(-0.2, 1)
    # axis[0].legend()
    # axis[1].legend()

    # # set thresholds:
    # axis[1].axhline(-0.5, linestyle="--", color="red", alpha=0.5)
    # axis[1].axhline(0.5, linestyle="--", color="red", alpha=0.5)
    # axis[1].axvline(-0.02, linestyle="--", color="red", alpha=0.5)
    # axis[1].axvline(0.02, linestyle="--", color="red", alpha=0.5)

    # sel_cc = cc.loc[
    #     (cc["mean"].abs() > 0.02) &
    #     (cc["fc"].abs() > 0.5)]

    # fig, axis = plt.subplots(1, 2, figsize=(8, 4))
    # for i, phase in enumerate(phases):
    #     axis[0].scatter(
    #         sel_cc.loc[sel_cc["Phase"] == phase, "S.Score"],
    #         sel_cc.loc[sel_cc["Phase"] == phase, "G2M.Score"],
    #         c=plt.get_cmap("Paired")(i), label=phase, s=3, alpha=0.1)
    #     axis[1].scatter(
    #         sel_cc.loc[sel_cc["Phase"] == phase, "mean"],
    #         sel_cc.loc[sel_cc["Phase"] == phase, "fc"],
    #         c=plt.get_cmap("Paired")(i), label=phase, s=3, alpha=0.1)

    # axis[0].plot((-0.2, 1), (-0.2, 1), linestyle="--", color="black", alpha=0.5)
    # axis[1].axhline(0, linestyle="--", color="black", alpha=0.5)
    # axis[1].axvline(0, linestyle="--", color="black", alpha=0.5)
    # axis[0].set_xlabel("S score")
    # axis[0].set_ylabel("G2M score")
    # axis[1].set_xlabel("Mean")
    # axis[1].set_ylabel("log2 fold change (G2M / S)")
    # axis[0].set_xlim(-0.2, 1)
    # axis[0].set_ylim(-0.2, 1)
    # axis[0].legend()
    # axis[1].legend()

    # # summary stats
    # fractions = (
    #     cc.groupby(['patient_id', 'cellType', 'timepoint'])
    #     ['Phase']
    #     .apply(lambda x: x.value_counts() / x.count())
    #     .reset_index()
    #     .rename(columns={"level_3": "Phase", "Phase": "fraction"}))
    # fractions2 = fractions.loc[
    #     (fractions['cellType'].isin(["CD4", "CD8", "CLL", "Mono"])) &
    #     (~fractions['patient_id'].isin(["KI"]))]
    # g = sns.factorplot(
    #     data=fractions2, x="timepoint", y='fraction',
    #     row="patient_id", col="cellType", sharey=False, hue='Phase')
    # g.savefig(os.path.join("results", "cell_cycle_changes.all_cells.factorplot.svg"), dpi=300, bbox_inches="tight")
    # g = sns.factorplot(
    #     data=fractions2,#.loc[fractions2['timepoint'] <= 120, :],
    #     x="timepoint", y='fraction',
    #     col="cellType", sharey=False, hue='Phase')
    # g.savefig(os.path.join("results", "cell_cycle_changes.all_cells.cross_patients.factorplot.svg"), dpi=300, bbox_inches="tight")

    # fractions = (
    #     sel_cc.groupby(['patient_id', 'cellType', 'timepoint'])
    #     ['Phase']
    #     .apply(lambda x: x.value_counts() / x.count())
    #     .reset_index()
    #     .rename(columns={"level_3": "Phase", "Phase": "fraction"}))
    # fractions2 = fractions.loc[
    #     (fractions['cellType'].isin(["CD4", "CD8", "CLL", "Mono"])) &
    #     (~fractions['patient_id'].isin(["KI"]))]
    # g = sns.factorplot(
    #     data=fractions2, x="timepoint", y='fraction',
    #     row="patient_id", col="cellType", sharey=False, hue='Phase')
    # g.savefig(os.path.join("results", "cell_cycle_changes.filtered_cells.factorplot.svg"), dpi=300, bbox_inches="tight")
    # # analysis.results_dir


def offtarget_on_atac(analysis):

    # get chromatin accessibility at gene-level
    analysis.matrix_norm = pd.read_csv(os.path.abspath(os.path.join(
            "results", analysis.name + ".accessibility.annotated_metadata.csv")),
        header=list(range(9)), index_col=0)
    check = analysis.gene_annotation.index.str.contains(":").all() and analysis.gene_annotation.index.str.contains("-").all()
    if check:
        analysis.gene_annotation = analysis.gene_annotation.reset_index()
        analysis.gene_annotation.index = (
            analysis.gene_annotation['chrom'] + ":" +
            analysis.gene_annotation['start'].astype(str) + "-" +
            analysis.gene_annotation['end'].astype(str))
    mean_acc = analysis.get_gene_level_accessibility()

    # get signature
    common = pd.read_csv(os.path.join("results", "offtarget_signature.csv"), index_col=0, header=None, squeeze=True)
    g = common.abs().argmin()
    common = (common - common.min()) / (common.max() - common.min())
    common -= common.loc[g]

    sig = common.index.to_series()

    sig_acc = mean_acc.reindex(sig).dropna()
    c = analysis.get_level_colors(index=sig_acc.columns, levels=['patient_id', 'timepoint', 'cell_type'])
    sample_color = pd.DataFrame(c, index=['patient_id', 'timepoint', 'cell_type'], columns=sig_acc.columns)
    direction_color = plt.get_cmap("RdBu_r")(common.reindex(sig_acc.index))

    g = sns.clustermap(
        sig_acc.T,
        metric="correlation",
        col_colors=direction_color, row_colors=sample_color.T,
        rasterized=True, robust=True,
        xticklabels=True, yticklabels=sig_acc.columns.get_level_values("sample_name"))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=4)
    g.savefig(os.path.join("results", "offtarget_signature.atacseq_level.svg"), dpi=300)

    g = sns.clustermap(
        sig_acc.T,
        metric="correlation",
        z_score=1, cmap="RdBu_r", center=0,
        col_colors=direction_color, row_colors=sample_color.T, rasterized=True, robust=True,
        xticklabels=True, yticklabels=sig_acc.columns.get_level_values("sample_name"))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=4)
    g.savefig(os.path.join("results", "offtarget_signature.atacseq_level.z_score.svg"), dpi=300)

    g = sns.clustermap(
        sig_acc.T.groupby(['cell_type', 'timepoint']).mean(),
        metric="correlation",
        z_score=1, cmap="RdBu_r", center=0,
        col_colors=direction_color,rasterized=True, robust=True,
        xticklabels=True, yticklabels=sig_acc.columns.get_level_values("sample_name"))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join("results", "offtarget_signature.atacseq_level.timepoint_mean.z_score.svg"), dpi=300)

    g = sns.clustermap(
        sig_acc.T.groupby(['cell_type', 'timepoint']).mean(),
        metric="correlation",
        row_cluster=False,
        z_score=1, cmap="RdBu_r", center=0,
        col_colors=direction_color, rasterized=True, robust=True,
        xticklabels=True, yticklabels=True)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join("results", "offtarget_signature.atacseq_level.timepoint_mean.z_score.ordered.svg"), dpi=300)

    for cell_type in sig_acc.columns.levels[sig_acc.columns.names.index('cell_type')]:
        print(cell_type)
        p = sig_acc.loc[:, sig_acc.columns.get_level_values("cell_type") == cell_type]
        p = p.sort_index(axis=1, level=['timepoint'])

        g = sns.clustermap(
            p.T,
            metric="correlation", row_cluster=False,
            z_score=1, cmap="RdBu_r", center=0,
            col_colors=direction_color, row_colors=sample_color.T, rasterized=True, robust=True, figsize=(12, 5),
            xticklabels=True, yticklabels=p.columns.get_level_values("sample_name"))
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=4)
        g.savefig(os.path.join("results", "offtarget_signature.atacseq_level.{}.z_score.ordered.svg".format(cell_type)), dpi=300)

        g = sns.clustermap(
            p.T,
            metric="correlation",
            z_score=1, cmap="RdBu_r", center=0,
            col_colors=direction_color, row_colors=sample_color.T, rasterized=True, robust=True, figsize=(12, 5),
            xticklabels=True, yticklabels=p.columns.get_level_values("sample_name"))
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=4)
        g.savefig(os.path.join("results", "offtarget_signature.atacseq_level.{}.z_score.svg".format(cell_type)), dpi=300)

        g = sns.clustermap(
            p.T.drop(["003d", "008d", "150d"], level="timepoint", axis=0),
            metric="correlation",
            z_score=1, cmap="RdBu_r", center=0,
            col_colors=direction_color,
            row_colors=sample_color.T.drop(["003d", "008d", "150d"], level="timepoint", axis=0),
            row_cluster=False,
            rasterized=True, robust=True, figsize=(12, 5),
            xticklabels=True,
            yticklabels=p.T.drop(["003d", "008d", "150d"], level="timepoint", axis=0).index.get_level_values("sample_name"))
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=4)
        g.savefig(os.path.join("results", "offtarget_signature.atacseq_level.{}.z_score.selection.ordered.svg".format(cell_type)), dpi=300)

        g = sns.clustermap(
            p.T.groupby("timepoint").mean(),
            metric="correlation",
            z_score=1, cmap="RdBu_r", center=0,
            col_colors=direction_color, rasterized=True, robust=True, figsize=(12, 2.5),
            xticklabels=True, yticklabels=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=4)
        g.savefig(os.path.join("results", "offtarget_signature.atacseq_level.{}.timepoint_mean.z_score.svg".format(cell_type)), dpi=300)

        g = sns.clustermap(
            p.T.groupby("timepoint").mean(), row_cluster=False,
            metric="correlation",
            z_score=1, cmap="RdBu_r", center=0,
            col_colors=direction_color, rasterized=True, robust=True, figsize=(12, 2.5),
            xticklabels=True, yticklabels=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.savefig(os.path.join("results", "offtarget_signature.atacseq_level.{}.timepoint_mean.z_score.ordered.svg".format(cell_type)), dpi=300)

        g = sns.clustermap(
            p.T.groupby("timepoint").mean().drop(["003d", "008d", "150d"], axis=0), row_cluster=False,
            metric="correlation",
            z_score=1, cmap="RdBu_r", center=0,
            col_colors=direction_color, rasterized=True, robust=True, figsize=(12, 2.5),
            xticklabels=True, yticklabels=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.savefig(os.path.join("results", "offtarget_signature.atacseq_level.{}.timepoint_mean.z_score.ordered.selection.svg".format(cell_type)), dpi=300)


def investigate_unclustered_regions():

    root_dir = os.path.join("/home/afr/Downloads", "Ibrutinib_time_course", "tables_and_supplement", "previous")
    df = pd.read_excel(os.path.join(root_dir, "Supplementary Table 5.xls"))

    for cell_type in df['cell_type'].unique():
        if not os.path.exists(os.path.join(root_dir, cell_type)):
            os.makedirs(os.path.join(root_dir, cell_type))
        regions = df.loc[(df['cell_type'] == cell_type) & (df['assigned'] == False), "region"]
        df2 = pd.DataFrame([
            regions.str.split(":").apply(lambda x: x[0]),
            regions.str.split(":").apply(lambda x: x[1].split("-")[0]),
            regions.str.split(":").apply(lambda x: x[1].split("-")[1])]).T
        df2.to_csv(
            os.path.join(root_dir, cell_type, "unclustered_regions.{}.bed".format(cell_type)),
            index=False, sep="\t", header=False)

        specific_universe = os.path.join(root_dir, cell_type, "dynamic_regions.{}.bed".format(cell_type))
        regions = df.loc[(df['cell_type'] == cell_type), "region"]
        df2 = pd.DataFrame([
            regions.str.split(":").apply(lambda x: x[0]),
            regions.str.split(":").apply(lambda x: x[1].split("-")[0]),
            regions.str.split(":").apply(lambda x: x[1].split("-")[1])]).T
        df2.to_csv(specific_universe, index=False, sep="\t", header=False)

    prefix = "unclustered_regions"
    output_dir = root_dir = os.path.join("results", "unclustered_regions")
    genome = "hg19"

    universe_bed = os.path.abspath(os.path.join("results", "cll-time_course_peak_set.bed"))
    for cell_type in ["Bcell", "CD4", "CD8", "CLL", "Mono", "NK"]:

        background_type = "universe"
        name = ".".join([prefix, cell_type, background_type])
        log = os.path.join(output_dir, "log", name + ".log")
        job = " ".join([
            "Rscript /home/arendeiro/jobs/run_LOLA.R",
            " {}".format(os.path.join(root_dir, cell_type, "unclustered_regions.{}.bed".format(cell_type))),
            " {}".format(universe_bed),
            " {}".format(genome)])
        cmd = " ".join([
            "sbatch",
            "-J {}".format(name),
            "-o {}".format(log),
            "-p {}".format("shortq"),
            "-c {}".format(4),
            "--mem {}".format(20000),
            "--wrap",
            "'{}'".format(job)])
        os.system(cmd)

        background_type = "specific_universe"
        name = ".".join([prefix, cell_type, background_type])
        log = os.path.join(output_dir, "log", name + ".log")
        job = " ".join([
            "Rscript /home/arendeiro/jobs/run_LOLA.R",
            " {}".format(os.path.join(root_dir, cell_type, "dynamic_regions.{}.bed".format(cell_type))),
            " {}".format(os.path.join(root_dir, cell_type, "dynamic_regions.{}.bed".format(cell_type))),
            " {}".format(genome)])
        cmd = " ".join([
            "sbatch",
            "-J {}".format(name),
            "-o {}".format(log),
            "-p {}".format("shortq"),
            "-c {}".format(4),
            "--mem {}".format(20000),
            "--wrap",
            "'{}'".format(job)])
        os.system(cmd)


def temporal_ibrutinib_signature():
    from ngs_toolkit.general import standard_score

    # Gene expression
    expr = pd.read_csv("results/single_cell_RNA/MeanMatrix.csv", index_col=0)
    c = pd.Series(expr.columns).str.split("_").apply(pd.Series)
    c[1] = c[1].str.replace('d', '').astype(int)
    c = c[[1, 2, 0]]
    expr.columns = pd.MultiIndex.from_arrays(c.values.T, names=['timepoint', 'cell_type', 'patient_id'])
    expr = expr.sort_index(level=['timepoint', 'cell_type'], axis=1)

    # Z-score values within each cell type
    expr_z = expr.copy()
    for cell_type in expr.columns.levels[1]:
        expr_z.loc[:, expr_z.columns.get_level_values("cell_type") == cell_type] = scipy.stats.zscore(
            expr_z.loc[:, expr_z.columns.get_level_values("cell_type") == cell_type], axis=1)

    # Gene expression diff
    expr_diff = pd.read_csv(os.path.join(
        "results", "single_cell_RNA", "13_4_Overtime_nUMI_Cutoff", "SigGenes_overTime.tsv"), sep="\t").rename(
            columns={"cellType": "cell_type", "logFC": "log_fold_change"}).set_index("gene")
    expr_diff['intercept'] = 1
    expr_diff['patient_id'] = expr_diff['patient'].str.split("_").apply(lambda x: x[0])
    expr_diff['timepoint'] = expr_diff['patient'].str.split("_").apply(lambda x: x[1]).str.replace("d", "").astype(int)
    expr_diff['direction'] = (expr_diff['log_fold_change'] > 0).astype(int).replace(0, -1)

    expr_diff2 = expr_diff.loc[expr_diff['qvalue'] < 0.05]

    # Get signature
    common = pd.read_csv(os.path.join("results", "offtarget_signature.csv"), index_col=0, header=None, squeeze=True)

    # Get signature score for each patient, timepoint, cell type
    e_z = (expr_z.loc[common.index, expr_z.columns.get_level_values("cell_type").isin(["CLL", "CD4", "CD8", "Mono"])]
        .dropna().T.sort_index(axis=0, level=['timepoint', 'cell_type']))

    t_sigs = e_z.loc[:, common[common > 0].index].mean(axis=1) - e_z.loc[:, common[common < 0].index].mean(axis=1)

    # Aggregate across cell_types for a combined patient sig
    p_sig = t_sigs[t_sigs.index.get_level_values('cell_type') == "CLL"].groupby(level=['timepoint', 'patient_id']).mean()
    # p_sig = t_sigs.groupby(level=['timepoint', 'patient_id']).mean()

    # Get FACS data
    facs = pd.read_csv(os.path.join("metadata", "facs.cell_type.quantification.csv"))
    facs = facs.loc[facs['value_type'] == 'percentage', :]
    facs = pd.melt(facs[["timepoint", "patient_code", "CD5+ CLL"]], id_vars=["timepoint", "patient_code"], value_name="CLL").dropna()
    facs['CLL_scale'] = standard_score(facs['CLL'])
    facs['timepoint'] = facs['timepoint'].astype(int)
    facs = facs.drop("variable", axis=1).set_index(['timepoint', 'patient_code'])


    # Patient
    p = p_sig.index.levels[1]
    fig, axis = plt.subplots(len(p), 1, figsize=(6 * 1, 3 * len(p)), sharex=True)
    for i, patient in enumerate(p):
        p_sig_p = p_sig.loc[p_sig.index.get_level_values("patient_id") == patient]
        axis[i].plot(
            p_sig_p.index.get_level_values("timepoint") + 1, p_sig_p,
            color="blue", linestyle='--', marker='o', label="Ibrutinib response signature")

        ax2 = axis[i].twinx()
        facs_p = facs.loc[facs.index.get_level_values("patient_code") == patient, "CLL"]
        ax2.plot(
            facs_p.index.get_level_values("timepoint") + 1, facs_p,
            color="orange", linestyle='--', marker='o', label="CLL cells per uL")

        axis[i].set_ylim((-1.5, 1.5))
        # ax2.set_ylim((-1, 1))

        # axis[i].set_xscale("log", basex=2)
        # ax2.set_xscale("log", basex=2)

        # axis[i].legend()
        # ax2.legend()

        axis[i].axhline(0, linestyle="--", alpha=0.5, color="black", zorder=1000)
        ax2.axhline(0, linestyle="--", alpha=0.5, color="black", zorder=1000)

        axis[i].set_ylabel("Signature score", color="blue")
        ax2.set_ylabel("CLL cells per uL", color="orange")

        axis[i].set_title(patient)

    axis[i].set_xlabel("Days since ibrutinib start")
    ax2.set_xlabel("Days since ibrutinib start")

    sns.despine(fig, top=True, left=False, right=False, bottom=False)
    fig.savefig(os.path.join("results", "temporal_ibrutinib_signature.per_patient2.svg"), dpi=300, bbox_inches="tight")

    # Cell type
    # get non-zscored values
    # Aggregate across patients for a combined cell type sig
    e = (expr.loc[common.index, expr.columns.get_level_values("cell_type").isin(["CLL", "CD4", "CD8", "Mono"])]
         .dropna().T.sort_index(axis=0, level=['timepoint', 'cell_type']))
    t_sigs = e.loc[:, common[common > 0].index].mean(axis=1) - e.loc[:, common[common < 0].index].mean(axis=1)
    ct_sig = t_sigs.groupby(level=['timepoint', 'cell_type']).mean()

    ct = ct_sig.index.levels[1]
    fig, axis = plt.subplots(1, 1, figsize=(6 * 1, 3 * 1))
    for i, cell_type in enumerate(ct):
        ct_sig_ct = ct_sig.loc[ct_sig.index.get_level_values("cell_type") == cell_type]
        axis.plot(
            ct_sig_ct.index.get_level_values("timepoint") + 1, ct_sig_ct,
            linestyle='--', marker='o', label=cell_type)

    # axis.set_ylim((-1.5, 1.5))
    # axis.set_xscale("log")
    axis.axhline(0, linestyle="--", alpha=0.5, color="black", zorder=1000)
    axis.set_ylabel("Signature score", color="blue")
    axis.set_xlabel("Days since ibrutinib start")
    axis.legend()

    sns.despine(fig, top=True, left=False, right=True, bottom=False)
    fig.savefig(os.path.join("results", "temporal_ibrutinib_signature.per_cell_type.svg"), dpi=300, bbox_inches="tight")


def cd5_expression():
    df = pd.read_csv("results/single_cell_RNA/MeanMatrix.csv", index_col=0)
    df.columns = pd.MultiIndex.from_arrays(pd.Series(df.columns).str.split("_").apply(pd.Series).values.T, names=['patient_id', 'timepoint', 'cell_type'])

    m = df.mean(axis=1)

    fig, axis = plt.subplots(1, figsize=(3, 3))
    sns.distplot(m, kde=False, bins=30, ax=axis)
    axis.set_yscale("log")
    axis.axvline(m.loc['CD5'], linestyle="--", color="brown")
    axis.text(m.loc['CD5'], 10**4, s="CD5")
    axis.set_xlabel("Mean expression (log)")
    axis.set_ylabel("Genes (log)")
    sns.despine(fig)
    fig.savefig(os.path.join("results", "cd5_expression.distplot.svg"), bbox_inches="tight", dpi=300)


def pre_treatment_correlation_with_response(analysis):
    def test_pca_significance(
            matrix, iterations=None,
            standardize_matrix=True,
            adjust_pvalues=True,
            multipletest_method="fdr_bh"):
        """
        Empirical testing of Principal Component Analysis.

        Both principal components and feature's weights in components are tested.
        The approach is to simply shuffle the input data, perform PCA,
        and get empirical p-values for the PCA fitted on original data compared with the
        randomized background.

        :param matrix: Dataframe of shape (n_samples, n_variables)
        :type matrix: pandas.DataFrame
        :param iterations: Number of iterations to perform randomization.
                           Defaults to length of samples.
        :type iterations: int, optional
        :param standardize_matrix: Whether to standardize the input data prior to fitting.
                                   Defaults to False
        :type standardize_matrix: bool, optional
        :returns: Metrics per component, for real or random data (shape n_pcs,4);
                  Weights of every feature in every component, for real or random data (shape n_pcs,n_vars, 2);
                  P-values of every feature in every component (shape n_pcs,n_vars);
        :rtype: (pandas.DataFrame, pandas.DataFrame, pandas.DataFrame)
        """
        from sklearn.decomposition import PCA
        from statsmodels.distributions.empirical_distribution import ECDF
        from statsmodels.stats.multitest import multipletests
        from sklearn.preprocessing import StandardScaler

        def shuffle_array_rows_by_col(a):
            """
            Shuflle each of the 2D array's columns independently.
            Rows should be much smaller than columns
            :param np.array a: 2D array
            :returns np.array: Shuffled array
            """
            import itertools
            all_perm = np.array(list(itertools.permutations(list(range(a.shape[0])))))
            idx = all_perm[np.random.randint(0, all_perm.shape[0], size=a.shape[1])].T
            idx += a.shape[0] * np.arange(a.shape[1])
            assert a.shape == idx.shape
            b = a.flatten('F')[idx.flatten('C')].reshape(a.shape)
            assert np.allclose(a.sum(axis=0), b.sum(axis=0))
            return b

        def empirical_interpolated_p_value(a, b):
            from scipy.interpolate import interp1d
            ecdf = ECDF(b)
            slope_changes = sorted(set(a))
            inverted_edf = interp1d(ecdf(a), slope_changes, fill_value="extrapolate")
            return inverted_edf(a)

        if standardize_matrix:
            matrix = pd.DataFrame(
                StandardScaler().fit_transform(matrix),
                index=matrix.index, columns=matrix.columns)
        n, m = matrix.shape
        iterations = pd.Index(
            range(iterations if iterations is not None else n), name="iterations")
        var = pd.Index(matrix.columns, name="feature")

        pca = PCA(n_components=min(n, m) - 1)
        # First fit to whole data
        pca.fit(matrix)
        pcs = pd.Index(range(1, pca.n_components_ + 1), name="pc")
        v = pd.Series(pca.explained_variance_ratio_, index=pcs)
        w = pd.DataFrame(pca.components_.T, index=var, columns=pcs)

        # Now fit n times shuffling data across features
        random_v = pd.DataFrame(index=pcs, columns=iterations)
        random_w = pd.DataFrame(
            index=pd.MultiIndex.from_product([pcs.tolist(), var.tolist()], names=['pc', 'feature']),
            columns=iterations)
        for i in tqdm(iterations, total=len(iterations)):
            pca.fit(shuffle_array_rows_by_col(matrix.values.copy()))
            # pca.fit(pd.Series(matrix.values.flatten()).sample(frac=1).values.reshape(matrix.shape))
            random_v.loc[:, i] = pd.Series(pca.explained_variance_ratio_, index=pcs)
            random_w.loc[:, i] = pca.components_.flatten()

        # Test PCs
        ecdf = ECDF(-random_v.values.flatten())
        sig_pcs = pd.Series(ecdf(-v), index=pcs)
        if adjust_pvalues:
            adj_sig_pcs = pd.Series(multipletests(sig_pcs, method=multipletest_method)[1], index=pcs)
        ecdf = ECDF(-random_w.abs().values.flatten())
        sig_vars = pd.DataFrame(ecdf(-w.abs()), index=var, columns=pcs)
        if adjust_pvalues:
            adj_sig_vars = pd.DataFrame(index=var, columns=pcs)
            for pc in pcs:
                adj_sig_vars.loc[:, pc] = multipletests(
                    sig_vars.loc[:, pc].values.squeeze(), method=multipletest_method)[1]

        v = v.to_frame(name="real")
        v['randomized'] = random_v.mean(axis=1)
        v['log_ratio_over_random'] = np.log10(v['real'] / v['randomized'])
        v['p_value'] = sig_pcs
        v['q_value'] = adj_sig_pcs

        random_w = random_w.mean(axis=1).reset_index().pivot_table(
            index="feature", columns="pc", values=0)
        w = pd.concat([w, random_w], axis=1, keys=['real', 'random'])

        sig_vars = pd.concat([sig_vars, adj_sig_vars], axis=1, keys=['raw', 'adjusted'])

        return v, w, sig_vars

    # 1. First approach, find a axis correlating with response in unsupervised analysis
    # # whole data
    analysis.unsupervised_analysis(
        steps=['pca', 'pca_association'], quant_matrix="accessibility",
        attributes_to_plot=plotting_attributes, output_prefix="20190111.accessibility")

    for cell_type in analysis.matrix_norm.columns.levels[3]:
        # if cell_type == "CLL": continue
        # All samples from one cell type
        analysis.matrix_norm_ct = analysis.matrix_norm.loc[
            :,
            (analysis.matrix_norm.columns.get_level_values("cell_type") == cell_type)]
        exclude = ['cell_type', 'compartment', 'timepoint']
        analysis.unsupervised_analysis(
            steps=['pca', 'pca_association'], quant_matrix="accessibility_ct",
            attributes_to_plot=[x for x in plotting_attributes if x not in exclude],
            output_prefix="20190111.accessibility_{}".format(cell_type),
            standardize_matrix=True)

        # Only T0
        analysis.matrix_norm_t0_ct = analysis.matrix_norm.loc[
            :, (analysis.matrix_norm.columns.get_level_values("timepoint") == "000d") &
            (analysis.matrix_norm.columns.get_level_values("cell_type") == cell_type)]
        exclude = ['cell_type', 'compartment', 'timepoint']
        analysis.unsupervised_analysis(
            steps=['pca', 'pca_association'], quant_matrix="accessibility_t0_ct",
            attributes_to_plot=[x for x in plotting_attributes if x not in exclude],
            output_prefix="20190111.accessibility_t0_{}".format(cell_type),
            standardize_matrix=True)

        v, w, sig_vars = test_pca_significance(
            analysis.matrix_norm_t0_ct.T,
            iterations=None,
            standardize_matrix=True)

        # Assemble table of results
        v.to_csv(
            os.path.join("20190111.accessibility_t0_{}".format(cell_type) +
                         ".pc_significance_stats.csv"))

    # # # plot all association p-values jointly
    ps = pd.DataFrame()
    vs = pd.DataFrame()
    for cell_type in analysis.matrix_norm.columns.levels[3]:
        p = pd.read_csv(
            os.path.join(
                'results', 'unsupervised_analysis_ATAC-seq',
                'cll-time_course.20190111.accessibility_t0_{}.pca.variable_principle_components_association.csv'
                .format(cell_type)))
        p.loc[:, 'cell_type'] = cell_type
        ps = ps.append(p)

        r = pd.read_csv(
            os.path.join("20190111.accessibility_t0_{}".format(cell_type) +
                         ".pc_significance_stats.csv"))
        r.loc[:, 'cell_type'] = cell_type
        vs = vs.append(r)

    ps.loc[:, '-log_pvalue'] = -np.log10(ps.loc[:, 'p_value'])
    pss = ps.pivot_table(index=['cell_type', 'pc'], columns="attribute", values='-log_pvalue')

    fig, axis = plt.subplots(1, 1, figsize=(8, 8))
    sns.heatmap(
        pss.T,
        xticklabels=True, yticklabels=True,
        square=True, ax=axis)
    fig.savefig(os.path.join(
        "results", "unsupervised_analysis_ATAC-seq",
        "cll-time_course.20190111.accessibility_t0.all_cell_types.pca.variable_principle_components_association.svg"),
        bbox_inches="tight")

    vs['-log10(p)'] = (-np.log10(vs['p_value'])).replace(np.inf, 3)
    fig, axis = plt.subplots(3, 1, figsize=(3 * 8, 8))
    sns.heatmap(
        vs[['pc', 'cell_type', 'real', 'randomized']]
        .melt(id_vars=['cell_type', 'pc'])
        .pivot_table(index=['cell_type', "pc"], columns=['variable'], values='value').T,
        xticklabels=True, yticklabels=True,
        cmap="Reds",
        square=True, ax=axis[0])
    sns.heatmap(
        vs[['pc', 'cell_type', 'log_ratio_over_random']]
        .melt(id_vars=['cell_type', 'pc'])
        .pivot_table(index=['cell_type', "pc"], columns=['variable'], values='value').T,
        xticklabels=True, yticklabels=True,
        cmap="PuOr_r", center=0,
        square=True, ax=axis[1])
    sns.heatmap(
        vs[['pc', 'cell_type', '-log10(p)']]
        .melt(id_vars=['cell_type', 'pc'])
        .pivot_table(index=['cell_type', "pc"], columns=['variable'], values='value').T,
        xticklabels=True, yticklabels=True,
        cmap="Reds",
        square=True, ax=axis[2])
    fig.savefig(os.path.join(
        "results", "unsupervised_analysis_ATAC-seq",
        "cll-time_course.20190111.accessibility_t0.all_cell_types.pca.variance_explained.svg"),
        bbox_inches="tight")


    # look into CLL PC2
    output_dir = os.path.join(analysis.results_dir, "unsupervised_analysis_" + analysis.data_type)
    prefix = os.path.join(output_dir, analysis.name + ".20190111.accessibility_t0_cll")
    alpha = 0.1
    analysis.matrix_norm_t0_cll = analysis.matrix_norm.loc[
        :, (analysis.matrix_norm.columns.get_level_values("timepoint") == "000d") &
        (analysis.matrix_norm.columns.get_level_values("cell_type") == "CLL")]

    # # Test significance of PCs
    v, w, sig_vars = test_pca_significance(
        analysis.matrix_norm_t0_cll.T,
        iterations=None,
        standardize_matrix=True)

    # Assemble table of results
    v['log_pvalue'] = (-np.log10(v['p_value'])).replace(np.inf, 3)
    v.to_csv(os.path.join(prefix + ".pc_significance_stats.csv"))

    unsup_results = w['real'].reset_index().melt(id_vars='feature', value_name="real_weight").set_index(['feature', 'pc'])
    unsup_results = unsup_results.join(w['random'].reset_index().melt(id_vars='feature', value_name="random_weight").set_index(['feature', 'pc']))
    unsup_results = unsup_results.join((w['real'] - w['random']).reset_index().melt(id_vars='feature', value_name="weight_over_random").set_index(['feature', 'pc']))
    unsup_results = unsup_results.join(sig_vars['raw'].reset_index().melt(id_vars='feature', value_name="pvalue").set_index(['feature', 'pc']))
    unsup_results = unsup_results.join(sig_vars['adjusted'].reset_index().melt(id_vars='feature', value_name="padj").set_index(['feature', 'pc']))
    m = analysis.matrix_norm_t0_cll.mean(axis=1)
    m.name = 'mean'
    m.index.name = 'feature'
    unsup_results = unsup_results.join(m)
    unsup_results.to_csv(os.path.join(prefix + ".pc_feature_significance_stats.csv"))
    # unsup_results = pd.read_csv(os.path.join(prefix + ".pc_feature_significance_stats.csv"), index_col=[0, 1])

    # Plot
    # v['log_pvalue'] = v['log_pvalue'].replace(np.inf, v['log_pvalue'][v['log_pvalue'] != np.inf].max())

    # # PC stats
    fig, axis = plt.subplots(1, 3, figsize=(3 * 3, 3), tight_layout=True)
    v_melted = v.reset_index()[['pc', 'real', 'randomized']].melt(
            id_vars='pc', var_name='data_type', value_name="variance_explained_ratio")
    sns.barplot(
        data=v_melted,
        x="variance_explained_ratio", y="pc", hue="data_type", orient="h", ax=axis[0])
    axis[1].plot(v.index, v['real'], "o-", label="real")
    axis[1].plot(v.index, v['randomized'], "o-", label="randomized")
    axis[1].legend()
    axis[1].set_xlabel("PC")
    axis[1].set_ylabel("variance_explained_ratio")
    axis[2].scatter(v['log_ratio_over_random'], v['log_pvalue'], alpha=0.5)
    axis[2].axvline(0, linestyle="--", color="grey")
    axis[2].set_xlabel("log_ratio_over_random")
    axis[2].set_ylabel("-log10(p-value)")
    for pc in v.sort_values('log_ratio_over_random').tail(5).index:
        axis[2].text(
            v.loc[pc, 'log_ratio_over_random'], v.loc[pc, 'log_pvalue'],
            s="PC {}".format(pc), fontsize=4)
    fig.savefig(
        os.path.join(prefix + ".pc_significance_association.svg"), dpi=300, bbox_inches="tight")

    # # real vs random
    # # # scatter
    # # # distplot
    for label, (a, b) in [
            ("explained_variance_ratio", ("real_weight", "random_weight")),
            ("p_value", ("pvalue", "padj"))]:
        n = int(np.ceil(np.sqrt(unsup_results.index.levels[1].shape)))
        fig1, axis1 = plt.subplots(n, n, figsize=(n * 3, n * 3), sharex=True, sharey=True)
        fig2, axis2 = plt.subplots(n, n, figsize=(n * 3, n * 3), sharex=True, sharey=True)
        axes = [axis1.flatten(), axis2.flatten()]
        for i, pc in enumerate(unsup_results.index.levels[1]):
            for ax in axes:
                ax[i].set_title("PC {}".format(pc))
                ax[i].axvline(0, linestyle="--", color="grey")
                ax[i].axhline(0, linestyle="--", color="grey")
            p = unsup_results.loc[unsup_results.index.get_level_values("pc") == pc]
            axes[0][i].scatter(p[b], p[a], alpha=0.05, s=0.5, rasterized=True)
            sns.distplot(p[b], kde=False, label=b, ax=axes[1][i])
            sns.distplot(p[a], kde=False, label=a, ax=axes[1][i])
            axes[1][i].legend()
            if label != 'p_value':
                m = max(p[b].abs().max(), p[a].abs().max())
                axes[0][i].set_xlim((-m, m))
                axes[0][i].set_ylim((-m, m))
            else:
                axes[0][i].set_xlim((-0.1, 1.1))
                axes[0][i].set_ylim((-0.1, 1.1))
        for ax in axis1[:, 0]:
            ax.set_ylabel(a)
        for ax in axis1[-1, :]:
            ax.set_xlabel(b)
        for ax in axis2[:, 0]:
            ax.set_ylabel("Features")
        for ax in axis2[-1, :]:
            ax.set_xlabel("Weights")
        fig1.savefig(
            os.path.join(prefix + ".pc.{}.real_vs_random.scatter.svg".format(label)), dpi=300, bbox_inches="tight")
        fig2.savefig(
            os.path.join(prefix + ".pc.{}.real_vs_random.distplot.svg".format(label)), dpi=300, bbox_inches="tight")

    # # plot pc2 vs response
    fits = pd.read_csv(prefix + ".pca.fit.csv")
    fig, axis = plt.subplots(1, 1, figsize=(1 * 3, 1 * 3))
    axis.scatter(fits['2'], fits['response_at_120'])
    for i, r in fits.iterrows():
        axis.text(r['2'], r['response_at_120'], s=r['patient_id'])
    axis.set_xlabel("PC 2")
    axis.set_ylabel("Response at day 120")
    fig.savefig(
        os.path.join(prefix + ".pc_vs_response.scatter.svg"), dpi=300, bbox_inches="tight")

    # # # for significant PCs associated with response, get significant features
    from ngs_toolkit.general import plot_differential, differential_enrichment
    attrs = [x for x in plotting_attributes if x not in exclude]

    # # # now only PC2 of CLLs
    pc = 2
    alpha = 0.1
    cell_type = "CLL"
    p = "cll-time_course.20190111.accessibility_t0_cll"
    output_prefix = p + ".PC{}".format(pc)

    pp = unsup_results.loc[unsup_results.index.get_level_values('pc') == pc, :].reset_index(level=1)

    plot_differential(
        analysis, results=pp,
        output_dir=output_dir, output_prefix=output_prefix,
        mean_column='mean', log_fold_change_column='weight_over_random',
        p_value_column='pvalue', adjusted_p_value_column='padj', comparison_column="pc",
        matrix=analysis.matrix_norm, samples=[s for s in analysis.samples if s.cell_type == cell_type],
        robust=True, group_wise_colours=True, group_variables=attrs)

    # get significant
    # diff = pp.loc[pp['padj'] < alpha, :].copy()
    # diff['pc'] = diff['pc'].astype(str)
    # diff = diff.rename(columns={'real_weight': 'log2FoldChange', 'pc': 'comparison_name'})

    # since no significant, look at top/bottom 500
    pp['direction'] = (pp['weight_over_random'] > 0).replace(True, "up").replace(False, "down")
    diff = pp.groupby(['direction'])['pvalue'].nsmallest(500).reset_index().set_index('feature')
    diff = diff.rename(columns={'direction': 'comparison_name'})
    diff['comparison_name'] = "PC{} ".format(pc) + diff['comparison_name']
    diff['log2FoldChange'] = 0
    diff.to_csv(os.path.join(prefix + ".PC2.top_regions.csv"))
    # diff = pd.read_csv(os.path.join(prefix + ".PC2.top_regions.csv"), index_col=0)

    """
    # first let's change to hg38
    # # universe
    liftOver \
    results/cll-time_course_peak_set.bed \
    ~/resources/tools/liftOver/hg19ToHg38.over.chain \
    results/cll-time_course_peak_set.liftedover_hg38.bed \
    results/cll-time_course_peak_set.liftedover_hg38.unmapped.bed
    # # down regions
    liftOver \
    results/unsupervised_analysis_ATAC-seq/PC2\ down.all/cll-time_course.20190111.accessibility_t0_cll.PC2_regions.bed \
    ~/resources/tools/liftOver/hg19ToHg38.over.chain \
    results/unsupervised_analysis_ATAC-seq/PC2\ down.all/cll-time_course.20190111.accessibility_t0_cll.PC2_regions.liftedover_hg38.bed \
    results/unsupervised_analysis_ATAC-seq/PC2\ down.all/cll-time_course.20190111.accessibility_t0_cll.PC2_regions.liftedover_hg38.unmapped.bed
    # # up regions
    liftOver \
    results/unsupervised_analysis_ATAC-seq/PC2\ up.all/cll-time_course.20190111.accessibility_t0_cll.PC2_regions.bed \
    ~/resources/tools/liftOver/hg19ToHg38.over.chain \
    results/unsupervised_analysis_ATAC-seq/PC2\ up.all/cll-time_course.20190111.accessibility_t0_cll.PC2_regions.liftedover_hg38.bed \
    results/unsupervised_analysis_ATAC-seq/PC2\ up.all/cll-time_course.20190111.accessibility_t0_cll.PC2_regions.liftedover_hg38.unmapped.bed
    # # lola again

    sbatch \
    -p develop --mem 20000 -c 4 \
    -o results/unsupervised_analysis_ATAC-seq/PC2\ down.all/cll-time_course.20190111.accessibility_t0_cll.PC2_regions.liftedover_hg38.log \
    -J cll-time_course.20190111.accessibility_t0_cll.PC2_regions.liftedover_hg38.down \
    --wrap "Rscript ~/jobs/run_LOLA.R \
    results/unsupervised_analysis_ATAC-seq/PC2\ down.all/cll-time_course.20190111.accessibility_t0_cll.PC2_regions.liftedover_hg38.bed \
    results/cll-time_course_peak_set.liftedover_hg38.bed \
    hg38"
    sbatch \
    -p develop --mem 20000 -c 4 \
    -o results/unsupervised_analysis_ATAC-seq/PC2\ up.all/cll-time_course.20190111.accessibility_t0_cll.PC2_regions.liftedover_hg38.log \
    -J cll-time_course.20190111.accessibility_t0_cll.PC2_regions.liftedover_hg38.up \
    --wrap "Rscript ~/jobs/run_LOLA.R \
    results/unsupervised_analysis_ATAC-seq/PC2\ up.all/cll-time_course.20190111.accessibility_t0_cll.PC2_regions.liftedover_hg38.bed \
    results/cll-time_course_peak_set.liftedover_hg38.bed \
    hg38"
    """

    differential_enrichment(
        analysis, differential=diff, directional=False, genome=analysis.genome,
        output_dir=output_dir, output_prefix=output_prefix, as_jobs=1)

    collect_differential_enrichment(
        diff,
        directional=False,
        # steps=['homer_consensus'],
        data_type="ATAC-seq",
        output_dir=output_dir,
        input_prefix=output_prefix,
        output_prefix=output_prefix,
        permissive=True)

    enrichments_dir = os.path.join(output_dir, "enrichments_PC{}".format(pc))
    enrichment_table = pd.read_csv(os.path.join(
        output_dir, output_prefix + ".lola.csv"))
    plot_differential_enrichment(
        enrichment_table,
        "lola",
        data_type="ATAC-seq",
        direction_dependent=True,
        output_dir=enrichments_dir,
        comp_variable="comparison_name",
        output_prefix=p + ".PC3.top_8",
        top_n=8, z_score=1,
        rasterized=False,
        clustermap_metric="euclidean")
    enrichment_table = pd.read_csv(os.path.join(
        output_dir, output_prefix + ".enrichr.csv"))
    plot_differential_enrichment(
        enrichment_table,
        "enrichr",
        data_type="ATAC-seq",
        direction_dependent=True,
        output_dir=enrichments_dir,
        comp_variable="comparison_name",
        output_prefix=output_prefix,
        top_n=12, z_score=1)

    # observe accessibility of those regions
    g = sns.clustermap(
        analysis.matrix_norm.loc[
            diff.index,
            (analysis.matrix_norm.columns.get_level_values("cell_type") == "CLL") &
            (analysis.matrix_norm.columns.get_level_values("timepoint") == "000d")
        ].T,
        metric="correlation",
        robust=True,
        xticklabels=False, yticklabels=True,
        row_colors=pd.DataFrame(
            analysis.get_level_colors(analysis.matrix_norm.columns),
            index=analysis.matrix_norm.columns.names, columns=analysis.matrix_norm.columns).T,
        cbar_kws={"label": "Accessibility"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_ylabel("ATAC-seq samples")
    g.ax_heatmap.set_xlabel("PC2 regions (n={})".format(diff.shape[0]))
    g.ax_heatmap.get_children()[0].set_rasterized(True)
    g.savefig(
        os.path.join(prefix + ".pc_regions.t0_accessibility.clustermap.svg"), dpi=300, bbox_inches="tight")
    g = sns.clustermap(
        analysis.matrix_norm.loc[
            diff.index,
            (analysis.matrix_norm.columns.get_level_values("cell_type") == "CLL") &
            (analysis.matrix_norm.columns.get_level_values("timepoint") == "000d")
        ].T,
        metric="correlation", z_score=1, center=0,
        cmap="RdBu_r", robust=True,
        xticklabels=False, yticklabels=True,
        row_colors=pd.DataFrame(
            analysis.get_level_colors(analysis.matrix_norm.columns),
            index=analysis.matrix_norm.columns.names, columns=analysis.matrix_norm.columns).T,
        cbar_kws={"label": "Accessibility"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_ylabel("ATAC-seq samples")
    g.ax_heatmap.set_xlabel("PC2 regions (n={})".format(diff.shape[0]))
    g.ax_heatmap.get_children()[0].set_rasterized(True)
    g.savefig(
        os.path.join(prefix + ".pc_regions.t0_accessibility.clustermap.z_score.svg"), dpi=300, bbox_inches="tight")

    # # now for all samples
    g = sns.clustermap(
        analysis.matrix_norm.loc[
            diff.index,
            analysis.matrix_norm.columns.get_level_values("cell_type") == "CLL"
        ].T,
        metric="correlation",
        robust=True,
        xticklabels=False, yticklabels=True,
        row_colors=pd.DataFrame(
            analysis.get_level_colors(analysis.matrix_norm.columns),
            index=analysis.matrix_norm.columns.names, columns=analysis.matrix_norm.columns).T,
        cbar_kws={"label": "Accessibility"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_ylabel("ATAC-seq samples")
    g.ax_heatmap.set_xlabel("PC2 regions (n={})".format(diff.shape[0]))
    g.ax_heatmap.get_children()[0].set_rasterized(True)
    g.savefig(
        os.path.join(prefix + ".pc_regions.accessibility.clustermap.svg"), dpi=300, bbox_inches="tight")

    g = sns.clustermap(
        c,
        metric="correlation", z_score=1, center=0,
        cmap="RdBu_r", robust=True,
        xticklabels=False, yticklabels=True,
        row_colors=pd.DataFrame(
            analysis.get_level_colors(analysis.matrix_norm.columns),
            index=analysis.matrix_norm.columns.names, columns=analysis.matrix_norm.columns).T,
        cbar_kws={"label": "Accessibility (Z-score)"}, row_cluster=False)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_ylabel("ATAC-seq samples")
    g.ax_heatmap.set_xlabel("PC2 regions (n={})".format(diff.shape[0]))
    g.ax_heatmap.get_children()[0].set_rasterized(True)
    g.savefig(
        os.path.join(prefix + ".pc_regions.accessibility.clustermap.z_score.svg"), dpi=300, bbox_inches="tight")

    #
    #
    # 2. Second approach: regress the response directly
    from ngs_toolkit.general import least_squares_fit

    for cell_type in analysis.matrix_norm.columns.levels[3]:
        print(cell_type)
        # All samples from one cell type
        # analysis.matrix_norm_ct = analysis.matrix_norm.loc[
        #     :,
        #     (analysis.matrix_norm.columns.get_level_values("cell_type") == cell_type)]
        # exclude = ['cell_type', 'compartment', 'timepoint']

        # quant_matrix = analysis.matrix_norm_ct.T
        # design_matrix = analysis.matrix_norm_ct.columns.to_frame()
        # results = least_squares_fit(
        #     quant_matrix, design_matrix,
        #     test_model="~ response_at_120", null_model="~ 1", standardize_data=True,
        #     multiple_correction_method="fdr_bh")
        # results.to_csv(os.path.join("results", "response_association", "regression.{}.csv".format(cell_type)))

        # Only T0
        analysis.matrix_norm_t0_ct = analysis.matrix_norm.loc[
            :, (analysis.matrix_norm.columns.get_level_values("timepoint") == "000d") &
            (analysis.matrix_norm.columns.get_level_values("cell_type") == cell_type)]

        quant_matrix = analysis.matrix_norm_t0_ct.T
        design_matrix = analysis.matrix_norm_t0_ct.columns.to_frame()
        results = least_squares_fit(
            quant_matrix, design_matrix,
            test_model="~ response_at_120", null_model="~ 1", standardize_data=True,
            multiple_correction_method="fdr_bh")
        results.to_csv(os.path.join("results", "response_association", "regression.{}.t0_only.csv".format(cell_type)))


def response_day0_rna():
    analysis = ATACSeqAnalysis(
        from_pep=os.path.join("metadata", "project_config.yaml"))
    analysis.load_data()

    # look into CLL PC2
    output_dir = os.path.join(analysis.results_dir, "unsupervised_analysis_" + analysis.data_type)
    prefix = os.path.join(output_dir, analysis.name + ".20190111.accessibility_t0_cll")
    alpha = 0.1
    analysis.matrix_norm_t0_cll = analysis.matrix_norm.loc[
        :, (analysis.matrix_norm.columns.get_level_values("timepoint") == "000d") &
        (analysis.matrix_norm.columns.get_level_values("cell_type") == "CLL")]

    # # # now only PC2 of CLLs
    pc = 2
    alpha = 0.1
    cell_type = "CLL"
    p = "cll-time_course.20190111.accessibility_t0_cll"
    output_prefix = p + ".PC{}".format(pc)
    diff = pd.read_csv(os.path.join(prefix + ".PC2.top_regions.csv"), index_col=0)

    if not os.path.exists("gene.csv.gz"):
        gene = analysis.get_gene_level_matrix(analysis.matrix_norm.loc[diff.index])
        gene.to_csv("gene.csv.gz")
    else:
        gene = pd.read_csv("gene.csv.gz", index_col=0, header=list(range(18)))

    colors = analysis.get_level_colors(analysis.matrix_norm.columns, as_dataframe=True).T

    # plot
    # # for day 0 only
    g = sns.clustermap(
        gene.loc[:, (gene.columns.get_level_values("cell_type") == "CLL") & (gene.columns.get_level_values("timepoint") == "000d")].T,
        metric="correlation",
        robust=True,
        xticklabels=False, yticklabels=True,
        row_colors=colors.loc['response_at_120', :].T,
        cbar_kws={"label": "Accessibility"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_ylabel("ATAC-seq samples")
    g.ax_heatmap.set_xlabel("PC2 regions (n={})".format(diff.shape[0]))
    g.ax_heatmap.get_children()[0].set_rasterized(True)
    g.savefig(
        os.path.join(prefix + ".pc_regions.t0_accessibility.gene_level.clustermap.svg"), dpi=300, bbox_inches="tight")

    g = sns.clustermap(
        gene.loc[:, (gene.columns.get_level_values("cell_type") == "CLL") & (gene.columns.get_level_values("timepoint") == "000d")].T,
        metric="correlation", z_score=1, center=0,
        cmap="RdBu_r", robust=True,
        xticklabels=False, yticklabels=True,
        row_colors=colors.loc['response_at_120', :].T,
        cbar_kws={"label": "Accessibility (Z-score)"}, row_cluster=False)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_ylabel("ATAC-seq samples")
    g.ax_heatmap.set_xlabel("PC2 regions (n={})".format(diff.shape[0]))
    g.ax_heatmap.get_children()[0].set_rasterized(True)
    g.savefig(
        os.path.join(prefix + ".pc_regions.t0_accessibility.gene_level.clustermap.z_score.svg"), dpi=300, bbox_inches="tight")

    # # now for all samples
    g = sns.clustermap(
        gene.loc[:, gene.columns.get_level_values("cell_type") == "CLL"].T,
        metric="correlation",
        robust=True,
        # xticklabels=None, yticklabels=True,
        # row_colors=colors.loc['response_at_120', :].T,
        cbar_kws={"label": "Accessibility"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_ylabel("ATAC-seq samples")
    g.ax_heatmap.set_xlabel("PC2 regions (n={})".format(diff.shape[0]))
    g.ax_heatmap.get_children()[0].set_rasterized(True)
    g.savefig(
        os.path.join(prefix + ".pc_regions.accessibility.gene_level.clustermap.svg"), dpi=300, bbox_inches="tight")

    g = sns.clustermap(
        gene.loc[:, gene.columns.get_level_values("cell_type") == "CLL"].T,
        metric="correlation", z_score=1, center=0,
        cmap="RdBu_r", robust=True,
        yticklabels=True,
        row_colors=colors.loc['response_at_120', :].T,
        cbar_kws={"label": "Accessibility (Z-score)"}, row_cluster=False)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_ylabel("ATAC-seq samples")
    g.ax_heatmap.set_xlabel("PC2 regions (n={})".format(diff.shape[0]))
    g.ax_heatmap.get_children()[0].set_rasterized(True)
    g.savefig(
        os.path.join(prefix + ".pc_regions.accessibility.gene_level.clustermap.z_score.svg"), dpi=300, bbox_inches="tight")

    #
    # Now take those genes and see expression
    import scanpy as sc
    sc_prefix = "cll-time_course-scRNA-seq.all_samples.250-50_filter"
    adata = sc.read(sc_prefix + ".dca_denoised-zinb.processed.h5ad", cached=True)
    adata.obs = pd.read_csv(sc_prefix + ".dca_denoised-zinb.processed.obs.csv", index_col=0)

    cll = adata[adata.obs['cell_type'] == "CLL", :]

    # Get expression per timepoint per patient
    if not os.path.exists("single_cell.mean_expression.csv.gz"):
        d = pd.DataFrame(cll.X, index=cll.obs.index, columns=cll.var.index).join(cll.obs[['patient_id', 'timepoint']])
        dt = d.groupby(['patient_id', 'timepoint']).mean().T.sort_index()
        dt.to_csv("single_cell.mean_expression.csv.gz")
    else:
        dt = pd.read_csv("single_cell.mean_expression.csv.gz", index_col=0, header=[0, 1])

    # Read both
    gene = pd.read_csv("gene.csv.gz", index_col=0, header=list(range(18)))
    dt = pd.read_csv("single_cell.mean_expression.csv.gz", index_col=0, header=[0, 1])

    # get vector of strength of association
    diff.loc[:, 'log_p'] = -np.log10(diff['pvalue'])
    diff.loc[diff['comparison_name'].str.endswith("down"), 'log_p'] *= -1
    diff = diff.sort_values("log_p")

    diff_gene = analysis.get_gene_level_matrix(diff[['log_p', 'log2FoldChange']])
    gene_colors = pd.Series([tuple(x) for x in plt.get_cmap("RdBu_r")(diff_gene['log_p'])], diff_gene.index)

    # visualize expression
    g = sns.clustermap(
        dt.loc[gene.index, dt.columns.get_level_values("timepoint") == "000d"].dropna().T,
        metric="correlation",
        robust=True,
        yticklabels=True,
        # row_colors=colors.loc['response_at_120', :].T,
        col_colors=gene_colors,
        cbar_kws={"label": "Accessibility"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_ylabel("ATAC-seq samples")
    g.ax_heatmap.set_xlabel("Genes associated with PC2 regions (n={})".format(len(gene.index)))
    g.ax_heatmap.get_children()[0].set_rasterized(True)
    g.savefig(
        os.path.join(prefix + ".pc_regions.gene_expression.t0_expression.reduced_patient.clustermap.svg"), dpi=300, bbox_inches="tight")

    g = sns.clustermap(
        dt.loc[gene.index, :].dropna().T,
        metric="correlation", z_score=1, center=0,
        cmap="RdBu_r", robust=True,
        yticklabels=True,
        # row_colors=colors.loc['response_at_120', :].T,
        col_colors=gene_colors.loc[gene.index],
        # col_colors=gene_colors,
        cbar_kws={"label": "Expression (Z-score)"}, row_cluster=False)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_ylabel("ATAC-seq samples")
    g.ax_heatmap.set_xlabel("Genes associated with PC2 regions (n={})".format(len(gene.index)))
    g.ax_heatmap.get_children()[0].set_rasterized(True)
    g.savefig(
        os.path.join(prefix + ".pc_regions.gene_expression.reduced_patient.clustermap.z_score.svg"), dpi=300, bbox_inches="tight")

    g = sns.clustermap(
        dt.loc[gene.index, dt.columns.get_level_values("timepoint") == "000d"].dropna().T,
        metric="correlation", z_score=1, center=0,
        cmap="RdBu_r", robust=True,
        yticklabels=True,
        # row_colors=colors.loc['response_at_120', :].T,
        col_colors=gene_colors.loc[gene.index],
        # col_colors=gene_colors,
        cbar_kws={"label": "Expression (Z-score)"}, row_cluster=False)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_ylabel("ATAC-seq samples")
    g.ax_heatmap.set_xlabel("Genes associated with PC2 regions (n={})".format(len(gene.index)))
    g.ax_heatmap.get_children()[0].set_rasterized(True)
    g.savefig(
        os.path.join(prefix + ".pc_regions.gene_expression.t0_expression.reduced_patient.clustermap.z_score.svg"), dpi=300, bbox_inches="tight")


def inspect_coefficients():
    # export coefficients from R
    """
    library("nmslibR")
    loadRData <- function(fileName){
        load(fileName)
        get(ls()[ls() != "fileName"])
    }
    reg <- loadRData("results/single_cell_RNA/80_3_2_01_PredictTime6_CV/Reg.models.RData")

    writeMM(reg$CLL1$model$beta, "results/single_cell_RNA.80_3_2_01_PredictTime6_CV.Reg.models.CLL1.beta.mm")
    writeMM(reg$CLL5$model$beta, "results/single_cell_RNA.80_3_2_01_PredictTime6_CV.Reg.models.CLL5.beta.mm")
    writeMM(reg$CLL6$model$beta, "results/single_cell_RNA.80_3_2_01_PredictTime6_CV.Reg.models.CLL6.beta.mm")
    writeMM(reg$CLL7$model$beta, "results/single_cell_RNA.80_3_2_01_PredictTime6_CV.Reg.models.CLL7.beta.mm")

    write(colnames(reg$CLL7$model$beta), "results/single_cell_RNA.80_3_2_01_PredictTime6_CV.Reg.models.CLL7.beta.colnames.txt")
    write(rownames(reg$CLL7$model$beta), "results/single_cell_RNA.80_3_2_01_PredictTime6_CV.Reg.models.CLL7.beta.rownames.txt")
    """
    from scipy.io import mmread
    from scipy.stats import fisher_exact

    output_dir = os.path.join("results", "scrna_prediction")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    col = pd.read_csv("results/single_cell_RNA.80_3_2_01_PredictTime6_CV.Reg.models.CLL7.beta.colnames.txt", header=None, squeeze=True)
    row = pd.read_csv("results/single_cell_RNA.80_3_2_01_PredictTime6_CV.Reg.models.CLL7.beta.rownames.txt", header=None, squeeze=True)

    clls = ["CLL1", "CLL5", "CLL6", "CLL7"]
    lambdas = [20, 21, 18, 16]

    coefs = list()
    for i, (cll, lamb) in enumerate(zip(clls, lambdas)):
        d = pd.DataFrame(
            mmread("results/single_cell_RNA.80_3_2_01_PredictTime6_CV.Reg.models.{}.beta.mm".format(cll))
            .todense(), index=row, columns=col).iloc[:, lamb]
        coefs.append(d.rename(cll))
    coefs = pd.concat(coefs, 1)
    coefs.index.name = "gene_name"
    coefs = coefs.rename(columns={"CLL7": "CLL8"})

    # drop zeros everywhere
    coefs = coefs.loc[~(coefs == 0).all(1), :]

    # add info on how each gene is related with response
    annot = pd.read_csv(os.path.join("metadata", "annotation.csv"))
    response = annot.loc[:, ['patient_id', 'response_at_120']].drop_duplicates().dropna().set_index("patient_id").squeeze()
    corr = coefs.T.corrwith(response.loc[coefs.columns])

    # search for consistent coefficients
    mean = coefs.mean(1)
    sign = (mean > 0).astype(int).replace(0, -1)
    log_mean = mean.abs() ** (1 / 3) * sign
    std = coefs.std(1) / 2
    coefs = coefs.assign(
        corr=corr,
        mean=mean,
        median=coefs.median(1),
        log_mean=log_mean,
        abs_mean=mean.abs(),
        std=std * 2, v=std,
        diff=mean.abs() - std).sort_values("diff")

    coefs.to_csv(os.path.join(output_dir, "regression_coefficients.csv"))
    coefs = coefs.query("(diff < -1e-3) | (diff > 1e-3)")

    # get colors for the genes
    div = plt.get_cmap("RdBu_r")
    cont = plt.get_cmap("plasma")
    norm0 = matplotlib.colors.Normalize(vmin=0, vmax=2)
    norm1 = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    norm2 = matplotlib.colors.Normalize(vmin=-0.3, vmax=0.3)
    gene_color_dataframe = pd.DataFrame(
        [[tuple(x) for x in div(norm1(coefs['corr'].values))],
         [tuple(x) for x in div(norm2(coefs['log_mean'].values))],
         [tuple(x) for x in div(norm2(coefs['median'].values))],
         [tuple(x) for x in cont(norm0(coefs['std'].values))],
         [tuple(x) for x in div(norm1(coefs['diff'].values))]],
        index=['corr', 'log_mean', 'median', 'std', 'diff'], columns=coefs.index).T

    # Display all coefficients
    grid = sns.clustermap(
        coefs.loc[:, coefs.columns.str.startswith("CLL")],
        cmap="RdBu_r", center=0, robust=True, rasterized=True, metric="correlation",
        cbar_kws={"label": "Regression coefficient"},
        row_colors=gene_color_dataframe)
    grid.savefig(
        os.path.join("results", "scrna_prediction", "regression_coefficients.clustermap.svg"),
        bbox_inches="tight", dpi=300)

    # Display how patients relate to each other
    grid = sns.clustermap(
        coefs.loc[:, coefs.columns.str.startswith("CLL")].corr(),
        cmap="RdBu_r", center=0, robust=True, rasterized=True, metric="correlation",
        cbar_kws={"label": "Correlation of regression coefficient"})
    grid.savefig(
        os.path.join("results", "scrna_prediction", "regression_coefficients.patient_correlation.clustermap.svg"),
        bbox_inches="tight", dpi=300)

    # Display how genes relate to each other
    grid = sns.clustermap(
        coefs.loc[:, coefs.columns.str.startswith("CLL")].T.corr(),
        cmap="RdBu_r", center=0, robust=True, rasterized=True, metric="correlation",
        cbar_kws={"label": "Correlation of regression coefficient"},
        row_colors=gene_color_dataframe)
    grid.savefig(
        os.path.join("results", "scrna_prediction", "regression_coefficients.gene_correlation.clustermap.svg"),
        bbox_inches="tight", dpi=300)

    # Select strongest, consistent, variable, correlated genes across patients
    strength_thresold = 1
    consistence_threshold = 0.25
    variability_threshold = 1.5
    correlation_threshold = 0.90
    attrs = ["strong", "consistent", "variable", "correlated"]
    coefs = coefs.assign(
        strong=coefs['mean'].abs() > strength_thresold,
        consistent=coefs['diff'] > consistence_threshold,
        variable=coefs['std'] > variability_threshold,
        correlated=coefs['corr'].abs() > correlation_threshold)
    coefs.to_csv(os.path.join(output_dir, "regression_coefficients.filtered.classified.csv"))

    # Illustration of selection procedure
    fig, axis = plt.subplots(1, 7, figsize=(7 * 4, 4))
    v = coefs['log_mean'].abs().max()
    v += v * 0.1
    axis[0].scatter(
        coefs['log_mean'], np.log(coefs['v']),
        c=coefs['corr'], cmap="coolwarm",
        alpha=0.2, rasterized=True)
    axis[0].set_xlabel("Coefficient mean (log)")
    axis[0].set_ylabel("Coefficient standard deviation (log)")
    axis[0].set_xlim(-v, v)
    axis[0].axvline(0, linestyle="--", color="grey")
    axis[0].axhline(np.log(variability_threshold), linestyle="--", color="black")
    axis[0].axvline(-np.log(strength_thresold), linestyle="--", color="black")
    vv = coefs.query("strong == True")['log_mean'].abs().min()
    axis[0].axvline(vv, linestyle="--", color="black")
    axis[0].axvline(-vv, linestyle="--", color="black")

    axis[1].scatter(
        coefs['log_mean'], coefs['diff'],
        c=coefs['corr'], cmap="coolwarm",
        alpha=0.2, rasterized=True)
    axis[1].set_xlabel("Coefficient log_mean")
    axis[1].set_ylabel("Deviation from expected")
    axis[1].set_xlim(-v, v)
    axis[1].axhline(consistence_threshold, linestyle="--", color="black")
    axis[1].axhline(0, linestyle="--", color="grey")

    axis[2].scatter(
        coefs['log_mean'], coefs['corr'],
        c=coefs['corr'], cmap="coolwarm",
        alpha=0.2, rasterized=True)
    axis[2].set_xlabel("Coefficient log_mean")
    axis[2].set_ylabel("Correlation to response")
    axis[2].set_xlim(-v, v)
    axis[2].axhline(-correlation_threshold, linestyle="--", color="black")
    axis[2].axhline(correlation_threshold, linestyle="--", color="black")
    axis[2].axhline(0, linestyle="--", color="grey")

    for i, label in enumerate(attrs):
        p = coefs.query(f"{label} == True")
        axis[3 + i].set_title(f"Most {label} genes")
        axis[3 + i].scatter(
            p['corr'].rank(), p['corr'],
            c=p['corr'], cmap="coolwarm",
            alpha=0.9, s=5, rasterized=True)
        axis[3 + i].set_xlabel("Correlation to response (rank)")
        axis[3 + i].set_ylabel("Correlation to response")
        axis[3 + i].axhline(0, linestyle="--", color="grey")
        r = p['corr'].rank()
        for l in p.index:
            axis[3 + i].text(r.loc[l], p.loc[l, 'corr'], s=l, fontsize=5, rotation=90, va="bottom")

    fig.savefig(
        os.path.join("results", "scrna_prediction", "regression_coefficients.inspection.scatter.svg"),
        bbox_inches="tight", dpi=300)

    # Display top genes
    for i, label in enumerate(attrs):
        p = coefs.query(f"{label} == True")
        grid = sns.clustermap(
            p.loc[:, p.columns.str.startswith("CLL")],
            cmap="RdBu_r", center=0, robust=True, rasterized=True, metric="correlation",
            cbar_kws={"label": "Regression coefficient"},
            row_colors=gene_color_dataframe)
        grid.savefig(
            os.path.join("results", "scrna_prediction", f"regression_coefficients.{label}.top.clustermap.svg"),
            bbox_inches="tight", dpi=300)

        # # display relationship between genes
        grid = sns.clustermap(
            p.loc[:, p.columns.str.startswith("CLL")].T.corr(),
            cmap="RdBu_r", center=0, robust=True, rasterized=True, metric="correlation",
            cbar_kws={"label": "Correlation of regression coefficient"},
            row_colors=gene_color_dataframe,
            xticklabels=False)
        grid.savefig(
            os.path.join("results", "scrna_prediction", f"regression_coefficients.{label}.top.gene_correlation.clustermap.svg"),
            bbox_inches="tight", dpi=300)

    #
    # See what fraction of genes is differentially expressed

    # # load diff expression
    scrna = pd.read_csv(os.path.join(
        "results", "single_cell_RNA", "13_4_Overtime_nUMI_Cutoff", "SigGenes_overTime.tsv"), sep="\t")
    all_genes = scrna['gene'].drop_duplicates()
    scrna2 = scrna.loc[scrna['qvalue'] < 0.05]
    scrna2.groupby('cellType')['gene'].nunique()
    scrna2['log_pvalue'] = -np.log10(scrna2['pvalue'])
    scrna2 = scrna_diff2.rename(
        columns={"logFC": "log2FoldChange", "cellType": "comparison_name", "gene": "gene_name"})
    diff = scrna2.query("(comparison_name == 'CLL') and (qvalue < 0.05)")

    # # load offtarget signature
    atac_sig = pd.read_csv(os.path.join("results", "offtarget_signature.csv"), index_col=0, header=None, squeeze=True)

    # # load ATAC-seq day 0 signature
    atac_d0_genes = pd.read_csv("gene.csv.gz", index_col=0, header=list(range(18)))

    # # assemble all gene sets
    gene_sets = {
        "differential_expression": {
            "pos": set(diff['gene_name']),
            "neg": set(all_genes[~all_genes.isin(diff['gene_name'])])},
        "cross_cell_type_signature": {
            "pos": set(atac_sig.index),
            "neg": set(all_genes[~all_genes.isin(atac_sig.index)])},
        "atacseq_d0_signature": {
            "pos": set(atac_d0_genes.index),
            "neg": set(all_genes[~all_genes.isin(atac_d0_genes.index)])},
        "regression_coefficients": {
            "pos": set(coefs.index),
            "neg": set(all_genes[~all_genes.isin(coefs.index)])}}
    for label in attrs:
        p = coefs.query(f"{label} == True")
        gene_sets["regression_coefficients-" + label] = {
            "pos": set(p.index),
            "neg": set(all_genes[~all_genes.isin(p.index)])}

    # # test
    all_tests = list()
    for s1 in gene_sets:
        for s2 in gene_sets:
            s1_ = gene_sets[s1]
            s1_pos, s1_neg = s1_['pos'], s1_['neg']
            s2_ = gene_sets[s2]
            s2_pos, s2_neg = s2_['pos'], s2_['neg']
            # # diff expressed and in coefficients
            a = len(s1_pos.intersection(s2_pos))
            # # not diff expressed but in coefficients
            b = len(s2_neg.intersection(s1_pos))
            # # diff expressed and not in coefficients
            c = len(s2_pos.intersection(s1_neg))
            # # not diff expressed and not in coefficients
            # d = (~s2_neg.isin(s1_neg)).sum()
            d = len(set(s1_neg).union(set(s2_neg)))
            odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            s = pd.Series(
                [s1, s2, 'fisher_exact', a, b, c, d, odds, p],
                name=f'{s1}_vs_{s2}',
                index=['set1', 'set2', 'test', 'a', 'b', 'c', 'd', "odds_ratio", "p_value"]).to_frame().T
            all_tests.append(s)

    s = pd.concat(all_tests).infer_objects()

    # Correct p-values
    # from statsmodels.stats.multitest import multipletests
    # s = s.assign(p_adj=multipletests(s['p_value'], method="fdr_bh")[1])
    s = s.assign(p_adj=s['p_value'] * s.shape[0])  # for some reason statsmodels capped the p-values
    s.loc[s['p_adj'] > 1, 'p_adj'] = 1
    s.index.name = "comparison"
    s.to_csv(os.path.join("results", "scrna_prediction", "regression_coefficients.overlap_with_all_signatures.csv"))

    # # plot pairwise relationship of signature
    # s = s.loc[s['odds_ratio'] != np.inf]
    s = s.query("set1 != set2")
    odds = np.log2(s.pivot_table(index="set1", columns="set2", values="odds_ratio"))
    p = -np.log10(s.pivot_table(index="set1", columns="set2", values="p_adj"))

    fig, axis = plt.subplots(1, 2, figsize=(2 * 5, 5))
    v = odds.max().replace(np.inf, 0).max()
    v += v * 0.1
    sns.heatmap(
        odds, cmap="RdBu_r", vmin=-v, vmax=v,
        ax=axis[0], cbar_kws={"label": "log2(Odds ratio)"}, square=True, annot=True, fmt=".2f")
    v = p.max().max()
    sns.heatmap(
        p, vmax=v + v * 0.1,
        ax=axis[1], cbar_kws={"label": "-log10(Bonferroni p-value)"}, square=True, annot=True, fmt=".2f",
        yticklabels=False)
    fig.savefig(
        os.path.join("results", "scrna_prediction", "regression_coefficients.overlap_with_all_signatures.heatmap.svg"),
        bbox_inches="tight")

    #
    # Ilustrate differential expression between patients

    # # read single cell expression
    prefix = "cll-time_course-scRNA-seq.all_samples.250-50_filter"
    adata = sc.read(prefix + ".dca_denoised-zinb.processed.h5ad")
    adata.obs = pd.read_csv(prefix + ".dca_denoised-zinb.processed.obs.csv", index_col=0)
    adata.obs = adata.obs.assign(sample=adata.obs['patient_id'] + "_" + adata.obs['timepoint'])

    # inspect mean across cells
    gene = "CD24"
    x = pd.Series(adata[:, gene].X, index=adata.obs.index, name=gene)
    adata.obs.join(x).query("cell_type == 'CLL'").groupby(['patient_id', 'timepoint'])[gene].mean()

    genes = ["CD24", "B2M", "CXCR4", "BTG1", "TNF"]
    ngenes = len(genes)
    fig, axis = plt.subplots(1, ngenes, figsize=(3 * ngenes, 3))
    for i, gene in enumerate(genes):
        ax = sc.pl.violin(adata, gene, groupby="sample", use_raw=False, ax=axis[i], show=False, rasterized=True)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    fig.savefig(os.path.join("results", "scrna_prediction", "regression_coefficients.marker_illustration.svg"), dpi=300, bbox_inches="tight")


def revisit_normalization():
    import os

    from ngs_toolkit import ATACSeqAnalysis
    from ngs_toolkit.utils import normalize_quantiles_r
    from combat import combat

    analysis = ATACSeqAnalysis(from_pep=os.path.join("metadata", "project_config.yaml"))
    cell_types = list(set([sample.cell_type for sample in analysis.samples]))

    ss = [
        "ATAC-seq_CLL4_0d_CD8",
        "ATAC-seq_CLL4_0d_NK",
        "ATAC-seq_CLL2_0d_CD4",
        "ATAC-seq_CLL2_0d_CD8",
    ]
    for sn in ss:
        s = [s for s in analysis.samples if s.name == sn][0]
        s.good_batch = 'FALSE'

    for sn in ["ATAC-seq_CLL2_30d_NK"]:
        s = [s for s in analysis.samples if s.name == sn][0]
        s.good_batch = 'TRUE'

    for s in analysis.samples:
        if not hasattr(s, "good_batch"):
            s.good_batch = 'TRUE'
        if pd.isnull(s.good_batch):
            s.good_batch = 'TRUE'

    analysis.load_data(only_these_keys=['matrix_raw', 'matrix_norm'])

    analysis.v = analysis.matrix_norm
    analysis.v = analysis.v.loc[analysis.v.index.str.startswith("chr")]
    analysis.v.columns = analysis.v.columns.get_level_values("sample_name")
    analysis.matrix_norm = analysis.annotate_samples(matrix="v", save=False)

    # data normalization
    # # first quantile normalize as before
    no_norm = np.log2(1 + pd.DataFrame(
        normalize_quantiles_r(analysis.matrix_raw.values),
        index=analysis.matrix_raw.index,
        columns=analysis.matrix_raw.columns
    ))

    # # previous PCA-based method
    matrix_pca = analysis.matrix_norm.copy()

    # Combat on batch
    batch = pd.Series([s.batch for s in analysis.samples], index=[s.name for s in analysis.samples])
    matrix_combat = combat(data=no_norm, batch=batch)

    # Combat on latent vector directly
    batch_lv = pd.Series([s.good_batch for s in analysis.samples], index=[s.name for s in analysis.samples])
    matrix_combat_lv = combat(data=no_norm, batch=batch_lv)

    # Do unsupervised analysis with all
    output_dir = os.path.join("results", "revisit_normalization")
    kwargs = {
        "output_dir": output_dir,
        "steps": ["pca", "pca_association"],
        "attributes_to_plot": ['cell_type', 'patient_id', 'timepoint', 'batch', 'good_batch'],
        "plot_max_pcs": 4,
        "output_dir": output_dir
    }
    analysis.unsupervised_analysis(
        matrix=no_norm,
        output_prefix="No_norm-2",
        **kwargs)
    analysis.unsupervised_analysis(
        matrix=matrix_pca,
        output_prefix="PCA-2",
        **kwargs)
    analysis.unsupervised_analysis(
        matrix=matrix_combat,
        output_prefix="Combat-2",
        **kwargs)
    analysis.unsupervised_analysis(
        matrix=matrix_combat_lv,
        output_prefix="Combat_lv-2",
        **kwargs)

    # replot again with less factors just for better illustration
    kwargs = {
        "output_dir": output_dir,
        "steps": ["pca"],
        "attributes_to_plot": ['cell_type', 'batch', 'good_batch'],
        "plot_max_pcs": 4,
        "output_dir": output_dir
    }
    analysis.unsupervised_analysis(
        matrix=matrix_pca,
        output_prefix="PCA.viz",
        **kwargs)
    analysis.unsupervised_analysis(
        matrix=matrix_combat_lv,
        output_prefix="Combat.viz",
        **kwargs)

    associations = dict()
    for i, prefix in enumerate(
            ["No_norm-2", "PCA-2", "Combat_lv-2"],
            start=1):
        associations["{}. {}".format(i, prefix)] = pd.read_csv(os.path.join(
            output_dir,
            "cll-time_course.{}.pca.variable_principle_components_association.csv".format(prefix)))
    associations = pd.concat(associations)
    associations.loc[:, 'attribute'] = associations.loc[:, 'attribute'].replace("good_batch", "latent_factor")
    associations.to_csv(os.path.join(output_dir, "cross_method.associations.csv"))

    #
    fig, axis = plt.subplots(2, 2, figsize=(4 * 4, 2 * 4))
    for i, var in enumerate(['p_value', 'adj_pvalue']):
        p = associations.query("pc <= 12").reset_index().pivot_table(columns="pc", index=['level_0', 'attribute'], values=var, aggfunc=min)
        p_masked = (p >= 0.01).astype(int)

        sns.heatmap(
            -np.log10(p),
            ax=axis[i, 0], square=True, cbar_kws={"label": "-log10(adjusted p-value)"})
        sns.heatmap(
            p_masked,
            ax=axis[i, 1], square=True, cbar_kws={"label": "-log10(adjusted p-value)"})
    fig.savefig(os.path.join(output_dir, "cross_method.associations.heatmap.svg"), bbox_inches="tight", dpi=300)


if __name__ == '__main__':
    import sys
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)


# analysis = ATACSeqAnalysis(from_pep=os.path.join("metadata", "project_config.yaml"))
# cell_types = list(set([sample.cell_type for sample in analysis.samples]))

# analysis.load_data()

# analysis.set_matrix(
#     matrix_name="matrix_raw",
#     csv_file=os.path.join("results", "cll-time_course_peaks.raw_coverage.csv"),
#     index_col=0)
# analysis.set_matrix(
#     matrix_name="matrix_norm",
#     csv_file=os.path.join("results", "cll-time_course.accessibility.annotated_metadata.csv"),
#     header=list(range(len(analysis.sample_attributes))), index_col=0)
# analysis.set_matrix(
#     matrix_name="gene_annotation",
#     csv_file=os.path.join("results", "cll-time_course_peaks.gene_annotation.csv"),
#     index_col=None)
# from ngs_toolkit.utils import bed_to_index
# analysis.gene_annotation.index = bed_to_index(analysis.gene_annotation)
# gp_output_dir = os.path.join(analysis.results_dir, "gp_fits")
# mohgp_output_dir = os.path.join(analysis.results_dir, "mohgp_fits")
# prefix = "accessibility.qnorm_pcafix_cuberoot"
# analysis.fits = pd.read_csv(os.path.join(
#     gp_output_dir, prefix + ".GP_fits.all_cell_types.csv"), index_col=0)


# for a in [
# 'coverage',
# 'coverage_annotated',
# 'accessibility',
# 'fits',
# 'support',
# 'sites',
# 'coverage_qnorm',
# 'region_annotation',
# 'chrom_state_annotation_b',
# 'region_annotation_b',
# 'assignments',
# 'gene_annotation',
# 'chrom_state_annotation']:
#     setattr(analysis, a, getattr(analysis2, a))


# order = [
# "ATAC-seq_CLL1_120d_CLL",
# "ATAC-seq_CLL7_30d_CLL",
# "ATAC-seq_CLL7_1d_CLL",
# "ATAC-seq_CLL7_8d_CLL",
# "ATAC-seq_CLL7_3d_CLL",
# "ATAC-seq_CLL7_150d_CLL",
# "ATAC-seq_CLL7_0d_CLL",
# "ATAC-seq_CLL6_0d_CLL",
# "ATAC-seq_CLL2_0d_CLL",
# "ATAC-seq_CLL2_3d_CLL",
# "ATAC-seq_CLL6_280d_CLL",
# "ATAC-seq_CLL6_120d_CLL",
# "ATAC-seq_CLL6_030d_CLL",
# "ATAC-seq_CLL2_120d_CLL",
# "ATAC-seq_CLL1_030d_CLL",
# "ATAC-seq_CLL1_3d_CLL",
# "ATAC-seq_CLL2_30d_CLL",
# "ATAC-seq_CLL2_8d_CLL",
# "ATAC-seq_CLL8_030d_CLL",
# "ATAC-seq_CLL8_3d_CLL",
# "ATAC-seq_CLL5_240d_CLL",
# "ATAC-seq_CLL5_2d_CLL",
# "ATAC-seq_CLL5_1d_CLL",
# "ATAC-seq_CLL5_3d_CLL",
# "ATAC-seq_CLL5_0d_CLL",
# "ATAC-seq_CLL4_0d_CLL",
# "ATAC-seq_CLL8_0d_CLL",
# "ATAC-seq_CLL5_30d_CLL",
# "ATAC-seq_CLL1_0d_CLL",
# "ATAC-seq_CLL5_150d_CLL",
# "ATAC-seq_CLL4_3d_CLL",
# "ATAC-seq_CLL4_2d_CLL",
# "ATAC-seq_CLL5_8d_CLL",
# ]
