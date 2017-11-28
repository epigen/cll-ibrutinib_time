#!/usr/bin/env python

"""
cll-time_course
"""

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')

import multiprocessing
import os
import random
import string
import time

import GPy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import parmap
import pybedtools
import seaborn as sns
import tqdm
from scipy.cluster.hierarchy import fcluster
from scipy.stats import pearsonr, spearmanr, zscore
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder
from statsmodels.stats.multitest import multipletests

from looper.models import Project, Sample
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.general import (collect_differential_enrichment,
                                 differential_enrichment,
                                 normalize_quantiles_r,
                                 plot_differential_enrichment,
                                 subtract_principal_component)

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
    prj = Project(os.path.join("metadata", "project_config.yaml"))
    prj.samples = [sample for sample in prj.samples if sample.library == "ATAC-seq"]
    for sample in prj.samples:
        sample.filtered = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
        sample.peaks = os.path.join(sample.paths.sample_root, "peaks", sample.name + "_peaks.narrowPeak")

        for attr in sample.__dict__.keys():
            if type(getattr(sample, attr)) is str:
                if getattr(sample, attr) == "nan":
                    setattr(sample, attr, np.nan)

    analysis = ATACSeqAnalysis(name="cll-time_course", prj=prj, samples=prj.samples)

    # Sample's attributes
    sample_attributes = ["sample_name", "patient_id", "timepoint", "cell_type", "compartment", "response", "cell_number", "batch"]
    numerical_attributes = ["cell_number"]
    plotting_attributes = ["patient_id", "timepoint", "cell_type", "compartment", "batch"]
    cell_types = list(set([sample.cell_type for sample in analysis.samples]))


    # Get accessibility matrix with normalization
    # generate consensus peak set from all samples, cell types together
    analysis.get_consensus_sites()

    # annotate consensus peak set
    analysis.get_peak_gene_annotation()
    analysis.get_peak_genomic_location()
    analysis.get_peak_chromatin_state(chrom_state_file=os.path.join("data", "external", "E032_15_coreMarks_mnemonics.bed"))
    analysis.calculate_peak_support(region_type="summits")

    # get coverage values for each region in each sample
    analysis.measure_coverage()
    analysis.coverage = analysis.coverage[~analysis.coverage.index.str.contains("chrX|chrY")]

    # data normalization
    analysis.accessibility = data_normalization(analysis)

    # annotate matrix with sample metadata and save
    analysis.accessibility = analysis.annotate_with_sample_metadata(quant_matrix="accessibility", attributes=sample_attributes + ["good_batch"], save=True, assign=False)
    analysis.to_pickle()

    analysis.annotate(quant_matrix="coverage")
    analysis.to_pickle()


    # Unsupervised analysis
    unsupervised(analysis, quant_matrix="accessibility", samples=None, attributes_to_plot=plotting_attributes + ["good_batch"], plot_prefix="accessibility")
    for cell_type in cell_types:
        unsupervised(
            analysis, quant_matrix="accessibility", samples=[s for s in analysis.samples if (s.cell_type == cell_type) ],
            attributes_to_plot=['patient_id', 'timepoint', 'batch'] + ["good_batch"], plot_prefix="accessibility_{}_only".format(cell_type))


    # Time Series Analysis
    matrix_file = os.path.abspath(os.path.join("results", analysis.name + ".accessibility.annotated_metadata.csv"))
    cell_types = ['Bcell', 'CD4', 'CD8', 'CLL', 'NK', 'Mono']
    gp_output_dir = os.path.join(analysis.results_dir, "gp_fits")
    mohgp_output_dir = os.path.join(analysis.results_dir, "mohgp_fits")

    # fit GPs with varying and constant kernel to detect variable regulatory elements
    prefix = "accessibility.qnorm_pcafix_cuberoot"
    fit_gaussian_processes(
        analysis.accessibility,
        cell_types=cell_types,
        matrix_file=matrix_file,
        prefix=prefix)  # wait for jobs to complete

    analysis.fits = gather_gaussian_processes(
        analysis.accessibility,
        matrix_file=matrix_file,
        prefix=prefix)
    analysis.fits.to_csv(os.path.join(gp_output_dir, prefix + ".GP_fits.all_cell_types.csv"), index=True)
    analysis.fits = pd.read_csv(os.path.join(gp_output_dir, prefix + ".GP_fits.all_cell_types.csv"), index_col=0)

    visualize_gaussian_process_fits(analysis, fits, output_dir=gp_output_dir, prefix=prefix)

    # cluster variable regulatory elements with hierarchical mixtures of GPs (MOHGP)
    for alpha in [0.01, 0.05]:
        clust_prefix = prefix + ".p={}".format(alpha)
        fit_MOHGP(
            analysis.accessibility,
            matrix_file=os.path.abspath(os.path.join(
                "results", analysis.name + ".accessibility.annotated_metadata.csv")),
            fits_file=os.path.join(gp_output_dir, prefix + ".GP_fits.all_cell_types.csv"),
            cell_types=cell_types,
            n_clust=[3, 5, 4, 4, 4, 3],
            prefix=clust_prefix, output_dir=mohgp_output_dir, alpha=alpha)
        # wait for jobs to finish

        analysis.assignments = gather_MOHGP(
            analysis.accessibility, fits,
            cell_types=cell_types,
            n_clust=[3, 5, 4, 4, 4, 3],
            prefix=clust_prefix, fits_dir=mohgp_output_dir, alpha=alpha, posterior_threshold=0.8)
        analysis.assignments.to_csv(os.path.join(mohgp_output_dir, clust_prefix + ".GP_fits.mohgp_clusters.csv"), index=True)
        analysis.assignments = pd.read_csv(os.path.join(mohgp_output_dir, clust_prefix + ".GP_fits.mohgp_clusters.csv"), index_col=0)

        visualize_clustering_assignments(
            analysis.accessibility, analysis.assignments, prefix=clust_prefix,
            output_dir=mohgp_output_dir)

        # export clusters at gene level
        clusters_to_signatures(assignments)


    # Linear time
    gp_linear(analysis)
    # Randomized samples
    gp_random(analysis)


    # Plot distribution of clusters per cell type dependent on their dynamic pattern
    gp_output_dir = os.path.join(analysis.results_dir, "gp_fits")
    mohgp_output_dir = os.path.join(analysis.results_dir, "mohgp_fits")
    alpha = 0.05
    prefix = "accessibility.qnorm_pcafix_cuberoot" + ".p={}".format(alpha)
    assignments = pd.read_csv(os.path.join(mohgp_output_dir, prefix + ".GP_fits.mohgp_clusters.csv"))
    l_assignments = cluster_stats(assignments, prefix=prefix)


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

    collect_differential_enrichment(
        l_assignments,
        directional=False,
        data_type="ATAC-seq",
        output_dir=enrichments_dir,
        output_prefix=prefix,
        permissive=True)
    
    for simple, label in [(False, ""), (True, "-simple")]:
        enrichment_table = pd.read_csv(os.path.join(
            enrichments_dir, prefix + ".meme_ame.csv"))
        if simple:
            enrichment_table = enrichment_table[enrichment_table['comparison_name'].str.contains("_")]
        plot_differential_enrichment(
            enrichment_table,
            "motif",
            data_type="ATAC-seq",
            direction_dependent=False,
            output_dir=enrichments_dir,
            comp_variable="comparison_name",
            output_prefix=prefix + label,
            top_n=5)
        enrichment_table = pd.read_csv(os.path.join(
            enrichments_dir, prefix + ".lola.csv"))
        if simple:
            enrichment_table = enrichment_table[enrichment_table['comparison_name'].str.contains("_")]
        plot_differential_enrichment(
            enrichment_table,
            "lola",
            data_type="ATAC-seq",
            direction_dependent=False,
            output_dir=enrichments_dir,
            comp_variable="comparison_name",
            output_prefix=prefix + label,
            top_n=5)
        enrichment_table = pd.read_csv(os.path.join(
            enrichments_dir, prefix + ".enrichr.csv"))
        if simple:
            enrichment_table = enrichment_table[enrichment_table['comparison_name'].str.contains("_")]
        plot_differential_enrichment(
            enrichment_table,
            "enrichr",
            data_type="ATAC-seq",
            direction_dependent=False,
            output_dir=enrichments_dir,
            comp_variable="comparison_name",
            output_prefix=prefix + label,
            top_n=5)

    # each cell type independently
    for cell_type in cell_types:
        enrichment_table = pd.read_csv(os.path.join(
            enrichments_dir, prefix + ".lola.csv"))
        enrichment_table = enrichment_table[enrichment_table['comparison_name'].str.contains("{}_".format(cell_type))]
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
        enrichment_table = enrichment_table[enrichment_table['comparison_name'].str.contains("{}_".format(cell_type))]
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

    # Indentity analysis
    correlation_to_bcell(analysis)
    differentiation_assessment(analysis)

    # Microenvironment analysis
    cytokine_receptor_repertoire(analysis)
    cytokine_interplay(analysis)


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
    to_norm = analysis.coverage.drop(['chrom', 'start', 'end'], axis=1)
    counts_qnorm = pd.DataFrame(
        normalize_quantiles_r(to_norm.values),
        index=to_norm.index,
        columns=to_norm.columns
    )
    # normalize batch effect
    counts_qnorm_pcafix = subtract_principal_component_by_attribute(counts_qnorm.T, attributes=cell_types[:-1], pcs=[1] * 5).T
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

    for cell_type in tqdm.tqdm(cell_types, desc="cell_type", dynamic_ncols=True):
        for start, end in tqdm.tqdm(zip(r, r[1:]) + [(r[-1], matrix.shape[0])], desc="chunk", dynamic_ncols=True):
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
    for cell_type in tqdm.tqdm(cell_types, desc="cell_type"):
        for start, end in tqdm.tqdm(zip(r, r[1:]) + [(r[-1], matrix.shape[0])], desc="chunk"):
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
    example_acc = analysis.accessibility.loc[examples, :]
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
            g.savefig(os.path.join(output_dir, prefix + "." + cell_type +
                                ".mohgp.fitted_model.clustermap.cluster_labels.{}.with_posterior_probs.svg".format(label)), dpi=300, bbox_inches="tight")

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
            g.savefig(os.path.join(output_dir, prefix + "." + cell_type +
                                ".mohgp.fitted_model.clustermap.cluster_labels.{}.only_posterior_above_threshold.svg".format(label)), dpi=300, bbox_inches="tight")

            # only variable and with assignments above threshold: mean per timepoint
            matrix_mean = matrix2.loc[regions_assigned.index,
                                    tp.index].T.groupby(level="timepoint").mean()
            background_mean = background.loc[:,
                                    tp.index].T.groupby(level="timepoint").mean()
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
            g.savefig(os.path.join(output_dir, prefix + "." + cell_type +
                                ".mohgp.fitted_model.mean_acc.clustermap.cluster_labels.{}.only_posterior_above_threshold.svg".format(label)), dpi=300, bbox_inches="tight")

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
                g.savefig(os.path.join(output_dir, prefix + "." + cell_type +
                                    ".mohgp.fitted_model.mean_acc.clustermap.cluster_labels.{}.{}only_posterior_above_threshold.svg".format(label, label2)), dpi=300, bbox_inches="tight")

                g = sns.clustermap(
                    cluster_matrix_mean_norm,
                    row_cluster=False, col_cluster=col_cluster, z_score=z_score, xticklabels=True, yticklabels=True,
                    rasterized=True, square=True, metric="correlation", robust=True, cbar_kws={"label": cbar_label})
                g.ax_heatmap.set_xticklabels(
                    g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
                g.ax_heatmap.set_yticklabels(
                    g.ax_heatmap.get_yticklabels(), rotation=0)
                g.ax_col_dendrogram.set_rasterized(True)
                g.savefig(os.path.join(output_dir, prefix + "." + cell_type +
                                    ".mohgp.fitted_model.mean_acc.clustermap.cluster_labels.total_norm{}.{}only_posterior_above_threshold.svg".format(label, label2)), dpi=300, bbox_inches="tight")

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
            all_g_coding = all_g[(~all_g.str.startswith("LINC")) & (~all_g.str.startswith("LOC")) & (~all_g.str.startswith("MIR"))].str.replace("-AS1", "")
            top_g_coding = top_g[(~top_g.str.startswith("LINC")) & (~top_g.str.startswith("LOC")) & (~top_g.str.startswith("MIR"))].str.replace("-AS1", "")

            all_g.to_csv(os.path.join("results", "single_cell_RNA", "ATAC-seq_signatures", "atac-seq.{}.cluster_{}.gene_level.csv".format(cell_type, cluster)), index=False)
            top_g.to_csv(os.path.join("results", "single_cell_RNA", "ATAC-seq_signatures", "atac-seq.{}.cluster_{}.gene_level.top200.csv".format(cell_type, cluster)), index=False)
            all_g_coding.to_csv(os.path.join("results", "single_cell_RNA", "ATAC-seq_signatures", "atac-seq.{}.cluster_{}.gene_level.only_coding.csv".format(cell_type, cluster)), index=False)
            top_g_coding.to_csv(os.path.join("results", "single_cell_RNA", "ATAC-seq_signatures", "atac-seq.{}.cluster_{}.gene_leveltop200.only_coding.csv".format(cell_type, cluster)), index=False)


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
        analysis.accessibility,
        cell_types=cell_types, 
        matrix_file=matrix_file,
        linear_time=True,
        prefix=prefix)  # wait for jobs to complete
    fits = gather_gaussian_processes(
        analysis.accessibility,
        matrix_file=matrix_file,
        prefix=prefix)
    fits.to_csv(os.path.join(gp_output_dir, prefix + ".GP_fits.all_cell_types.csv"), index=True)
    fits = pd.read_csv(os.path.join(gp_output_dir, prefix + ".GP_fits.all_cell_types.csv"), index_col=0)

    visualize_gaussian_process_fits(analysis, fits, output_dir=gp_output_dir, prefix=prefix)

    # cluster variable regulatory elements with hierarchical mixtures of GPs (MOHGP)
    
    for alpha in [0.01, 0.05]:
        prefix = prefix + ".p={}".format(alpha)
        fit_MOHGP(
            analysis.accessibility,
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
            analysis.accessibility, fits,
            cell_types=cell_types,
            n_clust=[3, 4, 4, 4, 4, 3],
            prefix=prefix, fits_dir=mohgp_output_dir, alpha=alpha, posterior_threshold=0.8)
        assignments.to_csv(os.path.join(mohgp_output_dir, prefix + ".GP_fits.linear_time.mohgp_clusters.csv"), index=True)

        visualize_clustering_assignments(
            analysis.accessibility, assignments, prefix=prefix,
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
            analysis.accessibility,
            cell_types=['CLL'],
            matrix_file=os.path.abspath(os.path.join(
                "results", analysis.name + ".accessibility.annotated_metadata.csv")),
            prefix=prefix + ".random_{}".format(i),
            randomize=True)  # wait for jobs to complete
        fits = gather_gaussian_processes(
            analysis.accessibility,
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
    example_acc = analysis.accessibility.loc[examples, :]
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
        analysis, output_dir="{results_dir}/tf_accessibility",
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
        transcription_factor_a = analysis.accessibility.loc[transcription_factor_r.index].dropna()

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
    mean = analysis.accessibility.mean(axis=1)

    res = pd.DataFrame()
    for l_quantile, u_quantile in zip(r, r[1:]):
        i = mean[(mean.quantile(l_quantile) > mean) & (mean < mean.quantile(u_quantile))].index

        m = analysis.accessibility.loc[i, :].mean()
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
    d["binding_sites"] = analysis.accessibility.shape[0]
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
            all_res.loc[(all_res['cell_type'] == cell_type) & (all_res['transcription_factor'] == tf), 'norm_accessibility'] = ((s - b) / b.std()).values

    all_res.to_csv(os.path.join(output_dir, "all_factor_binding.normalized.csv"), index=False)
    all_res = pd.read_csv(os.path.join(output_dir, "all_factor_binding.normalized.csv"))

    # Plot
    g = sns.factorplot(data=all_res[all_res['transcription_factor'] == "background"], x="timepoint", y="accessibility", hue="cell_type")
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.cell_type_hue.lineplots.background_only.svg"), dpi=300, bbox_inches="tight")

    g = sns.factorplot(data=all_res[all_res['transcription_factor'] != "background"], x="timepoint", y="norm_accessibility", hue="cell_type", col="transcription_factor", col_wrap=5, sharey=False)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.cell_type_hue.lineplots.free_y.svg"), dpi=300, bbox_inches="tight")

    g = sns.factorplot(data=all_res[all_res['transcription_factor'] != "background"], x="timepoint", y="norm_accessibility", hue="cell_type", col="transcription_factor", col_wrap=5)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_patients.cell_type_hue.lineplots.svg"), dpi=300, bbox_inches="tight")


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
        p_z.loc[:, p_z.columns.get_level_values("cell_type") == cell_type] = scipy.stats.zscore(
            p_z.loc[:, p_z.columns.get_level_values("cell_type") == cell_type], axis=1)

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



def correlation_to_bcell(analysis):
    cll = analysis.accessibility.loc[:, (analysis.accessibility.columns.get_level_values("cell_type") == "CLL")]
    bcell = analysis.accessibility.loc[:, (analysis.accessibility.columns.get_level_values("cell_type") == "Bcell")]

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
    accessibility = analysis.accessibility
    accessibility.columns = accessibility.columns.get_level_values("sample_name")

    for sample_set in ["fzhao", "dlara"]:
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
        tmp_analysis.coverage = tmp_analysis.coverage.loc[:, ~tmp_analysis.coverage.columns.str.contains("D199")]
        # Normalize cell types jointly (quantile normalization)
        tmp_analysis.normalize_coverage_quantiles(samples=[s for s in tmp_analysis.samples if "D199" not in s.name])
        tmp_analysis.to_pickle()

        # Join matrix with CLL samples from major tmp_analysis
        q = accessibility.join(tmp_analysis.coverage_qnorm).drop(['chrom', 'start', 'end'], axis=1)
        # Compute correlation

        c = q.corr()

        g = sns.clustermap(analysis.coverage_qnorm.drop(['chrom', 'start', 'end'], axis=1).corr(), xticklabels=False, figsize=(8 ,8), rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.savefig(os.path.join("results", "differentiation_assessment.{}_only.map.svg".format(sample_set)), dpi=200)

        g = sns.clustermap(c, xticklabels=False, figsize=(8 ,8), rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.savefig(os.path.join("results", "differentiation_assessment.{}_together.map.svg".format(sample_set)), dpi=200)

        mask = c.index.str.contains("d_CLL") | (c.index.str.contains("CB") & c.index.str.contains("NaiveBcell|chedBcell|Tcell"))
        g = sns.clustermap(c.loc[mask, mask], xticklabels=False, figsize=(8 ,8), rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.savefig(os.path.join("results", "differentiation_assessment.{}_together.relevant.map.svg".format(sample_set)), dpi=200)

        c2 = pd.melt(c.loc[mask, mask].reset_index(), "index")
        annot = pd.DataFrame(map(pd.Series, c2['index'].str.split("_"))).iloc[:, [1, 2, 3]]
        annot[2] = annot[2].str.replace("d|ND", "").astype(int)
        c2.index = pd.MultiIndex.from_arrays(annot.T.values, names=["patient_id", "timepoint", "cell_type"])

        c3 = c2[c2['variable'].str.contains("ND") & ~c2.index.get_level_values("cell_type").str.contains("ND") & ~c2['index'].str.contains("ND")]
        c3['healty_subset'] = [x[3] for x in c3['variable'].str.split("_")]
        g = sns.factorplot(data=c3.reset_index(), x="timepoint", y="value", col="patient_id", hue="healty_subset", sharey=False)
        g.savefig(os.path.join("results", "differentiation_assessment.{}_correlation_to_normal.time_dependent.svg".format(sample_set)), dpi=300)

        g = sns.factorplot(data=c3.reset_index(), x="timepoint", y="value", col="patient_id", hue="healty_subset", sharey=False, sharex=False)
        g.savefig(os.path.join("results", "differentiation_assessment.{}_correlation_to_normal.time_dependent.svg".format(sample_set)), dpi=300)

        data = pd.pivot_table(data=c3.reset_index(), index=['patient_id', 'timepoint', 'cell_type'], columns='healty_subset', values="value")
        data_z = pd.DataFrame(zscore(data, axis=0), index=data.index, columns=data.columns)
        data_z = (data - data.mean(0)) / data.std(0)

        g = sns.factorplot(
            data=pd.melt(data_z.reset_index(), id_vars=['patient_id', 'timepoint', 'cell_type'], var_name="healty_subset"),
            x="timepoint", y="value", col="patient_id", hue="healty_subset", sharey=False, sharex=False)
        g.savefig(os.path.join("results", "differentiation_assessment.{}_correlation_to_normal.time_dependent.zscore.svg".format(sample_set)), dpi=300)


        from ngs_toolkit.graphics import radar_plot

        f = radar_plot(
            data[data['patient_id'] != "KI"].reset_index(),
            subplot_var="patient_id", group_var="timepoint",
            radial_vars=["NaiveBcell", "SwitchedBcell", "UnswitchedBcell"],
            cmap="inferno", scale_to_max=False)
        f.savefig(os.path.join("results", "differentiation_assessment.{}_correlation_to_normal.time_dependent.radar_plot.svg".format(sample_set)), dpi=300)

        f = radar_plot(
            data[data['patient_id'] != "KI"].reset_index(),
            subplot_var="patient_id", group_var="timepoint",
            radial_vars=["NaiveBcell", "SwitchedBcell", "UnswitchedBcell"],
            cmap="inferno", scale_to_max=True)
        f.savefig(os.path.join("results", "differentiation_assessment.{}_correlation_to_normal.time_dependent.radar_plot.scale_to_max.svg".format(sample_set)), dpi=300)


def get_gene_level_accessibility(analysis):
    """
    Get gene-level measurements of chromatin accessibility.
    """
    assert hasattr(analysis, "gene_annotation")
    acc = analysis.accessibility.copy()

    g = analysis.gene_annotation['gene_name'].str.split(",").apply(pd.Series).stack()
    g.index = g.index.droplevel(1)
    g.name = "gene_name"
    acc2 = analysis.accessibility.join(g).drop("gene_name", axis=1)
    acc2.index = analysis.accessibility.join(g).reset_index().set_index(['index', 'gene_name']).index
    acc2.columns = analysis.accessibility.columns
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
    acc = analysis.accessibility.copy()
    acc.index = analysis.accessibility.join(analysis.gene_annotation).reset_index().set_index(['index', 'gene_name']).index

    # Get all reg. elements which are associated to each gene
    full_acc = acc.loc[acc.index.get_level_values("gene_name").isin(lr_genes), :]
    red_acc = full_acc.groupby(level="gene_name").mean()

    full_acc_time = full_acc.T.groupby(level=["cell_type", "timepoint"]).mean().T
    red_acc_time = red_acc.T.groupby(level=["cell_type", "timepoint"]).mean().T

    for variable, label1 in [("full", "region_level"), ("red", "gene_level")]:
        for ext, label2 in [("", ""), ("_time", ".timepoint_mean")]:
            for z, label3 in [(None, ""), (1, "_zscore")]:
                m = eval(variable + "_" + "acc" + ext)
                w, h = m.shape[0] * 0.01, m.shape[1] * 0.12
                g = sns.clustermap(
                    m.T, col_cluster=True, row_cluster=True, z_score=z,
                    figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
                g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.sig_only{}.clustermap.svg".format(label1, label2, label3)), dpi=300, bbox_inches="tight")

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
                m = eval(variable + "_" + "acc" + ext + "_sig")

                g = sns.clustermap(
                    m.T, col_cluster=True, row_cluster=True, z_score=z,
                    robust=True, xticklabels=False, rasterized=True)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
                g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.sig_only{}.clustermap{}.svg".format(label1, label2, label3)), dpi=300, bbox_inches="tight")

    # for each cell type
    for cell_type in acc.columns.levels[3]:
        variable_regions = analysis.assignments[(analysis.assignments['cell_type'] == cell_type)].index.drop_duplicates().tolist()
        for variable, label1 in [("full", "region_level"), ("red", "gene_level")]:
            for ext, label2 in [("", ""), ("_time", ".timepoint_mean")]:
                for row_cluster, label3 in [(False, "sorted"), (True, "clustered")]:
                    for z, label4 in [(None, ""), (1, ".zscore")]:
                        m = eval(variable + "_" + "acc" + ext + "_sig").sort_index(axis=1, level=['timepoint'])

                        if variable == "full":
                            m = m.loc[variable_regions, m.columns.get_level_values("cell_type") == cell_type]
                        else:
                            variable_genes = full_acc.index[full_acc.index.get_level_values("index").isin(variable_regions)].get_level_values("gene_name").drop_duplicates().tolist()
                            m = m.loc[variable_genes, m.columns.get_level_values("cell_type") == cell_type]

                        w, h = m.shape[0] * 0.1, m.shape[1] * 1.2

                        if variable == "red":
                            square = True
                            xticklabels = True
                        else:
                            square = False
                            xticklabels = False

                        g = sns.clustermap(
                            m.T, col_cluster=True, row_cluster=row_cluster, z_score=z, square=square,
                            figsize=(h, w), robust=True, xticklabels=xticklabels, rasterized=True)
                        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize=h / 3.)
                        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="center", fontsize=w / 3.)
                        g.savefig(os.path.join(output_dir, "ligand-receptor_repertoire.{}.{}.sig_only{}.clustermap.{}{}.svg".format(cell_type, label1, label2, label3, label4)), dpi=300, bbox_inches="tight")

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
    acc = analysis.accessibility.copy()
    g = analysis.gene_annotation['gene_name'].str.split(",").apply(pd.Series).stack()
    g.index = g.index.droplevel(1)
    g.name = "gene_name"
    variable_genes = list(set(g.loc[g.index.isin(variable)]))

    # Get gene-level accessibility averaged per cell type per timepoint
    acc = analysis.accessibility.loc[variable]
    acc.index = analysis.accessibility.join(analysis.gene_annotation).reset_index().set_index(['index', 'gene_name']).index

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

    mean_signatures = pd.DataFrame()
    for i, cell_type in enumerate(assignments['cell_type'].drop_duplicates()):
        for cluster in range(5):
            ass = assignments.loc[
                (assignments['cell_type'] == cell_type) &
                (assignments['cluster_assignment'] == cluster), :].sort_values("p_value")
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

    # Plot signature change with time within cell types for each ATAC-seq cluster
    piv = pd.pivot_table(data=mean_signatures, index=['cell_type', 'timepoint'], columns=['ATAC-seq_cell_type', 'cluster'])

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
    df2 = df.loc[:,
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

            if p.empty or p.shape[1] < 2:
                continue

            g = sns.clustermap(p, xticklabels=True, metric="correlation", rasterized=True, row_cluster=False, figsize=(p.shape[1] * 0.05, p.shape[0] * 0.12))
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=3)
            g.savefig(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.{}.cluster_{}.diff_rna.clustermap.svg".format(cell_type, cluster)), bbox_inches="tight", dpi=300)

            g = sns.clustermap(p, z_score=1, xticklabels=True, metric="correlation", rasterized=True, row_cluster=False, figsize=(p.shape[1] * 0.05, p.shape[0] * 0.12), robust=True)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=3)
            g.savefig(os.path.join(output_dir, "scRNA-seq_signal_on_ATAC-seq_clusters.{}.cluster_{}.diff_rna.clustermap.zscore.svg".format(cell_type, cluster)), bbox_inches="tight", dpi=300)



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
            fig, axis = plt.subplots(1, 2, figsize=(4, 3))
            fig.savefig(os.path.join(output_dir, "lola_vs_encode_enrichments.CLL_down.scatter.svg"), dpi=300, bbox_inches="tight")
            
            # barplots
            l = lola_red.loc[('CLL_down', 'all')].sort_values('log_p_value', ascending=False).head(20).reset_index()
            e = encode_red.loc[('CLL', 'down')].sort_values('log_p_value', ascending=False).head(20).reset_index()
            fig, axis = plt.subplots(1, 2, figsize=(4, 3))
            sns.barplot(data=l, x="log_p_value", y="antibody", orient="horizontal", color=sns.color_palette("colorblind")[0], ax=axis[0])
            sns.barplot(data=e, x="log_p_value", y="tf", orient="horizontal", color=sns.color_palette("colorblind")[0], ax=axis[1])
            sns.despine(fig)
            fig.savefig(os.path.join(output_dir, "lola_vs_encode_enrichments.CLL_down.barplot.svg"), dpi=300, bbox_inches="tight")

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
        expr = scrna_diff2.loc[~scrna_diff2['gene_name'].str.contains("RPL|RP-|RPS|MT-|HLA"), :]


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






if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
