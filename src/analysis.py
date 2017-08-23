
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
from looper.models import Project, Sample
from scipy.cluster.hierarchy import fcluster
from scipy.stats import pearsonr, spearmanr, zscore
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder
from statsmodels.stats.multitest import multipletests

from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.general import (collect_differential_enrichment,
                                 differential_enrichment,
                                 normalize_quantiles_r,
                                 plot_differential_enrichment,
                                 subtract_principal_component)

sns.set_style("white")
plt.rcParams['svg.fonttype'] = 'none'

# random seed
SEED = int("".join(
    LabelEncoder()
    .fit(list(string.ascii_uppercase))
    .transform(list("BOCKLAB")).astype(str)))
random.seed(SEED)
np.random.seed(SEED)



def subtract_principal_component_by_attribute(df, pc=1, attributes=["CLL"]):
    """
    Given a matrix (n_samples, n_variables), remove `pc` (1-based) from matrix.
    """
    from sklearn.decomposition import PCA

    pc -= 1

    X2 = pd.DataFrame(index=df.index, columns=df.columns)
    for attr in attributes:
        print(attr)
        sel = df.index[df.index.str.contains(attr)]
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


def fit_gaussian_processes(matrix, matrix_name="sorted"):
    """
    Estimate temporal variability of regulatory elements by comparing
    the fit of a Gaussian Process (GP) regression model with a
    variable kernel with and another with a static kernel.
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

    # Fit GPs in parallel jobs per cell type and in chunks
    chunks = 2000
    total_job_lim = 800
    refresh_time = 10
    in_between_time = 0.01
    output_prefix = "gp_fit_job"
    output_dir = "/scratch/users/arendeiro/gp_fit_job"
    library = "GPy"

    r = np.arange(0, matrix.shape[0], chunks)

    for cell_type in tqdm.tqdm(matrix.columns.get_level_values("cell_type").drop_duplicates(), desc="cell_type"):
        for start, end in tqdm.tqdm(zip(r, r[1:]) + [(r[-1], matrix.shape[0])], desc="chunk"):
            range_name = "{}-{}".format(start, end)
            name = ".".join([output_prefix, matrix_name, cell_type, range_name, library])
            log = os.path.join(output_dir, "log", name + ".log")
            job = """python ~/jobs/gp_fit_job.py --data-range {} --range-delim - --matrix-type {} --cell-type {} --library {} --output-prefix {} --output-dir {}""".format(
                range_name, matrix_name, cell_type, library, name, os.path.join(output_dir, "output"))
            cmd = """sbatch -J {} -o {} -p shortq -c 1 --mem 8000 --wrap "{}" """.format(
                name, log, job)

            if not os.path.exists(os.path.join(output_dir, "output", name + ".csv")):
                submit_job_if_possible(cmd, total_job_lim=total_job_lim, refresh_time=refresh_time, in_between_time=in_between_time)


def gather_gaussian_processes(matrix, matrix_name="sorted"):
    """
    Estimate temporal variability of regulatory elements by comparing
    the fit of a Gaussian Process (GP) regression model with a
    variable kernel with and another with a static kernel.
    """

    # Collect output of parallel jobs
    fits = pd.DataFrame()
    r = np.arange(0, matrix.shape[0], chunks)
    for cell_type in tqdm.tqdm(matrix.columns.get_level_values("cell_type").drop_duplicates(), desc="cell_type"):
        for start, end in tqdm.tqdm(zip(r, r[1:]) + [(r[-1], matrix.shape[0])], desc="chunk"):
            range_name = "{}-{}".format(start, end)
            name = ".".join([output_prefix, matrix_name, cell_type, range_name, library])
            df = pd.read_csv(os.path.join(output_dir, "output", name + ".csv"), index_col=0)
            df['cell_type'] = cell_type

            fits = fits.append(df)

    # correct p-values
    fits['q_value'] = np.concatenate(fits.groupby("cell_type", sort=False)['p_value'].apply(lambda x: multipletests(x, method="fdr_bh")[1]))

    # save
    fits.to_csv(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "all_fits.csv"])), index=True)
    # fits = pd.read_csv(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "all_fits.csv"])), index_col=0)

    # Visualize the relationship between the fits and parameters
    g = sns.PairGrid(fits.drop("cell_type", axis=1).sample(n=2000))
    g.map(plt.scatter, alpha=0.2, s=2, rasterized=True)
    g.savefig(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "fits.parameters.pairwise.all_cell_types.svg"])), dpi=300, bbox_inches="tight")

    # Plot likelihood relationships
    g = sns.FacetGrid(data=fits, col="cell_type", col_wrap=2)
    g.map(plt.scatter, "RBF", "White", alpha=0.1, s=2, rasterized=True)
    g.savefig(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "fits.RBF_vs_White.cell_types.svg"])), dpi=300, bbox_inches="tight")

    g = sns.FacetGrid(data=fits, col="cell_type", col_wrap=2)
    g.map(plt.scatter, "White", "D", alpha=0.1, s=2, rasterized=True)
    g.savefig(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "fits.D_vs_White.cell_types.svg"])), dpi=300, bbox_inches="tight")

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
    fig.savefig(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "fits.D_vs_White.mean_posterior_std.cell_types.svg"])), dpi=300, bbox_inches="tight")


    # Let's rank regions
    n_top = 6
    e = fits[~fits['cell_type'].str.contains("NK")].sort_values("p_value").head(n_top)
    examples = e.index
    example_ct = e['cell_type']
    example_acc = matrix.loc[examples]
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
    fig.savefig(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "top_variable.scatter.all_samples.svg"])), dpi=300, bbox_inches="tight")

    return fits

def fit_MOHGP(matrix, matrix_name="sorted"):
    """
    Cluster temporaly variable regulatory elements with
    a Hierarchical Mixture of Gaussian Processes (MOHGP).
    """

    # Fit MOHGPs in parallel jobs per cell type
    output_prefix = "gp_fit_job"
    output_dir = "/scratch/users/arendeiro/gp_fit_job"
    library = "GPy"

    for cell_type in ["Bcell", "Bulk", "CLL", "CD4", "CD8", "Mono", "NK"]: 
    # for cell_type in matrix.columns.get_level_values("cell_type").drop_duplicates():
        name = ".".join([output_prefix, matrix_name, cell_type, "GPclust"])
        log = os.path.join("results_deconvolve", "log", name + ".log")
        job = """python -u ~/jobs/gpclust_job.py --cell-type {} --output-prefix {} --output-dir {}""".format(
            cell_type, name, "/home/arendeiro/cll-time_course/results_deconvolve")
        cmd = """sbatch -J {} -o {} -p shortq -c 12 --mem 80000 --wrap "{}" """.format(
            name, log, job)
        os.system(cmd)


def gather_MOHGP(matrix, matrix_name="sorted", posterior_threshold=0.8):
    """
    Cluster temporaly variable regulatory elements with
    a Hierarchical Mixture of Gaussian Processes (MOHGP).
    """
    fits = pd.read_csv(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "all_fits.csv"])), index_col=0)

    assignments = pd.DataFrame()
    for cell_type in matrix.columns.get_level_values("cell_type").drop_duplicates():
        print(cell_type)

        # Get variable regions for cell type
        variable = fits[(fits['p_value'] < alpha) & (fits['cell_type'] == cell_type)].index

        # Read in their posterior probabilities matrix (Phi) of cluster assignments
        name = ".".join([output_prefix, matrix_name, cell_type, "GPclust"])

        phi = pd.DataFrame(
            np.fromfile(os.path.join("results_deconvolve", name + "." + cell_type + ".MOHCP.posterior_probs_phi.np")).reshape((len(variable), 4)),
            index=variable)

        # Threshold probabilities to filter out some regions
        assigned = (phi > posterior_threshold).any(axis=1)

        # Append
        cluster_labels = phi.loc[assigned].apply(np.argmax, axis=1).to_frame(name="cluster")
        cluster_labels["cell_type"] = cell_type
        assignments = assignments.append(cluster_labels)

        # Plot
        matrix2 = matrix.loc[variable, matrix.columns.get_level_values("cell_type") == cell_type]
        tp = pd.Series(matrix2.columns.get_level_values("timepoint").str.replace("d", "").astype(int), index=matrix2.columns).sort_values()

        # all variable regions with assignments and threshold mask
        col_colors = [[plt.get_cmap("Paired")(i) for i in phi.apply(np.argmax, axis=1)], [plt.get_cmap("binary")(i) for i in assigned.astype(float)]]
        g = sns.clustermap(
            matrix2.loc[:, tp.index].T,
            col_colors=col_colors,
            row_cluster=False, col_cluster=True, z_score=1, xticklabels=False, yticklabels=matrix2.loc[:, tp.index].columns.get_level_values("sample_name"),
            rasterized=True, figsize=(8, 0.2 * matrix2.shape[1]), metric="correlation", robust=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.ax_col_dendrogram.set_rasterized(True)
        g.savefig(os.path.join("results_deconvolve", output_prefix + "." + cell_type + ".MOHCP.fitted_model.clustermap.cluster_labels.with_posterior_probs.svg"), dpi=300, bbox_inches="tight")

        # only variable and with assignments above threshold
        g = sns.clustermap(
            matrix2.loc[assigned, tp.index].T,
            col_colors=[plt.get_cmap("Paired")(i) for i in phi.loc[assigned].apply(np.argmax, axis=1)],
            row_cluster=False, col_cluster=True, z_score=1, xticklabels=False, yticklabels=matrix2.loc[:, tp.index].columns.get_level_values("sample_name"),
            rasterized=True, figsize=(8, 0.2 * matrix2.shape[1]), metric="correlation", robust=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.ax_col_dendrogram.set_rasterized(True)
        g.savefig(os.path.join("results_deconvolve", output_prefix + "." + cell_type + ".MOHCP.fitted_model.clustermap.cluster_labels.only_posterior_above_threshold.svg"), dpi=300, bbox_inches="tight")

        # only variable and with assignments above threshold: mean per timepoint
        matrix_mean = matrix2.loc[assigned, tp.index].T.groupby(level="timepoint").mean()
        g = sns.clustermap(
            matrix_mean,
            col_colors=[plt.get_cmap("Paired")(i) for i in phi.loc[assigned].apply(np.argmax, axis=1)],
            row_cluster=False, col_cluster=True, z_score=1, xticklabels=False, yticklabels=True,
            rasterized=True, figsize=(8, 0.2 * matrix_mean.shape[0]), metric="correlation", robust=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.ax_col_dendrogram.set_rasterized(True)
        g.savefig(os.path.join("results_deconvolve", output_prefix + "." + cell_type + ".MOHCP.fitted_model.mean_acc.clustermap.cluster_labels.only_posterior_above_threshold.svg"), dpi=300, bbox_inches="tight")

        # Cluster patterns
        # this is a mock of the MOHGP underlying posterior.
        cs = cluster_labels['cluster'].drop_duplicates().dropna().shape[0]

        fig, axis = plt.subplots(1, cs, figsize=(cs * 4, 1 * 4), sharex=True, sharey=True)
        for i, cluster in enumerate(cluster_labels['cluster'].drop_duplicates().sort_values()):
            regions = cluster_labels[cluster_labels['cluster'] == cluster].index

            X = matrix2.loc[regions, tp.index].mean(axis=0).T.reset_index()
            X['time'] = np.log2(1 + X["timepoint"].str.replace("d", "").astype(int).values)
            # d = X.groupby('time').mean().squeeze().to_frame(name="mean")
            # d['upper_q'] = X.groupby('time').quantile(.975)
            # d['lower_q'] = X.groupby('time').quantile(.025)

            kernel = GPy.kern.RBF(1.0, variance=0.5) + GPy.kern.Bias(1.0, variance=0.05)
            m = GPy.models.GPRegression(X=X[['time']], Y=X[[0]], kernel=kernel)
            m.optimize()
            m.plot([0 - 0.5, max(X['time']) + 0.5], ax=axis[i], legend=None)

            # axis[i].set_ylim((1.5, 3.5))
            axis[i].set_title("Cluster {}\n(n = {})".format(cluster, regions.shape[0]))
        axis[0].set_ylabel("Chromatin accessibility")
        for ax in axis:
            ax.set_xlabel("Time (log2)")
        fig.savefig(os.path.join("results_deconvolve", output_prefix + "." + cell_type + ".MOHCP.fitted_model.mean_acc.clustermap.cluster_labels.only_posterior_above_threshold.variable.cluster_means.svg"), dpi=300, bbox_inches="tight")

    return assignments


def cluster_dynamics(assignments):
    cluster_labels = pd.DataFrame({
        "Bcell": {0: "down", 1: "up", 2: "other", 3: "other"},
        # "Bulk": {0: "down", 1: "up2", 2: "up1", 3: "other"},
        "Bulk": {0: "down", 1: "up", 2: "up", 3: "other"},
        "CD4": {0: "down", 1: "up", 2: "middle", 3: "extremes"},
        "CD8": {0: "down", 1: "up", 2: "middle", 3: "extremes"},
        "CLL": {0: "down", 1: "up", 2: "extremes", 3: "middle"},
        # "NK": {0: "down", 1: "up", 2: "extremes2", 3: "extremes1"},
        "NK": {0: "down", 1: "up", 2: "extremes", 3: "extremes"},
        "Mono": {0: "extremes", 1: "up", 2: "down", 3: "other"}
    })
    cluster_labels.index.name = "cluster"

    cluster_counts = assignments.groupby(['cell_type', 'cluster'])['index'].count().to_frame(name="count")

    counts = pd.merge(cluster_counts.reset_index(), pd.melt(cluster_labels.reset_index(), id_vars="cluster", var_name="cell_type", value_name='cluster_name'))

    counts2 = counts.groupby(['cell_type', 'cluster_name'])['count'].sum().to_frame(name="count")

    counts_fraction = (counts2 / counts2.groupby(level='cell_type').sum()) * 100.

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
    fig.savefig(os.path.join("results_deconvolve", output_prefix + "." + cell_type + ".MOHCP.cluster_name.counts.barplot.svg"), bbox_inches="tight")
    

# Start project and analysis objects
prj = Project("metadata/project_config.yaml")
prj.samples = [sample for sample in prj.samples if sample.library == "ATAC-seq"]
analysis = ATACSeqAnalysis(name="cll-time_course", prj=prj, samples=prj.samples)

# Project's attributes
sample_attributes = ["sample_name", "patient_id", "timepoint", "cell_type", "compartment", "response", "cell_number", "batch"]
numerical_attributes = ["CD38_cells_percentage", "cll_cells_%", "cell_number"]
cell_types = ["Bulk", "Bcell", "CLL", "CD4", "CD8", "NK", "Mono"]


# GET CONSENSUS PEAK SET, ANNOTATE IT, PLOT FEATURES
# Get consensus peak set from all samples
analysis.get_consensus_sites(analysis.samples)

# GET CHROMATIN OPENNESS MEASUREMENTS
# Get coverage values for each peak in each sample of ATAC-seq
analysis.measure_coverage(analysis.samples)
analysis.coverage = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.raw_coverage.csv"), index_col=0)

analysis.samples = [s for s in analysis.samples if s.final_pass != "0" and s.include != "0"]
analysis.coverage = analysis.coverage[[s.name for s in analysis.samples] + ['chrom', 'start', 'end']]
analysis.coverage = analysis.coverage[~analysis.coverage.index.str.contains("chrX|chrY")]

# Normalize cell types jointly (quantile normalization)
analysis.normalize_coverage_quantiles()

# Fix batches for cell types independently
to_norm = analysis.coverage_qnorm.drop(["chrom", "start", "end"], axis=1)
analysis.coverage_qnorm_batchfix = subtract_principal_component_by_attribute(to_norm.T, pc=1, attributes=cell_types[:-1]).T

# Transform values (cube root to handle simetrically negative, zero and positive values - pseudo-log scale)
sign = (analysis.coverage_qnorm_batchfix >= 0).astype(int).replace(0, -1)
analysis.accessibility = sign * np.absolute(analysis.coverage_qnorm_batchfix) ** (1 / 3.)
analysis.accessibility = analysis.annotate_with_sample_metadata(quant_matrix="accessibility", attributes=sample_attributes, save=False, assign=False)

analysis.accessibility.to_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.coverage.joint_qnorm.pca_fix.power.csv"), index=True)
analysis.accessibility = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.coverage.joint_qnorm.pca_fix.power.csv"), index_col=0, header=range(8))


# Vizualize all transformations
analysis.coverage_qnorm = analysis.annotate_with_sample_metadata(quant_matrix="coverage_qnorm", attributes=sample_attributes, save=False, assign=False)
analysis.unsupervised(quant_matrix="coverage_qnorm", samples=None, attributes_to_plot=sample_attributes, plot_prefix="coverage_qnorm")

analysis.coverage_qnorm_batchfix = analysis.annotate_with_sample_metadata(quant_matrix="coverage_qnorm_batchfix", attributes=sample_attributes, save=False, assign=False)
analysis.unsupervised(quant_matrix="coverage_qnorm_batchfix", samples=None, attributes_to_plot=sample_attributes, plot_prefix="coverage_qnorm_batchfix")


analysis.unsupervised(quant_matrix="accessibility", samples=None, attributes_to_plot=sample_attributes, plot_prefix="accessibility")


# Annotate peaks
analysis.get_peak_gene_annotation()
analysis.get_peak_genomic_location()
analysis.get_peak_chromatin_state()
analysis.calculate_peak_support(analysis.samples, "summits")

analysis.gene_annotation = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.gene_annotation.csv"))
analysis.region_annotation = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.region_annotation.csv"))
analysis.chrom_state_annotation = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.chromatin_state.csv"))
analysis.support = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.support.csv"))

analysis.to_annot = analysis.accessibility.copy()
analysis.to_annot.columns = analysis.to_annot.columns.get_level_values("sample_name")
analysis.to_annot = analysis.to_annot.join(analysis.coverage[["chrom", "start", "end"]])
analysis.annotate(quant_matrix="to_annot")


# # to annotate with genes
# analysis.gene_annotation.index = (
#     analysis.gene_annotation["chrom"] + ":" +
#     analysis.gene_annotation["start"].astype(str) + "-" +
#     analysis.gene_annotation["end"].astype(str))


# TIME SERIES ANALYSIS

# filter some samples out
samples_to_exclude = [
    "PBGY_1d_NK",
    "240d_CD8", "280d_CD8",
    "240d_CD4", "280d_CD4",
    "PBGY_1d_Mono", "PBGY_8d_Mono", "PBGY_150d_Mono",
    "KZ_240d_Bulk",
    "FE_3d_Bulk",
    "VZS_2d_Bulk",
]
matrix = analysis.accessibility.loc[:, ~analysis.accessibility.columns.get_level_values("sample_name").str.contains("|".join(samples_to_exclude))]

# Fit GPs with varying and constant kernel to detect variable regulatory elements
fit_gaussian_processes(matrix, matrix_name="sorted")
fits = gather_gaussian_processes(matrix, matrix_name="sorted")

# Cluster variable regulatory elements with hierarchical mixtures of GPs (MOHGP)
fit_MOHGP(analysis.accessibility, matrix_name="sorted")
assignments = gather_MOHGP(matrix, matrix_name="sorted", posterior_threshold=0.8)

assignments = pd.merge(assignments.reset_index(), fits.reset_index(), on=['index', 'cell_type'], how='left')
assignments.to_csv(os.path.join("results_deconvolve", output_prefix + ".GP_fits.MOHCP_clusters.csv"), index=False)


# Plot distribution of clusters per cell type dependent on their dynamic pattern
cluster_dynamics(assignments)


# Get enrichments of region clusters
assignments['comparison_name'] = assignments['cell_type'] + "_" + assignments['cluster'].astype(str)


differential_enrichment(
    analysis,
    assignments.set_index("index"),
    data_type="ATAC-seq",
    output_dir="results_deconvolve/GPclust",
    output_prefix=".".join([output_prefix, matrix_name, "GPclust"]),
    genome="hg19",
    directional=False,
    max_diff=10000,
    sort_var="D",
    as_jobs=True
)

diff = collect_differential_enrichment(
    assignments,
    directional=False,
    data_type="ATAC-seq",
    output_dir="results_deconvolve/GPclust",
    output_prefix=".".join([output_prefix, matrix_name, "GPclust"]),
    permissive=True)
enrichment_table = pd.read_csv(os.path.join("results_deconvolve", "GPclust", "gp_fit_job.sorted.GPclust.meme_ame.csv"))
plot_differential_enrichment(
    enrichment_table,
    "motif",
    data_type="ATAC-seq",
    direction_dependent=False,
    output_dir="results_deconvolve/GPclust",
    comp_variable="comparison_name",
    output_prefix=".".join([output_prefix, matrix_name, "GPclust"]),
    top_n=5)
enrichment_table = pd.read_csv(os.path.join("results_deconvolve", "GPclust", "gp_fit_job.sorted.GPclust.lola.csv"))
plot_differential_enrichment(
    enrichment_table,
    "lola",
    data_type="ATAC-seq",
    direction_dependent=False,
    output_dir="results_deconvolve/GPclust",
    comp_variable="comparison_name",
    output_prefix=".".join([output_prefix, matrix_name, "GPclust"]),
    top_n=5)
enrichment_table = pd.read_csv(os.path.join("results_deconvolve", "GPclust", "gp_fit_job.sorted.GPclust.enrichr.csv"))
plot_differential_enrichment(
    enrichment_table,
    "enrichr",
    data_type="ATAC-seq",
    direction_dependent=False,
    output_dir="results_deconvolve/GPclust",
    comp_variable="comparison_name",
    output_prefix=".".join([output_prefix, matrix_name, "GPclust"]),
    top_n=5)






# INTEGRATION with scRNA-seq
# Comparison with scRNA-seq diff genes

# read in 
diff_genes = pd.read_csv(os.path.join("metadata", "scRNA_diff_genes_over_time.tsv"), sep="\t")
diff_genes = diff_genes.loc[diff_genes['qvalue'] < 0.05, :]
all_diff_genes = diff_genes.loc[diff_genes.index, "gene"].drop_duplicates()

diff_peaks = analysis.coverage_annotated[analysis.coverage_annotated['gene_name'].isin(all_diff_genes)].drop_duplicates().index

g = sns.clustermap(
    analysis.accessibility.loc[diff_peaks, :],
    z_score=0, metric="correlation", robust=True,
    yticklabels=False, xticklabels=analysis.accessibility.columns.get_level_values("sample_name"), rasterized=True,
    figsize=(18, 6)
)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
g.ax_row_dendrogram.set_rasterized(True)
g.savefig(os.path.join("results", "scRNA_diff_genes_over_time.clustermap.svg"), dpi=300, bbox_inches="tight")


for cell_type in diff_genes['cellType'].drop_duplicates():
    print(cell_type)
    if cell_type in ["NurseLikeCells", "Tcells1"]:
        mask = analysis.accessibility.columns.get_level_values("cell_type").str.contains("CD4|CD8")
    elif cell_type == "Bcells":
        mask = analysis.accessibility.columns.get_level_values("cell_type").str.contains("Bcell|Bulk|CLL")
    elif cell_type == "Monos":
        mask = analysis.accessibility.columns.get_level_values("cell_type").str.contains("Mono")
    else:
        mask = analysis.accessibility.columns.get_level_values("cell_type").str.contains(cell_type.replace("cells", ""))

    df = analysis.accessibility[analysis.accessibility.columns[mask]]
    df = df[df.columns[df.columns.get_level_values("timepoint") != 240]]
    df = df.sort_index(axis=1, level=['cell_type', 'patient_id', 'timepoint'])

    dg = diff_genes.loc[(diff_genes['qvalue'] < 1e-4) & (diff_genes['cellType'] == cell_type), "gene"].drop_duplicates()
    diff_peaks = analysis.coverage_annotated[analysis.coverage_annotated['gene_name'].isin(dg)].drop_duplicates().index

    g = sns.clustermap(
        df.loc[diff_peaks, :],
        z_score=0, metric="correlation", robust=True, col_cluster=False,
        yticklabels=False, xticklabels=df.columns.get_level_values("sample_name"), rasterized=True,
        figsize=(max(0.12 * df.shape[1], 6) , 6)
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_row_dendrogram.set_rasterized(True)
    g.savefig(os.path.join("results", "scRNA_diff_genes_over_time.{}_specific.clustermap.svg".format(cell_type)), dpi=300, bbox_inches="tight")

    for patient_id in diff_genes['patient'].drop_duplicates():

        df2 = df[df.columns[df.columns.get_level_values("patient_id") == patient_id]]
        df2 = df2[df2.columns[df2.columns.get_level_values("timepoint") != 240]]
        df2 = df2.sort_index(axis=1, level=['cell_type', 'patient_id', 'timepoint'])

        if df2.shape[1] <= 2:
            continue
        g = sns.clustermap(
            df2.loc[diff_peaks, :],
            z_score=0, metric="correlation", robust=True, col_cluster=False,
            yticklabels=False, xticklabels=df2.columns.get_level_values("sample_name"), rasterized=True,
            figsize=(max(0.12 * df2.shape[1], 6) , 6)
        )
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        g.ax_row_dendrogram.set_rasterized(True)
        g.savefig(os.path.join("results", "scRNA_diff_genes_over_time.{}-{}_specific.clustermap.svg".format(cell_type, patient_id)), dpi=300, bbox_inches="tight")
