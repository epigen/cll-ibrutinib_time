
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


def bulk_pca_plots(df):
    """
    Plots for figure 1C.
    """
    from sklearn.decomposition import PCA

    # PCA
    pca = PCA()
    X_hat = pca.fit_transform(df.T)

    X_hat = pd.DataFrame(X_hat, index=df.columns)


    # Get cell type centroids
    cell_type_X_hat = X_hat.groupby(level=['cell_type']).mean()


    fig, axis = plt.subplots(1, 3, figsize=(4 * 3, 1 * 4), tight_layout=True)

    # Plot cell type centroids in all subplots
    for ax in axis:
        for ct in cell_type_X_hat.index:
            ax.scatter(
                cell_type_X_hat.loc[ct, 0], cell_type_X_hat.loc[ct, 1], label=ct,
                edgecolors='black', linewidths=0.4, alpha=1.0)
    axis[0].legend()

    bulk_X = X_hat.loc[X_hat.index.get_level_values("cell_type") == "Bulk", :]

    # Plot bulk samples in first subplot
    axis[0].scatter(
        bulk_X.loc[:, 0], bulk_X.loc[:, 1],
        linewidths=0.0, alpha=0.5)

    # Label patient's bulk samples in second subplot
    pat_enc = LabelEncoder().fit(bulk_X.index.get_level_values("patient_id"))

    for i in bulk_X.index:
        axis[1].scatter(
            bulk_X.loc[i, 0], bulk_X.loc[i, 1],
            label=i[1],
            c=plt.get_cmap("Paired")(pat_enc.transform([i[1]])),
            linewidths=0.0, alpha=0.5)
    axis[1].legend()

    # Label time in third subplot
    time_colors = plt.get_cmap("inferno")(bulk_X.index.get_level_values("timepoint"))

    axis[2].scatter(
        bulk_X.loc[:, 0], bulk_X.loc[:, 1],
        c=time_colors,
        linewidths=0.0, alpha=0.5)

    fig.savefig(os.path.join("results", "fig1c.svg"), bbox_inches="tight")


    # Only PBGY bulk
    pbgy_X = X_hat.loc[(X_hat.index.get_level_values("cell_type") == "Bulk") & (X_hat.index.get_level_values("patient_id") == "PBGY"), :]

    fig, axis = plt.subplots(1, 1, figsize=(4 * 1, 1 * 4), tight_layout=True)
    for ct in cell_type_X_hat.index:
        axis.scatter(
            cell_type_X_hat.loc[ct, 0], cell_type_X_hat.loc[ct, 1], label=ct,
            edgecolors='black', linewidths=0.4, alpha=1.0)
    # Label time in third subplot
    time_colors = plt.get_cmap("inferno")(pbgy_X.index.get_level_values("timepoint"))
    axis.scatter(
        pbgy_X.loc[:, 0], pbgy_X.loc[:, 1],
        c=time_colors,
        linewidths=0.0, alpha=0.5)
    fig.savefig(os.path.join("results", "fig1c.pbgy.svg"), bbox_inches="tight")



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
    chunks = 2000
    output_dir = "/scratch/users/arendeiro/gp_fit_job"
    library = "GPy"

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
    import itertools
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
    fig.savefig(os.path.join("results", output_prefix + ".MOHCP.cluster_name.counts.barplot.svg"), bbox_inches="tight")

    # Physical overlap
    def overlap((a, b), func=max):
        return (
            func(
                len(set(a).intersection(set(b))),
                len(set(b).intersection(set(a))))
            /
            float(func(len(a), len(b)))
        ) * 100

    a = assignments.sort_values(['cell_type', 'cluster']).set_index("index")
    g = a.groupby(['cell_type', 'cluster']).groups.values()
    n = map(lambda x: "_".join([str(i) for i in x]), a.groupby(['cell_type', 'cluster']).groups.keys())
    comb = pd.DataFrame(np.array(map(overlap, itertools.product(g, repeat=2))).reshape((28, 28)))
    ns = pd.DataFrame(np.array(map(lambda x: ":".join(x), itertools.product(n, repeat=2))).reshape((28, 28)))
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
    fig.savefig(os.path.join("results", output_prefix + ".MOHCP.cluster_name.overlap.heatmap.svg"), dpi=300)


def plot_lola_enrichments():
    import scipy
    comp_variable ='comparison_name'
    enrichment_table = pd.read_csv(os.path.join("results_deconvolve", "GPclust", "gp_fit_job.sorted.GPclust.lola.csv"))

    top_n = 20

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
        g.fig.savefig(os.path.join("results", "fig2c.{}".format(cell_type) + ".lola.cluster_specific.svg"), bbox_inches="tight", dpi=300)


    # All cell types and clusters together
    top_n = 20
    all_samples = enrichment_table.index
    no_mix = enrichment_table[enrichment_table[comp_variable].str.contains("_")].index
    no_bulk = enrichment_table[(enrichment_table[comp_variable].str.contains("_")) & (~enrichment_table[comp_variable].str.contains("ulk"))].index

    for mask, label in [
        (all_samples, "all_samples"),
        (all_samples, "no_mix"),
        (all_samples, "no_bulk")
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
        g.fig.savefig(os.path.join("results", "fig2c.{}".format(label) + ".lola.cluster_specific.svg"), bbox_inches="tight", dpi=300)

        # plot sorted heatmap
        g = sns.clustermap(
            lola_pivot.loc[top_terms, :], figsize=(max(6, 0.12 * shape[1]), max(6, 0.12 * shape[0])), square=True,
            cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"}, metric="correlation", col_cluster=False)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join("results", "fig2c.{}.sorted".format(label) + ".lola.cluster_specific.svg"), bbox_inches="tight", dpi=300)


        # plot correlation
        g = sns.clustermap(
            lola_pivot.loc[top_terms, :].corr(), figsize=(max(6, 0.12 * shape[1]), max(6, 0.12 * shape[1])), square=True,
            cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join("results", "fig2c.{}.corr".format(label) + ".lola.cluster_specific.svg"), bbox_inches="tight", dpi=300)



def plot_enrichments(top_n=10):
    import scipy

    comp_variable ='comparison_name'
    enrichment_table = pd.read_csv(os.path.join("results_deconvolve", "GPclust", "gp_fit_job.sorted.GPclust.enrichr.csv"))
    enrichment_table = enrichment_table[(enrichment_table[comp_variable].str.contains("_")) & (~enrichment_table[comp_variable].str.contains("ulk"))]

    # enrichment_table["description"] = enrichment_table["description"].str.decode("utf-8")
    enrichment_table["log_p_value"] = (-np.log10(enrichment_table["p_value"])).replace({np.inf: 300})

    for gene_set_library in enrichment_table["gene_set_library"].unique():
        print(gene_set_library)
        if gene_set_library == "Epigenomics_Roadmap_HM_ChIP-seq":
            continue

        enr = enrichment_table[enrichment_table['gene_set_library'] == gene_set_library]

        # Normalize enrichments per dataset with Z-score prior to comparing various region sets
        for comparison_name in enr['comparison_name'].drop_duplicates():
            mask = enr['comparison_name'] == comparison_name
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

        # plot clustered heatmap
        shape = enrichr_pivot[list(set(top_terms))].shape
        g = sns.clustermap(enrichr_pivot[list(set(top_terms))].T, figsize=(max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
            cbar_kws={"label": "-log10(p-value) of enrichment\nof differential genes"}, metric="correlation", square=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join("results", "fig2d.all_cell_types_clusters" + ".enrichr.{}.cluster_specific.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)

        # plot sorted heatmap
        g = sns.clustermap(enrichr_pivot[list(set(top_terms))].T, figsize=(max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
            cbar_kws={"label": "-log10(p-value) of enrichment\nof differential genes"}, metric="correlation", square=True, col_cluster=False)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join("results", "fig2d.all_cell_types_clusters" + ".sorted.enrichr.{}.cluster_specific.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)


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
    shape = enrichr_pivot[list(set(top_terms))].shape
    g = sns.clustermap(
        enrichr_pivot[list(set(top_terms))].T.corr(), square=True,
        cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"})
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.fig.savefig(os.path.join("results", "fig2d.all_cell_types_clusters.corr" + ".enrichr.svg"), bbox_inches="tight", dpi=300)


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


def fig2e():
    cll = analysis.accessibility.loc[:, 
        (analysis.accessibility.columns.get_level_values("patient_id") != "KI") &
        (analysis.accessibility.columns.get_level_values("cell_type") == "CLL")]

    bcell = analysis.accessibility.loc[:, 
        (analysis.accessibility.columns.get_level_values("patient_id") != "KI") &
        (analysis.accessibility.columns.get_level_values("cell_type") == "Bcell")]

    # fc = np.log2(bcell.mean(axis=1) / cll.mean(axis=1)).sort_values()
    # sns.distplot(fc)

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


    diff = res.groupby(['patient_id']).apply(lambda x: x.loc[x['timepoint'].argmax(), 'ratio'].squeeze() / x.loc[x['timepoint'] == 0, 'ratio'].squeeze())
    g = sns.barplot(data=diff.reset_index(), x='patient_id', y=0)
    

    g = sns.clustermap(cll.loc[fc.abs().sort_values().tail(2000).index, :], yticklabels=False, metric="correlation")
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


def transcription_factor_accessibility():
    from glob import glob
    import pybedtools

    all_res = pd.DataFrame()
    for factor_name in ["NFKB", "PU1", "ATF2", "BATF", "NFIC", "IRF4", "RUNX3", "BCL11A", "POL2", "GATA1", "NANOG", "POU5F", "CTCF"]:
        print(factor_name)
        # get consensus NFKB regions from LOLA database
        transcription_factor_regions_sets = glob("/home/arendeiro/resources/regions/LOLACore/hg19/encode_tfbs/regions/*{}*".format(factor_name.capitalize()))[:5]
        bt = pybedtools.BedTool(transcription_factor_regions_sets[0]).sort()
        for fn in transcription_factor_regions_sets[1:]:
            bt = bt.cat(pybedtools.BedTool(fn).sort().merge())
        bt = bt.merge()
        bt.saveas(os.path.join("data", "external", "TF_binding_sites" + factor_name + ".bed"))

        # get regions overlapping with NFKB sites
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
        all_res = all_res.append(d, ignore_index=True)


        g = sns.clustermap(res.dropna().T, col_cluster=False, z_score=1, rasterized=True, figsize=(res.shape[0] * 0.12, res.shape[1] * 0.12), row_cluster=False)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.savefig(os.path.join("results", "fig2X.{}_binding.per_quantile.sorted.svg".format(factor_name)), dpi=300, bbox_inches="tight")

        res_mean = res.dropna().T.groupby(level=['cell_type', 'timepoint']).mean()
        g = sns.clustermap(res_mean.dropna(), col_cluster=False, z_score=1, rasterized=False, square=True, row_cluster=False)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.savefig(os.path.join("results", "fig2X.{}_binding.per_quantile.mean_patients.sorted.svg".format(factor_name)), dpi=300, bbox_inches="tight")

        # d = pd.melt(res.dropna().T.groupby(level=['cell_type', 'timepoint']).mean().reset_index(), id_vars=['cell_type', 'timepoint'], var_name="quantile", value_name='accessibility')
        # g = sns.factorplot(data=d, x="timepoint", y="accessibility", hue="quantile", col="cell_type")
        # g.savefig(os.path.join("results", "fig2X.{}_binding.per_quantile.mean_patients_quantiles.col.lineplots.svg".format(factor_name)), dpi=300, bbox_inches="tight")

        d = res.dropna().T.groupby(level=['cell_type', 'timepoint']).mean().mean(1).reset_index()
        # g = sns.factorplot(data=d, x="timepoint", y=0, col="cell_type")
        # g.savefig(os.path.join("results", "fig2X.{}_binding.per_quantile.mean_patients_all_sites.col.lineplots.svg".format(factor_name)), dpi=300, bbox_inches="tight")
        g = sns.factorplot(data=d, x="timepoint", y=0, hue="cell_type")
        g.savefig(os.path.join("results", "fig2X.{}_binding.per_quantile.mean_patients_all_sites.hue.lineplots.svg".format(factor_name)), dpi=300, bbox_inches="tight")

        d = res.dropna().T.groupby(level=['cell_type', 'timepoint']).mean().mean(1).reset_index()
        d["transcription_factor"] = factor_name
        all_res = all_res.append(d, ignore_index=True)

    all_res = all_res.rename(columns={0: "accessibility"})
    g = sns.factorplot(data=all_res, x="timepoint", y="accessibility", hue="cell_type", col="transcription_factor", col_wrap=3)
    g.savefig(os.path.join("results", "fig2X.all_factor_binding.mean_patients.cell_type_hue.lineplots.svg".format(factor_name)), dpi=300, bbox_inches="tight")

    d = all_res[all_res['cell_type'] == "CLL"]
    g = sns.factorplot(data=d, x="timepoint", y="accessibility", hue="transcription_factor")
    g.savefig(os.path.join("results", "fig2X.all_factor_binding.mean_patients.CLL_only.lineplots.svg".format(factor_name)), dpi=300, bbox_inches="tight")


    d_h = d.groupby('transcription_factor')['accessibility'].apply(scipy.stats.zscore).apply(pd.Series)
    d_h.columns = d['timepoint'].unique()

    g = sns.clustermap(d_h, col_cluster=False, row_cluster=True, vmin=-5, vmax=5, square=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join("results", "fig2X.all_factor_binding.mean_patients.CLL_only.clustermap.svg".format(factor_name)), dpi=300, bbox_inches="tight")


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


def cytokine_receptor_repertoire():
    """
    """
    assert hasattr(analysis, "gene_annotation")

    # Get ligand-receptor info 
    lr = pd.read_excel("https://images.nature.com/original/nature-assets/ncomms/2015/150722/ncomms8866/extref/ncomms8866-s3.xlsx", 1)
    lr = lr.loc[
        ~lr['Pair.Evidence'].str.contains("EXCLUDED"),
        ['Pair.Name', 'Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol']]
    lr.columns = ['pair', 'ligand', 'receptor']
    lr_genes = lr['ligand'].unique().tolist() + lr['receptor'].unique().tolist()

    # Get cell adhesion molecules (CAMs)
    from bioservices.kegg import KEGG
    import requests

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
    full_acc = acc.loc[
        acc.index.get_level_values("gene_name").isin(lr_genes),
        (acc.columns.get_level_values("patient_id") != "KI") &
        (acc.columns.get_level_values("timepoint") <= 150) &
        (acc.columns.get_level_values("cell_type") != "Bulk")]
    red_acc = full_acc.groupby(level="gene_name").mean()

    full_acc_time = full_acc.T.groupby(level=["cell_type", "timepoint"]).mean().T
    red_acc_time = red_acc.T.groupby(level=["cell_type", "timepoint"]).mean().T

    w, h = red_acc.shape[0] * 0.01, red_acc.shape[1] * 0.12
    g = sns.clustermap(
        red_acc.T, col_cluster=True, row_cluster=True, z_score=1,
        figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.gene_level.clustermap.svg"), dpi=300, bbox_inches="tight")

    w, h = full_acc.shape[0] * 0.01, full_acc.shape[1] * 0.12
    g = sns.clustermap(
        full_acc.T, col_cluster=True, row_cluster=True, z_score=1,
        figsize=(w, h), robust=True, xticklabels=False, row_linkage=g.dendrogram_row.linkage, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.region_level.clustermap.svg"), dpi=300, bbox_inches="tight")

    w, h = red_acc_time.shape[0] * 0.01, red_acc_time.shape[1] * 0.12
    g = sns.clustermap(
        red_acc_time.T, col_cluster=True, row_cluster=False, z_score=1,
        figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.gene_level.timepoint_mean.clustermap.svg"), dpi=300, bbox_inches="tight")

    w, h = full_acc_time.shape[0] * 0.001, max(6, full_acc_time.shape[1] * 0.12)
    g = sns.clustermap(
        full_acc_time.T, col_cluster=True, row_cluster=False, z_score=1,
        figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.region_level.timepoint_mean.clustermap.svg"), dpi=300, bbox_inches="tight")

    # Differential regulatory elements with changing ligand or receptor
    output_dir = "/scratch/users/arendeiro/gp_fit_job"
    output_prefix = "gp_fit_job"
    library = "GPy"
    matrix_name="sorted"
    alpha = 0.01
    fits = pd.read_csv(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "all_fits.csv"])), index_col=0)
    # Get variable regions for all cell types
    variable = fits[(fits['p_value'] < alpha)].index.drop_duplicates()

    full_acc_sig = full_acc.loc[full_acc.index.get_level_values("index").isin(variable.tolist())]
    red_acc_sig = full_acc_sig.groupby(level="gene_name").mean()
    full_acc_time_sig = full_acc_sig.T.groupby(level=["cell_type", "timepoint"]).mean().T
    red_acc_time_sig = red_acc_sig.T.groupby(level=["cell_type", "timepoint"]).mean().T

    w, h = red_acc_sig.shape[0] * 0.01, red_acc_sig.shape[1] * 0.12
    g = sns.clustermap(
        red_acc_sig.T, col_cluster=True, row_cluster=True, z_score=1,
        figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.gene_level.sig_only.clustermap.svg"), dpi=300, bbox_inches="tight")

    w, h = full_acc_sig.shape[0] * 0.01, full_acc_sig.shape[1] * 0.12
    g = sns.clustermap(
        full_acc_sig.T, col_cluster=True, row_cluster=True, z_score=1,
        figsize=(w, h), robust=True, xticklabels=False, row_linkage=g.dendrogram_row.linkage, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.region_level.sig_only.clustermap.svg"), dpi=300, bbox_inches="tight")

    w, h = red_acc_time_sig.shape[0] * 0.01, red_acc_time_sig.shape[1] * 0.12
    g = sns.clustermap(
        red_acc_time_sig.T, col_cluster=True, row_cluster=False, z_score=1,
        figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.gene_level.sig_only.timepoint_mean.clustermap.svg"), dpi=300, bbox_inches="tight")

    w, h = full_acc_time_sig.shape[0] * 0.01, max(6, full_acc_time_sig.shape[1] * 0.12)
    g = sns.clustermap(
        full_acc_time_sig.T, col_cluster=True, row_cluster=False, z_score=1,
        figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.region_level.sig_only.timepoint_mean.clustermap.svg"), dpi=300, bbox_inches="tight")


    # Get variable regions for each cell type
    for cell_type in acc.columns.levels[3]:
        if cell_type == "Bulk":
            continue
        variable = fits[(fits['p_value'] < alpha) & (fits['cell_type'] == cell_type)].index.drop_duplicates()

        full_acc_sig = full_acc.loc[
            full_acc.index.get_level_values("index").isin(variable.tolist()),
            full_acc.columns.get_level_values("cell_type") == cell_type]
        red_acc_sig = full_acc_sig.groupby(level="gene_name").mean()
        full_acc_time_sig = full_acc_sig.T.groupby(level=["cell_type", "timepoint"]).mean().T
        red_acc_time_sig = red_acc_sig.T.groupby(level=["cell_type", "timepoint"]).mean().T

        w, h = red_acc_sig.shape[0] * 0.12, red_acc_sig.shape[1] * 0.12
        g = sns.clustermap(
            red_acc_sig.T, col_cluster=True, row_cluster=True, z_score=1,
            figsize=(w, h), square=True, robust=True, xticklabels=True, rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="left", fontsize="xx-small")
        g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.{}.gene_level.sig_only.clustermap.svg".format(cell_type)), dpi=300, bbox_inches="tight")

        # same as above, by same order but average per timepoint
        g = sns.clustermap(
            red_acc_time_sig.T.iloc[:, g.dendrogram_col.reordered_ind], col_cluster=False, row_cluster=False, z_score=1,
            figsize=(w, h), square=True, robust=True, xticklabels=True, rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="left", fontsize="xx-small")
        g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.{}.gene_level.sig_only.timepoint_mean.clustermap.sorted.svg".format(cell_type)), dpi=300, bbox_inches="tight")

        w, h = full_acc_sig.shape[0] * 0.05, full_acc_sig.shape[1] * 0.12
        g = sns.clustermap(
            full_acc_sig.T, col_cluster=True, row_cluster=True, z_score=1,
            figsize=(w, h), robust=True, xticklabels=False, row_linkage=g.dendrogram_row.linkage, rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
        g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.{}.region_level.sig_only.clustermap.svg".format(cell_type)), dpi=300, bbox_inches="tight")

        w, h = red_acc_time_sig.shape[0] * 0.12, red_acc_time_sig.shape[1] * 0.12
        g = sns.clustermap(
            red_acc_time_sig.T, col_cluster=True, row_cluster=False, z_score=1,
            figsize=(w, h), robust=True, xticklabels=True, rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="left", fontsize="xx-small")
        g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.{}.gene_level.sig_only.timepoint_mean.clustermap.svg".format(cell_type)), dpi=300, bbox_inches="tight")

        w, h = full_acc_time_sig.shape[0] * 0.05, max(6, full_acc_time_sig.shape[1] * 0.12)
        g = sns.clustermap(
            full_acc_time_sig.T, col_cluster=True, row_cluster=False, z_score=1,
            figsize=(w, h), robust=True, xticklabels=False, rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
        g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.{}.region_level.sig_only.timepoint_mean.clustermap.svg".format(cell_type)), dpi=300, bbox_inches="tight")

        red_acc_time_sig.to_csv(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.{}.gene_level.sig_only.timepoint_mean.clustermap.csv".format(cell_type)), index=True)

        # Only genes belonging to CAM pathway
        if cell_type == "CLL":
            red_acc_time_sig_specific = red_acc_time_sig.loc[
                (red_acc_time_sig.index.isin(cam_genes)) |
                (red_acc_time_sig.index.str.startswith("VEG"))]
            w, h = red_acc_time_sig_specific.shape[0] * 0.12, red_acc_time_sig_specific.shape[1] * 0.12
            g = sns.clustermap(
                red_acc_time_sig_specific.T, col_cluster=True, row_cluster=False, z_score=1,
                figsize=(w, h), robust=True, xticklabels=True, rasterized=True)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="left", fontsize="xx-small")
            g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.{}.gene_level.sig_only.cam_pathway.timepoint_mean.clustermap.svg".format(cell_type)), dpi=300, bbox_inches="tight")

            red_acc_time_sig_specific = red_acc_time_sig.loc[
                (red_acc_time_sig.index.str.startswith("TGF")) |
                (red_acc_time_sig.index.str.startswith("WNT")) |
                (red_acc_time_sig.index.str.startswith("NOTCH")) |
                (red_acc_time_sig.index.str.startswith("NOTCH")) |
                (red_acc_time_sig.index.str.startswith("TNF")) |
                (red_acc_time_sig.index.str.startswith("TLR"))]
            w, h = red_acc_time_sig_specific.shape[0] * 0.12, red_acc_time_sig_specific.shape[1] * 0.12
            g = sns.clustermap(
                red_acc_time_sig_specific.T, col_cluster=True, row_cluster=False, z_score=1,
                figsize=(w, h), robust=True, xticklabels=True, rasterized=True)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="left", fontsize="xx-small")
            g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.{}.gene_level.sig_only.pathways.timepoint_mean.clustermap.svg".format(cell_type)), dpi=300, bbox_inches="tight")

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
            g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.{}.gene_level.sig_only.cytokine-receptors.timepoint_mean.clustermap.svg".format(cell_type)), dpi=300, bbox_inches="tight")



    variable = fits[(fits['p_value'] < alpha)].index.drop_duplicates()

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
    g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.all_cell_types.gene_level.sig_only.clustermap.svg"), dpi=300, bbox_inches="tight")

    # same as above, by same order but average per timepoint
    g = sns.clustermap(
        red_acc_time_sig.T.iloc[:, g.dendrogram_col.reordered_ind], col_cluster=False, row_cluster=False, z_score=1,
        figsize=(w, h), square=True, robust=True, xticklabels=True, rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="left", fontsize="xx-small")
    g.savefig(os.path.join("results", analysis.name + ".ligand-receptor_repertoire.all_cell_types.gene_level.sig_only.timepoint_mean.clustermap.sorted.svg".format(cell_type)), dpi=300, bbox_inches="tight")


def cytokine_interplay():
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
    output_dir = "/scratch/users/arendeiro/gp_fit_job"
    output_prefix = "gp_fit_job"
    library = "GPy"
    matrix_name="sorted"
    alpha = 0.01
    fits = pd.read_csv(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "all_fits.csv"])), index_col=0)
    # Get variable regions for all cell types
    variable = fits[(fits['p_value'] < alpha)].index.drop_duplicates()

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
    acc_mean_gene = acc_mean_gene.loc[:, (acc_mean_gene.columns.get_level_values("patient_id") != "KI")]
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


# BULK ONLY ANALYSIS (Fig 1)
analysis.unsupervised(
    quant_matrix="accessibility",
    samples=[s for s in analysis.samples if s.cell_type == "Bulk" and s.patient_id != "KI"],
    attributes_to_plot=["patient_id", "timepoint", "batch"], plot_prefix="bulk_only")

bulk_pca_plots(analysis.accessibility.loc[:, analysis.accessibility.columns.get_level_values("patient_id") != "KI"])



# TIME SERIES ANALYSIS (Fig 2, 3)

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

output_prefix = 'gp_fit_job'
assignments = pd.read_csv(os.path.join("results_deconvolve", output_prefix + ".GP_fits.MOHCP_clusters.csv"))


# Plot distribution of clusters per cell type dependent on their dynamic pattern
cluster_dynamics(assignments)


# Get enrichments of region clusters
assignments['comparison_name'] = assignments['cell_type'] + "_" + assignments['cluster'].astype(str)


# In addition, get enrichments for all of changing regions for each cell type
_a = assignments.copy()
_a['comparison_name'] = _a['comparison_name'].str.replace("_.*", "")
assignments = assignments.append(_a)

output_prefix = 'gp_fit_job'
matrix_name = "sorted"

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

collect_differential_enrichment(
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



# Plot specially for CLL the enrichments
enrichment_table = pd.read_csv(os.path.join("results_deconvolve", "GPclust", "gp_fit_job.sorted.GPclust.meme_ame.csv"))

enrichment_table = enrichment_table[enrichment_table['comparison_name'].str.contains("CLL")]

plot_differential_enrichment(
    enrichment_table,
    "motif",
    data_type="ATAC-seq",
    direction_dependent=False,
    output_dir="results_deconvolve/GPclust",
    comp_variable="comparison_name",
    output_prefix=".".join([output_prefix, matrix_name, "GPclust", "CLL_only"]),
    top_n=5)
enrichment_table = pd.read_csv(os.path.join("results_deconvolve", "GPclust", "gp_fit_job.sorted.GPclust.lola.csv"))
enrichment_table = enrichment_table[enrichment_table['comparison_name'].str.contains("CLL")]

plot_differential_enrichment(
    enrichment_table,
    "lola",
    data_type="ATAC-seq",
    direction_dependent=False,
    output_dir="results_deconvolve/GPclust",
    comp_variable="comparison_name",
    output_prefix=".".join([output_prefix, matrix_name, "GPclust", "CLL_only"]),
    top_n=25)
enrichment_table = pd.read_csv(os.path.join("results_deconvolve", "GPclust", "gp_fit_job.sorted.GPclust.enrichr.csv"))
enrichment_table = enrichment_table[enrichment_table['comparison_name'].str.contains("CLL")]

plot_differential_enrichment(
    enrichment_table,
    "enrichr",
    data_type="ATAC-seq",
    direction_dependent=False,
    output_dir="results_deconvolve/GPclust",
    comp_variable="comparison_name",
    output_prefix=".".join([output_prefix, matrix_name, "GPclust", "CLL_only"]),
    top_n=5)


# Figure 2c, 3c
plot_lola_enrichments()


# Figure 2d, 3d-e
plot_enrichments()


# Figure 2e
transcription_factor_accessibility()



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
