#!/usr/bin/env python

"""
"""

import matplotlib
matplotlib.use('Agg')

import argparse
import os
import random
import string
import sys

import GPclust
import GPy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import LabelEncoder

sns.set_style("white")


# random seed
SEED = int("".join(
    LabelEncoder()
    .fit(list(string.ascii_uppercase))
    .transform(list("BOCKLAB")).astype(str)))
random.seed(SEED)
np.random.seed(SEED)


def parse_arguments():
    """
    Argument Parsing.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-type", dest="cell_type", default="CLL", type=str)
    parser.add_argument("--alpha", dest="alpha", default=0.01, type=float)
    parser.add_argument("--output-prefix", dest="output_prefix", default="gp_fit_job.MOHGP", type=str)
    parser.add_argument("--output-dir", dest="output_dir", default="/scratch/users/arendeiro/gp_fit_job/combat/output", type=str)

    args = parser.parse_args()

    return args


def read_matrix(cell_type):
    """
    """
    print("Reading matrix for cell type {}.".format(cell_type))
    matrix = pd.read_csv(os.path.join(
        "/", "home", "arendeiro", "cll-time_course", "results",
        "cll-time_course" + "_peaks.coverage.log2.z_score.qnorm.pca_fix.csv"), index_col=0, header=range(9))
    X = matrix.loc[:, matrix.columns.get_level_values("cell_type").isin([cell_type])]

    return X


def read_fits(cell_type, alpha=0.05):
    """
    """
    print("Reading fits from GPs.")
    fits = pd.read_csv(os.path.join(
        "/", "home", "arendeiro", "cll-time_course", "results",
        ".".join(["gp_fit_job", "combat", "GPy", "all_fits.csv"])), index_col=0)

    return fits[(fits['p_value'] < alpha) & (fits['cell_type'] == cell_type)].index


def fit_gaussian_process(matrix, output_dir, output_prefix, cell_type):
    """
    """

    # Sort by timepoint
    matrix = matrix.sort_index(axis=1, level="timepoint")

    # z score row-wise
    matrix = matrix.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

    # Get timepoints
    times = np.log2(1 + matrix.columns.get_level_values("timepoint").str.replace("d", "").astype(int).values)

    # Make kernels representing underlying process and mean and deviation from it
    k_underlying = GPy.kern.Matern52(input_dim=1, variance=1.0, lengthscale=times.max() / 3.)
    k_corruption = GPy.kern.Matern52(input_dim=1, variance=0.5, lengthscale=times.max() / 3.) + GPy.kern.White(1, variance=0.01)

    print("Fitting.")
    model = GPclust.MOHGP(
        X=times.reshape(-1, 1), Y=matrix.values,
        kernF=k_underlying, kernY=k_corruption,
        K=4, alpha=1., prior_Z="DP")
    model.hyperparam_opt_interval = 1000
    model.hyperparam_opt_args['messages'] = True

    model.optimize(verbose=True)

    print("Finished optimization.")

    model.reorder()

    # Save optimized model parameters
    model.save(os.path.join(output_dir, output_prefix + "." + cell_type + ".MOHCP.fitted_model.hd5"))

    # Save posterior probability matrix Phi
    model.phi.tofile(os.path.join(output_dir, output_prefix + "." + cell_type + ".MOHCP.posterior_probs_phi.np"))

    return model


def plot(model, matrix, output_dir, output_prefix, cell_type):

    print("Plotting cluster posteriors.")
    # Plot clusters
    fig = plt.figure()
    model.plot(newfig=False, on_subplots=True, colour=True, in_a_row=False, joined=False, errorbars=False)
    for ax in fig.axes:
        ax.set_rasterized(True)
        ax.set_ylabel("Chromatin accessibility")
        ax.set_xlabel("Time (log2)")
    fig.savefig(os.path.join(output_dir, output_prefix + "." + cell_type + ".MOHCP.fitted_model.clusters.svg"), dpi=300, bbox_inches="tight")

    print("Plotting parameters/probabilities.")
    # Posterior parameters
    fig, axis = plt.subplots(2, 1,
        gridspec_kw={'height_ratios':[12, 1]},
        figsize=(3 * 4, 1 * 4),
        tight_layout=True)
    mat = axis[0].imshow(model.phi.T, cmap=plt.get_cmap("hot"), vmin=0, vmax=1, aspect='auto')
    axis[0].set_xlabel('Region index')
    axis[0].set_ylabel('Cluster index')
    axis[1].set_aspect(0.1)
    plt.colorbar(mat, cax=axis[1], label="Posterior probability", orientation="horizontal")
    fig.savefig(os.path.join(output_dir, output_prefix + "." + cell_type + ".MOHCP.fitted_model.posterior_probs.svg"), dpi=300, bbox_inches="tight")

    # Assignment probabilities
    g = sns.clustermap(
        model.phi.T,
        cmap=plt.get_cmap("hot"), vmin=0, vmax=1, xticklabels=False, rasterized=True,
        figsize=(3, 0.2 * model.phi.T.shape[0]), cbar_kws={"label": "Posterior probability"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + "." + cell_type + ".MOHCP.fitted_model.posterior_probs.clustermap.svg"), dpi=300, bbox_inches="tight")

    # Clustermap with cluster assignments
    print("Plotting clusters.")
    tp = pd.Series(matrix.columns.get_level_values("timepoint").str.replace("d", "").astype(int), index=matrix.columns).sort_values()

    g2 = sns.clustermap(
        matrix.loc[:, tp.index].T,
        col_colors=[plt.get_cmap("Paired")(i) for i in np.argmax(model.phi,1)],
        row_cluster=False, col_cluster=True, z_score=1, xticklabels=False, yticklabels=matrix.loc[:, tp.index].columns.get_level_values("sample_name"),
        rasterized=True, figsize=(8, 0.2 * matrix.shape[1]), metric="correlation", robust=True)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_col_dendrogram.set_rasterized(True)
    g2.savefig(os.path.join(output_dir, output_prefix + "." + cell_type + ".MOHCP.fitted_model.clustermap.cluster_labels.svg"), dpi=300, bbox_inches="tight")

    matrix_mean = matrix.loc[:, tp.index].T.groupby(level="timepoint").mean()
    g3 = sns.clustermap(
        matrix_mean,
        col_colors=[plt.get_cmap("Paired")(i) for i in np.argmax(model.phi,1)],
        row_cluster=False, col_cluster=True, z_score=1, xticklabels=False, yticklabels=True,
        rasterized=True, figsize=(8, 0.2 * matrix_mean.shape[0]), metric="correlation", robust=True)
    g3.ax_heatmap.set_xticklabels(g3.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g3.ax_heatmap.set_yticklabels(g3.ax_heatmap.get_yticklabels(), rotation=0)
    g3.ax_col_dendrogram.set_rasterized(True)
    g3.savefig(os.path.join(output_dir, output_prefix + "." + cell_type + ".MOHCP.fitted_model.mean_acc.clustermap.cluster_labels.svg"), dpi=300, bbox_inches="tight")


    # # Filter more stringently for cluster assignments
    # # Threshold probabilities to filter out some regions
    # posterior_threshold = 0.8
    # phi = pd.DataFrame(model.phi.T)
    # assigned = (phi > posterior_threshold).any(axis=1)

    # g2 = sns.clustermap(
    #     matrix.loc[:, tp.index].T,
    #     col_colors=[plt.get_cmap("Paired")(i) for i in np.argmax(model.phi,1)],
    #     row_cluster=False, col_cluster=True, z_score=1, xticklabels=False, yticklabels=matrix.loc[:, tp.index].columns.get_level_values("sample_name"),
    #     rasterized=True, figsize=(8, 0.2 * matrix.shape[1]), metric="correlation", robust=True)
    # g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    # g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    # g2.ax_col_dendrogram.set_rasterized(True)
    # g2.savefig(os.path.join(output_dir, output_prefix + "." + cell_type + ".MOHCP.fitted_model.clustermap.cluster_labels.svg"), dpi=300, bbox_inches="tight")


def main():
    """
    Program's main entry point.
    """
    # Parse command-line arguments.
    args = parse_arguments()

    print("Starting.")
    print("Arguments: {}.".format(args))

    # Do it.
    matrix = read_matrix(args.cell_type)

    varying = read_fits(args.cell_type, alpha=args.alpha)

    model = fit_gaussian_process(matrix.loc[varying, :], args.output_dir, args.output_prefix, args.cell_type)

    plot(model, matrix.loc[varying, :], args.output_dir, args.output_prefix, args.cell_type)

    print("Done.")


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
