#!/usr/bin/env python

"""
Fit hierarchical gaussian processes to cluster chromatin accessibility
features on their temporal dynamic pattern.
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
    """
    Script's entry point.
    """
    # Parse command-line arguments.
    args = parse_arguments()

    print("Starting.")
    print("Arguments: {}.".format(args))

    # Read matrix
    matrix = read_matrix(args.matrix_file, args.cell_type, args.matrix_header_range)

    varying = read_fits(args.fits_file, args.cell_type, alpha=args.alpha)

    model = fit_gaussian_process(
        matrix.loc[varying, :], args.output_dir, args.output_prefix, args.cell_type, args.n_clust, not args.linear)

    plot(model, matrix.loc[varying, :], args.output_dir,
         args.output_prefix, args.cell_type)

    print("Done.")


def parse_arguments():
    """
    Argument Parsing.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix-file", dest="matrix_file", type=str)
    parser.add_argument("--fits-file", dest="fits_file", type=str)
    parser.add_argument("--matrix-header-range",
                        dest="matrix_header_range", default=9, type=int)
    parser.add_argument("--n_clust",
                        dest="n_clust", default=4, type=int)
    parser.add_argument("--cell-type", dest="cell_type",
                        default="CLL", type=str)
    parser.add_argument("--output-prefix", dest="output_prefix",
                        default="gp_fit_job", type=str)
    parser.add_argument("--output-dir", dest="output_dir",
                        default="/scratch/users/arendeiro/cll-time_course/mohgp_fit_job/fits", type=str)
    parser.add_argument("--alpha", dest="alpha", default=0.05, type=float)
    parser.add_argument("--linear", dest="linear", action="store_true")

    args = parser.parse_args()

    return args


def read_matrix(matrix_file, cell_type="CLL", header_range=9):
    """
    Read matrix of (n_features, n_samples) and return samples from respective cell type
    (sample annotation as column metadata) with length `header_range`.
    """
    print("Reading matrix {} for cell type {}.".format(matrix_file, cell_type))
    matrix = pd.read_csv(matrix_file, index_col=0, header=range(header_range))
    X = matrix.loc[:, matrix.columns.get_level_values("cell_type").isin([cell_type])]
    t = ["240d", "280d"]
    if cell_type in ["CD4", "CD8"]:
        t += ["1d"]
    elif cell_type in ["Bcell"]:
        t += ["150d"]
    X = X.loc[:, ~X.columns.get_level_values("timepoint").isin(t)]
    X = X.astype(float).T.groupby(['patient_id', 'timepoint']).mean().T

    print("Finished reading matrix with {} features and {} samples of '{}' cell type.".format(
        X.shape[0], X.shape[1], cell_type))

    return X


def read_fits(fits_file, cell_type, alpha=0.05):
    """
    Read Gaussian Process fits from `gp_fit_job.py` in and return variable regions for given cell type.
    """
    print("Reading fits from GPs.")
    fits = pd.read_csv(fits_file, index_col=0)

    return fits[(fits['p_value'] < alpha) & (fits['cell_type'] == cell_type)].index


def fit_gaussian_process(matrix, output_dir, output_prefix, cell_type, n_clust_guess, x_log_transform=True):
    """
    Fit Mixture of Hierarchical Gaussian Processes for clustering of features.
    """
    # Sort by timepoint
    matrix = matrix.sort_index(axis=1, level="timepoint")

    # z score row-wise
    matrix = matrix.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

    # Get timepoints
    times = pd.Series(
        matrix
        .columns
        .get_level_values('timepoint')
        .str.replace("d", "")
        .astype(int))
    # Make sure each timepoint is represented at least twice
    times = times[times.isin(times.value_counts()[times.value_counts() > 1].index)]
    if x_log_transform:
        times = np.log2(1 + times)
    matrix = matrix.iloc[:, times.index]

    # Make kernels representing underlying process and mean and deviation from it
    k_underlying = GPy.kern.Matern52(input_dim=1, variance=1.0, lengthscale=times.max() / 3.)
    k_corruption = GPy.kern.Matern52(input_dim=1, variance=0.5, lengthscale=times.max() / 3.) + GPy.kern.White(1, variance=0.01)

    print("Fitting.")
    model = GPclust.MOHGP(
        X=times.reshape(-1, 1), Y=matrix.values,
        kernF=k_underlying, kernY=k_corruption,
        K=n_clust_guess, alpha=1., prior_Z="DP")
    model.hyperparam_opt_interval = 1000
    model.hyperparam_opt_args['messages'] = True

    model.optimize(verbose=True)

    print("Finished optimization.")

    # Order clusters by their size
    model.reorder()

    # Save optimized model parameters
    model.save(os.path.join(output_dir, output_prefix + "." + cell_type + ".mohgp.fitted_model.hd5"))

    # Save posterior probability matrix Phi
    model.phi.tofile(os.path.join(output_dir, output_prefix + "." + cell_type + ".mohgp.posterior_probs_phi.npy"))

    return model


def plot(model, matrix, output_dir, output_prefix, cell_type):
    """
    Fitted model's parameters and probabilities.
    """

    print("Plotting cluster posteriors.")
    # Plot clusters
    fig = plt.figure()
    model.plot(newfig=False, on_subplots=True, colour=True, in_a_row=False, joined=False, errorbars=False)
    for ax in fig.axes:
        ax.set_rasterized(True)
        ax.set_ylabel("Chromatin accessibility")
        ax.set_xlabel("Time (log2)")
    fig.savefig(os.path.join(output_dir, output_prefix + "." + cell_type + ".mohgp.fitted_model.clusters.svg"), dpi=300, bbox_inches="tight")

    print("Plotting parameters/probabilities.")
    # Posterior parameters
    fig, axis = plt.subplots(
        2, 1,
        gridspec_kw={'height_ratios':[12, 1]},
        figsize=(3 * 4, 1 * 4),
        tight_layout=True)
    mat = axis[0].imshow(model.phi.T, cmap=plt.get_cmap("hot"), vmin=0, vmax=1, aspect='auto')
    axis[0].set_xlabel('Region index')
    axis[0].set_ylabel('Cluster index')
    axis[1].set_aspect(0.1)
    plt.colorbar(mat, cax=axis[1], label="Posterior probability", orientation="horizontal")
    fig.savefig(os.path.join(output_dir, output_prefix + "." + cell_type + ".mohgp.fitted_model.posterior_probs.svg"), dpi=300, bbox_inches="tight")

    # Assignment probabilities
    g = sns.clustermap(
        model.phi.T,
        cmap=plt.get_cmap("hot"), vmin=0, vmax=1, xticklabels=False, rasterized=True,
        figsize=(3, 0.2 * model.phi.T.shape[0]), cbar_kws={"label": "Posterior probability"})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + "." + cell_type + ".mohgp.fitted_model.posterior_probs.clustermap.svg"), dpi=300, bbox_inches="tight")

    # Clustermap with cluster assignments
    print("Plotting clusters.")
    tp = pd.Series(matrix.columns.get_level_values("timepoint").str.replace("d", "").astype(int), index=matrix.columns).sort_values()

    g2 = sns.clustermap(
        matrix.loc[:, tp.index].T,
        col_colors=[plt.get_cmap("Paired")(i) for i in np.argmax(model.phi,1)],
        row_cluster=False, col_cluster=True, z_score=1, xticklabels=False,
        rasterized=True, figsize=(8, 0.2 * matrix.shape[1]), metric="correlation", robust=True)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_col_dendrogram.set_rasterized(True)
    g2.savefig(os.path.join(output_dir, output_prefix + "." + cell_type + ".mohgp.fitted_model.clustermap.cluster_labels.svg"), dpi=300, bbox_inches="tight")

    matrix_mean = matrix.loc[:, tp.index].T.groupby(level="timepoint").mean()
    g3 = sns.clustermap(
        matrix_mean,
        col_colors=[plt.get_cmap("Paired")(i) for i in np.argmax(model.phi,1)],
        row_cluster=False, col_cluster=True, z_score=1, xticklabels=False, yticklabels=True,
        rasterized=True, figsize=(8, 0.2 * matrix_mean.shape[0]), metric="correlation", robust=True)
    g3.ax_heatmap.set_xticklabels(g3.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g3.ax_heatmap.set_yticklabels(g3.ax_heatmap.get_yticklabels(), rotation=0)
    g3.ax_col_dendrogram.set_rasterized(True)
    g3.savefig(os.path.join(output_dir, output_prefix + "." + cell_type + ".mohgp.fitted_model.mean_acc.clustermap.cluster_labels.svg"), dpi=300, bbox_inches="tight")


    # # Filter more stringently for cluster assignments
    # # Threshold probabilities to filter out some regions
    # posterior_threshold = 0.8
    # phi = pd.DataFrame(model.phi.T)
    # assigned = (phi > posterior_threshold).any(axis=1)

    # g2 = sns.clustermap(
    #     matrix.loc[:, tp.index].T,
    #     col_colors=[plt.get_cmap("Paired")(i) for i in np.argmax(model.phi,1)],
    #     row_cluster=False, col_cluster=True, z_score=1, xticklabels=False,
    #     rasterized=True, figsize=(8, 0.2 * matrix.shape[1]), metric="correlation", robust=True)
    # g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    # g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    # g2.ax_col_dendrogram.set_rasterized(True)
    # g2.savefig(os.path.join(output_dir, output_prefix + "." + cell_type + ".mohgp.fitted_model.clustermap.cluster_labels.svg"), dpi=300, bbox_inches="tight")


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
