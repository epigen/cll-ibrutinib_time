#!/usr/bin/env python

"""
"""

import argparse
import os
import numpy as np
import pandas as pd
import GPy
import sys


def parse_arguments():
    """
    Argument Parsing.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-range", dest="data_range", type=str)
    parser.add_argument("--range-delim", dest="range_delim",default="-", type=str)
    parser.add_argument("--cell-type", dest="cell_type", default="CLL", type=str)
    parser.add_argument("--output-prefix", dest="output_prefix", default="gp_fit_job", type=str)
    parser.add_argument("--output-dir", dest="output_dir", default="/scratch/users/arendeiro/gp_fit_job/output", type=str)

    args = parser.parse_args()

    return args


def read_matrix(cell_type):
    """
    """
    deconv_norm = pd.read_csv(os.path.join("/", "home", "arendeiro", "cll-time_course", "results_deconvolve", "coverage.cell_type_deconvoluted.qnorm.csv"), index_col=0, header=range(4))
    X = deconv_norm.loc[:, ~deconv_norm.columns.get_level_values("patient_id").isin(["KI"])]
    X = X.loc[:, X.columns.get_level_values("cell_type").isin([cell_type])]
    X = X.loc[:, ~X.columns.get_level_values("timepoint").isin(["240d"])]

    return X


def parse_range(range_string, delim="___"):
    """
    """
    return tuple([int(x) for x in range_string.split(delim)])


def fit_gaussian_process(matrix):
    """
    """
    def fit_optimize(index, matrix):
        y = np.log2(matrix.loc[index, :].values[:, None])
        x = np.log2(1 + matrix.columns.get_level_values('timepoint').str.replace("d", "").astype(int).values.reshape((y.shape[0], 1)))

        kernel = GPy.kern.RBF(1) + GPy.kern.Bias(1)
        white_kernel = GPy.kern.Bias(1)

        m = GPy.models.gp_regression.GPRegression(x, y, kernel)  # GPHeteroscedasticRegression
        try:
            m.optimize()
        except RuntimeWarning:
            return [(np.nan, np.nan)]
        w_m = GPy.models.gp_regression.GPRegression(x, y, white_kernel)  # GPHeteroscedasticRegression
        try:
            w_m.optimize()
        except RuntimeWarning:
            return [(np.nan, np.nan)]
        return m.log_likelihood(), w_m.log_likelihood()

    ll = [fit_optimize(i, matrix) for i in matrix.index]

    return pd.DataFrame(ll, index=matrix.index, columns=['RBF', "White"])


def main():
    """
    Program's main entry point.
    """
    # Parse command-line arguments.
    args = parse_arguments()

    # Do it.
    matrix = read_matrix(args.cell_type)

    index_range = parse_range(args.data_range, delim=args.range_delim)

    fits = fit_gaussian_process(matrix.iloc[range(*index_range), :])

    fits.to_csv(os.path.join(args.output_dir, args.output_prefix + args.data_range + ".csv"))


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
