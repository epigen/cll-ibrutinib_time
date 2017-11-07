#!/usr/bin/env python

"""
Fit gaussian processes to a matrix of chromatin accessibility
with samples collected at different timepoints.
"""

import argparse
import os
import random
import string
import sys

import numpy as np
from sklearn.preprocessing import LabelEncoder


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
    matrix = read_matrix(args.matrix_file, args.cell_type,
                         args.matrix_header_range)
    # Fit GPs for requested range
    index_range = parse_range(args.data_range, delim=args.range_delim)
    fits = fit_gaussian_process(
        matrix.iloc[range(*index_range), :], library=args.library)
    # Save fits
    fits.to_csv(os.path.join(args.output_dir,
                             args.output_prefix + ".csv"), index=True)

    print("Done.")


def parse_arguments():
    """
    Argument Parsing.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-range", dest="data_range", type=str)
    parser.add_argument("--matrix-file", dest="matrix_file", type=str)
    parser.add_argument("--range-delim", dest="range_delim",
                        default="-", type=str)
    parser.add_argument("--matrix-header-range",
                        dest="matrix_header_range", default=9, type=int)
    parser.add_argument("--library", dest="library",
                        default="GPy", choices=['GPy', 'sklearn'], type=str)
    parser.add_argument("--cell-type", dest="cell_type",
                        default="CLL", type=str)
    parser.add_argument("--output-prefix", dest="output_prefix",
                        default="mohgp_fit_job", type=str)
    parser.add_argument("--output-dir", dest="output_dir",
                        default="/scratch/users/arendeiro/cll-time_course/gp_fit_job/output", type=str)

    args = parser.parse_args()

    return args


def read_matrix(matrix_file, cell_type="CLL", header_range=9):
    """
    Read matrix of (n_features, n_samples) and return samples from respective cell type
    (sample annotation as column metadata) with length `header_range`.
    """
    print("Reading matrix {} for cell type {}.".format(matrix_file, cell_type))
    matrix = pd.read_csv(matrix_file, index_col=0, header=range(header_range))
    X = matrix.loc[:, matrix.columns.get_level_values(
        "cell_type").isin([cell_type])]
    X = X.loc[:, ~X.columns.get_level_values("timepoint").isin(["240d"])]

    print("Finished reading matrix with {} features and {} samples of '{}' cell type.".format(
        X.shape[0], X.shape[1], cell_type))

    return X


def parse_range(range_string, delim="-"):
    """
    Given a string representing a numeric range and a delimiter,
    return the start and end of range.
    """
    return tuple([int(x) for x in range_string.split(delim)])


def fit_gaussian_process(matrix, library="GPy"):
    """
    Given a matrix of (n_features, n_samples), fit Gaussian Process with
    'timepoint' (from sample/column metadata) as x values and measurements as y
    for each feature.
    Will calculate log likelihoods and further statistics for two models.
    """
    def gpy_fit_optimize(index, matrix):
        import GPy
        import scipy

        y = matrix.loc[index, :].values[:, None]
        x = np.log2(1 +
                    matrix
                    .columns
                    .get_level_values('timepoint')
                    .str.replace("d", "")
                    .astype(int).values.reshape((y.shape[0], 1)))

        kernel = GPy.kern.RBF(input_dim=1) + GPy.kern.Bias(input_dim=1)
        white_kernel = GPy.kern.Bias(input_dim=1)

        m = GPy.models.GPRegression(x, y, kernel)
        try:
            m.optimize()
        except RuntimeWarning:
            return [np.nan] * 11
        w_m = GPy.models.GPRegression(x, y, white_kernel)
        try:
            w_m.optimize()
        except RuntimeWarning:
            return [np.nan] * 11

        # D statistic
        d = 2 * (m.log_likelihood() - w_m.log_likelihood())

        # p-value
        # the RBF + Bias kernel has 4 parameters and the Bias 2, so the chisquare has 2 degrees of freedom is 2
        p = scipy.stats.chi2.sf(d, df=2)

        # Let's calculate the STD of the posterior mean
        # because we have several y inputs for each x value
        # the posterior mean values retrieved will also be duplicated
        # let's make sure our STD is computed on the unique values only
        mean_posterior_std = (
            pd.DataFrame([x.squeeze(), m.posterior.mean.squeeze()],
                         index=['x', 'posterior'])
            .T
            .groupby('x')['posterior']
            .apply(lambda i: i.unique()[0])
            .std())

        return [d, p, mean_posterior_std, m.log_likelihood()] + m.param_array.tolist() + [w_m.log_likelihood()] + w_m.param_array.tolist()

    def sklearn_fit_optimize(index, matrix):
        raise NotImplementedError

        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

        y = matrix.loc[index, :].values
        x = np.log2(1 +
                    matrix
                    .columns
                    .get_level_values('timepoint')
                    .str.replace("d", "")
                    .astype(int).values.reshape((y.shape[0], 1)))
        kernel = 1 * (
            RBF(length_scale=1.0, length_scale_bounds=(1e-05, 100000.0)))  # +
        # WhiteKernel(noise_level=1.0, noise_level_bounds=(1e-05, 100000.0)))# +
        # ConstantKernel(constant_value=1.0, constant_value_bounds=(1e-05, 100000.0)))
        white_kernel = 1 * (  # WhiteKernel(noise_level=1.0, noise_level_bounds=(1e-05, 100000.0)) )# +
            ConstantKernel(constant_value=0., constant_value_bounds=(1e-05, 100000.0)))
        # white_kernel = ConstantKernel(constant_value=1.0, constant_value_bounds=(1e-05, 100000.0))

        m = GaussianProcessRegressor(kernel, normalize_y=False)
        try:
            m.fit(x, y)
        except:
            return np.nan, np.nan
        w_m = GaussianProcessRegressor(white_kernel, normalize_y=False)
        try:
            w_m.fit(x, y)
        except:
            return np.nan, np.nan
        return m.log_marginal_likelihood(), w_m.log_marginal_likelihood()

    def gpy_hierarchical_fit_optimize(index, matrix):
        import GPy
        import scipy

        y = matrix.loc[index, :].values[:, None]
        x = np.log2(1 +
                    matrix
                    .columns
                    .get_level_values('timepoint')
                    .str.replace("d", "")
                    .astype(int).values.reshape((y.shape[0], 1)))

        # construct a hierarchical GPy kernel.
        kern_upper = GPy.kern.RBF(input_dim=1, variance=1.0, name='upper')
        kern_lower = GPy.kern.RBF(input_dim=1, variance=0.1, name='lower')
        kernel = GPy.kern.Hierarchical(kernels=[kern_upper, kern_lower])

        # construct a 'grid' on which to examine samples
        # the first column contains time (T). the second columns contains
        # zeros, ones and twos to represent three replications.
        from sklearn.preprocessing import LabelEncoder
        t = np.log2(
            1 + matrix.columns.get_level_values("timepoint").str.replace("d", "").astype(int))
        patient_encoder = LabelEncoder()
        patient_encoder.fit(matrix.columns.get_level_values("patient_id"))
        p = patient_encoder.transform(
            matrix.columns.get_level_values("patient_id"))
        X = pd.DataFrame(np.concatenate(
            [t.values.reshape(-1, 1), p.reshape(-1, 1)], 1), columns=['T', 'patient_id'])
        X['y'] = y
        X = X.sort_values('T')

        m = GPy.models.GPRegression(
            X[['T', 'patient_id']].values,
            X[['y']].values,
            kernel)
        m.likelihood.variance = 0.01
        m.optimize('bfgs', messages=1)

        fig, axis = plt.subplots(1, len(patient_encoder.classes_), figsize=(
            20, 2), sharex=True, sharey=True)
        # axis[0].scatter(t, y)
        Xplot = X[['T']].values
        mu, var = m.predict(Xplot, kern=kern_upper)
        GPy.plotting.matplot_dep.base_plots.gpplot(
            Xplot, mu, mu - 2 * np.sqrt(var), mu + 2 * np.sqrt(var), ax=axis[0], edgecol='r', fillcol='r')
        # mu, var = m.predict(Xplot, kern=kern_lower)
        # GPy.plotting.matplot_dep.base_plots.gpplot(Xplot, mu, mu - 2 * np.sqrt(var), mu + 2 * np.sqrt(var), ax=axis[0], edgecol='b', fillcol='b')
        axis[0].set_title('Underlying\nfunction $f_{nr}(t)$')
        axis[0].set_ylabel('Chromatin accessibility')

        # plot each of the functions f_{nr}(t)
        for patient in range(1, len(patient_encoder.classes_)):
            m.plot(fixed_inputs=[(1, patient)], ax=axis[patient], which_data_rows=(
                X.loc[:, "patient_id"] == patient).values, legend=None)
            axis[patient].set_title("Patient {}".format(
                patient_encoder.inverse_transform(patient)))
            axis[patient].plot(Xplot, mu, 'r--', linewidth=1)
        for ax in axis:
            ax.set_xlabel('Time (log2)')
        fig.savefig(os.path.join("results_deconvolve", ".".join(
            [output_prefix, matrix_name, library, "hierarchy.top_variable.example.svg"])), dpi=300, bbox_inches="tight")

        white_kernel = GPy.kern.Bias(input_dim=1)
        try:
            m.optimize()
        except RuntimeWarning:
            return [np.nan] * 11
        w_m = GPy.models.GPRegression(x, y, white_kernel)
        try:
            w_m.optimize()
        except RuntimeWarning:
            return [np.nan] * 11

        # D statistic
        d = 2 * (m.log_likelihood() - w_m.log_likelihood())

        # p-value
        # the RBF + Bias kernel has 4 parameters and the Bias 2, so the chisquare has 2 degrees of freedom is 2
        p = scipy.stats.chi2.sf(d, df=2)

        # Let's calculate the STD of the posterior mean
        # because we have several y inputs for each x value
        # the posterior mean values retrieved will also be duplicated
        # let's make sure our STD is computed on the unique values only
        mean_posterior_std = (
            pd.DataFrame([x.squeeze(), m.posterior.mean.squeeze()],
                         index=['x', 'posterior'])
            .T
            .groupby('x')['posterior']
            .apply(lambda i: i.unique()[0])
            .std())

        return [d, p, mean_posterior_std, m.log_likelihood()] + m.param_array.tolist() + [w_m.log_likelihood()] + w_m.param_array.tolist()

    print("Fitting with library {}.".format(library))

    if library == "GPy":
        ll = [gpy_fit_optimize(i, matrix) for i in matrix.index]
        print("Finished.")
        return pd.DataFrame(ll, index=matrix.index, columns=["D", "p_value", "mean_posterior_std"] +
                            ['RBF'] + ['sum.rbf.variance', 'sum.rbf.lengthscale', 'sum.bias.variance', 'rbf.Gaussian_noise.variance'] +
                            ['White'] + ['bias.variance', 'bias.Gaussian_noise.variance'])
    elif library == "sklearn":
        ll = [sklearn_fit_optimize(i, matrix) for i in matrix.index]
        print("Finished.")
        return pd.DataFrame(ll, index=matrix.index, columns=['RBF'] + ['White'])

    elif library == "GPy_hierarchical":
        ll = [gpy_hierarchical_fit_optimize(i, matrix) for i in matrix.index]
        print("Finished.")
        return pd.DataFrame(ll, index=matrix.index, columns=["D", "p_value", "mean_posterior_std"] +
                            ['RBF'] + ['sum.rbf.variance', 'sum.rbf.lengthscale', 'sum.bias.variance', 'rbf.Gaussian_noise.variance'] +
                            ['White'] + ['bias.variance', 'bias.Gaussian_noise.variance'])


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
