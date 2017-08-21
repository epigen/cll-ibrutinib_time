
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
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.manifold import SpectralEmbedding
from sklearn.metrics import (explained_variance_score, mean_absolute_error,
                             mean_squared_error, r2_score)
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.preprocessing import LabelEncoder
from statsmodels.stats.multitest import multipletests

from mapalign.embed import DiffusionMapEmbedding
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.general import (collect_differential_enrichment,
                                 differential_enrichment,
                                 normalize_quantiles_r,
                                 plot_differential_enrichment,
                                 subtract_principal_component)

sns.set_style("white")

# random seed
SEED = int("".join(
    LabelEncoder()
    .fit(list(string.ascii_uppercase))
    .transform(list("BOCKLAB")).astype(str)))
random.seed(SEED)
np.random.seed(SEED)


def z_score(df, axis=1):
    from scipy.stats import zscore
    return pd.DataFrame(zscore(df, axis=axis), index=df.index, columns=df.columns)


def deseq_analysis(
        count_matrix, experiment_matrix, formula,
        output_dir, output_prefix,
        overwrite=True, alpha=0.05, independent_filtering=False):
    """
    Perform differential comparisons with DESeq2.
    """
    import pandas as pd
    from rpy2.robjects import numpy2ri, pandas2ri
    import rpy2.robjects as robjects
    numpy2ri.activate()
    pandas2ri.activate()

    def r2pandas_df(r_df):
        import numpy as np
        df = pd.DataFrame(np.asarray(r_df)).T
        df.columns = [str(x) for x in r_df.colnames]
        df.index = [str(x) for x in r_df.rownames]
        return df

    robjects.r('require("DESeq2")')
    _as_formula = robjects.r('as.formula')
    _DESeqDataSetFromMatrix = robjects.r('DESeqDataSetFromMatrix')
    _DESeq = robjects.r('DESeq')
    _results = robjects.r('results')
    _as_data_frame = robjects.r('as.data.frame')
    _resultsNames = robjects.r('resultsNames')

    # order experiment and count matrices in same way
    if experiment_matrix.index.name != "sample_name":
        try:
            experiment_matrix = experiment_matrix.set_index("sample_name")
        except KeyError:
            pass
    experiment_matrix = experiment_matrix.loc[count_matrix.columns, :]

    # save the matrices just in case
    count_matrix.to_csv(os.path.join(output_dir, output_prefix + ".count_matrix.tsv"), sep="\t")
    experiment_matrix.to_csv(os.path.join(output_dir, output_prefix + ".experiment_matrix.tsv"), sep="\t")

    # Run DESeq analysis
    dds = _DESeqDataSetFromMatrix(
        countData=count_matrix.astype(int),
        colData=experiment_matrix,
        design=_as_formula(formula))
    dds = _DESeq(dds, parallel=True)
    # _save(dds, file=os.path.join(output_dir, output_prefix + ".deseq_dds_object.Rdata"))

    comps = [str(x) for x in _resultsNames(dds)][1:]
    results = pd.DataFrame()
    for comp in comps:
        out_file = os.path.join(output_dir, output_prefix + ".deseq_result.{}.csv".format(comp))
        if not overwrite and os.path.exists(out_file):
            continue
        print("Doing comparison '{}'".format(comp))

        res = _as_data_frame(
            _results(dds, contrast=[comp], alpha=alpha, independentFiltering=independent_filtering, parallel=True))

        # convert to pandas dataframe
        res2 = r2pandas_df(res)
        res2.loc[:, "comparison_name"] = comp

        # save
        res2.to_csv(out_file)
        # append
        results = results.append(res2.reset_index(), ignore_index=True)

    # save all
    results.to_csv(os.path.join(output_dir, output_prefix + ".deseq_result.all_comparisons.csv"), index=False)


# Start project and analysis objects
prj = Project("metadata/project_config.yaml")
prj.samples = [sample for sample in prj.samples if sample.library == "ATAC-seq"]

# add samples from subprojects
for subproject in ["cll-chromatin", "stanford_atacseq"]:
    prj.samples += Project("metadata/project_config.yaml", subproject=subproject).samples

for s in prj.samples:
    s.filtered = os.path.join(s.paths.sample_root, "mapped", s.name + ".trimmed.bowtie2.filtered.bam")
    s.peaks = os.path.join(s.paths.sample_root, "peaks", s.name + "_peaks.narrowPeak")

analysis = ATACSeqAnalysis(name="cll-time_course", prj=prj, samples=prj.samples, results_dir="results_deconvolve")

# Project's attributes
sample_attributes = ["sample_name", "patient_id", "timepoint", "cell_type", "compartment", "response", "cell_number", "batch"]
numerical_attributes = ["CD38_cells_percentage", "cll_cells_%", "cell_number"]
cell_types = ["Bulk", "Bcell", "CLL", "CD4", "CD8", "NK", "Mono"]


# GET CONSENSUS PEAK SET, ANNOTATE IT, PLOT FEATURES
# Get consensus peak set from all samples
# analysis.get_consensus_sites([s for s in analysis.samples if s.cell_type == "Bulk"])

# GET CHROMATIN OPENNESS MEASUREMENTS
# Get coverage values for each peak in each sample of ATAC-seq
analysis.sites = pybedtools.BedTool(os.path.join(analysis.results_dir, analysis.name + "_peak_set.bed"))
analysis.measure_coverage(
    [s for s in analysis.samples if not hasattr(s, "pass_counts")] +
    [s for s in analysis.samples if s.cell_type == "Bulk"])


# Normalize Stanford samples + CLL merged
to_norm = analysis.coverage.loc[:, ~analysis.coverage.columns.str.contains("Bulk")].drop(['chrom', 'start', 'end'], axis=1)
ref_norm = pd.DataFrame(
    normalize_quantiles_r(
        to_norm
        .values),
        index=to_norm.index,
        columns=to_norm.columns).T
ref_norm.to_csv(os.path.join(analysis.results_dir, "coverage.cell_type_reference.qnorm.csv"))


# Normalize Bulk samples independently
to_norm = analysis.coverage.loc[:, analysis.coverage.columns.str.contains("Bulk")]
to_deconv = pd.DataFrame(
    normalize_quantiles_r(
        to_norm
        .values),
        index=to_norm.index,
        columns=to_norm.columns)

cols = pd.Series(to_deconv.columns.str.split('_')).apply(pd.Series)
to_deconv.columns = pd.MultiIndex.from_arrays(cols.T.values, names=["library", "patient_id", "timepoint", "cell_type", "replicate"])


to_deconv = to_deconv.loc[:, to_deconv.columns.get_level_values("replicate").isnull()]

to_deconv.to_csv(os.path.join(analysis.results_dir, "coverage.Bulk.qnorm.csv"))


# Get mean of cell type and fractional contributions
ref_norm['cell_type'] = ["CLL"] + map(lambda x: x[1], ref_norm.index[1:].str.split("_"))
ref_norm = ref_norm.groupby(['cell_type']).mean().T

contrib = (ref_norm.T / ref_norm.sum(axis=1)).T

# Let's make synthetic samples from Bulk

# get the amount of each cell type in each bulk
facs = pd.read_csv(os.path.join("metadata", "facs_quantification.annotated.csv"))

deconv = pd.DataFrame()

patient_ids = list()
timepoints = list()
cell_types = list()
i = 0
for patient_id in to_deconv.columns.levels[1]:
    for timepoint in to_deconv.columns.levels[2]:
        for cell_type in contrib.columns:
            atac = to_deconv.loc[:, (to_deconv.columns.get_level_values("patient_id") == patient_id) & (to_deconv.columns.get_level_values("timepoint") == timepoint)].squeeze()
            if len(atac.shape) != 1:
                continue
            c = contrib.loc[:, cell_type]
            f = facs.loc[(facs["patient_id"] == patient_id) & (facs["timepoint"] == timepoint), cell_type].squeeze()

            deconv.loc[:, i] = atac * c * f
            i += 1
            patient_ids.append(patient_id)
            timepoints.append(timepoint)
            cell_types.append(cell_type)
names = pd.DataFrame([patient_ids, timepoints, cell_types]).apply(lambda x: "_".join(x))
deconv.columns = pd.MultiIndex.from_arrays([names, patient_ids, timepoints, cell_types], names=["sample_name", "patient_id", "timepoint", "cell_type"])

deconv.to_csv(os.path.join(analysis.results_dir, "coverage.cell_type_deconvoluted.raw.csv"))

# quantile normalize
to_norm = deconv.loc[:, ~deconv.isnull().all()]
to_norm = to_norm.loc[:, ~(deconv == 0).all()]
deconv_norm = pd.DataFrame(
    normalize_quantiles_r(
        to_norm
        .values),
        index=to_norm.index,
        columns=to_norm.columns)

deconv_norm.to_csv(os.path.join(analysis.results_dir, "coverage.cell_type_deconvoluted.qnorm.csv"))

deconv_norm = pd.read_csv(os.path.join(analysis.results_dir, "coverage.cell_type_deconvoluted.qnorm.csv"), index_col=0, header=range(4))

# Do some unsupervised analysis
analysis.samples = deconv_norm.columns.to_frame().apply(Sample, axis=1).tolist()
analysis.deconvolved_data = deconv_norm
analysis.unsupervised(quant_matrix="deconvolved_data", samples=None, attributes_to_plot=["patient_id", "timepoint", "cell_type"], plot_prefix="deconvolved_data")


c = np.log2(1 + deconv_norm).corr()
g = sns.clustermap(c, xticklabels=False, figsize=(0.12 * c.shape[0], 0.12 * c.shape[0]), cbar_kws={"label": "Pearson correlation"})
g.ax_row_dendrogram.set_rasterized(True)
g.ax_col_dendrogram.set_rasterized(True)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize='xx-small', rasterized=True)
g.savefig(os.path.join(analysis.results_dir, "cll-time_course.deconvolved_data.qnorm.log2.corr.clustermap.svg"), bbox_inches="tight", dpi=300)


# Check patient characteristics are still preserved
corr = deconv_norm.corr()
np.fill_diagonal(corr.values, np.nan)
corr.index = corr.index.get_level_values("sample_name")
corr.columns = corr.columns.get_level_values("sample_name")

# for each cell type, each timepoint, check if within-patient correlation is higher than between
c = corr.stack(level=0)
c.index.names = ['s1', 's2']
c = c.reset_index()

res = pd.DataFrame()
for patient_id in deconv_norm.columns.levels[1]:
    for cell_type in deconv_norm.columns.levels[3]:
        print(patient_id, cell_type)
        for timepoint in deconv_norm.columns.levels[2]:
            current = c.loc[
                (c["s1"].str.contains(patient_id)) &
                (c["s1"].str.contains(cell_type)) &
                (c["s1"].str.contains(timepoint))
            ]
            if current.shape[0] == 0:
                continue
            # within patient
            within = current.loc[
                (current["s2"].str.contains(patient_id)) &
                (current["s2"].str.contains(cell_type))
            ].mean().squeeze()
            # between patients
            between = current.loc[
                (~current["s2"].str.contains(patient_id)) &
                (current["s2"].str.contains(cell_type))
            ].mean().squeeze()
            res = res.append(
                pd.Series([patient_id, cell_type, timepoint, within, between]),
                ignore_index=True
            )
res.columns = ["patient_id", "cell_type", "timepoint", "within", "between"]
res.to_csv(os.path.join(analysis.results_dir, "coverage.cell_type_deconvoluted.qnorm.correlation_comparison.csv"), index=False)


fig, axis = plt.subplots(1, 3,
    gridspec_kw={'width_ratios':[3, 1, 1]},
    figsize=(3 * 6, 1 * 6),
    tight_layout=True)
res.groupby(['cell_type', 'timepoint']).mean().plot(kind='bar', ax=axis[0])
res.groupby(['timepoint']).mean().plot(kind='bar', ax=axis[1])
res.groupby(['cell_type']).mean().plot(kind='bar', ax=axis[2])
for ax in axis:
    ax.set_ylim((0.5, 1))
sns.despine(fig)
fig.savefig(
    os.path.join(analysis.results_dir, "coverage.cell_type_deconvoluted.qnorm.correlation_comparison.barplots.svg"),
    bbox_inches="tight")


# Now compare with sorted samples
deconv_samples = analysis.samples
pickle.dump(deconv_samples, open(os.path.join(analysis.results_dir, "deconv_samples.pickle"), "wb"))

analysis.samples = prj.samples
analysis.measure_coverage(
    [s for s in analysis.samples if hasattr(s, "pass_counts")])


# Normalize Stanford samples + CLL merged
to_norm = analysis.coverage.drop(['chrom', 'start', 'end'], axis=1)
sorted_norm = pd.DataFrame(
    normalize_quantiles_r(
        to_norm
        .values),
        index=to_norm.index,
        columns=to_norm.columns)
sorted_norm.to_csv(os.path.join(analysis.results_dir, "coverage.sorted_samples.qnorm.csv"))

sorted_norm = pd.read_csv(os.path.join(analysis.results_dir, "coverage.sorted_samples.qnorm.csv"), index_col=0)

cols = pd.DataFrame(sorted_norm.columns.str.split("_").tolist(), columns=["library", "patient_id", "timepoint", "cell_type", "replicate"]).drop("library", axis=1)
cols['cell_type'] = cols['cell_type'].str.replace("CD4", "CD4Tcell").str.replace("CD8", 'CD8Tcell')
cols['sample_name'] = sorted_norm.columns
sorted_norm.columns = pd.MultiIndex.from_arrays(cols.values.T, names=["patient_id", "timepoint", "cell_type", "replicate", "sample_name"])

# now, for each cell type, each timepoint, check if within-patient correlation is higher than between
# across data types (sorted vs deconvolved)
sorted_norm = sorted_norm.loc[:, sorted_norm.columns.get_level_values("replicate").isnull()]

metrics = [pearsonr, spearmanr, mean_squared_error, mean_absolute_error, explained_variance_score, r2_score]
res = pd.DataFrame()
for patient_id in sorted_norm.columns.levels[0]:
    for cell_type in sorted_norm.columns.levels[2]:
        if cell_type == "Bulk":
            continue
        print(patient_id, cell_type)
        for timepoint in sorted_norm.columns.levels[1]:
            sorted_current = sorted_norm.loc[:,
                (sorted_norm.columns.get_level_values("patient_id") == patient_id) &
                (sorted_norm.columns.get_level_values("cell_type") == cell_type) &
                (sorted_norm.columns.get_level_values("timepoint") == timepoint)
            ]
            # within patient
            within = deconv_norm.loc[:,
                (deconv_norm.columns.get_level_values("patient_id") == patient_id) &
                (deconv_norm.columns.get_level_values("cell_type") == cell_type) &
                (deconv_norm.columns.get_level_values("timepoint") == timepoint)
            ]
            # between patients
            between = deconv_norm.loc[:,
                (deconv_norm.columns.get_level_values("patient_id") != patient_id) &
                (deconv_norm.columns.get_level_values("cell_type") == cell_type) &
                (deconv_norm.columns.get_level_values("timepoint") == timepoint)
            ]
            if sorted_current.shape[1] == 0:
                continue
            if within.shape[1] == 0:
                continue

            # calculate metrics
            for metric in metrics:
                w = metric(sorted_current.squeeze(), within.squeeze())
                
                if between.shape[1] > 0:
                    bs = list()
                    for other in between.columns:
                        m = metric(sorted_current.squeeze(), between[other].squeeze())
                        if hasattr(m, "correlation"):
                            m = m.correlation
                        if type(m) is tuple:
                            m = m[0]
                        bs.append(m)
                    bmean = np.mean(bs)
                    bmax = np.max(bs)
                else:
                    bmean = np.nan
                    bmax = np.nan

                if hasattr(w, "correlation"):
                    w = w.correlation
                if type(w) is tuple:
                    w = w[0]

                res = res.append(
                    pd.Series([patient_id, cell_type, timepoint, metric.__name__, w, bmean, bmax]),
                    ignore_index=True
                )
res.columns = ["patient_id", "cell_type", "timepoint", "metric", "within", "between_mean", "between_max"]
res.to_csv(os.path.join(analysis.results_dir, "coverage.sorted_vs_deconvoluted.correlation_comparison.csv"), index=False)


fig, axis = plt.subplots(6, 4,
    gridspec_kw={'width_ratios':[3, 1, 1, 1]},
    figsize=(4 * 6, 6 * 6), sharey='row',
    tight_layout=True)
for i, metric in enumerate(metrics):
    axis[i, 0].set_title(metric.__name__)
    res.dropna().drop('between_max', axis=1)[res['metric'] == metric.__name__].groupby(['cell_type', 'timepoint']).mean().plot(kind='bar', ax=axis[i, 0])
    res.dropna().drop('between_max', axis=1)[res['metric'] == metric.__name__].groupby(['timepoint']).mean().plot(kind='bar', ax=axis[i, 1])
    res.dropna().drop('between_max', axis=1)[res['metric'] == metric.__name__].groupby(['cell_type']).mean().plot(kind='bar', ax=axis[i, 2])
    res.dropna().drop('between_max', axis=1)[res['metric'] == metric.__name__].groupby(['patient_id']).mean().plot(kind='bar', ax=axis[i, 3])
sns.despine(fig)
fig.savefig(
    os.path.join(analysis.results_dir, "coverage.sorted_vs_deconvoluted.correlation_comparison.metric_mean.barplots.svg"),
    bbox_inches="tight")



fig, axis = plt.subplots(6, 4,
    gridspec_kw={'width_ratios':[3, 1, 1, 1]},
    figsize=(4 * 6, 6 * 6), sharey='row',
    tight_layout=True)
for i, metric in enumerate(metrics):
    axis[i, 0].set_title(metric.__name__)
    res.dropna().drop('between_mean', axis=1)[res['metric'] == metric.__name__].groupby(['cell_type', 'timepoint']).mean().plot(kind='bar', ax=axis[i, 0])
    res.dropna().drop('between_mean', axis=1)[res['metric'] == metric.__name__].groupby(['timepoint']).mean().plot(kind='bar', ax=axis[i, 1])
    res.dropna().drop('between_mean', axis=1)[res['metric'] == metric.__name__].groupby(['cell_type']).mean().plot(kind='bar', ax=axis[i, 2])
    res.dropna().drop('between_mean', axis=1)[res['metric'] == metric.__name__].groupby(['patient_id']).mean().plot(kind='bar', ax=axis[i, 3])
sns.despine(fig)
fig.savefig(
    os.path.join(analysis.results_dir, "coverage.sorted_vs_deconvoluted.correlation_comparison.metric_max.barplots.svg"),
    bbox_inches="tight")



# Get differential peaks across patients between timpoints within cell types
count_matrix = deconv.loc[:, (~deconv.isnull().all()) & (~(deconv == 0).all())].multiply(100).astype(int)
experiment_matrix = pd.DataFrame([s.as_series() for s in analysis.samples])[['sample_name', 'patient_id', 'cell_type', 'timepoint']]
experiment_matrix = experiment_matrix[experiment_matrix['sample_name'].isin(count_matrix.columns.get_level_values("sample_name"))]

for cell_type in experiment_matrix["cell_type"].drop_duplicates():
    c = count_matrix.loc[:, count_matrix.columns.get_level_values("cell_type") == cell_type]
    c.columns = c.columns.get_level_values("sample_name")
    deseq_analysis(
        count_matrix=c,
        experiment_matrix=experiment_matrix[experiment_matrix["cell_type"] == cell_type],
        formula="~timepoint",
        output_dir=os.path.join(analysis.results_dir, "deseq"),
        output_prefix=cell_type,
        overwrite=True,
        alpha=0.05,
        independent_filtering=False
    )

sorted_norm = pd.read_csv(os.path.join(analysis.results_dir, "coverage.sorted_samples.qnorm.csv"), index_col=0)

quantity = "accessibility"
var_name = "regions"

for cell_type in experiment_matrix["cell_type"].drop_duplicates():
    df = pd.read_csv(os.path.join(analysis.results_dir, "deseq", cell_type + ".deseq_result.all_comparisons.csv"), index_col=0)
    print(cell_type, df[(df['padj'] < 0.05)].shape)

    all_diff =  df[(df['padj'] < 0.05) & (df['comparison_name'].isin('timepoint' + pd.Series(['0d', '30d', '120d'])))].index.drop_duplicates()

    matrix_deconv = deconv_norm.loc[
        all_diff,
        (deconv_norm.columns.get_level_values("patient_id") != "KI") &
        (deconv_norm.columns.get_level_values("cell_type") == cell_type) &
        (deconv_norm.columns.get_level_values("timepoint") != "3d")]

    group_matrix = matrix_deconv.T.reset_index().groupby("timepoint").mean().T.loc[all_diff]

    if type(matrix_deconv.columns) is pd.core.indexes.multi.MultiIndex:
        matrix_deconv.columns = matrix_deconv.columns.get_level_values("sample_name")

    # Deconv data
    figsize = (max(5, 0.12 * matrix_deconv.shape[1]), 5)
    g = sns.clustermap(matrix_deconv,
        yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
        xticklabels=True,
        metric="correlation", rasterized=True, figsize=figsize)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.fig.savefig(os.path.join(analysis.results_dir, cell_type + ".diff_{}.samples.deconv_samples.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=300)

    # Deconv data by group
    figsize = (max(5, 0.12 * group_matrix.shape[1]), 5)
    g = sns.clustermap(group_matrix,
        yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
        xticklabels=True,
        metric="correlation", rasterized=True, figsize=figsize)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.fig.savefig(os.path.join(analysis.results_dir, cell_type + ".diff_{}.groups.deconv_samples.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=300)

    # Now sorted data
    matrix_sorted = sorted_norm.loc[
        all_diff,
        (~sorted_norm.columns.str.contains("KI")) &
        (sorted_norm.columns.str.contains(cell_type.replace("Tcell", "").replace("cell", ""))) &
        (~sorted_norm.columns.str.contains("3d"))]
    if type(matrix_sorted.columns) is pd.core.indexes.multi.MultiIndex:
        matrix_sorted.columns = matrix_sorted.columns.get_level_values("sample_name")
    figsize = (max(5, 0.12 * matrix_sorted.shape[1]), 5)

    g = sns.clustermap(matrix_sorted,
        yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
        xticklabels=True,
        metric="correlation", rasterized=True, figsize=figsize)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.fig.savefig(os.path.join(analysis.results_dir, cell_type + ".diff_{}.samples.sorted_samples.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=300)

    # Now Bulk data
    matrix_bulk = sorted_norm.loc[
        all_diff,
        (~sorted_norm.columns.str.contains("KI")) &
        (sorted_norm.columns.str.contains("Bulk")) &
        (~sorted_norm.columns.str.contains("3d"))]
    if type(matrix_bulk.columns) is pd.core.indexes.multi.MultiIndex:
        matrix_bulk.columns = matrix_bulk.columns.get_level_values("sample_name")
    figsize = (max(5, 0.12 * matrix_bulk.shape[1]), 5)

    g = sns.clustermap(matrix_bulk,
        yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
        xticklabels=True,
        metric="correlation", rasterized=True, figsize=figsize)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.fig.savefig(os.path.join(analysis.results_dir, cell_type + ".diff_{}.samples.bulk_samples.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=300)

    # Now joint synthetic and sorted
    matrix_joint = z_score(matrix_sorted).join(z_score(matrix_deconv))
    figsize = (max(5, 0.12 * matrix_joint.shape[1]), 5)

    g = sns.clustermap(matrix_joint.sort_index(axis=1),
        yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
        col_cluster=True,
        xticklabels=True,
        metric="correlation", rasterized=True, figsize=figsize)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.fig.savefig(os.path.join(analysis.results_dir, cell_type + ".diff_{}.samples.joint_synthetic_sorted.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=300)


# Get enrichment of regions
all_diff = pd.DataFrame()
for cell_type in experiment_matrix["cell_type"].drop_duplicates():
    df = pd.read_csv(os.path.join(analysis.results_dir, "deseq", cell_type + ".deseq_result.all_comparisons.csv"), index_col=0)
    print(cell_type, df[(df['padj'] < 0.05)].shape)

    df2 =  df.loc[(df['padj'] < 0.05) & (df['comparison_name'].isin('timepoint' + pd.Series(['0d', '30d', '120d']))), :]

    df2.loc[:, 'comparison_name'] = cell_type + "_" + df2['comparison_name'].str.replace("timepoint", "")

    all_diff = all_diff.append(df2.reset_index(), ignore_index=True)
all_diff.to_csv(os.path.join(analysis.results_dir, "deseq", "all_cell_types" + ".deseq_result.all_comparisons.csv"), index=False)


all_diff = pd.read_csv(os.path.join(analysis.results_dir, "deseq", "all_cell_types" + ".deseq_result.all_comparisons.csv"))


# Get enrichments of regions found in the differential comparisons
output_dir = analysis.results_dir
output_prefix = "deseq"

differential_enrichment(
    analysis, differential=all_diff.set_index("index"),
    data_type="ATAC-seq", genome="hg19",
    directional=False,
    output_dir=output_dir, output_prefix=output_prefix)

collect_differential_enrichment(
    differential=all_diff,
    directional=False,
    data_type="ATAC-seq",
    permissive=True,
    output_dir=output_dir, output_prefix=output_prefix)

# Visualize enrichments
enrichment_table = pd.read_csv(os.path.join(analysis.results_dir, output_prefix + ".lola.csv"))
plot_differential_enrichment(
    enrichment_table,
    enrichment_type="lola",
    data_type="ATAC-seq", top_n=25,
    output_dir=output_dir, output_prefix=output_prefix)



# Annotate peaks
analysis.samples = prj.samples
analysis.get_peak_gene_annotation()
analysis.get_peak_genomic_location()
analysis.get_peak_chromatin_state(os.path.join("data", "external", "E032_15_coreMarks_mnemonics.bed"))
analysis.calculate_peak_support([s for s in analysis.samples if s.cell_type == "Bulk"], "summits")
analysis.annotate(quant_matrix="coverage", samples=[s for s in prj.samples if s.name in analysis.coverage.columns])

analysis.gene_annotation = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.gene_annotation.csv"))
analysis.region_annotation = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.region_annotation.csv"))
analysis.chrom_state_annotation = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.chromatin_state.csv"))
analysis.support = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.support.csv"))











# on a subset of regions (e.g. differential)
all_diff = pd.read_csv(os.path.join(analysis.results_dir, "deseq", "all_cell_types" + ".deseq_result.all_comparisons.csv"))
r = all_diff['index'].drop_duplicates()


d = pd.read_csv("results/cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.all.diff_timepoint.limma.csv")
r = d.loc[d['q_value'] < 0.05, "region"].drop_duplicates()



deconv_norm = pd.read_csv(os.path.join(analysis.results_dir, "coverage.cell_type_deconvoluted.qnorm.csv"), index_col=0, header=range(4))
# to_deconv = pd.read_csv(os.path.join(analysis.results_dir, "coverage.Bulk.qnorm.csv"), index_col=0, header=range(5))
# to_deconv = to_deconv.loc[:, ~to_deconv.columns.get_level_values("patient_id").isin(["KI", "FE", "KZ"])]
# X = deconv_norm.loc[r, deconv_norm.columns.get_level_values("cell_type") == "CLL"]
X = np.log2(deconv_norm.loc[r, deconv_norm.columns.get_level_values("cell_type") != "Q"])

X = X.loc[:, ~X.columns.get_level_values("patient_id").isin(["KI"])]

# group_mean = X.T.groupby(level=['timepoint']).mean()
t = 100
de = DiffusionMapEmbedding(alpha=0.5, diffusion_time=t, affinity='markov',
                           n_components=10).fit_transform(X.T)
ed = (de - de[0, :])
ed = np.sqrt(np.sum(ed * ed , axis=1))
ed /= max(ed)

rank = pd.Series(ed).rank().astype(int) - 1
g = sns.clustermap(X.dropna(), z_score=0, row_cluster=True, metric="correlation", figsize=(0.12 * X.shape[1], 8), rasterized=True)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, visible=False, fontsize="xx-small")
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
g.savefig("/home/arendeiro/m.png", dpi=300, bbox_inches="tight")


q = SpectralEmbedding().fit_transform(X.T)
plt.scatter(q[:, 0], q[:, 1], c=ed)
[plt.text(q[i, 0], q[i, 1], X.columns[i]) for i in range(len(X.columns))]






deconv_norm = pd.read_csv(os.path.join(analysis.results_dir, "coverage.cell_type_deconvoluted.qnorm.csv"), index_col=0, header=range(4))
X = deconv_norm.loc[:, ~deconv_norm.columns.get_level_values("patient_id").isin(["KI"])]
X = X.loc[:, X.columns.get_level_values("cell_type").isin(["CLL"])]
X = X.loc[:, ~X.columns.get_level_values("timepoint").isin(["240d"])]



pca = PCA()
pcs = pd.DataFrame(pca.fit_transform(np.log2(1 + X.T)), index=X.columns)

t = 100
de = DiffusionMapEmbedding(alpha=0.5, diffusion_time=t, affinity='markov',
                           n_components=2).fit_transform(pcs.iloc[:, :10])

times = pcs.index.get_level_values("timepoint").str.replace("d", "").astype(int)
plt.scatter(zscore(de[:, 0]), zscore(de[:, 1]), c=np.log10(1 + times), cmap="inferno")

plt.scatter(zscore(pcs.loc[:, 0]), zscore(pcs.loc[:, 1]), c=np.log10(1 + times), cmap="inferno")



# Gaussian Process based




white_kernel = 1.0 * WhiteKernel()
kernel = 1.0 * RBF() + WhiteKernel()
white_gp = GaussianProcessRegressor(kernel=white_kernel, alpha=0.0)
gp = GaussianProcessRegressor(kernel=kernel, alpha=0.0)

white_lls = list()
lls = list()
for region in X.index:
    print(region)
    y = X.loc[region, :]
    x = np.log2(1 + X.columns.get_level_values('timepoint').str.replace("d", "").astype(int).values.reshape((y.shape[0], 1)))
    y = y.values

    white_lls.append(white_gp.fit(x, y).log_marginal_likelihood())
    lls.append(gp.fit(x, y).log_marginal_likelihood())

res = pd.DataFrame([white_lls, lls], columns=X.index[:len(lls)], index=["constant", "varying"]).T


res["lldiff"] = res["varying"] - res["constant"]
res["2lldiff"] = 2 * res["varying"] - 2 * res["constant"]
plt.scatter(res["constant"], res["2lldiff"])


# Visualize variable regions
n_regions = 4
top_regions = res["2lldiff"].sort_values().tail(n_regions).index

fig, axis = plt.subplots(n_regions, 2, figsize=(2 * 4, n_regions * 4))
for i, region in enumerate(top_regions):
    y = X.loc[region, :].groupby(level="patient_id").apply(zscore).apply(pd.Series).stack()
    x = np.log2(1 + X.columns.get_level_values('timepoint').str.replace("d", "").astype(int).values.reshape((y.shape[0], 1)))
    y = y.values
    X_ = np.linspace(-2, max(x) + 2)

    white_gp = white_gp.fit(x, y)
    y_mean, y_cov = white_gp.predict(X_[:, np.newaxis], return_cov=True)  # , return_std=True)
    axis[i, 0].plot(X_, y_mean, 'k', lw=3, zorder=9)
    axis[i, 0].fill_between(X_, y_mean - np.sqrt(np.diag(y_cov)),
                    y_mean + np.sqrt(np.diag(y_cov)),
                    alpha=0.5, color='k')
    axis[i, 0].scatter(x[:, 0], y, c='r', s=50, zorder=10)
    gp = gp.fit(x, y)
    y_mean, y_cov = gp.predict(X_[:, np.newaxis], return_cov=True)
    axis[i, 1].plot(X_, y_mean, 'k', lw=3, zorder=9)
    axis[i, 1].fill_between(X_, y_mean - np.sqrt(np.diag(y_cov)),
                    y_mean + np.sqrt(np.diag(y_cov)),
                    alpha=0.5, color='k')
    axis[i, 1].scatter(x[:, 0], y, c='r', s=50, zorder=10)
fig.tight_layout()


# With GPy




def count_jobs_running(cmd="squeue", sep="\n"):
    """
    Count running jobs on a cluster by invoquing a command that lists the jobs.
    """
    import subprocess
    return subprocess.check_output(cmd).split(sep).__len__()


def submit_job_if_possible(cmd, total_job_lim=800, refresh_time=10, in_between_time=5):
    import time
    import os

    submit = count_jobs_running() < total_job_lim
    while not submit:
        time.sleep(refresh_time)
        submit = count_jobs_running() < total_job_lim
    os.system(cmd)
    time.sleep(in_between_time)

# deconv_norm = pd.read_csv(os.path.join(analysis.results_dir, "coverage.cell_type_deconvoluted.qnorm.csv"), index_col=0, header=range(4))
sorted_norm = pd.read_csv(os.path.join("results", "cll-time_course" + "_peaks.coverage.joint_qnorm.pca_fix.power.csv"), index_col=0, header=range(8))
deconv_norm = pd.read_csv(os.path.join("results_deconvolve", "coverage.cell_type_deconvoluted.qnorm.csv"), index_col=0, header=range(4))

chunks = 2000
total_job_lim = 800
refresh_time = 10
in_between_time = 0.01
output_prefix = "gp_fit_job"
output_dir = "/scratch/users/arendeiro/gp_fit_job"

library = "GPy"

for matrix_name, matrix in tqdm.tqdm([("sorted", sorted_norm), ("deconv", deconv_norm)], desc="matrix"):
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
    from statsmodels.stats.multitest import multipletests
    # TODO: check this holds the order
    fits['q_value'] = np.concatenate(fits.groupby("cell_type", sort=False)['p_value'].apply(lambda x: multipletests(x, method="fdr_bh")[1]))

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
    fits['mean_posterior_std_rank'] = fits['mean_posterior_std'].rank(ascending=False)
    fits['D_rank'] = fits['D'].rank(ascending=False)
    fits['RBF_rank'] = fits['RBF'].rank(ascending=False)
    fits['White_rank'] = fits['White'].rank(ascending=True)

    fits['rank_max'] = fits[["mean_posterior_std_rank", "D_rank", "RBF_rank", "White_rank"]].max(axis=1)
    fits['rank_min'] = fits[["mean_posterior_std_rank", "D_rank", "RBF_rank", "White_rank"]].min(axis=1)
    fits['rank_mean'] = fits[["mean_posterior_std_rank", "D_rank", "RBF_rank", "White_rank"]].mean(axis=1)

    n_top = 6
    e = fits[~fits['cell_type'].str.contains("NK")].sort_values("rank_mean").head(n_top)
    examples = e.index
    example_ct = e['cell_type']
    if matrix_name == "deconv":
        example_acc = np.log2(matrix.loc[examples])
    else:
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
            white_kernel = GPy.kern.Bias(input_dim=1)

            m = GPy.models.GPRegression(x.reshape(-1, 1), cur.values.reshape(-1, 1), kernel)
            m.optimize()
            m.plot_f(ax=axis[i, j], lower=2.5, upper=97.5, legend=None, plot_density=True, plot_data=False, color="red")

            # w_m = GPy.models.GPRegression(x.reshape(-1, 1), cur.values.reshape(-1, 1), white_kernel)
            # w_m.optimize()

    for i, ax in enumerate(axis[:, 0]):
        ax.set_ylabel(examples[i])
        # ax.set_title(example_acc.iloc[i]['cell_type'].squeeze())
    for i, ax in enumerate(axis[0, :]):
        ax.set_title(cell_types[i])
    fig.savefig(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "top_variable.scatter.all_samples.svg"])), dpi=300, bbox_inches="tight")


    # # Hierarchycal model
    # from sklearn.preprocessing import LabelEncoder
    # patient_encoder = LabelEncoder()
    # patient_encoder.fit(matrix.columns.get_level_values("patient_id"))

    # # construct a hierarchical GPy kernel. 
    # kern_upper = GPy.kern.RBF(input_dim=1, lengthscale=max(x) / 3., variance=1.0, name='upper')
    # kern_lower = GPy.kern.RBF(input_dim=1, lengthscale=max(x) / 3., variance=0.1, name='lower')
    # kernel = GPy.kern.Hierarchical(kernels=[kern_upper, kern_lower])

    # size = (len(examples), len(patient_encoder.classes_))
    # fig, axis = plt.subplots(*size, figsize=tuple(s * 2 for s in size[::-1]), sharex=True, sharey=True)
    # for i, region in enumerate(examples):
    #     print(region)
    #     samples = matrix.columns.get_level_values("cell_type") == cell_type
    #     cur = matrix.loc[
    #         region,
    #         samples
    #     ].drop_duplicates().squeeze()

    #     # construct a 'grid' on which to examine samples
    #     # the first column contains time (T). the second columns contains
    #     # zeros, ones and twos to represent three replications.
    #     t = np.log2(1 + cur.index.get_level_values("timepoint").str.replace("d", "").astype(int))
    #     p = patient_encoder.transform(cur.index.get_level_values("patient_id"))
    #     X = pd.DataFrame(np.concatenate([t.values.reshape(-1, 1), p.reshape(-1, 1)], 1), columns=['T', 'patient_id'])
    #     X['y'] = scipy.stats.zscore(cur.values)
    #     X = X.sort_values('T')

    #     m = GPy.models.GPRegression(
    #         X[['T', 'patient_id']].values,
    #         X[['y']].values,
    #         kernel)
    #     m.likelihood.variance = 0.01
    #     m.optimize('bfgs', messages=1)

    #     # axis[i, 0].scatter(t, y)
    #     Xplot = X[['T']].sort_values('T').values
    #     mu, var = m.predict(Xplot, kern=kern_upper)
    #     GPy.plotting.matplot_dep.base_plots.gpplot(Xplot, mu, mu - 2 * np.sqrt(var), mu + 2 * np.sqrt(var), ax=axis[i, 0], edgecol='r', fillcol='r')
    #     # mu, var = m.predict(Xplot, kern=kern_lower)
    #     # GPy.plotting.matplot_dep.base_plots.gpplot(Xplot, mu, mu - 2 * np.sqrt(var), mu + 2 * np.sqrt(var), ax=axis[i, 0], edgecol='b', fillcol='b')

    #     # plot each of the functions f_{nr}(t)
    #     for patient in range(1, len(patient_encoder.classes_)):
    #         m.plot(fixed_inputs=[(1, patient)], ax=axis[i, patient], which_data_rows=(X.loc[:, "patient_id"]==patient).values, legend=None)
    #         axis[i, patient].plot(Xplot, mu, 'r--', linewidth=1)

    # l = ['Underlying\nfunction $f_{nr}(t)$'] + ["Patient {}".format(p) for p in patient_encoder.classes_]
    # for i, ax in zip(l, axis[0, :]):
    #     ax.set_title(i)
    # for ax in axis[:, 0]:
    #     ax.set_ylabel('Chromatin accessibility')
    # for ax in axis[-1, :]:
    #     ax.set_xlabel('Time (log2)')
    # fig.savefig(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "hierarchy.top_variable.examples.svg"])), dpi=300, bbox_inches="tight")


# Visualize regions, cluster them and get enrichments
from scipy.cluster.hierarchy import fcluster

samples_to_exclude = [
    "PBGY_1d_NK",
    "240d_CD8", "280d_CD8",
    "240d_CD4", "280d_CD4",
    "PBGY_1d_Mono", "PBGY_8d_Mono", "PBGY_150d_Mono",
    "KZ_240d_Bulk",
    "FE_3d_Bulk",
    "VZS_2d_Bulk",
]

cluster_params = {
    "Bcell": {
        "n_clust": 3,
        "pass_clusters": list(range(4)),
        "cluster_labels": {
            1: "up",
            2: "down",
            3: "middle",
            np.nan: np.nan}
    },
    "Bulk": {
        "n_clust": 2,
        "pass_clusters": list(range(3)),
        "cluster_labels": {
            1: "down",
            2: "up",
            np.nan: np.nan}
    },
    "CLL": {
        "n_clust": 6,
        "pass_clusters": [1, 2, 3, 6],
        "cluster_labels": {
            6: "down",
            3: "up",
            2: "extremes",
            1: "middle",
            np.nan: np.nan}
    },
    "CD4": {
        "n_clust": 3,
        "pass_clusters": list(range(4)),
        "cluster_labels": {
            1: "down",
            2: "middle",
            3: "up",
            np.nan: np.nan}
    },
    "CD8": {
        "n_clust": 2,
        "pass_clusters": list(range(3)),
        "cluster_labels": {
            1: "up",
            2: "down",
            np.nan: np.nan}
    },
    "NK": {
        "n_clust": 2,
        "pass_clusters": list(range(3)),
        "cluster_labels": {
            1: "up",
            2: "down",
            np.nan: np.nan}
    },
    "Mono": {
        "n_clust": 2,
        "pass_clusters": list(range(3)),
        "cluster_labels": {
            1: "up",
            2: "down",
            np.nan: np.nan}
    },
}

diff_regions = pd.DataFrame()

# Visualize accessibility of changing regions
# params
alpha = 0.01
varying = fits[(fits['p_value'] < alpha)].index.drop_duplicates()

if matrix_name == "deconv":
    matrix = np.log2(matrix)

to_plot = matrix.loc[:, ~matrix.columns.get_level_values("sample_name").str.contains("|".join(samples_to_exclude))]

g = sns.clustermap(to_plot.loc[varying].T, z_score=1, xticklabels=False, rasterized=True, figsize=(8, 0.2 * matrix.shape[1]), metric="correlation", robust=True)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
g.ax_col_dendrogram.set_rasterized(True)
g.savefig(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "all_D10_variable.clustermap.all_samples.svg"])), dpi=300, bbox_inches="tight")

for cell_type in fits['cell_type'].drop_duplicates():
    print(cell_type)
    varying = fits[(fits["cell_type"] == cell_type) & (fits['p_value'] < alpha)].index.drop_duplicates()

    matrix2 = to_plot[to_plot.columns[
        (to_plot.columns.get_level_values("cell_type") == cell_type) &
        (to_plot.columns.get_level_values("patient_id") != "KI")
    ]]

    g = sns.clustermap(matrix2.loc[varying].T, z_score=1, xticklabels=False, rasterized=True, figsize=(8, 0.2 * matrix2.shape[1]), metric="correlation", robust=True)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_col_dendrogram.set_rasterized(True)
    g.savefig(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, cell_type, "variable.clustermap.svg"])), dpi=300, bbox_inches="tight")

    tp = pd.Series(matrix2.columns.get_level_values("timepoint").str.replace("d", "").astype(int), index=matrix2.columns).sort_values()

    g = sns.clustermap(matrix2.loc[varying, tp.index].T, row_cluster=False, col_cluster=True, z_score=1, xticklabels=False, rasterized=True, figsize=(8, 0.2 * matrix2.shape[1]), metric="correlation", robust=True)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_col_dendrogram.set_rasterized(True)
    g.savefig(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, cell_type, "variable.clustermap.sorted.svg"])), dpi=300, bbox_inches="tight")

    # Flat clusters
    labels = fcluster(Z=g.dendrogram_col.linkage, t=cluster_params[cell_type]["n_clust"], criterion="maxclust")
    labels = [l if l in cluster_params[cell_type]["pass_clusters"] else np.nan for l in labels]

    g2 = sns.clustermap(
        matrix2.loc[varying, tp.index].T,
        col_colors=[plt.get_cmap("Paired")(i) for i in labels],
        row_cluster=False, col_cluster=True, z_score=1, xticklabels=False, yticklabels=matrix2.loc[:, tp.index].columns.get_level_values("sample_name"),
        rasterized=True, figsize=(8, 0.2 * matrix2.shape[1]), metric="correlation", robust=True)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.ax_col_dendrogram.set_rasterized(True)
    g2.savefig(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, cell_type, "variable.clustermap.sorted.cluster_labels.svg"])), dpi=300, bbox_inches="tight")

    cluster_labels = pd.Series(labels, index=varying).to_frame(name="cluster")
    cluster_labels['cluster_name'] = cluster_labels['cluster'].apply(lambda x: cluster_params[cell_type]["cluster_labels"][x] if not pd.isnull(x) else np.nan)
    cluster_labels = cluster_labels.join(fits[fits["cell_type"] == cell_type])

    diff_regions = diff_regions.append(cluster_labels)


    # Cluster patterns
    # this is a mock of the MOHGP underlying posterior.
    cs = cluster_labels['cluster_name'].drop_duplicates().dropna().shape[0]

    fig, axis = plt.subplots(1, cs, figsize=(cs * 4, 1 * 4), sharex=True, sharey=True)
    for i, cluster in enumerate(cluster_labels['cluster_name'].drop_duplicates().dropna()):
        regions = cluster_labels[cluster_labels['cluster_name'] == cluster].index

        X = matrix2.loc[regions, tp.index].mean(axis=0).T.reset_index()
        X['time'] = np.log2(1 + X["timepoint"].str.replace("d", "").astype(int).values)
        # d = X.groupby('time').mean().squeeze().to_frame(name="mean")
        # d['upper_q'] = X.groupby('time').quantile(.975)
        # d['lower_q'] = X.groupby('time').quantile(.025)

        kernel = GPy.kern.RBF(1.0, variance=0.5) + GPy.kern.Bias(1.0, variance=0.05)
        m = GPy.models.GPRegression(X=X[['time']], Y=X[[0]], kernel=kernel)
        m.optimize()
        m.plot([0 - 0.5, max(x) + 0.5], ax=axis[i], legend=None)

        # axis[i].set_ylim((1.5, 3.5))
        axis[i].set_title("Cluster {}\n(n = {})".format(cluster, regions.shape[0]))
    axis[0].set_ylabel("Chromatin accessibility")
    for ax in axis:
        ax.set_xlabel("Time (log2)")
    fig.savefig(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, cell_type, "variable.cluster_means.svg"])), dpi=300, bbox_inches="tight")




diff_regions['comparison_name'] = diff_regions['cell_type'] + "_" + diff_regions['cluster_name'].astype(str)
diff_regions.index.name = "region"
diff_regions.to_csv(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, "variable.csv"])), index=True)


differential_enrichment(
    analysis,
    diff_regions.dropna(),
    data_type="ATAC-seq",
    output_dir="results_deconvolve",
    output_prefix=".".join([output_prefix, matrix_name, library]),
    genome="hg19",
    directional=False,
    max_diff=10000,
    sort_var="D",
    as_jobs=True
)

diff = collect_differential_enrichment(
    diff_regions.dropna(),
    directional=False,
    data_type="ATAC-seq",
    output_dir="results_deconvolve",
    output_prefix=".".join([output_prefix, matrix_name, library]),
    permissive=True)
enrichment_table = pd.read_csv(os.path.join("results_deconvolve", "gp_fit_job.sorted.GPy.meme_ame.csv"))
plot_differential_enrichment(
    enrichment_table,
    "motif",
    data_type="ATAC-seq",
    direction_dependent=False,
    output_dir="results_deconvolve",
    comp_variable="comparison_name",
    output_prefix=".".join([output_prefix, matrix_name, library]),
    top_n=5)
enrichment_table = pd.read_csv(os.path.join("results_deconvolve", "gp_fit_job.sorted.GPy.lola.csv"))
plot_differential_enrichment(
    enrichment_table,
    "lola",
    data_type="ATAC-seq",
    direction_dependent=False,
    output_dir="results_deconvolve",
    comp_variable="comparison_name",
    output_prefix=".".join([output_prefix, matrix_name, library]),
    top_n=5)
enr = pd.read_csv(os.path.join("results_deconvolve", "gp_fit_job.sorted.GPy.enrichr.csv"))
plot_differential_enrichment(
    enrichment_table,
    "enrichr",
    data_type="ATAC-seq",
    direction_dependent=False,
    output_dir="results_deconvolve",
    comp_variable="comparison_name",
    output_prefix=".".join([output_prefix, matrix_name, library]),
    top_n=5)




#         for cluster in cluster_params[cell_type]["pass_clusters"]:
#             var_regions = s_labels[s_labels == cluster].index
#             cluster_name = cluster_params[cell_type]["cluster_labels"][cluster]

#             output_dir = os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, cell_type, "cluster_{}".format(cluster_name)]))
#             if not os.path.exists(output_dir):
#                 os.makedirs(output_dir)

#             characterize_regions_function(
#                 analysis,
#                 df=analysis.coverage_annotated.loc[var_regions, :],
#                 output_dir=output_dir,
#                 prefix="cluster_{}".format(cluster_name),
#                 run=False
#             )



#     # Get enrichments
#     from gprofiler import Gprofiler
#     gp = Gprofiler()
#     res = gp.gprofile(query)
#     res = pd.DataFrame(res, columns=["#", "signf", "p-value", "T", "Q", "Q&T", "Q&T/Q", "Q&T/T", "term ID", "t type", "t group", "t name", "t depth", "Q&T list"])


#     # Cluster expression patterns
#     import GPy
#     import GPclust

#     # Get GP for mean and deviation from it
#     times = np.log2(1 + matrix.columns.get_level_values("timepoint").str.replace("d", "").astype(int).values)
#     k_underlying = GPy.kern.Matern52(input_dim=1, variance=1.0, lengthscale=times.max() / 3.)
#     k_corruption = GPy.kern.Matern52(input_dim=1, variance=0.5, lengthscale=times.max() / 3.) + GPy.kern.White(1, variance=0.01)

#     for cell_type in matrix.columns.get_level_values("cell_type").drop_duplicates():

#         varying = fits[(fits["cell_type"] == cell_type) & (fits['p_value'] < alpha)].index.drop_duplicates()

#         # Get acc matrix for cell type
#         if matrix_name == "sorted":
#             matrix2 = matrix.loc[varying, matrix.columns.get_level_values("cell_type") == cell_type].sort_index(1, level="timepoint")
#         else:
#             matrix2 = np.log2(matrix.loc[varying, matrix.columns.get_level_values("cell_type") == cell_type]).sort_index(1, level="timepoint")
        
#         # Get timepoints
#         times = np.log2(1 + matrix2.columns.get_level_values("timepoint").str.replace("d", "").astype(int).values)

#         # z score row-wise
#         matrix2 = matrix2.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

#         # mean per timepoint
#         # matrix3 = matrix2.T.groupby(level='timepoint').mean().T
#         # times = np.log2(1 + matrix3.columns.str.replace("d", "").astype(int).values)

#         m = GPclust.MOHGP(
#             X=times.reshape(-1, 1), Y=matrix2.values,
#             kernF=k_underlying, kernY=k_corruption,
#             K=4, alpha=1., prior_Z="DP")
#         m.hyperparam_opt_interval = 1000 # how often to optimize the hyperparameters

#         m.hyperparam_opt_args['messages'] = True # turn off the printing of the optimization
#         m.optimize(verbose=True)
#         m.systematic_splits(verbose=True)
#         m.reorder()

#         # Save optimized model
#         m.save(os.path.join("results_deconvolve", output_prefix + "." + cell_type + ".MOHCP.fitted_model.hd5"))

#         # Plot clusters
#         fig = plt.figure()
#         m.plot(newfig=False, on_subplots=True, colour=True, in_a_row=False, joined=False, errorbars=False)
#         [ax.set_rasterized(True) for ax in fig.axes]
#         fig.savefig(os.path.join("results_deconvolve", output_prefix + "." + cell_type + ".MOHCP.fitted_model.clusters.svg"), dpi=300, bbox_inches="tight")

#         # Posterior parameters
#         fig, axis = plt.subplots(2, 1,
#             gridspec_kw={'height_ratios':[12, 1]},
#             figsize=(3 * 4, 1 * 4),
#             tight_layout=True)
#         mat = axis[0].imshow(m.phi.T, cmap=plt.cm.hot, vmin=0, vmax=1, aspect='auto')
#         axis[0].set_xlabel('data index')
#         axis[0].set_ylabel('cluster index')
#         axis[1].set_aspect(0.1)
#         plt.colorbar(mat, cax=axis[1], label="Posterior probability", orientation="horizontal")
#         fig.savefig(os.path.join("results_deconvolve", output_prefix + "." + cell_type + ".MOHCP.fitted_model.posterior_probs.svg"), dpi=300, bbox_inches="tight")


#         g = sns.clustermap(m.phi.T, cmap=plt.cm.hot, vmin=0, vmax=1, xticklabels=False, rasterized=True, cbar_kws={"label": "Posterior probability"})
#         g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
#         g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
#         g.savefig(os.path.join("results_deconvolve", output_prefix + "." + cell_type + ".MOHCP.fitted_model.posterior_probs.clustermap.svg"), dpi=300, bbox_inches="tight")



# def dtw_distance(ts_a, ts_b, d=lambda x, y: abs(x - y), max_warping_window=10000):
#     """Returns the DTW similarity distance between two 2-D
#     timeseries numpy arrays.

#     Arguments
#     ---------
#     ts_a, ts_b : array of shape [n_samples, n_timepoints]
#         Two arrays containing n_samples of timeseries data
#         whose DTW distance between each sample of A and B
#         will be compared
    
#     d : DistanceMetric object (default = abs(x-y))
#         the distance measure used for A_i - B_j in the
#         DTW dynamic programming function
    
#     Returns
#     -------
#     DTW distance between A and B
#     """

#     # Create cost matrix via broadcasting with large int
#     ts_a, ts_b = np.array(ts_a), np.array(ts_b)
#     M, N = len(ts_a), len(ts_b)
#     cost = sys.maxint * np.ones((M, N))

#     # Initialize the first row and column
#     cost[0, 0] = d(ts_a[0], ts_b[0])
#     for i in xrange(1, M):
#         cost[i, 0] = cost[i-1, 0] + d(ts_a[i], ts_b[0])

#     for j in xrange(1, N):
#         cost[0, j] = cost[0, j-1] + d(ts_a[0], ts_b[j])

#     # Populate rest of cost matrix within window
#     for i in xrange(1, M):
#         for j in xrange(max(1, i - max_warping_window),
#                         min(N, i + max_warping_window)):
#             choices = cost[i - 1, j - 1], cost[i, j-1], cost[i-1, j]
#             cost[i, j] = min(choices) + d(ts_a[i], ts_b[j])

#     # Return DTW distance given window 
#     return cost[-1, -1]


# matrix2.loc['time', :] = matrix2.columns.get_level_values("timepoint").str.replace("d", "").astype(int)
# matrix2 = matrix2.T.sort_values("time").T.drop('time')

# dtw_distance(matrix2.loc[varying[0]].squeeze(), matrix2.loc[varying[1]].squeeze()

# scipy.spatial.distance.pdist(matrix2.loc[varying], metric=dtw_distance)

# d = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(matrix2.loc[varying], metric=dtw_distance))

# d = pairwise_distances(matrix2.loc[varying], metric=dtw_distance, n_jobs=-1)

# d.tofile(os.path.join("results_deconvolve", ".".join([output_prefix, matrix_name, library, cell_type, "variable.pairwise_distance_DTW.numpy"])))

# clusterer = hdbscan.HDBSCAN(metric="precomputed", core_dist_n_jobs=-1, algorithm='best')
# labels = clusterer.fit_predict(d)

# # Kshape clustering
# l = kshape(zscore(matrix2.loc[varying].head(), axis=1), cluster_num)
# labels = pd.concat([pd.Series(l[c][1], index=[c] * len(l[c][1])) for c in range(len(l))]).sort_values().index





# tp = pd.Series(matrix2.columns.get_level_values("timepoint").str.replace("d", "").astype(int), index=matrix2.columns).sort_values()

# g = sns.clustermap(
#     matrix2.loc[varying, tp.index].T,
#     col_colors=[sns.color_palette("colorblind")[i] for i in labels],
#     row_cluster=False, col_cluster=True, z_score=1, xticklabels=False,
#     rasterized=True, figsize=(8, 0.2 * matrix2.shape[1]), metric="correlation", robust=True)
# g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
# g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
# g.ax_col_dendrogram.set_rasterized(True)





#         # import hdbscan
#         # from sklearn.cluster import KMeans, AffinityPropagation, MeanShift, SpectralClustering

#         # clusterer = Kmeans(n_clusters=4)
#         # clusterer = AffinityPropagation(preference=-5.0, damping=0.95)
#         # clusterer = MeanShift(n_jobs=1)
#         # clusterer = hdbscan.HDBSCAN(metric="correlation", core_dist_n_jobs=-1, algorithm='generic')
#         # labels = clusterer.fit_predict(matrix2.loc[varying])
