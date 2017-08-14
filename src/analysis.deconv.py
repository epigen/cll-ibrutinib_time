
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools
import seaborn as sns
from looper.models import Project, Sample
from scipy.cluster.hierarchy import fcluster
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import (explained_variance_score, mean_absolute_error,
                             mean_squared_error, r2_score)
from statsmodels.stats.multitest import multipletests

from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.general import (collect_differential_enrichment,
                                 differential_enrichment,
                                 normalize_quantiles_r,
                                 plot_differential_enrichment,
                                 subtract_principal_component)

sns.set_style("white")


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










#  Try pseudotime ordering
from mapalign.embed import DiffusionMapEmbedding
from sklearn.manifold import SpectralEmbedding

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


from sklearn.decomposition import PCA
from scipy.stats import zscore

pca = PCA()
pcs = pd.DataFrame(pca.fit_transform(np.log2(1 + X.T)), index=X.columns)

t = 100
de = DiffusionMapEmbedding(alpha=0.5, diffusion_time=t, affinity='markov',
                           n_components=2).fit_transform(pcs.iloc[:, :10])

times = pcs.index.get_level_values("timepoint").str.replace("d", "").astype(int)
plt.scatter(zscore(de[:, 0]), zscore(de[:, 1]), c=np.log10(1 + times), cmap="inferno")

plt.scatter(zscore(pcs.loc[:, 0]), zscore(pcs.loc[:, 1]), c=np.log10(1 + times), cmap="inferno")



# Gaussian Process based

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel

import multiprocessing
import parmap


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

import GPy
import multiprocessing
import parmap
import matplotlib.pyplot as plt

# deconv_norm = pd.read_csv(os.path.join(analysis.results_dir, "coverage.cell_type_deconvoluted.qnorm.csv"), index_col=0, header=range(4))
deconv_norm = pd.read_csv(os.path.join("results_deconvolve", "coverage.cell_type_deconvoluted.qnorm.csv"), index_col=0, header=range(4))
X = deconv_norm.loc[:, ~deconv_norm.columns.get_level_values("patient_id").isin(["KI"])]
X = X.loc[:, X.columns.get_level_values("cell_type").isin(["CLL"])]
X = X.loc[:, ~X.columns.get_level_values("timepoint").isin(["240d"])]

chunks = 200
cell_type = "CLL"
output_prefix = "gp_fit_job"
output_dir = "/scratch/users/arendeiro/gp_fit_job"
r = np.arange(0, X.shape[0], chunks)

for start, end in zip(r, r[1:]) + [(r[-1], X.shape[0])]:
    range_name = "{}-{}".format(start, end)
    name = "-".join([output_prefix, range_name])
    log = os.path.join(output_dir, "log", name + ".log")
    job = """python ~/jobs/gp_fit_job.py --data-range {} --range-delim - --cell-type {} --output-prefix {} --output-dir {}""".format(
        range_name, cell_type, output_prefix, os.path.join(output_dir, "output"))
    cmd = """sbatch -J {} -o {} -p shortq -c 1 --mem 8000 --wrap "{}" """.format(
        name, log, job)

    if not os.path.exists(os.path.join(output_dir, "output", output_prefix + range_name + ".csv")):
        os.system(cmd)
