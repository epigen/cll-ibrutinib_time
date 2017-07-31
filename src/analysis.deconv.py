
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools
import seaborn as sns
from looper.models import Project, Sample
from scipy.cluster.hierarchy import fcluster
from statsmodels.stats.multitest import multipletests

from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.general import (normalize_quantiles_r,
                                 subtract_principal_component)
from sklearn.metrics import mean_squared_error, mean_absolute_error, explained_variance_score, r2_score
from scipy.stats import pearsonr, spearmanr


sns.set_style("white")


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

# Do some unsupervised analysis
analysis.samples = deconv_norm.columns.to_frame().apply(Sample, axis=1).tolist()
analysis.deconvolved_data = deconv_norm
analysis.unsupervised(quant_matrix="deconvolved_data", samples=None, attributes_to_plot=["patient_id", "timepoint", "cell_type"], plot_prefix="deconvolved_data")


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

    matrix = deconv_norm.loc[
        all_diff,
        (deconv_norm.columns.get_level_values("patient_id") != "KI") &
        (deconv_norm.columns.get_level_values("cell_type") == cell_type) &
        (deconv_norm.columns.get_level_values("timepoint") != "3d")]

    group_matrix = matrix.T.reset_index().groupby("timepoint").mean().T.loc[all_diff]

    if type(matrix.columns) is pd.core.indexes.multi.MultiIndex:
        matrix.columns = matrix.columns.get_level_values("sample_name")

    # Deconv data
    figsize = (max(5, 0.12 * matrix.shape[1]), 5)
    g = sns.clustermap(matrix,
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
    matrix = sorted_norm.loc[
        all_diff,
        (~sorted_norm.columns.str.contains("KI")) &
        (sorted_norm.columns.str.contains(cell_type.replace("Tcell", "").replace("cell", ""))) &
        (~sorted_norm.columns.str.contains("3d"))]
    if type(matrix.columns) is pd.core.indexes.multi.MultiIndex:
        matrix.columns = matrix.columns.get_level_values("sample_name")
    figsize = (max(5, 0.12 * matrix.shape[1]), 5)

    g = sns.clustermap(matrix,
        yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
        xticklabels=True,
        metric="correlation", rasterized=True, figsize=figsize)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.fig.savefig(os.path.join(analysis.results_dir, cell_type + ".diff_{}.samples.sorted_samples.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=300)

    # Now Bulk data
    matrix = sorted_norm.loc[
        all_diff,
        (~sorted_norm.columns.str.contains("KI")) &
        (sorted_norm.columns.str.contains("Bulk")) &
        (~sorted_norm.columns.str.contains("3d"))]
    if type(matrix.columns) is pd.core.indexes.multi.MultiIndex:
        matrix.columns = matrix.columns.get_level_values("sample_name")
    figsize = (max(5, 0.12 * matrix.shape[1]), 5)

    g = sns.clustermap(matrix,
        yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
        xticklabels=True,
        metric="correlation", rasterized=True, figsize=figsize)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.fig.savefig(os.path.join(analysis.results_dir, cell_type + ".diff_{}.samples.bulk_samples.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=300)












# Annotate peaks
analysis.get_peak_gene_annotation()
analysis.get_peak_genomic_location()
analysis.get_peak_chromatin_state()
analysis.calculate_peak_support([s for s in analysis.samples if s.cell_type == "Bulk"], "summits")

analysis.gene_annotation = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.gene_annotation.csv"))
analysis.region_annotation = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.region_annotation.csv"))
analysis.chrom_state_annotation = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.chromatin_state.csv"))
analysis.support = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.support.csv"))

analysis.to_annot = analysis.accessibility.copy()
analysis.to_annot.columns = analysis.to_annot.columns.get_level_values("sample_name")
analysis.to_annot = analysis.to_annot.join(analysis.coverage[["chrom", "start", "end"]])
analysis.annotate(quant_matrix="to_annot")


# # Let's try a naive differential accessibility
# design_matrix = pd.DataFrame([s.as_series() for s in analysis.samples])
# design_matrix.loc[design_matrix['timepoint'] == "150d", "timepoint"] = "120d"

# results = pd.DataFrame()
# for cell_type in design_matrix['cell_type'].drop_duplicates():

#     des = design_matrix[design_matrix['cell_type'] == cell_type]
#     acc = analysis.accessibility.loc[:, analysis.accessibility.columns.get_level_values("sample_name").isin(des["sample_name"])]

#     acc = acc.loc[acc.mean(axis=1) >= 1.5]

#     res = least_squares_fit(
#         quant_matrix=acc.T, design_matrix=des,
#         standardize_data=True,
#         test_model="~ timepoint", null_model="~ 1", multiple_correction_method="fdr_bh")

#     res['mean'] = acc.mean(1)
#     res['max_beta'] = res[res.columns[res.columns.str.contains("timepoint")]].apply(lambda x: max(abs(x)), axis=1)
#     res['cell_type'] = cell_type
#     res.index.name = "region"
#     res.sort_values("p_value").to_csv("~/res.csv")
