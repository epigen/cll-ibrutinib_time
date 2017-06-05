

import matplotlib.pyplot as plt
import pybedtools
import seaborn as sns
from looper.models import Project
from statsmodels.stats.multitest import multipletests

from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.general import subtract_principal_component

sns.set_style("white")


def subtract_principal_component_by_attribute(df, pc=1, attributes=["CLL"]):
    """
    Given a matrix (n_samples, n_variables), remove `pc` (1-based) from matrix.
    """
    import numpy as np
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


# Start project and analysis objects
prj = Project("metadata/project_config.yaml")
prj.samples = [sample for sample in prj.samples if sample.library == "ATAC-seq"]
analysis = ATACSeqAnalysis(name="cll-time_course", prj=prj, samples=prj.samples)

# Project's attributes
sample_attributes = ["sample_name", "patient_id", "timepoint", "cell_type", "compartment", "response", "cell_number", "batch"]
numerical_attributes = ["CD38_cells_percentage", "cll_cells_%", "cell_number"]
cell_types = ["Bulk", "Bcell", "CLL", "CD4", "CD8", "NK", "Mono"]


# # GET CONSENSUS PEAK SET, ANNOTATE IT, PLOT FEATURES
# # Get consensus peak set from all samples
# analysis.get_consensus_sites(analysis.samples)
# analysis.calculate_peak_support(analysis.samples, "summits")

# # Annotate peaks
# analysis.get_peak_gene_annotation()
# analysis.get_peak_genomic_location()
# analysis.get_peak_chromatin_state()

# # GET CHROMATIN OPENNESS MEASUREMENTS
# # Get coverage values for each peak in each sample of ATAC-seq
# analysis.measure_coverage(analysis.samples)
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
analysis.accessibility = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.coverage.joint_qnorm.pca_fix.power.csv"), index_col=0, header=range(len(sample_attributes)))


# Vizualize all transformations
analysis.coverage_qnorm = analysis.annotate_with_sample_metadata(quant_matrix="coverage_qnorm", attributes=sample_attributes, save=False, assign=False)
unsupervised(analysis, quant_matrix="coverage_qnorm", samples=None, attributes_to_plot=sample_attributes, plot_prefix="coverage_qnorm")

analysis.coverage_qnorm_batchfix = analysis.annotate_with_sample_metadata(quant_matrix="coverage_qnorm_batchfix", attributes=sample_attributes, save=False, assign=False)
unsupervised(analysis, quant_matrix="coverage_qnorm_batchfix", samples=None, attributes_to_plot=sample_attributes, plot_prefix="coverage_qnorm_batchfix")


unsupervised(analysis, quant_matrix="accessibility", samples=None, attributes_to_plot=sample_attributes, plot_prefix="accessibility")


"""
require("limma")
require("data.table")

# read accessibility
df = fread("results/cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.tsv", sep="\t")
acc = df[!c(1:8), !'sample_name']
colnames(acc) = gsub("-", ".", colnames(acc))
regions = df[["sample_name"]][9:nrow(df)]

# get design matrix
design = data.frame(df[c(1:5), ], row.names=1)
design[design == "150d"] <- "120d"

for (cell_type in c("Bulk", "Bcell", "CLL", "CD4", "CD8", "NK", "Mono")){
    print(cell_type)
    # subset, transform
    design2 = design[, design["cell_type", ] == cell_type]
    design3 = model.matrix(~timepoint -1, data.frame(t(design2)))
    # design3 = model.matrix(~timepoint + response -1, data.frame(t(design2)))

    acc2 = acc[, colnames(design[, design["cell_type", ] == cell_type]), with=FALSE]
    acc3 = acc2[, lapply(.SD, as.numeric)]

    # fit model
    model = lmFit(acc3, design3)

    # make contrasts
    if (cell_type == "Mono"){
        contrast.matrix=makeContrasts(
            timepoint001d-timepoint000d,
            timepoint003d-timepoint001d,
            timepoint008d-timepoint003d,
            timepoint120d-timepoint008d,
            timepoint240d-timepoint120d,
            timepoint003d-timepoint000d,
            timepoint008d-timepoint000d,
            timepoint120d-timepoint000d,
            timepoint240d-timepoint000d,
            levels=design3)
    } else if (cell_type == "NK"){
        contrast.matrix=makeContrasts(
            timepoint001d-timepoint000d,
            timepoint003d-timepoint001d,
            timepoint008d-timepoint003d,
            timepoint030d-timepoint008d,
            timepoint120d-timepoint030d,
            timepoint240d-timepoint120d,
            timepoint003d-timepoint000d,
            timepoint008d-timepoint000d,
            timepoint030d-timepoint000d,
            timepoint120d-timepoint000d,
            timepoint240d-timepoint000d,
            levels=design3)
    } else {
        contrast.matrix=makeContrasts(
            timepoint001d-timepoint000d,
            timepoint002d-timepoint001d,
            timepoint003d-timepoint002d,
            timepoint008d-timepoint003d,
            timepoint030d-timepoint008d,
            timepoint120d-timepoint030d,
            timepoint240d-timepoint120d,
            timepoint002d-timepoint000d,
            timepoint003d-timepoint000d,
            timepoint008d-timepoint000d,
            timepoint030d-timepoint000d,
            timepoint120d-timepoint000d,
            timepoint240d-timepoint000d,
            levels=design3)
    }

    # fit constrasts
    fit.contrast = contrasts.fit(model, contrast.matrix)
    # compute p-values
    efit.contrast = eBayes(fit.contrast)

    write.table(
        file=paste0("results/cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.", cell_type, ".diff_timepoint.limma.csv"),
        cbind(regions, model$coefficients, efit.contrast$coefficients, efit.contrast$p.value),
        row.names=F,
        col.names=c(
            "region",
            paste0(colnames(model$coefficients), "_mean"),
            paste0(colnames(efit.contrast$coefficients), "_foldchange"),
            paste0(colnames(efit.contrast$p.value), "_p_value")),
        sep=",")

}

"""

# recover regions
for i, cell_type in enumerate(cell_types):
    diff = pd.read_csv(os.path.join(
        "results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.{}.diff_timepoint.limma.csv".format(cell_type)))
    diff["cell_type"] = cell_type

    if i == 0:
        all_diff = diff
    else:
        for c in all_diff.columns:
            if c not in diff.columns.tolist():
                diff.loc[:, c] = pd.np.nan
        for c in diff.columns:
            if c not in all_diff.columns.tolist():
                all_diff.loc[:, c] = pd.np.nan
        diff = diff.loc[:, all_diff.columns]
        all_diff = all_diff.append(diff, ignore_index=True)

timepoints = all_diff.columns[all_diff.columns.str.contains("_mean")].str.replace("_mean", "")
comparisons = all_diff.columns[all_diff.columns.str.contains("_p_value")].str.replace("_p_value", "")

means = all_diff[["region", "cell_type"] + (timepoints + "_mean").tolist()].drop_duplicates()
means.to_csv(os.path.join(
        "results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.all.diff_timepoint.limma.mean.csv"), index=False)

# arrange in long format
# and correct p-values
all_diff2 = pd.DataFrame()
for cell_type in cell_types:
    for comparison in comparisons:
        print(cell_type, comparison)
        d = all_diff.loc[
            all_diff["cell_type"] == cell_type,
            ["region", "{}_foldchange".format(comparison), "{}_p_value".format(comparison), "cell_type"]
        ]
        d.columns = ["region", "fold_change", "p_value", "cell_type"]
        d["comparison"] = comparison
        d["q_value"] = multipletests(d["p_value"], method="fdr_bh")[1]
        all_diff2 = all_diff2.append(d, ignore_index=True)
all_diff2 = all_diff2.dropna()

all_diff2.to_csv(os.path.join(
        "results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.all.diff_timepoint.limma.csv"), index=False)

all_diff2 = pd.read_csv(os.path.join(
        "results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.all.diff_timepoint.limma.csv"))


# filter for differential
alpha = 5e-2
min_fold = 0.5
diff = all_diff2.loc[
    (all_diff2.loc[:, "q_value"] < alpha) &
    (np.absolute(all_diff2.loc[:, "fold_change"]) > min_fold),
    :
]

# reorder levels to have samples ordered by cell type, timepoint, patient
analysis.accessibility = analysis.accessibility[
    analysis.accessibility.columns.reorder_levels(
        ["cell_type", "timepoint", "patient_id", "sample_name", "compartment", "response", "cell_number", "batch"]).sort_values().get_level_values("sample_name")]

# make timepoints numeric
df = analysis.accessibility.columns.to_series().reset_index().drop(0, axis=1)
df["timepoint"] = df["timepoint"].str.replace("d", "").astype(int)
analysis.accessibility.columns = pd.MultiIndex.from_arrays(df.T.values, names=df.columns)

# get color dataframe
color_df = pd.DataFrame(
    analysis.get_level_colors(index=analysis.accessibility.columns),
    index=analysis.accessibility.columns.names,
    columns=analysis.accessibility.columns)

# see all changing region across all samples
g = sns.clustermap(
    analysis.accessibility.loc[diff['region'].unique(), :].T.astype(float),
    xticklabels=False,
    metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
    row_colors=color_df.iloc[1:, :].values,
    figsize=(20, 20), rasterized=True
)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
g.savefig(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.all_diff.svg"), dpi=300)

g = sns.clustermap(
    analysis.accessibility.loc[diff['region'].unique(), :].T.astype(float),
    row_cluster=False,
    xticklabels=False,
    metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
    row_colors=color_df.iloc[1:, :].values,
    figsize=(20, 20), rasterized=True
)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
g.savefig(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.all_diff.sorted.svg"), dpi=300)


alpha = 0.1
min_fold = 0.25
for cell_type in cell_types:
    print(cell_type)
    diff2 = all_diff2.loc[
        (all_diff2.loc[:, "cell_type"] == cell_type) &
        (all_diff2.loc[:, "q_value"] < alpha) &
        (np.absolute(all_diff2.loc[:, "fold_change"]) > min_fold),
        ['region']
    ]
    if diff2.shape[0] < 5:
        continue
    diff2 = diff2.squeeze().drop_duplicates()

    sample_mask = analysis.accessibility.columns.get_level_values("cell_type") == cell_type
    # accessibility of samples from cell type
    g = sns.clustermap(
        analysis.accessibility.loc[diff2, sample_mask].T.astype(float),
        xticklabels=False,
        metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
        row_colors=color_df.iloc[1:, sample_mask].values,
        figsize=(20, 20), rasterized=True
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.savefig(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.{}_diff.svg".format(cell_type)), dpi=300)
    # sorted
    g = sns.clustermap(
        analysis.accessibility.loc[diff2, sample_mask].T.astype(float),
        row_cluster=False,
        xticklabels=False,
        metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
        row_colors=color_df.iloc[1:, sample_mask].values,
        figsize=(12, 12), rasterized=True,
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.savefig(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.{}_diff.sorted.svg".format(cell_type)), dpi=300)

    # mean accessibility of the cell types
    g = sns.clustermap(
        means.loc[
            (means['region'].isin(diff2)) & (means['cell_type'] == cell_type),
            (timepoints + "_mean")
        ].T.dropna().astype(float),
        xticklabels=False,
        metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
        figsize=(8, 8)
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.savefig(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.{}_diff.mean.svg".format(cell_type)), dpi=300)

    g = sns.clustermap(
        means.loc[
            (means['region'].isin(diff2)) & (means['cell_type'] == cell_type),
            (timepoints + "_mean")
        ].T.dropna().astype(float),
        xticklabels=False, row_cluster=False,
        metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
        figsize=(8, 8)
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.savefig(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.{}_diff.mean.sorted.svg".format(cell_type)), dpi=300)
