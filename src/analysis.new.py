
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools
import seaborn as sns
from looper.models import Project
from scipy.cluster.hierarchy import fcluster
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
analysis.accessibility = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.coverage.joint_qnorm.pca_fix.power.csv"), index_col=0, header=range(len(sample_attributes)))


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

means = pd.read_csv(os.path.join(
        "results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.all.diff_timepoint.limma.mean.csv"))
means.index = means["region"]

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


timepoints = pd.Series(sorted(list(set(all_diff2["comparison"].drop_duplicates().str.split(" - ").sum()))))
comparisons = all_diff2["comparison"].drop_duplicates()

# to censor potentially:
samples_to_exclude = [
    "KAF_30d_Bulk",
    "PBGY_1d_NK",
    "PBGY_8d_Mono",
    "PBGY_8d_Bcell",
    "KAF_30d_Bcell"]

# to annotate with genes
analysis.gene_annotation.index = (
    analysis.gene_annotation["chrom"] + ":" +
    analysis.gene_annotation["start"].astype(str) + "-" +
    analysis.gene_annotation["end"].astype(str))

alpha = 0.1
min_fold = 0.25
n_top = 500
for cell_type in cell_types:
    print(cell_type)

    # Filter for cell type
    if cell_type != "CD8":
        diff2 = all_diff2.loc[(all_diff2.loc[:, "cell_type"] == cell_type) & (~all_diff2.loc[:, "comparison"].str.contains("240d|280d")), :]
    elif cell_type == "Bcell":
        diff2 = all_diff2.loc[(all_diff2.loc[:, "cell_type"] == cell_type) & (~all_diff2.loc[:, "comparison"].str.contains("8d|30d|240d|280d")), :]
    elif cell_type == "NK":
        diff2 = all_diff2.loc[(all_diff2.loc[:, "cell_type"] == cell_type) & (~all_diff2.loc[:, "comparison"].str.contains("1d|240d|280d")), :]
    elif cell_type == "Mono":
        diff2 = all_diff2.loc[(all_diff2.loc[:, "cell_type"] == cell_type) & (~all_diff2.loc[:, "comparison"].str.contains("8d|240d|280d")), :]
    else:
        diff2 = all_diff2.loc[(all_diff2.loc[:, "cell_type"] == cell_type) & (~all_diff2.loc[:, "comparison"].str.contains("2d|240d|280d")), :]

    # Filter for thresholds
    diff_regions = diff2.loc[
        (diff2.loc[:, "q_value"] < alpha) &
        (np.absolute(diff2.loc[:, "fold_change"]) > min_fold),
        ['region']
    ]
    # If there aren't any significant, take top N to visualize
    if diff_regions.shape[0] < 100:
        diff_regions = diff2.loc[(diff2.loc[:, "cell_type"] == cell_type), :].sort_values("p_value").head(n_top)['region']
    diff_regions = diff_regions.squeeze().drop_duplicates()

    sample_mask = (
        (~analysis.accessibility.columns.get_level_values("sample_name").str.contains("|".join(samples_to_exclude))) &
        (analysis.accessibility.columns.get_level_values("cell_type") == cell_type) &
        (analysis.accessibility.columns.get_level_values("timepoint") < 240))
    # accessibility of samples from cell type
    g = sns.clustermap(
        analysis.accessibility.loc[diff_regions, sample_mask].T.astype(float),
        xticklabels=False,
        metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
        row_colors=color_df.iloc[1:, sample_mask].values,
        figsize=(20, 20), rasterized=True
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.savefig(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.{}_diff.svg".format(cell_type)), dpi=300)

    # sorted
    g = sns.clustermap(
        analysis.accessibility.loc[diff_regions, sample_mask].T.astype(float),
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
            (means['region'].isin(diff_regions)) & (means['cell_type'] == cell_type),
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
            (means['region'].isin(diff_regions)) & (means['cell_type'] == cell_type),
            (timepoints + "_mean")
        ].T.dropna().astype(float),
        xticklabels=False, row_cluster=False,
        metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
        figsize=(8, 8)
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.savefig(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.{}_diff.mean.sorted.svg".format(cell_type)), dpi=300)

    # Annotate differential regions with genes
    diff_genes = pd.Series(analysis.gene_annotation.loc[diff_regions, "gene_name"].str.split(",").sum()).drop_duplicates().sort_values()
    diff_genes.to_csv(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.{}_diff.gene_signature.csv".format(cell_type)), index=False)

    # With weights too (change across all samples from earliest to latest timepoint)
    mean_diff = means.loc[
            (means['region'].isin(diff_regions)) & (means['cell_type'] == cell_type),
            (timepoints + "_mean")
    ].T.dropna().astype(float)

    mean_fold = np.log2(mean_diff.iloc[-1, :] / mean_diff.iloc[0, :])
    mean_fold.name = "fold_change"

    g = analysis.gene_annotation.loc[diff_regions, "gene_name"].str.split(",").apply(pd.Series).stack().reset_index(level=1, drop=True).to_frame(name="gene_name")

    diff_genes = g.join(mean_fold).dropna().groupby("gene_name")['fold_change'].mean().sort_values()
    diff_genes.to_csv(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.{}_diff.gene_signature.weighted.csv".format(cell_type)), index=True)

    # Get enrichments

    # for all regions
    out_dir = os.path.join("results", "enrichment." + cell_type)
    if cell_type in ["Mono", "NK"]:
        out_dir2 = os.path.join(out_dir, "1_cluster")

        for path in [out_dir, out_dir2]:
            if not os.path.exists(path):
                os.makedirs(path)

            characterize_regions_function(
                analysis,
                df=analysis.coverage_annotated.loc[diff_regions, :],
                output_dir=out_dir2,
                prefix="cluster_{}".format(1)
            )
    # for two clusters
    if cell_type in ["CLL", "CD4", "CD8"]:
        # cluster again
        g = sns.clustermap(
            analysis.accessibility.loc[diff_regions, sample_mask].T.astype(float),
            row_cluster=False, metric="correlation", z_score=1, vmin=-1, vmax=4
        )
        clusters = fcluster(g.dendrogram_col.linkage, t=2, criterion='maxclust')

        # plot again with cluster identities
        g = sns.clustermap(
            analysis.accessibility.loc[diff_regions, sample_mask].T.astype(float),
            row_cluster=False,
            xticklabels=False,
            metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
            row_colors=color_df.iloc[1:, sample_mask].values,
            col_colors=[sns.color_palette("colorblind", 2)[i - 1] for i in clusters],
            figsize=(12, 12), rasterized=True,
        )
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.savefig(os.path.join(out_dir2, "clustermap.sorted.labeled.svg"), dpi=300)

        for cluster in np.unique(clusters):
            out_dir2 = os.path.join(out_dir, "2_clusters.cluster_{}".format(cluster))
            if not os.path.exists(out_dir2):
                os.makedirs(out_dir2)

            index = diff_regions[clusters == cluster]

            characterize_regions_function(
                analysis,
                df=analysis.coverage_annotated.loc[index, :],
                output_dir=out_dir2,
                prefix="cluster_{}".format(cluster)
            )
    # for four clusters
    if cell_type in ["CD8", "CD4"]:
        out_dir2 = os.path.join(out_dir, "4_clusters")
        if not os.path.exists(out_dir2):
            os.makedirs(out_dir2)
        clusters = fcluster(g.dendrogram_col.linkage, t=4, criterion='maxclust')

        # plot again with cluster identities
        g = sns.clustermap(
            analysis.accessibility.loc[diff_regions, sample_mask].T.astype(float),
            row_cluster=False,
            xticklabels=False,
            metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
            row_colors=color_df.iloc[1:, sample_mask].values,
            col_colors=[sns.color_palette("colorblind", 4)[i - 1] for i in clusters],
            figsize=(12, 12), rasterized=True,
        )
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.savefig(os.path.join(out_dir2, "clustermap.sorted.labeled.svg"), dpi=300)

        for cluster in np.unique(clusters):
            out_dir2 = os.path.join(out_dir, "4_clusters.cluster_{}".format(cluster))
            if not os.path.exists(out_dir2):
                os.makedirs(out_dir2)

            index = diff_regions[clusters == cluster]

            characterize_regions_function(
                analysis,
                df=analysis.coverage_annotated.loc[index, :],
                output_dir=out_dir2,
                prefix="cluster_{}".format(cluster)
            )

    # for seven clusters
    if cell_type in ["CD4"]:
        out_dir2 = os.path.join(out_dir, "7_clusters")
        if not os.path.exists(out_dir2):
            os.makedirs(out_dir2)
        clusters = fcluster(g.dendrogram_col.linkage, t=7, criterion='maxclust')

        # plot again with cluster identities
        g = sns.clustermap(
            analysis.accessibility.loc[diff_regions, sample_mask].T.astype(float),
            row_cluster=False,
            xticklabels=False,
            metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
            row_colors=color_df.iloc[1:, sample_mask].values,
            col_colors=[sns.color_palette("colorblind", 7)[i - 1] for i in clusters],
            figsize=(12, 12), rasterized=True,
        )
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.savefig(os.path.join(out_dir2, "clustermap.sorted.labeled.svg"), dpi=300)

        for cluster in np.unique(clusters):
            out_dir2 = os.path.join(out_dir, "7_clusters.cluster_{}".format(cluster))
            if not os.path.exists(out_dir2):
                os.makedirs(out_dir2)

            index = diff_regions[clusters == cluster]

            characterize_regions_function(
                analysis,
                df=analysis.coverage_annotated.loc[index, :],
                output_dir=out_dir2,
                prefix="cluster_{}".format(cluster)
            )


# Run Enrichr
find results -name "*.symbols.txt" \
-exec sbatch -J {}.enrichr -o {}.enrichr.log -p shortq -c 1 --mem 4000 \
--wrap "python ~/jobs/run_Enrichr.py --input-file {} --output-file {}.enrichr.csv" \;

# Run LOLA
find results -name "*_regions.bed" \
-exec sbatch -J {}.lola -o {}.lola.log -p shortq -c 8 --mem 24000 \
--wrap "Rscript ~/jobs/run_LOLA.R {} ~/cll-time_course/results/cll-time_course_peak_set.bed hg19" \;

# Run AME
for F in `find results -name "*_regions.fa"`
do
DIR=`dirname $F`
echo $F $DIR
sbatch -J "${F}.meme_ame" -o "${F}.meme_ame.log" -p shortq -c 1 --mem 4000 \
--wrap \
"fasta-dinucleotide-shuffle -c 1 -f "$F" > "$F".shuffled.fa; \
ame --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 \
--control "$F".shuffled.fa -o "$DIR" "$F" ~/resources/motifs/motif_databases/HUMAN/HOCOMOCOv10.meme"
done


# Collect enrichments and plot
enrichments = pd.DataFrame()

for cell_type in cell_types:
    for n_clust in [1, 2, 4, 7]:
        out_dir = os.path.join("results", "enrichment." + cell_type, "{}_clusters".format(n_clust))
        for cluster in range(1, n_clust + 1):
            if cluster == 1:
                f = os.path.join(out_dir.replace("clusters", "cluster"), "cluster_{}_genes.symbols.txt.enrichr.csv".format(cluster))
            else:
                f = os.path.join(out_dir + ".cluster_{}".format(cluster), "cluster_{}_genes.symbols.txt.enrichr.csv".format(cluster))

            try:
                enr = pd.read_csv(f)
            except IOError:
                continue

            enr["cell_type"] = cell_type
            enr["cluster"] = cluster
            enrichments = enrichments.append(enr, ignore_index=True)

for n_top in [10, 20]:
    for gene_set_library in enrichments['gene_set_library'].drop_duplicates():
        enr = enrichments[enrichments['gene_set_library'] == gene_set_library]

        # transform p-values
        enr.loc[:, "q_value"] = -np.log10(enr["p_value"])

        # get top N enriched per set
        enr.loc[:, "label"] = enr["cell_type"] + " " + enr["cluster"].astype(str)
        t = enr.groupby("label")["q_value"].nlargest(n_top)
        t = t[~t.index.get_level_values(0).str.contains("NK|Mono")]
        terms = enr.loc[t.index.levels[1], "description"].drop_duplicates()

        # make pivot table
        enr_pivot = pd.pivot_table(enr, index="description", columns="label", values="q_value", fill_value=0).loc[terms, :]

        g = sns.clustermap(
            enr_pivot,
            col_cluster=False,
            metric="correlation", cbar_kws={"label": "-log10(P-value)"},
            figsize=(4, 20), rasterized=True,
        )
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.savefig(os.path.join("results", "enrichments.all_cell_types.all_clusters.{}.top{}.svg".format(gene_set_library, n_top)), dpi=300)

        # reduce to maximum per cell type (across clusters)
        enr_pivot.columns = enr_pivot.columns.str.replace(" .*", "")
        enr_pivot = enr_pivot.T.groupby(level=0).max().T
        g = sns.clustermap(
            enr_pivot,
            col_cluster=False,
            metric="correlation", cbar_kws={"label": "max(-log10(P-value))"},
            figsize=(4, 20), rasterized=True,
        )
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.savefig(os.path.join("results", "enrichments.all_cell_types.all_clusters.{}.top{}.reduced_max.svg".format(gene_set_library, n_top)), dpi=300)


# LOLA
lola = pd.DataFrame()

for cell_type in cell_types:
    for n_clust in [1, 2, 4, 7]:
        out_dir = os.path.join("results", "enrichment." + cell_type, "{}_clusters".format(n_clust))
        for cluster in range(1, n_clust + 1):
            if cluster == 1:
                f = os.path.join(out_dir.replace("clusters", "cluster"), "allEnrichments.txt")
            else:
                f = os.path.join(out_dir + ".cluster_{}".format(cluster), "allEnrichments.txt")

            try:
                enr = pd.read_csv(f, sep="\t")
            except IOError:
                continue

            enr["cell_type"] = cell_type
            enr["cluster"] = cluster
            lola = lola.append(enr, ignore_index=True)


for n_top in [10, 20]:
    #
    enr = lola.copy()
    enr["set_label"] = (
        enr[['description', 'cellType', 'tissue', 'antibody', 'treatment']]
        .astype(str)
        .apply(lambda x: " ".join([y for y in set(x) if not y in ['None', 'nan']]), axis=1)
    )

    # get top N enriched per set
    enr.loc[:, "label"] = enr["cell_type"] + " " + enr["cluster"].astype(str)
    t = enr.groupby("label")["pValueLog"].nlargest(n_top)
    t = t[~t.index.get_level_values(0).str.contains("NK|Mono")]
    terms = enr.loc[t.index.levels[1], "set_label"].drop_duplicates()

    # make pivot table
    enr_pivot = pd.pivot_table(enr, index="set_label", columns="label", values="pValueLog", fill_value=0).loc[terms, :]

    g = sns.clustermap(
        enr_pivot,
        col_cluster=False,
        metric="correlation", cbar_kws={"label": "-log10(P-value)"},
        figsize=(4, 20), rasterized=True,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.savefig(os.path.join("results", "lola.all_cell_types.all_clusters.top{}.svg".format(n_top)), dpi=300)
    # with z-score
    g = sns.clustermap(
        enr_pivot,
        col_cluster=False,
        metric="correlation", cbar_kws={"label": "-log10(P-value)"},
        figsize=(4, 20), rasterized=True, z_score=1
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.savefig(os.path.join("results", "lola.all_cell_types.all_clusters.top{}.z-score.svg".format(n_top)), dpi=300)

    # reduce to maximum per cell type (across clusters)
    enr_pivot.columns = enr_pivot.columns.str.replace(" .*", "")
    enr_pivot = enr_pivot.T.groupby(level=0).max().T
    g = sns.clustermap(
        enr_pivot,
        col_cluster=False,
        metric="correlation", cbar_kws={"label": "max(-log10(P-value))"},
        figsize=(4, 20), rasterized=True,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.savefig(os.path.join("results", "lola.all_cell_types.all_clusters.top{}.reduced_max.svg".format(n_top)), dpi=300)

