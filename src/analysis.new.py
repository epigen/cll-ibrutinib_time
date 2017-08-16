
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



def get_principal_component_by_attribute(df, pc=1, attributes=["CLL"]):
    """
    Given a matrix (n_samples, n_variables), remove `pc` (1-based) from matrix.
    """
    import numpy as np
    from sklearn.decomposition import PCA
    from ngs_toolkit.general import z_score

    pc -= 1

    X2 = pd.DataFrame(index=df.index)
    for attr in attributes:
        print(attr)
        sel = df.index[df.index.str.contains(attr)]
        X = df.loc[sel, :]

        # PCA
        pca = PCA()
        X_hat = pca.fit_transform(X)

        # Remove PC
        X2.loc[sel, "pc"] = X_hat[:, pc]
        X2.loc[sel, "pc_zscore"] = z_score(X_hat[:, pc])
    return X2


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


def remove_batch_effect(quant_matrix, design_matrix, test_model, null_model="~ 1", standardize_data=True):
    """
    Fit a least squares model with only categorical predictors.
    Gets p-values by getting the log likelihood ration compared to a `null_model`.
    
    `quant_matrix` is a (samples, variables) matrix.
    `design_matrix` is a (samples, variables) dataframe with all the variables in `test_model`.
    """
    from sklearn.preprocessing import StandardScaler
    import patsy
    from scipy.linalg import lstsq
    from scipy import stats
    from statsmodels.sandbox.stats.multicomp import multipletests

    # # to test
    # quant_matrix = np.random.random(10000000).reshape(100, 100000)
    # P = np.concatenate([[0] * 50, [1] * 50])
    # Q = np.concatenate([[0] * 25, [1] * 25] + [[0] * 25, [1] * 25])
    # design_matrix = pd.DataFrame([P, Q], index=["P", "Q"]).T
    # quant_matrix = quant_matrix.T * (1 + (design_matrix.sum(axis=1) * 4).values)
    # quant_matrix = pd.DataFrame(quant_matrix.T)
    # test_model = "~ Q + P"
    # null_model = "~ 1"

    if standardize_data:
        norm = StandardScaler()
        quant_matrix = pd.DataFrame(
            norm.fit_transform(quant_matrix),
            index=quant_matrix.index, columns=quant_matrix.columns)

    A1 = patsy.dmatrix(test_model, design_matrix)
    betas1, residuals1, _, _ = lstsq(A1, quant_matrix)

    A0 = patsy.dmatrix(null_model, design_matrix)
    betas0, residuals0, _, _ = lstsq(A0, quant_matrix)

    results = pd.DataFrame(betas1.T, columns=A1.design_info.column_names, index=quant_matrix.columns)

    # Calculate the log-likelihood ratios
    n = float(quant_matrix.shape[0])
    results['model_residuals'] = residuals1
    results['null_residuals'] = residuals0
    results['model_log_likelihood'] = (-n / 2.) * np.log(2 * np.pi) - n / 2. * np.log(results['model_residuals'] / n) - n / 2.
    results['null_log_likelihood'] = (-n / 2.) * np.log(2 * np.pi) - n / 2. * np.log(results['null_residuals'] / n) - n / 2.

    results['log_likelihood_ratio'] = results['model_log_likelihood'] - results['null_log_likelihood']
    results['p_value'] = stats.chi2.sf(2 * results['log_likelihood_ratio'], df=betas1.shape[0] - betas0.shape[0])
    results['q_value'] = multipletests(results['p_value'], method=multiple_correction_method)[1]

    if not standardize_data:
        results["mean"] = quant_matrix.mean(axis=0)

    return results


def deseq_analysis(
        count_matrix, experiment_matrix, formula,
        output_dir, output_prefix,
        overwrite=True, alpha=0.05, independent_filtering=False,):
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
    experiment_matrix = experiment_matrix.loc[count_matrix.columns, :]
    if experiment_matrix.index.name != "sample_name":
        try:
            experiment_matrix = experiment_matrix.set_index("sample_name")
        except KeyError:
            pass

    # save the matrices just in case
    count_matrix.to_csv(os.path.join(output_dir, output_prefix + ".count_matrix.tsv"), sep="\t")
    experiment_matrix.to_csv(os.path.join(output_dir, output_prefix + ".experiment_matrix.tsv"), sep="\t")
    comparison_table.to_csv(os.path.join(output_dir, output_prefix + ".comparison_table.tsv"), sep="\t")

    # Run DESeq analysis
    dds = _DESeqDataSetFromMatrix(
        countData=count_matrix.astype(int),
        colData=experiment_matrix,
        design=_as_formula(formula))
    dds = _DESeq(dds, parallel=True)
    # _save(dds, file=os.path.join(output_dir, output_prefix + ".deseq_dds_object.Rdata"))

    comps = [str(x) for x in _resultsNames(dds)][3:]
    results = pd.DataFrame()
    for comp in comps:
        out_file = os.path.join(output_dir, output_prefix + ".deseq_result.{}_vs_000d.csv".format(comp))
        if not overwrite and os.path.exists(out_file):
            continue
        print("Doing comparison '{}'".format(comp))

        res = _as_data_frame(
            _results(dds, contrast=np.array(["timepoint", comp.replace("timepoint", ""), "000d"]), alpha=alpha, independentFiltering=independent_filtering, parallel=True))

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



# # Let's try to add the PC1 location of each sample per cell type as covariate to DESeq2
# from ngs_toolkit.general import plot_differential
# cov = get_principal_component_by_attribute(
#     analysis.coverage.T.drop(["chrom", "start", "end"]),
#     pc=1, attributes=cell_types)
# cov.to_csv(os.path.join("results", "per_cell_type_pc1_position.csv"))

# comparison_table = pd.DataFrame(
#         [s.as_series() for s in analysis.samples]
#     ).set_index("sample_name")[['cell_type', 'timepoint', 'batch']]
# comparison_table = comparison_table.join(cov)

# comparison_table.loc[comparison_table['timepoint'] == "150d", 'timepoint'] = '120d'


# output_dir = os.path.join("results", "differential_analysis_ATAC-seq.deseq")
# output_prefix = "differential_analysis"

# analysis.support.index = (
#     analysis.support["chrom"] + ":" +
#     analysis.support["start"].astype(str) + "-" +
#     analysis.support["end"].astype(str))

# for cell_type in comparison_table['cell_type'].drop_duplicates():

#     c = comparison_table[comparison_table['cell_type'] == cell_type]
#     cov = analysis.coverage.loc[:, c.index]

#     deseq_analysis(
#         count_matrix=cov;
#         experiment_matrix=c;
#         formula="~ pc_zscore + timepoint";
#         output_dir=output_dir;
#         output_prefix=output_prefix + "." + cell_type;
#         overwrite=True; alpha=0.05; independent_filtering=True
#         )

#     results = pd.read_csv(os.path.join(output_dir, output_prefix + "." + cell_type + ".deseq_result.all_comparisons.csv"), index_col=0)

#     # exclude peaks not called in ~10% of samples
#     mask = analysis.support.columns.str.contains(cell_type)
#     support = (analysis.support.loc[:, mask] >= 1).sum(axis=1) / mask.sum()
#     pass_support = support[support > 0.5].index
#     # exclude sex chroms
#     pass_support = pass_support[~pass_support.str.contains("chrX|chrY|chrUn_|_random")]

#     results = results.loc[pass_support, :]

#     for comp in results["comparison_name"].drop_duplicates():
#         results.loc[results['comparison_name'] == comp, "padj"] = _padjust(
#             results.loc[results['comparison_name'] == comp, "pvalue"], method="BH")


#     # Visualize regions and samples found in the differential comparisons
#     plot_differential(
#         analysis,
#         results[~results['comparison_name'].str.contains("pc")],
#         None,
#         samples=[s for s in analysis.samples if s.name in c.index],
#         data_type="ATAC-seq",
#         alpha=0.05,
#         corrected_p_value=True,
#         fold_change=None,
#         output_dir=output_dir,
#         output_prefix=output_prefix + "." + cell_type)




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
analysis.support.index = (
    analysis.support["chrom"] + ":" +
    analysis.support["start"].astype(str) + "-" +
    analysis.support["end"].astype(str))

for i, cell_type in enumerate(cell_types):
    print(cell_type)
    diff = pd.read_csv(os.path.join(
        "results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.{}.diff_timepoint.limma.csv".format(cell_type)),
        index_col=0)
    diff["cell_type"] = cell_type

    # filter out regions never called as a peak
    # in the respective cell type
    sample_names = [s.name for s in analysis.samples if s.cell_type == cell_type]
    cell_type_support = (analysis.support[sample_names] >= 1).sum(1) / float(len(sample_names))
    specific_regions = cell_type_support[cell_type_support >= 0.5].index
    specific_regions = specific_regions[~specific_regions.str.contains("chrX|chrY|_random|Un_")]
    diff = diff.loc[specific_regions, :].reset_index()

    # Append
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

means = all_diff[["index", "cell_type"] + (timepoints + "_mean").tolist()].drop_duplicates()
means.to_csv(os.path.join(
        "results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.all.diff_timepoint.limma.mean.csv"), index=False)

means = pd.read_csv(os.path.join(
        "results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.all.diff_timepoint.limma.mean.csv"))
means.index = means["index"]


# arrange in long format
# and correct p-values
import pandas as pd
from rpy2.robjects import numpy2ri, pandas2ri
import rpy2.robjects as robjects
numpy2ri.activate()
pandas2ri.activate()
_padjust = robjects.r('p.adjust')


all_diff2 = pd.DataFrame()
for cell_type in cell_types:
    for comparison in comparisons:
        print(cell_type, comparison)
        d = all_diff.loc[
            all_diff["cell_type"] == cell_type,
            ["index", "{}_foldchange".format(comparison), "{}_p_value".format(comparison), "cell_type"]
        ]
        d = d.dropna()
        if d.shape[0] == 0:
            continue
        d.columns = ["index", "fold_change", "p_value", "cell_type"]
        d["comparison"] = comparison
        # d["q_value"] = multipletests(d["p_value"], method="fdr_bh")[1]
        d["q_value"] = _padjust(d["p_value"], method="BH")
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
# df["timepoint"] = df["timepoint"].str.replace("d", "").astype(int)
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
    figsize=(20, max(6, analysis.accessibility.shape[1] * 0.12)), rasterized=True
)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
g.ax_col_dendrogram.set_rasterized(True)
g.savefig(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.all_diff.svg"), dpi=300)

g = sns.clustermap(
    analysis.accessibility.loc[diff['region'].unique(), :].T.astype(float),
    row_cluster=False,
    xticklabels=False,
    metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
    row_colors=color_df.iloc[1:, :].values,
    figsize=(20, max(6, analysis.accessibility.shape[1] * 0.12)), rasterized=True
)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
g.ax_col_dendrogram.set_rasterized(True)
g.savefig(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.all_diff.sorted.svg"), dpi=300)


timepoints = pd.Series(sorted(list(set(all_diff2["comparison"].drop_duplicates().str.split(" - ").sum()))))
comparisons = all_diff2["comparison"].drop_duplicates()

# to censor potentially:
samples_to_exclude = [
    "KAF_30d_Bulk",
    "KAF_30d_CD4",
    "KAF_30d_CD8",
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
for cell_type in cell_types[1:]:
    print(cell_type)

    # # Filter for cell type
    # if cell_type != "CD8":
    #     diff2 = all_diff2.loc[(all_diff2.loc[:, "cell_type"] == cell_type) & (~all_diff2.loc[:, "comparison"].str.contains("240d|280d")), :]
    # elif cell_type == "Bcell":
    #     diff2 = all_diff2.loc[(all_diff2.loc[:, "cell_type"] == cell_type) & (~all_diff2.loc[:, "comparison"].str.contains("8d|30d|240d|280d")), :]
    # elif cell_type == "NK":
    #     diff2 = all_diff2.loc[(all_diff2.loc[:, "cell_type"] == cell_type) & (~all_diff2.loc[:, "comparison"].str.contains("1d|240d|280d")), :]
    # elif cell_type == "Mono":
    #     diff2 = all_diff2.loc[(all_diff2.loc[:, "cell_type"] == cell_type) & (~all_diff2.loc[:, "comparison"].str.contains("8d|240d|280d")), :]
    # else:
    #     diff2 = all_diff2.loc[(all_diff2.loc[:, "cell_type"] == cell_type) & (~all_diff2.loc[:, "comparison"].str.contains("2d|240d|280d")), :]
    diff2 = all_diff2.loc[
        (all_diff2.loc[:, "cell_type"] == cell_type) &
        (
            (all_diff2.loc[:, "comparison"] == "timepoint120d - timepoint000d")  |
            (all_diff2.loc[:, "comparison"] == "timepoint030d - timepoint000d")
        ), :]

    # # Filter for thresholds
    diff_regions = diff2.loc[
        (diff2.loc[:, "q_value"] < alpha),#  &
        # (np.absolute(diff2.loc[:, "fold_change"]) > min_fold),
        ['region']
    ]
    # # If there aren't any significant, take top N to visualize
    if diff_regions.shape[0] < 100:
        print("No diff")
        diff_regions = diff2.loc[(diff2.loc[:, "cell_type"] == cell_type), :].sort_values("p_value").head(n_top)['region']
    diff_regions = diff_regions.squeeze().drop_duplicates()

    sample_mask = (
        (~analysis.accessibility.columns.get_level_values("sample_name").str.contains("|".join(samples_to_exclude))) &
        (analysis.accessibility.columns.get_level_values("cell_type") == cell_type) &
        (analysis.accessibility.columns.get_level_values("timepoint") < 240))
    # accessibility of samples from cell type
    # g = sns.clustermap(
    #     analysis.accessibility.loc[diff_regions, sample_mask].T.astype(float),
    #     xticklabels=False,
    #     metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
    #     row_colors=color_df.iloc[1:, sample_mask].values,
    #     figsize=(20, 20), rasterized=True
    # )
    # g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    # g.savefig(os.path.join("results", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.{}_diff.svg".format(cell_type)), dpi=300)

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
    # if cell_type in ["Mono", "NK"]:
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
    # if cell_type in ["CLL", "CD4", "CD8"]:
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
    # # for four clusters
    # if cell_type in ["CD8", "CD4"]:
    #     out_dir2 = os.path.join(out_dir, "4_clusters")
    #     if not os.path.exists(out_dir2):
    #         os.makedirs(out_dir2)
    #     clusters = fcluster(g.dendrogram_col.linkage, t=4, criterion='maxclust')

    #     # plot again with cluster identities
    #     g = sns.clustermap(
    #         analysis.accessibility.loc[diff_regions, sample_mask].T.astype(float),
    #         row_cluster=False,
    #         xticklabels=False,
    #         metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
    #         row_colors=color_df.iloc[1:, sample_mask].values,
    #         col_colors=[sns.color_palette("colorblind", 4)[i - 1] for i in clusters],
    #         figsize=(12, 12), rasterized=True,
    #     )
    #     g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    #     g.savefig(os.path.join(out_dir2, "clustermap.sorted.labeled.svg"), dpi=300)

    #     for cluster in np.unique(clusters):
    #         out_dir2 = os.path.join(out_dir, "4_clusters.cluster_{}".format(cluster))
    #         if not os.path.exists(out_dir2):
    #             os.makedirs(out_dir2)

    #         index = diff_regions[clusters == cluster]

    #         characterize_regions_function(
    #             analysis,
    #             df=analysis.coverage_annotated.loc[index, :],
    #             output_dir=out_dir2,
    #             prefix="cluster_{}".format(cluster)
    #         )

    # # for seven clusters
    # if cell_type in ["CD4"]:
    #     out_dir2 = os.path.join(out_dir, "7_clusters")
    #     if not os.path.exists(out_dir2):
    #         os.makedirs(out_dir2)
    #     clusters = fcluster(g.dendrogram_col.linkage, t=7, criterion='maxclust')

    #     # plot again with cluster identities
    #     g = sns.clustermap(
    #         analysis.accessibility.loc[diff_regions, sample_mask].T.astype(float),
    #         row_cluster=False,
    #         xticklabels=False,
    #         metric="correlation", z_score=1, vmin=-1, vmax=4, cbar_kws={"label": "Accesibility (Z-score)"},
    #         row_colors=color_df.iloc[1:, sample_mask].values,
    #         col_colors=[sns.color_palette("colorblind", 7)[i - 1] for i in clusters],
    #         figsize=(12, 12), rasterized=True,
    #     )
    #     g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    #     g.savefig(os.path.join(out_dir2, "clustermap.sorted.labeled.svg"), dpi=300)

    #     for cluster in np.unique(clusters):
    #         out_dir2 = os.path.join(out_dir, "7_clusters.cluster_{}".format(cluster))
    #         if not os.path.exists(out_dir2):
    #             os.makedirs(out_dir2)

    #         index = diff_regions[clusters == cluster]

    #         characterize_regions_function(
    #             analysis,
    #             df=analysis.coverage_annotated.loc[index, :],
    #             output_dir=out_dir2,
    #             prefix="cluster_{}".format(cluster)
    #         )


# # Run Enrichr
# find results -not -path "results/_old*" -name "*.symbols.txt" \
# -exec sbatch -J enrichr.{} -o {}.enrichr.log -p shortq -c 1 --mem 4000 \
# --wrap "python ~/jobs/run_Enrichr.py --input-file {} --output-file {}.enrichr.csv" \;

# # Run LOLA
# find results -not -path "results/_old*" -name "*_regions.bed" \
# -exec sbatch -J lola.{} -o {}.lola.log -p shortq -c 8 --mem 24000 \
# --wrap "Rscript ~/jobs/run_LOLA.R {} ~/cll-time_course/results/cll-time_course_peak_set.bed hg19" \;

# # Run AME
# for F in `find results -not -path "results/_old*" -name "*_regions.fa"`
# do
# DIR=`dirname $F`
# echo $F $DIR
# sbatch -J "meme_ame.${F}" -o "${F}.meme_ame.log" -p shortq -c 1 --mem 4000 \
# --wrap \
# "fasta-dinucleotide-shuffle -c 1 -f "$F" > "$F".shuffled.fa; \
# ame --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 \
# --control "$F".shuffled.fa -o "$DIR" "$F" ~/resources/motifs/motif_databases/HUMAN/HOCOMOCOv10.meme"
# done


# Collect enrichments and plot
enrichments = pd.DataFrame()

for cell_type in cell_types:
    # for n_clust in [1, 2, 4, 7]:
    for n_clust in [1, 2]:
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
            enr["n_clust"] = n_clust
            enr["cluster"] = cluster
            enrichments = enrichments.append(enr, ignore_index=True)
enrichments.to_csv(os.path.join("results", "enrichment.all_cell_types.csv"), index=False)

for n_top in [10, 20]:
    for gene_set_library in enrichments['gene_set_library'].drop_duplicates():
        enr = enrichments[enrichments['gene_set_library'] == gene_set_library]

        # transform p-values
        enr.loc[:, "q_value"] = -np.log10(enr["p_value"])

        # get top N enriched per set
        enr.loc[:, "label"] = enr["cell_type"] + " " + enr["n_clust"].astype(str) + " " + enr["cluster"].astype(str)
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
    # for n_clust in [1, 2, 4, 7]:
    for n_clust in [1, 2]:
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
            enr["n_clust"] = n_clust
            enr["cluster"] = cluster
            lola = lola.append(enr, ignore_index=True)
lola.to_csv(os.path.join("results", "lola.all_cell_types.csv"), index=False)

for n_top in [10, 20]:
    #
    enr = lola.copy()
    enr["set_label"] = (
        enr[['description', 'cellType', 'tissue', 'antibody', 'treatment']]
        .astype(str)
        .apply(lambda x: " ".join([y for y in set(x) if not y in ['None', 'nan']]), axis=1)
    )

    # get top N enriched per set
    enr.loc[:, "label"] = enr["cell_type"] + " " + enr["n_clust"].astype(str) + " " + enr["cluster"].astype(str)
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



# # One step back:
# # just take the differential regions per cluster and measure overlaps
# import pyupset

# filenames = glob.glob(os.path.join("results", "enrichment.*", "*", "cluster_*_regions.bed"))
# diffs = dict()
# for filename in filenames:
#     df = pd.read_csv(filename, header=None, sep="\t")
#     diffs[filename.replace("results/enrichment.", "")] = df

# pyupset.plot(diffs)



# # One step back:
# # just take the differential genes per cluster and measure overlaps
# import pyupset

# filenames = glob.glob(os.path.join("results", "enrichment.*", "*", "cluster_*_genes.symbols.txt"))
# diffs = dict()
# for filename in filenames:
#     df = pd.read_csv(filename, header=None, names=['gene'])
#     diffs[filename.replace("results/enrichment.", "")] = df

# pyupset.plot(diffs)





# filenames = glob.glob(os.path.join("results", "enrichment.*", "*", "cluster_*_regions.bed"))
# diffs = pd.DataFrame()
# for filename in filenames:
#     df = pd.read_csv(filename, header=None, names=['region'])
#     df["cluster"] = filename.replace("results/enrichment.", "")
#     diffs = diffs.append(df, ignore_index=True)

# diffs["intercept"] = 1

# p = pd.pivot_table(diffs, index="region", columns="cluster", values="intercept", fill_value=0)

# q = pd.DataFrame()
# for col1 in p.columns:
#     for col2 in p.columns:
#         q.loc[col1, col2] = ((p[col1] == 1) & (p[col2] == 1)).sum()

# g = sns.clustermap(np.log(1 + q), figsize=(12, 12), cbar_kws={"label": "log(shared diff regions)"})
# g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
# g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)




# filenames = glob.glob(os.path.join("results", "enrichment.*", "*", "cluster_*_genes.symbols.txt"))
# diffs = pd.DataFrame()
# for filename in filenames:
#     df = pd.read_csv(filename, header=None, names=['gene'])
#     df["cluster"] = filename.replace("results/enrichment.", "")
#     diffs = diffs.append(df, ignore_index=True)

# diffs["intercept"] = 1

# p = pd.pivot_table(diffs, index="gene", columns="cluster", values="intercept", fill_value=0)

# q = pd.DataFrame()
# for col1 in p.columns:
#     for col2 in p.columns:
#         q.loc[col1, col2] = ((p[col1] == 1) & (p[col2] == 1)).sum()

# g = sns.clustermap(np.log(1 + q), figsize=(12, 12), cbar_kws={"label": "log(shared diff genes)"})
# g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
# g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)



#

#

#

#


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
