import os
import random
import string

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from sklearn.preprocessing import LabelEncoder

sns.set_style("white")
plt.rcParams['svg.fonttype'] = 'none'


# random seed
SEED = int("".join(
    LabelEncoder()
    .fit(list(string.ascii_uppercase))
    .transform(list("BOCKLAB")).astype(str)))
random.seed(SEED)
np.random.seed(SEED)


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def seurat_rdata_to_pandas(rdata_file, object_name="pbmc", preloaded=False, data_type="normalized"):
    """
    Load an Rdata file storing an Seurat analysis object and
    extract the processed expression matrix and its metadata into pandas dataframes.

    `rdata_file`: path to RData
    `object_name`: name of Seurat object
    `preloaded`: Attempt to 
    `data_type`: one of "raw" or "normalized".
    """
    from rpy2.robjects import numpy2ri, pandas2ri
    import rpy2.robjects as robjects
    numpy2ri.activate()
    pandas2ri.activate()

    _load = robjects.r('load')
    if not preloaded:
        _load(rdata_file)
    pbmc = robjects.r[object_name]

    # Data
    if data_type == "raw":
        _data_type = "raw.data"
    elif data_type == "normalized":
        _data_type = "data"
    data = pbmc.slots[_data_type]
    x = robjects.r('as.matrix({}@{})'.format(object_name, _data_type))
    gene_index = pd.Series(robjects.r('rownames({}@{})'.format(object_name, _data_type)), name='gene_index')
    cell_index = pd.Series(robjects.r('colnames({}@{})'.format(object_name, _data_type)), name='cell_index')

    # Metadata
    if data_type == "raw":
        metadata = None
    elif data_type == "normalized":
        metadata = pd.DataFrame(
            [np.asarray(m) for m in pbmc.slots['data.info']],
            columns=cell_index, index=np.asarray(pbmc.slots['data.info'].names)).T

    return (
        pd.DataFrame(x, index=gene_index, columns=cell_index),
        metadata
    )


def seurat_tsne_to_pandas(rdata_file, object_name="pbmc", preloaded=True):
    """
    Load an Rdata file storing an Seurat analysis object and
    extract the processed expression matrix into a pandas dataframe.
    """
    from rpy2.robjects import numpy2ri, pandas2ri
    import rpy2.robjects as robjects
    numpy2ri.activate()
    pandas2ri.activate()

    if not preloaded:
        _load = robjects.r('load')
        _load(rdata_file)
    pbmc = robjects.r[object_name]

    # Data
    data = pbmc.slots['data']
    cell_index = pd.Series(np.asarray(data.slots['Dimnames'][1]), name='cell_index')
    return pd.DataFrame(np.asarray(pbmc.slots['tsne.rot']), index=[0, 1], columns=cell_index)


def seurat_variable_genes_to_pandas(rdata_file, object_name="pbmc", preloaded=True):
    """
    Load an Rdata file storing an Seurat analysis object and
    extract the processed expression matrix into a pandas dataframe.
    """
    from rpy2.robjects import numpy2ri, pandas2ri
    import rpy2.robjects as robjects
    numpy2ri.activate()
    pandas2ri.activate()

    if not preloaded:
        _load = robjects.r('load')
        _load(rdata_file)
    pbmc = robjects.r[object_name]

    return np.asarray(pbmc.slots['var.genes'])


def predict_cell_cycle_r(rdata_file, object_name="pbmc", preloaded=True):
    from rpy2.robjects import numpy2ri, pandas2ri
    import rpy2.robjects as robjects
    numpy2ri.activate()
    pandas2ri.activate()

    robjects.r("""
    require(scran)
    require(org.Hs.eg.db)

    ccMarkers <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
    anno <- select(org.Hs.eg.db, keys=rownames({object_name}@data), keytype="SYMBOL", column="ENSEMBL")
    ensembl <- anno$ENSEMBL[match(rownames({object_name}@data), anno$SYMBOL)]
    assignments <- cyclone(as.matrix({object_name}@data), ccMarkers, gene.names=ensembl, BPPARAM=MulticoreParam(workers=12), verbose=TRUE)

    write.csv(assignments$scores, "results/single_cell_RNA/10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.cell_cyle_scores.csv")
    write.csv(assignments$phases, "results/single_cell_RNA/10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.cell_cyle_assignments.csv")
    write.csv(assignments$normalized.scores, "results/single_cell_RNA/10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.cell_cyle_normscores.csv")

    """.format(object_name=object_name))


# Read enrichments
enrichment_table = os.path.join("results/single_cell_RNA/13_3_Overtime_Together_tobit/All_EnrichR.tsv")

diff_enrichment = pd.read_table(enrichment_table)
diff_enrichment['cell_type'] = diff_enrichment['grp'].str.split("_").apply(lambda x: x[0])
diff_enrichment['patient_id'] = diff_enrichment['grp'].str.split("_").apply(lambda x: x[1])
diff_enrichment['direction'] = diff_enrichment['grp'].str.split("_").apply(lambda x: x[-1])


# Combine p-values across patients with Fisher's method
combined_diff_enrichment = (
    diff_enrichment
    .groupby(['cell_type', 'direction', "database", "category"])
    ['pval'].apply(lambda x: scipy.stats.combine_pvalues(x, method="fisher")[1])
    .reset_index())


# plot as bars
for database in combined_diff_enrichment['database'].drop_duplicates():
    ct = combined_diff_enrichment.loc[
        (combined_diff_enrichment["database"] == database), "cell_type"].drop_duplicates()
    n_col = ct.shape[0]
    fig, axis = plt.subplots(1, n_col, figsize=(n_col * 6, 4), tight_layout=True)
    for i, cell_type in enumerate(ct):
        axis[i].set_title(cell_type)

        tmp = combined_diff_enrichment.loc[
            (combined_diff_enrichment["database"] == database) &
            (combined_diff_enrichment["cell_type"] == cell_type), :]
        if tmp.shape[0] == 0:
            continue

        # convert p-values into signed log p-values
        tmp.loc[:, 'log_pvalue'] = -np.log10(tmp.loc[:, 'pval'])
        tmp.loc[:, 'log_signed_pvalue'] = tmp.apply(lambda x: -x['log_pvalue'] if x['direction'] == 'down' else x['log_pvalue'], axis=1)
        tmp.loc[:, 'category'] = tmp.apply(lambda x: x['category'] + "_" if x['direction'] == 'down' else x['category'], axis=1)

        # sort
        tmp = tmp.sort_values("log_signed_pvalue", ascending=False)

        sns.barplot(y=tmp['log_signed_pvalue'], x=tmp['category'], orient="vertical", color=sns.color_palette("colorblind")[i], ax=axis[i])
        axis[i].set_xticklabels(axis[i].get_xticklabels(), rotation=90)

    sns.despine(fig)
    fig.savefig(os.path.join("results", "scrna.timepoint_comparison.enrichment.combined_pvalues.{}.bar.svg".format(database)), bbox_inches="tight")


# Get Seurat-normalized gene expression and metadata from R
rdata_file = os.path.join("results/single_cell_RNA/10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.RData")
expression, metadata = seurat_rdata_to_pandas(rdata_file, "pbmc")

metadata = metadata.rename(columns={"nGene": "genes_covered", "nUMI": "umis_detected", "cellType": "assigned_cell_type", "sample": "sample_id"})
metadata["patient_id"] = metadata['sample_id'].str.replace(r"\d", "").str.split("_").apply(lambda x: x[0])
metadata["timepoint"] = metadata['sample_id'].str.split("_").apply(lambda x: x[-1]).str.replace("d", "").astype(int)

metadata.to_csv(os.path.join("results", "single_cell_RNA", "10_Seurat", "allDataBest_NoDownSampling_noIGH", "allDataBest_NoDownSampling_noIGH.metadata.csv"), index=True)

# Additional computed info
tsne_position = seurat_tsne_to_pandas(rdata_file, "pbmc")
var_genes = seurat_variable_genes_to_pandas(rdata_file, "pbmc")


# Let's replot the T-SNE to see if all's well
cell_type_enc = LabelEncoder()
metadata['assigned_cell_type_enc'] = cell_type_enc.fit_transform(metadata['assigned_cell_type'])
patient_enc = LabelEncoder()
metadata['patient_id_enc'] = patient_enc.fit_transform(metadata['patient_id'])

n_factors = 4
fig, axis = plt.subplots(2, n_factors, figsize=(n_factors * 4, 4 * 2), tight_layout=True)
cbar = axis[0, 0].scatter(x=tsne_position.loc[0, :], y=tsne_position.loc[1, :], c=np.log10(metadata['umis_detected'].astype(float)), cmap="inferno", alpha=0.2, s=2, rasterized=True)
plt.colorbar(cbar, ax=axis[1, 0])
axis[0, 0].set_title("UMIs")
cbar = axis[0, 1].scatter(x=tsne_position.loc[0, :], y=tsne_position.loc[1, :], c=metadata['assigned_cell_type_enc'], cmap="Paired", alpha=0.2, s=2, rasterized=True)
plt.colorbar(cbar, ax=axis[1, 1])
axis[0, 1].set_title("Cell type")
cbar = axis[0, 2].scatter(x=tsne_position.loc[0, :], y=tsne_position.loc[1, :], c=metadata['patient_id_enc'], cmap="viridis", alpha=0.2, s=2, rasterized=True)
plt.colorbar(cbar, ax=axis[1, 2])
axis[0, 2].set_title("Patient")
cbar = axis[0, 3].scatter(x=tsne_position.loc[0, :], y=tsne_position.loc[1, :], c=metadata['timepoint'], cmap="inferno", alpha=0.2, s=2, rasterized=True)
plt.colorbar(cbar, ax=axis[1, 3])
axis[0, 3].set_title("Timepoint")

for ax in axis[0, :]:
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")

fig.savefig(os.path.join("results", "scrna.t-sne.4_factors.scatter.svg"), dpi=300, bbox_inches="tight")


# Marker expression
markers = pd.read_csv(os.path.join("metadata", "CellMarkers.csv"))
markers = markers[markers['Organism'] == 'human']

m = markers['GeneSymbol'].drop_duplicates()

n_row = n_col = int(np.ceil(np.sqrt(m.shape[0])))

fig, axis = plt.subplots(n_row, n_col, figsize=(n_col * 3, n_row * 3), tight_layout=True)
for i, marker in enumerate(m):
    axis.flat[i].set_title(markers.loc[markers['GeneSymbol'] == marker, "Marker"].squeeze())

    if marker not in expression.index:
        continue
    axis.flat[i].scatter(x=tsne_position.loc[0, :], y=tsne_position.loc[1, :], c=expression.loc[marker, :], cmap="inferno", alpha=0.2, s=2, rasterized=True)

for ax in fig.axes[-1, :]:
    ax.set_xlabel("t-SNE 1")
for ax in fig.axes[:, 0]:
    ax.set_ylabel("t-SNE 2")
fig.tight_layout()
fig.savefig(os.path.join("results", "scrna.t-sne.cell_type_markers.scatter.svg"), dpi=200, bbox_inches="tight")



# Let's now plot the pathways enriched in the changing genes as heatmaps
database = "NCI-Nature_2016"
top_pathways = 3

for cell_type in combined_diff_enrichment['cell_type'].drop_duplicates():
    print(cell_type)
    # Get top differential pathways for cell type
    paths = combined_diff_enrichment.loc[
        (combined_diff_enrichment['database'] == database) &
        (combined_diff_enrichment['cell_type'] == cell_type)
    ].sort_values("pval").head(top_pathways)["category"]

    if paths.shape[0] == 0:
        continue

    # Get diff genes from those pathways
    diff_genes = (
        diff_enrichment.loc[(diff_enrichment['cell_type'] == cell_type) & (diff_enrichment['category'].isin(paths.tolist())), 'genes']
        .str.split(",")
        .apply(pd.Series)
        .stack()
        .drop_duplicates().tolist())

    if cell_type == "Bcells":
        cell_type = "CLL"
    elif cell_type == "Tcell1":
        cell_type = "GDT"
    elif cell_type == "Monos":
        cell_type = "Mono"
    elif cell_type == "NKcells":
        cell_type = "NK"

    cell_mask = metadata['assigned_cell_type'] == cell_type

    norm = matplotlib.colors.Normalize(vmin=0, vmax=np.log10(metadata['umis_detected'].astype(float).loc[cell_mask.values]).max())
    c1 = plt.get_cmap("inferno")(norm(np.log10(metadata['umis_detected'].astype(float).loc[cell_mask.values])))
    norm = matplotlib.colors.Normalize(vmin=0, vmax=metadata['assigned_cell_type_enc'].loc[cell_mask.values].max())
    c2 = plt.get_cmap("Paired")(norm(metadata['patient_id_enc'].astype(float).loc[cell_mask.values]))
    norm = matplotlib.colors.Normalize(vmin=0, vmax=metadata['patient_id_enc'].loc[cell_mask.values].max())
    c3 = plt.get_cmap("viridis")(norm(metadata['patient_id_enc'].loc[cell_mask.values]))
    norm = matplotlib.colors.Normalize(vmin=0, vmax=metadata['timepoint'].loc[cell_mask.values].max())
    c4 = plt.get_cmap("inferno")(norm(metadata['timepoint'].loc[cell_mask.values]))

    g = sns.clustermap(
        expression.loc[diff_genes, cell_mask.values],
        z_score=0, metric="correlation",
        col_colors=[c1, c2, c3, c4], xticklabels=False, rasterized=True, figsize=(12, 6),
        cbar_kws={"label": "Expression (Z-score)"}
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join("results", "scrna.diff_pathways.{}.gene_expression.clustermap.svg".format(cell_type)), dpi=200, bbox_inches="tight")

    # Group cells by patient, timepoint

    # (maybe worth showing all cell types to see that it is specific/shared)

    mean_expression = expression.loc[diff_genes, cell_mask].T.join(metadata).groupby(['patient_id', 'timepoint']).mean().loc[:, diff_genes].T

    g = sns.clustermap(
        mean_expression,
        metric="correlation", z_score=1,
        xticklabels=True, rasterized=True, figsize=(12, 6),
        cbar_kws={"label": "Expression (Z-score)"}
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join("results", "scrna.diff_pathways.{}.gene_expression.mean_patient_timepoint.clustermap.svg".format(cell_type)), dpi=200, bbox_inches="tight")

    for square in [True, False]:
        g = sns.clustermap(
            mean_expression.sort_index(1, 'patient_id'), col_cluster=False, square=square,
            metric="correlation", z_score=0,
            xticklabels=True, rasterized=True, figsize=(12, 6),
            cbar_kws={"label": "Expression (Z-score)"}
        )
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.savefig(os.path.join("results", "scrna.diff_pathways.{}.gene_expression.mean_patient_timepoint.sorted_patient.clustermap{}.svg".format(cell_type, square)), dpi=200, bbox_inches="tight")

        g = sns.clustermap(
            mean_expression.sort_index(1, 'timepoint'), col_cluster=False, square=square,
            metric="correlation", z_score=0,
            xticklabels=True, rasterized=True, figsize=(12, 6),
            cbar_kws={"label": "Expression (Z-score)"}
        )
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.savefig(os.path.join("results", "scrna.diff_pathways.{}.gene_expression.mean_patient_timepoint.sorted_timepoint.clustermap{}.svg".format(cell_type, square)), dpi=200, bbox_inches="tight")

# Let's plot the expression of the individual genes changing the most
database = "NCI-Nature_2016"
top_pathways = 3





# Let's observe the distributions of cell cycle assignments

# cell cycle phase prediction probabilities from old analysis
old_data, old_metadata = seurat_rdata_to_pandas("results/single_cell_RNA/10_Seurat/allDataBest_NoDownSampling_noIGH_RmSample_cellCycle/allDataBest_NoDownSampling_noIGH_RmSample_cellCycle.RData")
cycle_assignment = pd.read_csv(os.path.join("results", "single_cell_RNA", "10_Seurat", "allDataBest_NoDownSampling_noIGH_RmSample_cellCycle", "Cluster.tsv"), sep="\t")
cycle_assignment.index = old_metadata.index
metadata = metadata.join(cycle_assignment[['G1', 'S', 'G2M']])
metadata.to_csv(os.path.join("results", "single_cell_RNA", "10_Seurat", "allDataBest_NoDownSampling_noIGH", "allDataBest_NoDownSampling_noIGH.metadata.csv"), index=True)


# let's try also with cells which have high separation between G1 and G2/M
sel_cells = metadata[(metadata["G2M"] - metadata['G1']).abs() >= 0.75].index

for met, prefix in [(metadata, "all_cells"), (metadata.loc[sel_cells], "filtered_cells")]:
    # all cells together
    fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4 * 1))
    c = axis[0].scatter(met["G1"], met["G2M"], c=met["S"], cmap="BrBG", alpha=0.2, s=2, rasterized=True)
    plt.colorbar(c, ax=axis[1], label="S phase probability")
    axis[0].set_xlabel("G1 phase probability")
    axis[0].set_ylabel("G2/M phase probability")
    fig.savefig(os.path.join("results", "scrna.cell_cycle_assignment.{}.all_cells.scatter.svg".format(prefix)), dpi=200, bbox_inches="tight")

    # per cell type, per timepoint
    times = met['timepoint'].drop_duplicates().sort_values()
    times = times[times < 180]
    cell_types = met['assigned_cell_type'].drop_duplicates().sort_values()
    cell_types = cell_types[cell_types != "NA"]

    fig, axis = plt.subplots(len(times), len(cell_types), figsize=(2 * len(cell_types), 2 * len(times)), sharex=True, sharey=True)
    for i, timepoint in enumerate(times):
        for j, cell_type in enumerate(cell_types):
            sel_met = met[(met['timepoint'] == timepoint) & (met['assigned_cell_type'] == cell_type)]
            axis[i, j].scatter(sel_met["G1"], sel_met["G2M"], c=sel_met["S"], cmap="BrBG", alpha=0.5, s=2, rasterized=True)
    for i, ax in enumerate(axis[:, 0]):
        ax.set_ylabel(times[i])
    for j, ax in enumerate(axis[-1, :]):
        ax.set_xlabel(cell_types[j])
    fig.savefig(os.path.join("results", "scrna.cell_cycle_assignment.{}.per_patient_timepoint.scatter.svg".format(prefix)), dpi=200, bbox_inches="tight")

    # Mean values
    grouped_cycle = pd.melt(
        met.groupby(['patient_id', 'assigned_cell_type', 'timepoint'])['G1', 'S', 'G2M'].mean().reset_index(),
        id_vars=['patient_id', 'assigned_cell_type', 'timepoint'], var_name='phase', value_name='probability')
    grouped_cycle = grouped_cycle[grouped_cycle['timepoint'] < 180]

    g = sns.factorplot(data=grouped_cycle.reset_index(), x='phase', y='probability', hue='timepoint', col='assigned_cell_type', kind='bar', size=3)
    g.savefig(os.path.join("results", "scrna.cell_cycle_assignment.{}.mean.joint_patients.bar.svg".format(prefix)), dpi=200, bbox_inches="tight")

    try:
        g = sns.factorplot(data=grouped_cycle.reset_index(), x='timepoint', y='probability', hue='phase', col='assigned_cell_type')
        g.savefig(os.path.join("results", "scrna.cell_cycle_assignment.{}.mean.joint_patients.time_line.svg".format(prefix)), dpi=200, bbox_inches="tight")
    except:
        pass

    g = sns.factorplot(data=grouped_cycle.reset_index(), x='phase', y='probability', hue='timepoint', col='assigned_cell_type')
    g.savefig(os.path.join("results", "scrna.cell_cycle_assignment.{}.mean.joint_patients.line.svg".format(prefix)), dpi=200, bbox_inches="tight")

    g = sns.factorplot(data=grouped_cycle.reset_index(), x='phase', y='probability', hue='timepoint', col='assigned_cell_type', row='patient_id')
    g.savefig(os.path.join("results", "scrna.cell_cycle_assignment.{}.mean.per_patient.line.svg".format(prefix)), dpi=200, bbox_inches="tight")


# Let's play with the cells from one patient where we have 3 timepoints
patient_id = metadata.groupby(['patient_id'])['sample_id'].unique().apply(len).argmax()

sel_cells = metadata.loc[
    (metadata['patient_id'] == patient_id) &
    (metadata['assigned_cell_type'] == "CLL")
].index



# t-SNE inset only with these cells
fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4 * 1), tight_layout=True)

norm = matplotlib.colors.Normalize(vmin=0, vmax=np.log10(metadata['umis_detected'].astype(float).loc[cell_mask.values]).max())
c1 = plt.get_cmap("inferno")(norm(np.log10(metadata.loc[sel_cells, 'umis_detected'].astype(float))))
cbar = axis[0].scatter(x=tsne_position.loc[0, sel_cells], y=tsne_position.loc[1, sel_cells], c=c1, alpha=0.5, s=2, rasterized=True)
# plt.colorbar(cbar, ax=axis[1, 0])
axis[0].set_title("UMIs")

norm = matplotlib.colors.Normalize(vmin=0, vmax=metadata['timepoint'].max())
c2 = plt.get_cmap("inferno")(norm(metadata.loc[sel_cells, 'timepoint'].astype(float)))
cbar = axis[1].scatter(x=tsne_position.loc[0, sel_cells], y=tsne_position.loc[1, sel_cells], c=c2, alpha=0.5, s=2, rasterized=True)
# plt.colorbar(cbar, ax=axis[1, 1])
axis[1].set_title("Timepoint")

for ax in axis[:]:
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
fig.savefig(os.path.join("results", "scrna.t-sne.patient_{}_only.scatter.svg".format(patient_id)), dpi=300, bbox_inches="tight")




# PCA transform

from sklearn.decomposition import PCA
pca = PCA()

fit = pca.fit_transform(expression.loc[:, sel_cells].T)
fit = pca.fit_transform(expression.loc[var_genes, sel_cells].T)


plt.plot(range(len(pca.explained_variance_ratio_)), pca.explained_variance_ratio_)

plt.scatter(
    fit[0], fit[1],
    alpha=0.2, s=2, c=plt.get_cmap("inferno")(metadata.loc[sel_cells, "timepoint"]), cmap="inferno"
)


###
# Let's try a different normalization
rdata_file = "results/single_cell_RNA/10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.RData"

# load up raw gene expression
expression, _ = seurat_rdata_to_pandas(rdata_file, object_name="pbmc", preloaded=False, data_type="raw")

# Get cell's metadata
norm_e, metadata = seurat_rdata_to_pandas(rdata_file, object_name="pbmc", preloaded=True, data_type="normalized")
# Filter out cells not present in the Seurat normalized data
expression = expression.loc[:, metadata.index]

# Filter out genes
# not expressed
expression = expression.loc[~(expression.sum(axis=1) == 0)]

# ribossomal/mitochondrial genes
expression = expression.loc[
    (~expression.index.str.startswith("RP-")) &
    (~expression.index.str.contains("^RP\d+-")) &
    (~expression.index.str.contains("^RPS")) &
    (~expression.index.str.contains("^RPL")) &
    (~expression.index.str.startswith("MT-")), :]

# filter lowly expressed genes
# (not covered in at least x% cells)
perc_cells = 5.0
expression = expression.loc[(expression != 0).sum(axis=1) > (expression.shape[1] / perc_cells), :]

# normalize to TPM
expression_norm = (expression / expression.sum()) * 1e3
expression_norm = np.log2(1 + expression_norm)


expression.to_csv("scrna-seq.expression_matrix.csv.gz")
metadata.to_csv("scrna-seq.metadata.csv.gz")


"""
library("data.table")
library("MAST")


df = data.table::fread("scrna-seq.expression_matrix.csv.gz", sep=',', header=TRUE, data.table=FALSE)
rownames(df) = df$gene_index
df <- df[, 2:dim(df)[2]]
expression = log2(((df / colSums(df)) * 1e3) + 1)

metadata = data.table::fread("scrna-seq.metadata.csv.gz", sep=',', header=TRUE, data.table=FALSE)
rownames(metadata) = metadata$cell_index

fdat = data.frame(rownames(df))
colnames(fdat) <- c("gene_name")

scaRaw <- FromMatrix(expression, metadata, fdat)

"""


# Plot example genes from distribution of expression
n = 9
w = h = int(np.sqrt(n))
fig, axis = plt.subplots(h, w + 1, figsize=((w + 1) * 4, h * 4))
axis = axis.flatten()
axis[0].set_xlabel("Mean expression (log(TPM))")
axis[0].set_ylabel("Cell count")
sns.distplot(expression_norm.mean(axis=0), kde=False, ax=axis[0])
axis[1].set_xlabel("Mean expression (log(TPM))")
axis[1].set_ylabel("Gene count")
sns.distplot(expression_norm.mean(axis=1), kde=False, ax=axis[1])
axis[2].set_xlabel("Mean expression (log(TPM))")
axis[2].set_ylabel("Coeffient of deviation (STD/MEAN)")
axis[2].scatter(x=expression_norm.mean(axis=1), y=(expression_norm.std(axis=1) / expression_norm.mean(axis=1)), alpha=0.1, s=5, rasterized=True)

# pick genes from top 5% expression, middle 50%
gene_mean = expression_norm.mean(axis=1).sort_values()
m = int(gene_mean.shape[0] / 2.)
example_genes = (
    gene_mean.head(3).index.tolist() +
    gene_mean.iloc[[m, m + 1, m + 2]].index.tolist() + 
    gene_mean.tail(3).index.tolist())

for i, gene in enumerate(example_genes):
    axis[3 + i].set_title(gene)
    sns.distplot(expression_norm.loc[gene], kde=False, ax=axis[3 + i])
    axis[3 + i].set_xlabel("Expression (log(TPM))")
    axis[3 + i].set_ylabel("Cell count")
sns.despine(fig)
fig.savefig(os.path.join("results", "cll-time_course.single_cell.all_samples.no_bad_genes.tpm_norm.dist_and_example_genes.svg"), dpi=300, bbox_inches="tight")


# Label matrix with metadata (Multiindex columns)
metadata = metadata.loc[expression_norm.columns]

meta = metadata[["nGene", "nUMI", "sample", "cellType"]]
meta.columns = ["n_genes", "n_umis", "sample", "cell_type"]
meta['patient_id'] = meta['sample'].str.split("_").apply(lambda x: x[0])
meta['patient_id'] = meta['patient_id'].str.replace("\d+", "")
meta['timepoint'] = meta['sample'].str.split("_").apply(lambda x: int(x[-1].replace("d", "")))


expression_norm.columns = pd.MultiIndex.from_arrays(meta.T.values, names=meta.columns)

# See expression of Cytokines/Receptors with time
cytokine_genes = pd.read_csv(os.path.join("results", "cll-time_course.ligand-receptor_repertoire.CLL.gene_level.sig_only.timepoint_mean.clustermap.csv"), index_col=0)

g = sns.clustermap(expression_norm.loc[cytokine_genes.index].dropna().T.groupby(level=['cell_type', 'timepoint']).mean(), row_cluster=False)#, z_score=1)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")



g = sns.clustermap(
    expression_norm.loc[
        (expression_norm.index.str.contains("^CD\d+")) |
        (expression_norm.index.str.contains("^IL\d+"))].T.groupby(level=['cell_type', 'timepoint']).mean(),
    row_cluster=False)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")

e = expression_norm.T.groupby(level=['cell_type', 'timepoint']).mean()
e = e.loc[e.index.get_level_values("cell_type") == "CLL"]

e = e.loc[:, e.sum() != 0]
g = sns.clustermap(
    e.T.loc[
        (e.columns.str.contains("^TNF")) |
        (e.columns.str.contains("^WNT")) |
        (e.columns.str.contains("^TGF")), :],
    col_cluster=False, metric="correlation", z_score=0)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")



sns.violinplot(data=expression_norm.loc[["CD44"], expression_norm.columns.get_level_values("cell_type") == 'CLL'].T.reset_index(), x='timepoint', y ='CD44')
sns.swarmplot(data=expression_norm.loc[["CD44"], expression_norm.columns.get_level_values("cell_type") == 'CLL'].T.reset_index(), x='timepoint', y ='CD44')


# Differential expression
from ngs_toolkit.general import least_squares_fit
import gseapy

gene_set_libraries = [
    'GO_Biological_Process_2017b',
    # 'GO_Cellular_Component_2017b',
    'GO_Molecular_Function_2017b',
    "ChEA_2016",
    "KEGG_2016",
    # "ESCAPE",
    # "Epigenomics_Roadmap_HM_ChIP-seq",
    "ENCODE_TF_ChIP-seq_2015",
    "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
    # "ENCODE_Histone_Modifications_2015",
    # "OMIM_Expanded",
    # "TF-LOF_Expression_from_GEO",
    # "Single_Gene_Perturbations_from_GEO_down",
    # "Single_Gene_Perturbations_from_GEO_up",
    # "Disease_Perturbations_from_GEO_down",
    # "Disease_Perturbations_from_GEO_up",
    # "Drug_Perturbations_from_GEO_down",
    # "Drug_Perturbations_from_GEO_up",
    "WikiPathways_2016",
    "Reactome_2016",
    "BioCarta_2016",
    "NCI-Nature_2016"
]

results = pd.DataFrame()
enrichments = pd.DataFrame()
for cell_type in expression_norm.columns.get_level_values("cell_type").unique():
    print(cell_type)
    res = least_squares_fit(
        quant_matrix=expression_norm.loc[:,
            (expression_norm.columns.get_level_values("cell_type") == cell_type) &
            (expression_norm.columns.get_level_values("timepoint") <= 120)].T,
        design_matrix=meta[
            (meta["cell_type"] == cell_type) &
            (meta["timepoint"] <= 120)
        ],
        test_model="~ patient_id + timepoint", null_model="~ patient_id", standardize_data=True, multiple_correction_method="fdr_bh")
    res["cell_type"] = cell_type
    results = results.append(res)

    all_genes = res[res['q_value'] < 0.01].index.tolist()
    up_genes = res[(res['q_value'] < 0.01) & (res['timepoint'] > 0)].index.tolist()
    down_genes = res[(res['q_value'] < 0.01) & (res['timepoint'] < 0)].index.tolist()

    for name, gene_set in [("all_genes", all_genes), ("up_genes", up_genes), ("down_genes", down_genes)]:
        if len(gene_set) == 0:
            continue
        print(name)
        for gene_set_library in gene_set_libraries:
            print(gene_set_library)
            enr = gseapy.enrichr(
                    gene_list=gene_set, description='test_name', gene_sets=gene_set_library,
                    outdir='enrichr_kegg', cutoff=0.5, no_plot=False).res2d
            enr["gene_set_library"] = gene_set_library
            enr["gene_set"] = name
            enr["cell_type"] = cell_type
            enrichments = enrichments.append(enr, ignore_index=True)

results.to_csv(os.path.join("results", "cll-time_course.single_cell.all_samples.no_bad_genes.tpm_norm.differential_expression.csv"), index=True)
enrichments.to_csv(os.path.join("results", "cll-time_course.single_cell.all_samples.no_bad_genes.tpm_norm.differential_expression.enrichments.csv"), index=True)


enrichments = pd.read_csv(os.path.join("results", "cll-time_course.single_cell.all_samples.no_bad_genes.tpm_norm.differential_expression.enrichments.csv"), index_col=0)


top_n = 20

for gene_set_library in gene_set_libraries:
    enr = enrichments[enrichments['gene_set_library'] == gene_set_library]
    enr = enr[enr["gene_set"] != "all_genes"]

    enr['label'] = enr["cell_type"] + " " + enr["gene_set"]
    piv = -np.log10(pd.pivot_table(data=enr, index="Term", columns="label", values="Adjusted P-value", fill_value=1))

    terms = enr.set_index("Term").groupby(["cell_type", "gene_set"])["Adjusted P-value"].nsmallest(top_n).reset_index(level='Term')['Term'].unique()
    print(gene_set_library, terms.shape[0])

    piv = piv.apply(scipy.stats.zscore, axis=0)

    g = sns.clustermap(piv.loc[terms], figsize=(piv.shape[1] * 0.2, terms.shape[0] * 0.12), rasterized=True, cbar_kws={"label": "-log10(p-value)"})
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.savefig(os.path.join("results", "cll-time_course.single_cell.all_samples.no_bad_genes.tpm_norm.differential_expression.enrichments.{}.svg".format(gene_set_library)), dpi=300, bbox_inches="tight")

