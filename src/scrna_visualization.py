import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

sns.set_style("white")
plt.rcParams['svg.fonttype'] = 'none'


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


# Read enrichments
diff_enrichment = pd.read_table("/scratch/lab_bock/shared/projects/cll-time_course/results/single_cell_RNA/13_3_Overtime_Together_tobit/All_EnrichR.tsv")
diff_enrichment['cell_type'] = diff_enrichment['grp'].str.split("_").apply(lambda x: x[0])
diff_enrichment['patient_id'] = diff_enrichment['grp'].str.split("_").apply(lambda x: x[1])
diff_enrichment['direction'] = diff_enrichment['grp'].str.split("_").apply(lambda x: x[-1])


# Combine p-values across patients with Fisher's method
combined_diff_enrichment = diff_enrichment.groupby(['cell_type', 'direction', "database", "category"])['pval'].apply(lambda x: scipy.stats.combine_pvalues(x, method="fisher")[1]).reset_index()


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


# Get Seurat-normalized gene expression from R
from rpy2.robjects import numpy2ri, pandas2ri
import rpy2.robjects as robjects
numpy2ri.activate()
pandas2ri.activate()

robjects.r('library("Seurat")')
_load = robjects.r('load')
_rownames = robjects.r('rownames')
_colnames = robjects.r('colnames')

_load("/home/arendeiro/cll-time_course/results/single_cell_RNA/10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.RData")
pbmc = robjects.r['pbmc']

# Get metadata
sample_id = np.asarray(pbmc.slots['data.info'][-2])
patient_id = pd.Series(sample_id).str.replace(r"\d", "").str.split("_").apply(lambda x: x[0]).values
timepoint = pd.Series(sample_id).str.split("_").apply(lambda x: x[-1]).str.replace("d", "").astype(int)

assigned_cell_type = pd.Series(np.asarray(pbmc.slots['data.info'][-1]), name='assigned_cell_type')
genes_expressed = pd.Series(np.asarray(pbmc.slots['data.info'][0]).astype(int), name='genes_expressed')
umis_detected = pd.Series(np.asarray(pbmc.slots['data.info'][1]).astype(int), name='umis_detected')

# Get expression matrix (processed)
data = pbmc.slots['data']
x = pdata = robjects.r('as.matrix(pbmc@data)')
indexes = data.slots['Dimnames']
gene_index = pd.Series(np.asarray(indexes[0]), name='gene_index')
cell_index = pd.Series(np.asarray(indexes[1]), name='cell_index')

expression = pd.DataFrame(x, index=gene_index, columns=cell_index)

# Computed info
var_genes = pd.Series(np.asarray(pbmc.slots['var.genes']), name='var_genes')
tsne_position = pd.DataFrame(np.asarray(pbmc.slots['tsne.rot']), index=[0, 1])


# Let's replot the T-SNE to see if all's well

from sklearn.preprocessing import LabelEncoder
cell_type_enc = LabelEncoder()
assigned_cell_type_enc = pd.Series(cell_type_enc.fit_transform(assigned_cell_type))
patient_enc = LabelEncoder()
patient_id_enc = pd.Series(patient_enc.fit_transform(patient_id))

n_factors = 4
fig, axis = plt.subplots(2, n_factors, figsize=(n_factors * 4, 4 * 2), tight_layout=True)
cbar = axis[0, 0].scatter(x=tsne_position[0, :], y=tsne_position[1, :], c=np.log10(umis_detected), cmap="inferno", alpha=0.2, s=2, rasterized=True)
plt.colorbar(cbar, ax=axis[1, 0])
axis[0, 0].set_title("UMIs")
cbar = axis[0, 1].scatter(x=tsne_position[0, :], y=tsne_position[1, :], c=assigned_cell_type_enc, cmap="Paired", alpha=0.2, s=2, rasterized=True)
plt.colorbar(cbar, ax=axis[1, 1])
axis[0, 1].set_title("Cell type")
cbar = axis[0, 2].scatter(x=tsne_position[0, :], y=tsne_position[1, :], c=patient_id_enc, cmap="viridis", alpha=0.2, s=2, rasterized=True)
plt.colorbar(cbar, ax=axis[1, 2])
axis[0, 2].set_title("Patient")
cbar = axis[0, 3].scatter(x=tsne_position[0, :], y=tsne_position[1, :], c=timepoint, cmap="inferno", alpha=0.2, s=2, rasterized=True)
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

    if marker not in gene_index:
        continue
    axis.flat[i].scatter(x=tsne_position[0, :], y=tsne_position[1, :], c=expression.loc[marker, :], cmap="inferno", alpha=0.2, s=2, rasterized=True)

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

    cell_mask = assigned_cell_type == cell_type
    

    norm = matplotlib.colors.Normalize(vmin=0, vmax=np.log10(umis_detected.loc[cell_mask.values]).max())
    c1 = plt.get_cmap("inferno")(norm(np.log10(umis_detected.loc[cell_mask.values])))
    norm = matplotlib.colors.Normalize(vmin=0, vmax=assigned_cell_type_enc.loc[cell_mask.values].max())
    c2 = plt.get_cmap("Paired")(norm(umis_detected.loc[cell_mask.values]))
    norm = matplotlib.colors.Normalize(vmin=0, vmax=patient_id_enc.loc[cell_mask.values].max())
    c3 = plt.get_cmap("viridis")(norm(patient_id_enc.loc[cell_mask.values]))
    norm = matplotlib.colors.Normalize(vmin=0, vmax=timepoint.loc[cell_mask.values].max())
    c4 = plt.get_cmap("inferno")(norm(timepoint.loc[cell_mask.values]))

    g = sns.clustermap(
        expression.loc[diff_genes, cell_mask.values],
        z_score=0, metric="correlation",
        col_colors=[c1, c2, c3, c4], xticklabels=False, rasterized=True, figsize=(12, 6),
        cbar_kws={"label": "Expression (Z-score)"}
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join("results", "scrna.diff_pathways.{}.gene_expression.clustermap.svg".format(cell_type)), dpi=200, bbox_inches="tight")

