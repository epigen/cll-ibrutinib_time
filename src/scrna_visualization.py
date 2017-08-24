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


def seurat_rdata_to_pandas(rdata_file, object_name="pbmc", preloaded=False):
    """
    Load an Rdata file storing an Seurat analysis object and
    extract the processed expression matrix and its metadata into pandas dataframes.
    """
    from rpy2.robjects import numpy2ri, pandas2ri
    import rpy2.robjects as robjects
    numpy2ri.activate()
    pandas2ri.activate()

    _load = robjects.r('load')
    _load(rdata_file)
    pbmc = robjects.r[object_name]

    # Data
    data = pbmc.slots['data']
    x = robjects.r('as.matrix(pbmc@data)')
    indexes = data.slots['Dimnames']
    gene_index = pd.Series(np.asarray(indexes[0]), name='gene_index')
    cell_index = pd.Series(np.asarray(indexes[1]), name='cell_index')

    # Metadata
    metadata_names = np.asarray(pbmc.slots['data.info'].names)
    metadata = [np.asarray(m) for m in pbmc.slots['data.info']]

    return (
        pd.DataFrame(x, index=gene_index, columns=cell_index),
        pd.DataFrame(metadata, columns=cell_index, index=metadata_names).T
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
    c2 = plt.get_cmap("Paired")(norm(metadata['umis_detected'].astype(float).loc[cell_mask.values]))
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
