#!/usr/bin/env python


from collections import OrderedDict
import os
import random
import string

import matplotlib
import matplotlib.pyplot as plt
from natsort import natsorted as sorted
import numpy as np
import pandas as pd
from peppy import Project
import scanpy as sc
import scipy
import seaborn as sns
from sklearn.preprocessing import LabelEncoder

from natsort import natsorted

# from sklearn.metrics import silhouette_score
# from ngs_toolkit.analysis import Analysis
# from dca.utils import plot_mean_dropout


sc.set_figure_params(format="svg", dpi_save=300, vector_friendly=True)
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


def clustermap_fix_label_orientation(grid, fontsize="xx-small", **kwargs):
    grid.ax_heatmap.set_xticklabels(
        grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize=fontsize, **kwargs)
    grid.ax_heatmap.set_yticklabels(
        grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize=fontsize, **kwargs)


def clustermap_rasterize_heatmap(grid):
    import matplotlib
    q = [x for x in grid.ax_heatmap.get_children()
         if isinstance(x, matplotlib.collections.QuadMesh)][0]
    q.set_rasterized(True)


def clustermap_rasterize_dendrogram(grid, axis=[0, 1]):
    import matplotlib

    if isinstance(axis, int):
        axis = [axis]

    for ax in axis:
        if ax == 1:
            d = grid.ax_col_dendrogram
        elif ax == 0:
            d = grid.ax_row_dendrogram
        q = [x for x in d.get_children()
             if isinstance(x, matplotlib.collections.LineCollection)]
        if len(q) == 1:
            q = q[0]
            q.set_rasterized(True)


def savefig(fig, file_name, **kwargs):
    if isinstance(fig, sns.axisgrid.Grid):
        fig = fig.fig
    default_kwargs = {"dpi": 300, "bbox_inches": "tight"}
    fig.savefig(file_name, **default_kwargs)


def add_cell_type(data):
    metadata = pd.read_csv(
        os.path.join(
            "results", "single_cell_RNA", "10_Seurat_raw",
            "inclDay30_noIGHLK_negbinom", "inclDay30_noIGHLK.metadata.csv"))
    metadata.loc[:, 'sample_id'] = (
        metadata['sample_id']
        .str.replace("FE", "CLL1")
        .str.replace("PBGY", "CLL5")
        .str.replace("PT", "CLL6")
        .str.replace("VZS", "CLL8"))

    for sample in data.obs['sample'].unique():
        # q = "CLL5_d150" if sample == "CLL5_d120" else sample
        q = prj.samples[sample].name.replace("10X-scRNA-seq_", "")
        q = "CLL5_d150" if q == "CLL5_d120" else q
        mt = metadata.loc[metadata['sample_id'] == q, :].drop("timepoint", axis=1)
        c = data.obs.loc[data.obs['sample'] == sample, :]
        c.loc[:, 'new_index'] = c.index.str.extract(r"(.*?)-.*").values
        mt['cell_index'] = mt['cell_index'].str.extract(r"(.*?)-.*")
        data.obs.loc[data.obs['sample'] == sample, 'cell_type'] = (
            c[['new_index']].reset_index().merge(
                mt[['cell_index', 'CellType']],
                how="left", left_on="new_index", right_on="cell_index")
            .set_index("index")
            .loc[c.index, "CellType"])


def get_colors(obs, factors, as_df=True):
    v_colors = dict()
    for f in factors:
        v = obs.loc[:, f]
        d = OrderedDict(zip(sorted(v.unique()), colors[f]))
        v_colors[f] = pd.Series(v.map(d).values, index=v.index)
    return v_colors if not as_df else pd.DataFrame(v_colors)


def natsort_dataframe(df, by):
    from natsort import index_natsorted, order_by_index
    if not isinstance(by, list):
        by = [by]
    return df.reindex(index=order_by_index(df.index, index_natsorted(zip(*[df[x] for x in by]))))


mark1 = ['CD3D', 'CD3G', 'CD4', 'CD8A', 'NKG7', 'CD14', "CST3", "CD79A", "TCL1A", "CD19", 'GZMB', 'CD68', 'CD27', 'CD20', 'CD24', "NCR1"]
mark2 = ["IL32", "IFNG", "IFNGR1", "IFNGR1", "IL4R", "IL4", "JUN", "JUNB", "JUND", "JAK1", "JAK2", "GATA1", "JARID2", "KRAS", "MYC"]
mark3 = ["BTK", "LCK", "E2F4", "CXCR4", "ITGA4", "HBA1", "PTPRC"]
red_mark = ["CST3", "TCL1A", "GZMB", "NKG7", "CD3D"]
pbmc_markers = mark1 + mark2 + mark3

colors = {
    "cell_type": ["#b7cee3", "#c4e6a7", "#f6c4c2", "#f7cc82", "#fdfe97", "#cdbbdc", "#cdbbdc"],
    "patient_id": ["#84190d", "#8b8c00", "#318b00", "#5d5d5d"],
    "timepoint": ["#101013", "#67257d", "#ef9d00", "#fafca5"],
    "batch": sns.color_palette()[:4]}

url = "https://raw.githubusercontent.com/broadinstitute/inferCNV/master/example/gencode_downsampled.txt"
gene_order = pd.read_csv(url, sep="\t", header=None)
gene_order.columns = ['gene', 'chr', 'start', 'end']

prj = Project(os.path.join("metadata", "project_config.scRNA-seq.yaml"))
data_type = "scRNA-seq"


# Read data
prefix = "cll-time_course-scRNA-seq.all_samples.250-50_filter"
adata = sc.read(prefix + ".dca_denoised-zinb.processed.h5ad")
adata.obs = pd.read_csv(prefix + ".dca_denoised-zinb.processed.obs.csv", index_col=0)


# CNV analysis
C = pd.read_csv("RdBu_r.custom_colormap.csv", header=None)
cmap = matplotlib.colors.ListedColormap(C.values / 255.0)

g = gene_order.loc[gene_order['gene'].isin(adata.var.index.tolist()), 'gene'].squeeze().tolist()
q = adata[:, g].X
if isinstance(q, scipy.sparse.csr_matrix):
    q = q.to_dense()
q = pd.DataFrame(q, index=adata.obs.index, columns=g).T

# cap expression?

# remove cell median
q2 = q - q.median()
# standardize matrix (where parameters are from non-leukemic cell types)
# normal = q2.loc[:, adata.obs.cell_type != "CLL"]
q2 = ((q2.T - q2.median(1)) / q2.std(1)).T

# rolling mean per chromosome
rm = list()
for chrom in gene_order['chr'].unique():
    print(chrom)
    chr_genes = gene_order.loc[gene_order['chr'] == chrom, 'gene'].tolist()
    rm.append(q2.reindex(chr_genes).dropna().rolling(window=50).mean())
rm = pd.concat(rm)
# Subtract by mean of whole matrix (try to center on "normality")
# normal_rm = rm.loc[:, adata.obs.cell_type != "CLL"]
p = rm.dropna().T - rm.mean().mean()

# filter noise out / smooth
p = p ** 3 * 3

# cll = adata.obs.loc[adata.obs['cell_type'] == "CLL"].index
# rm_cll = rm.loc[:, cll]
# chr_genes = gene_order.loc[gene_order['chr'].isin(['chr6', 'chr11', 'chr12', 'chr13', 'chr17']), 'gene']
# p = rm_cll.loc[:, np.random.choice(rm_cll.columns, 5000)].dropna()


p_cll = p.loc[adata.obs.cell_type == "CLL", :]

chosen = list()
for _p in adata.obs.patient_id.unique():
    chosen += np.random.choice(p.loc[adata.obs.patient_id == _p, :].index, 2000).tolist()
p_cll = p_cll.loc[chosen, :].dropna()

g = sns.clustermap(
    p_cll,
    col_cluster=False, xticklabels=False, yticklabels=False,
    cmap="RdBu_r", center=0, robust=True, vmin=-2, vmax=2,
    row_colors=get_colors(p_cll.obs, factors=["patient_id", "timepoint"]),
    cbar_kws={"label": "log difference from normal"},
    metric="correlation",
    rasterized=True)
g.ax_heatmap.set_ylabel("Cells")
g.ax_heatmap.set_xlabel("Chromosome positions")
clustermap_rasterize_heatmap(g)
clustermap_fix_label_orientation(g)
clustermap_rasterize_dendrogram(g)
savefig(g, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.whole_genome.2kcells_per_patient.RdBu_r.svg"))

p_cll = p_cll.loc[adata.obs.sort_values(['patient_id', 'timepoint']).index, :].dropna()
g = sns.clustermap(
    p_cll,
    row_cluster=False, col_cluster=False, xticklabels=False, yticklabels=False,
    cmap="RdBu_r", center=0, robust=True, vmin=-2, vmax=2,
    row_colors=get_colors(p_cll.obs, factors=["patient_id", "timepoint"]),
    cbar_kws={"label": "log difference from normal"},
    metric="correlation",
    rasterized=True)
g.ax_heatmap.set_ylabel("Cells")
g.ax_heatmap.set_xlabel("Chromosome positions")
clustermap_rasterize_heatmap(g)
clustermap_fix_label_orientation(g)
clustermap_rasterize_dendrogram(g)
savefig(g, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.whole_genome.2kcells_per_patient.ordered.RdBu_r.svg"))


v = ['patient_id', 'timepoint', 'cell_type']
p2 = p.join(adata.obs[v])
p2 = p2.groupby(v).mean().dropna()
p2 = p2.loc[~p2.index.get_level_values("cell_type").isin(["nan", "NurseLikeCell"]), :]

# p2_normal = p2.loc[p2.index.get_level_values("cell_type") != "CLL"]
p2 = ((p2.T - p2.mean(1)) / p2.std(1)).T
p2 = p2 - p2.mean()


g = sns.clustermap(
    p2 - p2.mean().mean() / 2,
    col_cluster=False,
    xticklabels=False, yticklabels=True,
    cmap="RdBu_r", center=0,  # robust=True,
    figsize=(4, 8),
    metric="correlation",
    cbar_kws={"label": "log difference from normal"},
    rasterized=True)
clustermap_rasterize_heatmap(g)
clustermap_fix_label_orientation(g)
clustermap_rasterize_dendrogram(g)
savefig(g, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.whole_genome.grouped.svg"))

cnv = sc.AnnData(p)
cnv.obs = adata.obs
sc.pp.pca(cnv)
sc.pp.neighbors(cnv)
sc.tl.umap(cnv)

cnv.obs.to_csv(prefix + ".dca_denoised-zinb.batch_combat.processed.cnv.processed.obs.csv")
cnv.obs = pd.DataFrame(index=cnv.obs.index)
sc.write(prefix + ".dca_denoised-zinb.batch_combat.processed.cnv.processed.h5ad", cnv)

prefix = "cll-time_course-scRNA-seq.all_samples.250-50_filter"
cnv = sc.read(prefix + ".dca_denoised-zinb.batch_combat.processed.cnv.processed.h5ad")
cnv.obs = pd.read_csv(prefix + ".dca_denoised-zinb.batch_combat.processed.cnv.processed.obs.csv", index_col=0)

c = (
    pd.DataFrame(cnv.X, index=cnv.obs.index, columns=cnv.var.index)
    .T
    .join(gene_order.set_index("gene")))

c = natsort_dataframe(c, ['chr', 'start'])

c_cll = c.loc[:, cnv.obs.loc[cnv.obs['cell_type'] == "CLL", :].index]

chosen_cells = list()
for patient_id in cnv.obs.patient_id.unique():
    chosen_cells += cnv.obs.loc[c_cll.columns].loc[cnv.obs['patient_id'] == patient_id].sample(n=2500).index.tolist()
# r_c_cll = c_cll.T.sample(n=5000)
r_c_cll = c_cll.T.loc[chosen_cells]
r_c_cll.columns.name = "Chromosome positions"

g = sns.clustermap(
    r_c_cll,
    col_cluster=False,
    xticklabels=False, yticklabels=False,
    cmap="RdBu_r", center=0, vmin=-2, vmax=2, robust=False,
    row_colors=get_colors(cnv.obs.loc[r_c_cll.index], factors=["patient_id", "timepoint"]),
    figsize=(6, 3),
    metric="correlation",
    cbar_kws={"label": "log difference from normal"},
    rasterized=True)
clustermap_rasterize_heatmap(g)
clustermap_fix_label_orientation(g)
clustermap_rasterize_dendrogram(g)
savefig(g, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.whole_genome.single_cells.2500_per_patient.CLL.svg"))


r_c_cll_o = r_c_cll.loc[cnv.obs.sort_values(['patient_id', 'timepoint']).index, :].dropna()

g = sns.clustermap(
    r_c_cll_o ** 3 * 3,  # add some smoothing
    row_cluster=False, col_cluster=False,
    xticklabels=False, yticklabels=False,
    cmap="RdBu_r", center=0, vmin=-2, vmax=2, robust=False,
    row_colors=get_colors(cnv.obs.loc[r_c_cll_o.index], factors=["patient_id", "timepoint"]),
    figsize=(6, 3),
    metric="correlation",
    cbar_kws={"label": "log difference from normal"},
    rasterized=True)
clustermap_rasterize_heatmap(g)
clustermap_fix_label_orientation(g)
clustermap_rasterize_dendrogram(g)
savefig(g, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.whole_genome.single_cells.2500_per_patient.CLL.sorted.svg"))


# Plot same heatmap but with chromosomes separately
import pybedtools
chrom_sizes = dict(pybedtools.chromsizes('hg38'))
total_size = sum([x[1] for x in chrom_sizes.values()])
for chrom in gene_order['chr'].unique():
    print(chrom)
    chrom_size = chrom_sizes[chrom][1] / total_size
    p = r_c_cll_o.reindex(gene_order.loc[gene_order['chr'] == chrom, "gene"], axis=1).T.dropna().T

    if p.empty:
        continue

    g = sns.clustermap(
        p ** 3 * 3,  # add some smoothing
        row_cluster=False, col_cluster=False,
        xticklabels=False, yticklabels=False,
        cmap="RdBu_r", center=0, vmin=-2, vmax=2, robust=False,
        # row_colors=get_colors(cnv.obs.loc[p.index], factors=["patient_id", "timepoint"]),
        figsize=(8 * chrom_size, 3),
        # metric="correlation",
        # cbar_kws={"label": "log difference from normal"},
        rasterized=True, cbar=False)
    g.ax_heatmap.set_title(chrom)
    # remove colormap and any other marks
    g.cax.set_visible(False)
    g.ax_heatmap.set_xlabel("", visible=False)
    col = g.ax_col_dendrogram.get_position()
    g.ax_col_dendrogram.set_position([col.x0, col.y0, col.width, col.height*0.01])
    # add black square around
    g.ax_heatmap.get_children()[-1].set_linewidth(0.5)
    clustermap_rasterize_heatmap(g)
    clustermap_fix_label_orientation(g)
    clustermap_rasterize_dendrogram(g)
    savefig(g, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.whole_genome.single_cells.2500_per_patient.CLL.sorted.{}.svg".format(chrom)))
    plt.close()


cm = c.drop(['start', 'end'], axis=1).groupby(['chr']).mean().T
cm = cm.reindex(natsorted(cm.columns), axis=1)

# g = sns.clustermap(
#     cm,  # - cm.mean().mean(),
#     col_cluster=False,
#     xticklabels=False, yticklabels=True,
#     cmap="RdBu_r", center=0,  # robust=True,
#     figsize=(4, 8),
#     metric="correlation",
#     cbar_kws={"label": "log difference from normal"},
#     rasterized=True)
# clustermap_rasterize_heatmap(g)
# clustermap_fix_label_orientation(g)
# clustermap_rasterize_dendrogram(g)
# savefig(g, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.chrs.single_cells.svg"))

cnv.obs = cnv.obs.join(cm)

fig, axis = plt.subplots(3, 9, figsize=(9 * 3, 3 * 3))
axis = axis.flatten()
for i, attr in enumerate(['cell_type', 'patient_id', 'timepoint']):
    sc.pl.pca(cnv, color=attr, ax=axis[i], show=False)
for i, attr in enumerate(list(range(1, 23)) + ['X']):
    # v = cnv.obs[f"chr{attr}"].abs().max()
    # sc.pl.pca(cnv, color=f"chr{attr}", cmap="RdBu_r", vmin=-v, vmax=v, ax=axis[i + 3], show=False)
    sc.pl.pca(cnv, color=f"chr{attr}", cmap="RdBu_r", vmin=-2, vmax=2, ax=axis[i + 3], show=False)
savefig(fig, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.pca.single_cells_colored_by_chr.fixed_scale.svg"))

attrs = ["cell_type", "patient_id", "timepoint"]
fig, axis = plt.subplots(9, 3, figsize=(3 * 3, 9 * 3))
axis = axis.flatten()
for i, attr in enumerate(['cell_type', 'patient_id', 'timepoint']):
    sc.pl.umap(cnv, color=attr, ax=axis[i], show=False, palette=colors[attr])
for i, attr in enumerate(list(range(1, 23)) + ['X']):
    # v = cnv.obs[f"chr{attr}"].abs().max()
    # sc.pl.umap(cnv, color=attr, cmap="RdBu_r", vmin=-v, vmax=v, ax=axis[i + 3], show=False)
    sc.pl.umap(cnv, color=f"chr{attr}", cmap="RdBu_r", vmin=-2, vmax=2, ax=axis[i + 3], show=False)
savefig(fig, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.umap.single_cells_colored_by_chr.fixed_scale.svg"))

# Correlate CNVs with response

# # Intra-timepoint variability in CNVs between single-cells (clonal variability)
cll_cnv = cnv[cnv.obs['cell_type'] == "CLL", :]
cll_cnv.obs.loc[:, "sample"] = cll_cnv.obs.loc[:, "patient_id"].astype(str) + " - " + cll_cnv.obs.loc[:, "timepoint"].astype(str)

dists = dict()
for sample in cll_cnv.obs['sample'].unique():
    print(sample)
    cells = cll_cnv.obs["sample"] == sample
    intra_sample_corr = pd.DataFrame(cll_cnv.X[cells, :], index=cells[cells].index).T.corr()
    # take half of the values (without diagonal)
    np.fill_diagonal(intra_sample_corr.values, np.nan)
    ss = np.tri(*intra_sample_corr.shape)
    ss[ss == 0] = np.nan
    intra_sample_corr *= ss
    iv = intra_sample_corr.values.flatten()
    iv = iv[~pd.isnull(iv)]
    dists[sample] = iv

fig, axis = plt.subplots(4, 4, figsize=(2 * 4, 2 * 4))
axis = axis.flatten()
stats = dict()
for i, sample in enumerate(cnv.obs['sample'].unique()):
    print(sample)
    q = -dists[sample]
    sns.distplot(q, ax=axis[i], hist=False)
    axis[i].set_title(sample + "(n = {}".format(len(q)))
    stats[sample] = q.mean(), q.std()
savefig(fig, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.grouped.CLL.correlation_change.distplot.svg"))
stats = pd.DataFrame(stats).T
stats.columns = ['mean', 'std']
stats.loc[:, "qv"] = stats['std'] / stats['mean']
stats.loc[:, "patient_id"] = stats.index.str.split(" - ").map(lambda x: x[0])
stats.loc[:, "timepoint"] = stats.index.str.split(" - ").map(lambda x: x[1])

change = adata.obs[['patient_id', 'response_at_120']].drop_duplicates().set_index("patient_id")
change.loc[:, 'clonal_diversity_at_day0(mean)'] = stats.groupby("patient_id").apply(
    lambda x: x[x['timepoint'] == "000d"]["mean"].values.squeeze())
change.loc[:, 'clonal_diversity_at_day0(qv)'] = stats.groupby("patient_id").apply(
    lambda x: x[x['timepoint'] == "000d"]["qv"].values.squeeze())
change.loc[:, 'change_clonal_diversity(mean)'] = stats.groupby("patient_id").apply(
    lambda x:
        (x[x['timepoint'] == "120d"]["mean"].values -
         x[x['timepoint'] == "000d"]["mean"].values).squeeze())
change.loc[:, 'change_clonal_diversity(qv)'] = stats.groupby("patient_id").apply(
    lambda x:
        (x[x['timepoint'] == "120d"]["qv"].values -
         x[x['timepoint'] == "000d"]["qv"].values).squeeze())
change.loc[:, 'min_clonal_diversity'] = stats.groupby("patient_id")['qv'].min()
change.loc[:, 'mean_clonal_diversity'] = stats.groupby("patient_id")['qv'].mean()
change.loc[:, 'max_clonal_diversity'] = stats.groupby("patient_id")['qv'].max()

g = sns.pairplot(change.astype(float), height=2, aspect=1)
savefig(g, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.grouped.CLL.correlation_change.scatter.svg"))


# # Get signature score for each patient, timepoint, cell type
common = pd.read_csv(os.path.join("results", "offtarget_signature.csv"), index_col=0, header=None, squeeze=True)
e_z = adata.to_df().loc[:, common.index]
t_sigs = e_z.loc[:, common[common > 0].index].mean(axis=1) - e_z.loc[:, common[common < 0].index].mean(axis=1)
t_sigs.name = "ibrutinib_score"
adata.obs = adata.obs.join(t_sigs)

# # # add also to cnvs
cnv.obs = cnv.obs.join(adata.obs['ibrutinib_score'])
cnv.obs.loc[:, "sample"] = cnv.obs.loc[:, "patient_id"].astype(str) + " - " + cnv.obs.loc[:, "timepoint"].astype(str)

# # plot a few examples of clonality and their response to ibrutinib
cll_cnv = cnv[cnv.obs['cell_type'] == "CLL", :]

fig, axis = plt.subplots(1, 4, figsize=(4 * 3, 1 * 2.8))
pp = cll_cnv.obs.loc[cll_cnv.obs['sample'] == "CLL1 - 000d"]
axis[0].scatter(pp['chr11'], pp['ibrutinib_score'], s=3, alpha=0.3, rasterized=True)
axis[0].set_title("Day 0")
axis[0].set_xlabel("chr11")
axis[0].set_ylabel("Ibrutinib response score")
pp = cll_cnv.obs.loc[cll_cnv.obs['sample'] == "CLL5 - 030d"]
axis[1].scatter(pp['chr12'], pp['ibrutinib_score'], s=3, alpha=0.3, rasterized=True)
axis[1].set_title("Day 30")
axis[1].set_xlabel("chr12")
axis[1].set_ylabel("Ibrutinib response score")
pp = cll_cnv.obs.loc[cll_cnv.obs['sample'] == "CLL6 - 030d"]
axis[2].scatter(pp['chr12'], pp['ibrutinib_score'], s=3, alpha=0.3, rasterized=True)
axis[2].set_title("Day 30")
axis[2].set_xlabel("chr12")
axis[2].set_ylabel("Ibrutinib response score")
pp = cll_cnv.obs.loc[cll_cnv.obs['sample'] == "CLL8 - 120d"]
axis[3].scatter(pp['chr17'], pp['ibrutinib_score'], s=3, alpha=0.3, rasterized=True)
axis[3].set_title("Day 120")
axis[3].set_xlabel("chr17")
axis[3].set_ylabel("Ibrutinib response score")
savefig(fig, os.path.join(
    "results", prefix +
    ".dca_denoised-zinb.CNVs.single-cells.clones.examples_vs_ibrutinib_response.scatter.svg"))


# umap without nan cells
attrs = ["cell_type", "patient_id", "timepoint", "batch", "log_counts", "ibrutinib_score"]

fig, axis = plt.subplots(1, 6, figsize=(2 * 6.5, 2 * 1))
for i, attr in enumerate(attrs):
    if i < 4:
        k = {"palette": colors[attr]}
    elif i == 5:
        k = {"cmap": "RdBu_r", "vmin": -2, "vmax": 2}
    else:
        k = {}
    sc.pl.pca(
        adata,
        color=attr,
        alpha=0.75, s=1, show=False, ax=axis[i], **k)
savefig(fig, os.path.join("results", prefix + ".ibrutinib_score.pca.svg"))
fig, axis = plt.subplots(1, 6, figsize=(2 * 6.5, 2 * 1))
for i, attr in enumerate(attrs):
    if i < 4:
        k = {"palette": colors[attr]}
    elif i == 5:
        k = {"cmap": "RdBu_r", "vmin": -2, "vmax": 2}
    else:
        k = {}
    sc.pl.umap(
        adata,
        color=attr,
        alpha=0.75, s=1, show=False, ax=axis[i], **k)
savefig(fig, os.path.join("results", prefix + ".ibrutinib_score.umap.svg"))


adata.obs.loc[:, 'time'] = adata.obs.loc[:, 'timepoint'].str.replace("d", "").astype(int)
fig, axis = plt.subplots(1, 1, figsize=(3, 3))
sns.violinplot(data=adata.obs, x='time', y='ibrutinib_score', hue="patient_id", ax=axis)
savefig(fig, os.path.join("results", prefix + f".dca_denoised-zinb.ibrutinib_score_per_cell.violinplot.svg"))

fig, axis = plt.subplots(1, 1, figsize=(3, 3))
sns.factorplot(data=adata.obs, x='time', y='ibrutinib_score', hue="patient_id", ax=axis)
savefig(fig, os.path.join("results", prefix + f".dca_denoised-zinb.ibrutinib_score_per_cell.factorplot.svg"))

fig, axis = plt.subplots(1, 1, figsize=(3, 3))
sns.swarmplot(data=adata.obs, x='time', y='ibrutinib_score', hue="patient_id", ax=axis, rasterized=True)
savefig(fig, os.path.join("results", prefix + f".dca_denoised-zinb.ibrutinib_score_per_cell.swarmplot.svg"))

fig, axis = plt.subplots(1, 1, figsize=(3, 3))
sns.scatterplot(data=adata.obs, x='time', y='ibrutinib_score', hue="patient_id", ax=axis, rasterized=True)
savefig(fig, os.path.join("results", prefix + f".dca_denoised-zinb.ibrutinib_score_per_cell.scatterplot.svg"))

fig, axis = plt.subplots(1, 1, figsize=(3, 3))
axis.scatter(adata.obs['log_counts'], adata.obs['ibrutinib_score'], rasterized=True, s=0.5, alpha=0.1)
savefig(fig, os.path.join("results", prefix + f".dca_denoised-zinb.ibrutinib_score_per_cell_vs_counts.scatter.svg"))


# # # Illustrate for each patient
if "ibrutinib_score" not in cll_cnv.obs.columns:
    cll_cnv.obs = cll_cnv.obs.join(adata.obs['ibrutinib_score'])

# add expression of CD5 and TCL1 to see if these are CLL
# cll_cnv.obs = pd.DataFrame(cll_cnv.obs.values, index=cll_cnv.obs.index.astype(str).values, columns=cll_cnv.obs.columns.astype(str).values)
s = pd.Series(adata[:, "CD5"][cll_cnv.obs_names, :].X.toarray(), index=cll_cnv.obs_names, name="CD5")
cll_cnv.obs.loc[:, "CD5"] = np.nan
cll_cnv.obs.loc[s.index, "CD5"] = s.values
s = pd.Series(adata[:, "TCL1A"][cll_cnv.obs_names, :].X.toarray(), index=cll_cnv.obs_names, name="TCL1A")
cll_cnv.obs.loc[:, "TCL1A"] = np.nan
cll_cnv.obs.loc[s.index, "TCL1A"] = s.values
for patient in ["CLL1", "CLL5", "CLL6", "CLL8"]:
    print(patient)
    cells = (cll_cnv.obs["patient_id"] == patient)
    p = cll_cnv[cells, :]

    datas = list()

    cells = (p.obs["timepoint"] == "000d")
    if cells.sum() > 0:
        p_d0 = p[cells, :]
        sc.pp.pca(p_d0)
        sc.pp.neighbors(p_d0)
        sc.tl.umap(p_d0)
        datas.append((0, p_d0))

    cells = (p.obs["timepoint"] == "030d")
    if cells.sum() > 0:
        p_d30 = p[cells, :]
        sc.pp.pca(p_d30)
        sc.pp.neighbors(p_d30)
        sc.tl.umap(p_d30)
        datas.append((30, p_d30))

    cells = (p.obs["timepoint"] == "120d")
    if cells.sum() > 0:
        p_d120 = p[cells, :]
        sc.pp.pca(p_d120)
        sc.pp.neighbors(p_d120)
        sc.tl.umap(p_d120)
        datas.append((120, p_d120))

    if len(datas) == 0:
        continue

    x = len(datas)
    fig, axis = plt.subplots(x, 9, figsize=(9 * 3, x * 3))
    if len(axis.shape) == 1:
        axis = axis.reshape((1, axis.shape[0]))
    for i, (t, data) in enumerate(datas):
        for j, attr in enumerate(["log_counts", "chr1", "chr11", "chr12", "chr13", "chr17", "CD5", "TCL1A", "ibrutinib_score"]):
            if j in [0, 6, 7]:
                k = {}
            elif j == 8:
                k = {"cmap": "RdBu_r", "vmin": -2, "vmax": 2}
            else:
                k = {"cmap": "RdBu_r", "vmin": -1, "vmax": 1}
            sc.pl.umap(data, color=attr, ax=axis[i, j], show=False, **k)
            axis[i, j].set_title(f"{attr} - {t}")
    savefig(fig, os.path.join("results", prefix + f".dca_denoised-zinb.CNVs.{patient}_per_timepoint.umap.scatter.svg"))


# # Between time change in CNV profile (clonal change)
p2_cor = p2.T.corr()

g = sns.clustermap(
    p2_cor, xticklabels=False, yticklabels=True, cmap="RdBu_r", center=0)
clustermap_rasterize_heatmap(g)
clustermap_fix_label_orientation(g)
clustermap_rasterize_dendrogram(g)
savefig(g, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.grouped.correlation.svg"))

p2_cor = p2.loc[p2.index.get_level_values("cell_type") == "CLL"].T.corr()
g = sns.clustermap(
    p2_cor, xticklabels=False, yticklabels=True, cmap="RdBu_r", center=0, figsize=(3, 3))
clustermap_rasterize_heatmap(g)
clustermap_fix_label_orientation(g)
clustermap_rasterize_dendrogram(g)
savefig(g, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.grouped.CLL.correlation.svg"))

p2_cor = p2.loc[
    (p2.index.get_level_values("cell_type") == "CLL") &
    (p2.index.get_level_values("timepoint").isin(['120d', '000d']))
    ].T.corr(method=scipy.spatial.distance.euclidean)

change = adata.obs[['patient_id', 'response_at_120']].drop_duplicates().set_index("patient_id")
change.loc[:, 'clonal_change'] = pd.Series(np.diag(p2_cor, 1)[::2], p2_cor.index.levels[0])

fig, axis = plt.subplots(1, figsize=(3, 3))
axis.scatter(change['clonal_change'], change['response_at_120'])
for i in change.index:
    axis.text(change.loc[i, 'clonal_change'], change.loc[i, 'response_at_120'], s=i)
axis.set_xlabel("Relative clonal change (day 120 vs 0)")
axis.set_ylabel("Response to ibrutinib (day 120 vs 0)")
savefig(fig, os.path.join("results", prefix + ".dca_denoised-zinb.CNVs.grouped.CLL.euclidean_change.scatter.svg"))
