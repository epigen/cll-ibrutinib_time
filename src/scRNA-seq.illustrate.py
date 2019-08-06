#!/usr/bin/env python


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
# from sklearn.metrics import silhouette_score

from peppy import Project
# from ngs_toolkit.analysis import Analysis

import scanpy as sc
# from dca.utils import plot_mean_dropout

from natsort import natsorted as sorted


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


def violinplot_rasterize_jitter(ax):
    import matplotlib
    q = filter(
            lambda x: isinstance(x, matplotlib.collections.PathCollection),
            ax.get_children())
    for i in q:
        i.set_rasterized(True)


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


mark1 = ['CD3D', 'CD3G', 'CD4', 'CD8A', 'NKG7', 'CD14', "CST3", "CD79A", "TCL1A", "CD19", 'GZMB', 'CD68', 'CD27', 'CD20', 'CD24', "NCR1"]
mark2 = ["IL32", "IFNG", "IFNGR1", "IFNGR1", "IL4R", "IL4", "JUN", "JUNB", "JUND", "JAK1", "JAK2", "GATA1", "JARID2", "KRAS", "MYC"]
mark3 = ["BTK", "LCK", "E2F4", "CXCR4", "ITGA4", "HBA1", "PTPRC"]
red_mark = ["CST3", "TCL1A", "GZMB", "NKG7", "CD3D"]
pbmc_markers = mark1 + mark2 + mark3


prj = Project(os.path.join("metadata", "project_config.scRNA-seq.yaml"))
data_type = "scRNA-seq"


# Load data
prefix = "cll-time_course-scRNA-seq.all_samples.250-50_filter"
adata = sc.read(prefix + ".dca_denoised-zinb.processed.h5ad")
adata.obs = pd.read_csv(prefix + ".dca_denoised-zinb.processed.obs.csv", index_col=0)


# Fig1
# # Plotting markers

# umap without nan cells
attrs = ["cell_type", "patient_id", "timepoint", "batch", "log_counts"]
colors = {
    "cell_type": ["#b7cee3", "#c4e6a7", "#f6c4c2", "#f7cc82", "#fdfe97", "#cdbbdc", "#cdbbdc"],
    "patient_id": ["#84190d", "#8b8c00", "#318b00", "#5d5d5d"],
    "timepoint": ["#101013", "#67257d", "#ef9d00", "#fafca5"],
    "batch": sns.color_palette()[:4]}

fig, axis = plt.subplots(1, 5, figsize=(2 * 5, 2 * 1))
for i, attr in enumerate(attrs):
    sc.pl.pca(
        adata,
        color=attr,
        palette=colors[attr] if i != 4 else None,
        alpha=0.75, s=1, show=False, ax=axis[i])
savefig(fig, os.path.join("results", prefix + ".pca.paper_colors.svg"))
fig, axis = plt.subplots(1, 5, figsize=(2 * 5, 2 * 1))
for i, attr in enumerate(attrs):
    sc.pl.umap(
        adata,
        color=attr,
        palette=colors[attr] if i != 4 else None,
        alpha=0.75, s=1, show=False, ax=axis[i])
savefig(fig, os.path.join("results", prefix + ".umap.paper_colors.svg"))

markers = [x for x in pbmc_markers if x in adata.var_names]
n = int(np.ceil(np.sqrt(len(markers))))
fig, axis = plt.subplots(n, n, figsize=(2 * n, 2 * n))
axis = axis.flatten()
for i, mark in enumerate(markers):
    print(i)
    sc.pl.pca(
        adata,
        color=mark,
        alpha=0.75, s=1, show=False, ax=axis[i])
savefig(fig, os.path.join("results", prefix + ".pca.pbmc_markers.svg"))
fig, axis = plt.subplots(n, n, figsize=(2 * n, 2 * n))
axis = axis.flatten()
for i, mark in enumerate(markers):
    print(i)
    sc.pl.umap(
        adata,
        use_raw=False, color=mark,
        alpha=0.75, s=1, show=False, ax=axis[i])
savefig(fig, os.path.join("results", prefix + ".umap.pbmc_markers.svg"))


# # Differential analysis
# # # from Nikolaus' results
expr_diff = pd.read_csv(os.path.join(
    "results", "single_cell_RNA", "13_4_Overtime_nUMI_Cutoff", "SigGenes_overTime.tsv"), sep="\t").rename(
        columns={"cellType": "cell_type", "logFC": "log_fold_change"}).set_index("gene")
expr_diff.loc[:, 'intercept'] = 1
expr_diff.loc[:, 'patient_id'] = expr_diff['patient'].str.split("_").apply(lambda x: x[0])
expr_diff.loc[:, 'timepoint'] = expr_diff['patient'].str.split("_").apply(lambda x: x[1]).str.replace("d", "").astype(int)
expr_diff.loc[:, 'direction'] = (expr_diff['log_fold_change'] > 0).astype(int).replace(0, -1)
expr_diff.loc[:, '-log10_pvalue'] = (-np.log10(expr_diff['pvalue'])).replace(np.nan, 300)
# add mean expression
expr_diff = expr_diff.join(pd.Series(np.log10(1 + np.asarray(adata.raw.X.mean(0)).squeeze()), index=adata.var_names, name="raw_mean"))
expr_diff = expr_diff.join(pd.Series(adata.X.mean(0), index=adata.var_names, name="mean"))

# tracksplot of top diff genes
cll_markers = [
    "HMGN2", "CXCR4", "RPL22", "RP11−160E2.6", "H3F3A", "ZFP36L1", "IER2", "NFKBIA",
    "FOS", "LGALS1", "S100A4", "CRIP1", "RPL36A", "RPS28", "RPL17", "ZFP36",
    "DUSP1", "ZFP36L2", "CD27", "CCDC88A", "RPL41", "RPL39", "HIF1A", "CXCL8",
    "CCL3L3", "EEF1A1", "NBEAL1", "CCL3", "AC090498.1", "LAPTM5", "DDX5", "JUN"]

ordered_markers = [
    "FOS", "JUN", "CD27", "DDX5", "LAPTM5", "CD27", "ZFP36", "DUSP1",
    "CXCR4",
    "RPL17", "RPL36A"]

adata.obs['cell_type'] = adata.obs['cell_type'].replace("nan", pd.np.nan)
sc.pl.dotplot(
    adata, list(set([x for x in pbmc_markers if x in adata.var_names])),
    groupby='cell_type', use_raw=False, color_map="RdBu_r", vmin=-3, vmax=3, show=False)
savefig(plt.gca().figure, os.path.join("results", "scRNA-seq.pbmc_markers.dotplot.svg"))
sc.pl.dotplot(
    adata, list(set([x for x in cll_markers if x in adata.var_names])),
    groupby='cell_type', use_raw=False, color_map="RdBu_r", vmin=-3, vmax=3, show=False)
savefig(plt.gca().figure, os.path.join("results", "scRNA-seq.cll_changing_genes.dotplot.svg"))

sc.pl.dotplot(
    adata[adata.obs['cell_type'] == "CD19", :], list(set([x for x in cll_markers if x in adata.var_names])),
    groupby='timepoint', use_raw=False, color_map="RdBu_r", vmin=-3, vmax=3, show=False)
savefig(plt.gca().figure, os.path.join("results", "scRNA-seq.cll_changing_genes.cll_only.dotplot.svg"))

sc.pl.dotplot(
    adata[adata.obs['cell_type'] == "CD19", :], list(set([x for x in cll_markers if x in adata.var_names])),
    groupby='timepoint', show=False)
savefig(plt.gca().figure, os.path.join("results", "scRNA-seq.cll_changing_genes.cll_only.raw.dotplot.svg"))

sc.pl.stacked_violin(
    adata[adata.obs['cell_type'] == "CD19", :], list(set([x for x in cll_markers if x in adata.var_names])),
    groupby='timepoint', use_raw=False, color_map="RdBu_r", vmin=-3, vmax=3, show=False)
savefig(plt.gca().figure, os.path.join("results", "scRNA-seq.cll_changing_genes.cll_only.violin.svg"))

sc.pl.stacked_violin(
    adata[adata.obs['cell_type'] == "CD19", :], list(set([x for x in cll_markers if x in adata.var_names])),
    groupby='timepoint', show=False)
savefig(plt.gca().figure, os.path.join("results", "scRNA-seq.cll_changing_genes.cll_only.raw.violin.svg"))


# Plot a few usual suspects
cll_genes = [
    "CD5", "CD27", "CD38", "TCL1A", "CD24", "CD52", "CD19", "CD20",
    "JUN", "BTK", "PLCG2",
    "RELA", "NFKB1", "NFKB2", "NFKBIA", "ZAP70",
    "FAS", "FASLG", "CD28"]
for cell_type in ["CLL", "CD8", "CD4", "NK", "Mono"]:
    print(cell_type)
    data = adata[adata.obs['cell_type'] == cell_type, :]
    data = data[data.obs['timepoint'] != "280d", :]
    data.strings_to_categoricals()

    ax = sc.pl.violin(
        data,
        [x for x in cll_genes
            if x in data.var_names],
        "timepoint", use_raw=False, rasterize=False, show=False)
    [violinplot_rasterize_jitter(c) for c in ax]
    savefig(ax[0].figure, os.path.join("results", "scRNA-seq.differential.{}.cll_genes.norm_data.violinplot.svg"
                                       .format(cell_type)))


# For top genes per cell type plot violin plots and trackplots
n_top = 20
n_cells = 2000
for cell_type in ["CLL", "CD8", "CD4", "NK", "Mono"]:
    for direction, f in [
            ("down", expr_diff['log_fold_change'] < 0),
            ("up", expr_diff['log_fold_change'] > 0)]:
        plt.close("all")

        # get down genes
        sig_expr = expr_diff.loc[
            (expr_diff['qvalue'] < 0.05) &
            (expr_diff['cell_type'] == cell_type) &
            f]
        print(cell_type, direction, sig_expr.shape)
        # add ranks
        sig_expr.loc[:, "raw_mean_rank"] = sig_expr.loc[:, "raw_mean"].rank(ascending=False)
        sig_expr.loc[:, "mean_rank"] = sig_expr.loc[:, "mean"].rank(ascending=False)
        sig_expr.loc[:, "log_fold_change_rank"] = sig_expr.loc[:, "log_fold_change"].abs().rank(ascending=False)
        sig_expr.loc[:, "-log10_pvalue_rank"] = sig_expr.loc[:, "-log10_pvalue"].rank(ascending=False)
        # add max rank
        sig_expr.loc[:, "max_rank"] = sig_expr.loc[:, sig_expr.columns.str.endswith("_rank")].max(1)
        genes = list(set(sig_expr[["max_rank"]].drop_duplicates().sort_values('max_rank').head(n_top).index))  # sig.head(n_top).index.tolist() +
        data = adata[adata.obs['cell_type'] == cell_type, :]
        data = data[data.obs['timepoint'] != "280d", :]
        data.strings_to_categoricals()

        ax = sc.pl.violin(data, genes, "timepoint", use_raw=False, rasterize=False, show=False)
        [violinplot_rasterize_jitter(c) for c in ax]
        savefig(ax[0].figure, os.path.join("results", "scRNA-seq.differential.{}.top_{}_{}_genes.norm_data.violinplot.svg"
                                           .format(cell_type, n_top, direction)))
        # sample N random cells
        data = data[np.random.choice(data.obs_names, n_cells), :]
        data.strings_to_categoricals()

        ax = sc.pl.tracksplot(
            data,
            var_names=[x for x in genes if x in data.var_names],
            groupby="timepoint",
            # dendrogram=True,
            use_raw=True, log=False, show=False)
        ax[0].set_title(cell_type)
        savefig(ax[0].figure, os.path.join("results", "scRNA-seq.differential.{}.top_{}_{}_genes.zscore.only_{}_cells.tracksplot.svg"
                                           .format(cell_type, n_top, direction, n_cells)))
        ax = sc.pl.tracksplot(
            data,
            var_names=[x for x in genes if x in data.var_names],
            groupby="timepoint",
            # dendrogram=True,
            use_raw=False, log=False, show=False)
        ax[0].set_title(cell_type)
        savefig(ax[0].figure, os.path.join("results", "scRNA-seq.differential.{}.top_{}_{}_genes.zscore.only_{}_cells.tracksplot.norm.svg"
                                           .format(cell_type, n_top, direction, n_cells)))

# # plot extreme genes of each cell type

for n_cells in [2000, 500000]:
    for n_top in [25, 500, 100000]:
        for cell_type in ["CLL", "CD4", "CD8", "Mono", "NK"]:
            print(cell_type, n_top)
            sig_expr = expr_diff.loc[
                (expr_diff['qvalue'] < 0.05) &
                (expr_diff['cell_type'] == cell_type)]

            sig = sig_expr.groupby(level=0)['log_fold_change'].mean().sort_values()
            genes = list(set(sig.head(n_top).index.tolist() + sig.tail(n_top).index.tolist()))

            data = adata[adata.obs['cell_type'] == cell_type, :]

            # z_score
            ds = pd.DataFrame(
                scipy.stats.zscore(data.X, axis=0),
                index=data.obs.index, columns=data.var.index)
            # subsample
            ds = ds.sample(
                n=n_cells if n_cells < ds.shape[0] else ds.shape[0],
                replace=False).loc[:, genes]

            # sort by timepoint, patient
            d = data.obs.loc[ds.index].sort_values(['patient_id', 'timepoint'])
            ds = ds.loc[d.index, :]

            colors = list()
            for f, cmap in [
                    ('patient_id', 'tab20'),
                    ('timepoint', 'tab20'),
                    # ('cell_type', 'tab20')
            ]:
                q = data.obs.loc[ds.index, f]
                qq = sorted(q.unique())
                m = dict(zip(qq, range(len(qq))))
                cmap = plt.get_cmap(cmap)
                colors.append(cmap([m[x] for x in q]))

            g = sns.clustermap(
                ds.T.dropna().T,
                row_cluster=False,
                row_colors=colors,
                metric="euclidean",
                xticklabels=False, yticklabels=False, cbar_kws={"label": "Expression\n(Z-score)"},
                cmap="RdBu_r", center=0, robust=True)
            g.ax_heatmap.set_xlabel("Genes (n = {})".format(ds.shape[1]))
            g.ax_heatmap.set_ylabel("Cells (n = {})".format(ds.shape[0]))
            clustermap_rasterize_heatmap(g)
            clustermap_rasterize_dendrogram(g)
            clustermap_fix_label_orientation(g)
            savefig(g, os.path.join("results", "scRNA-seq.differential.{}.top_{}.zscore.only_{}_cells.clustermap.svg"
                                    .format(cell_type, n_top, ds.shape[0])))
        plt.close('all')


# Plot heatmaps
a = adata.to_df()

diff = os.path.join("results", "single_cell_RNA", "13_4_Overtime_nUMI_Cutoff", "SigGenes_overTime.tsv")
diff = pd.read_csv(diff, sep="\t")

d = diff.loc[
    # (diff['cellType'] == "CLL") &
    (diff['qvalue'] < 0.01) &
    (diff['logFC'].abs() > 0.5),
    "gene"].unique()

# subsample
pa = a.loc[np.random.choice(a.index, 20000), d]
# order
adata.obs.loc[:, 'cell_type'] = adata.obs['cell_type'].replace("CLL", "CD19")
pa = pa.loc[adata.obs.loc[pa.index].sort_values(['cell_type', 'patient_id', 'timepoint']).index]

colors = list()
for f, cmap in [
        ('patient_id', 'tab20'),
        ('timepoint', 'Set1'),
        ('cell_type', 'tab20')]:
    q = adata.obs.loc[pa.index, f]
    qq = sorted(q.unique())
    m = dict(zip(qq, range(len(qq))))
    cmap = plt.get_cmap(cmap)
    colors.append(cmap([m[x] for x in q]))

g = sns.clustermap(
    pa,
    row_colors=colors, rasterized=True, cbar_kws={"label": "Expression (Z-score)"},
    metric="correlation",
    xticklabels=False, yticklabels=False,
    robust=True, z_score=1, cmap="RdBu_r", center=0)
g.ax_heatmap.set_xlabel(f"Genes (n = {len(d)})")
g.ax_heatmap.set_ylabel("Cells")
g2 = sns.clustermap(
    pa,
    row_cluster=False,
    row_colors=colors, rasterized=True, cbar_kws={"label": "Expression (Z-score)"},
    metric="correlation",
    xticklabels=False, yticklabels=False,
    robust=True, z_score=1, cmap="RdBu_r", center=0)
g2.ax_heatmap.set_xlabel(f"Genes (n = {len(d)})")
g2.ax_heatmap.set_ylabel("Cells")

for q, label in [(g, ".z_score"), (g2, ".z_score.ordered")]:
    print(label)
    clustermap_rasterize_heatmap(q)
    clustermap_fix_label_orientation(q)
    clustermap_rasterize_dendrogram(q)
    savefig(q, os.path.join("results", prefix + ".dca_denoised-zinb.diff_genes_across_cell_types{}.svg".format(label)))


# heatmaps of diff per cell type
for cell_type in ["CLL", "CD8", "CD4", "NK", "Mono"]:
    if cell_type == "nan":
        continue
    print(cell_type)

    d = diff.loc[
        (diff['cellType'] == cell_type) &
        (diff['qvalue'] < 0.05) &
        (diff['logFC'].abs() > 0.5),
        "gene"].unique()

    c = adata.obs.loc[adata.obs['cell_type'] == cell_type, :].index
    d = d[pd.Series(d).isin(a.columns.tolist())]
    pa = a.loc[c, d].fillna(0)

    colors = list()
    for f, cmap in [
            ('patient_id', 'tab20'),
            ('timepoint', 'Set1')]:
        q = adata.obs.loc[pa.index, f]
        m = dict(zip(q.unique(), range(len(q.unique()))))
        cmap = plt.get_cmap(cmap)
        colors.append(cmap([m[x] for x in q]))

    g = sns.clustermap(
        pa,
        row_colors=colors, rasterized=True, cbar_kws={"label": "Expression (Z-score)"},
        metric="correlation",
        xticklabels=False, yticklabels=False,
        robust=True, z_score=1, cmap="RdBu_r", center=0)
    g.ax_heatmap.set_xlabel(f"Genes (n = {len(d)})")
    g.ax_heatmap.set_ylabel("Cells")
    clustermap_rasterize_heatmap(g)
    clustermap_fix_label_orientation(g)
    clustermap_rasterize_dendrogram(g)
    savefig(g, os.path.join("results", "scRNA-seq.{}.{}.specific_genes.zscore.svg".format("denoised", cell_type)))


# Cell cycle
url = "https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt"
cell_cycle_genes = pd.read_csv(url, header=None, squeeze=True)
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# cll = adata.obs.loc[adata.obs['cell_type'] == "CLL"]
sns.catplot(
    data=adata.obs,
    x="timepoint", y="S_score",
    hue="patient_id",
    kind="boxen")

s = adata.obs.groupby(['cell_type', 'timepoint'])['phase'].value_counts(normalize=True)
s.name = "cell_count"
sns.catplot(
    data=s.reset_index(),
    x="timepoint", y="cell_count",
    col="cell_type",  hue="phase",
    kind="boxen")

# # Try dpt ordering of CLL cells
# cll = adata[adata.obs['cell_type'] == "CLL", :]

# # # First regress out batch and patient
# cll_b = cll.copy()
# sc.pp.combat(cll_b, 'batch', inplace=True)
# cll_bp = cll_b.copy()
# sc.pp.combat(cll_bp, 'patient_id', inplace=True)


# # # regress out
# cll_reg = sc.pp.regress_out(cll[:, common.index], keys='batch', copy=True)
# sc.pp.neighbors(cll_reg, use_rep="X")
# sc.tl.diffmap(cll_reg)
# cll_reg2 = sc.pp.regress_out(cll_reg, keys='patient_id', copy=True)
# sc.pp.neighbors(cll_reg2, use_rep="X")
# sc.tl.diffmap(cll_reg2)

# sc.pl.diffmap(cll, color=['patient_id', 'timepoint', 'log_counts'], show=False)
# sc.pl.diffmap(cll_reg, color=['patient_id', 'timepoint', 'log_counts'], show=False)
# sc.pl.diffmap(cll_reg2, color=['patient_id', 'timepoint', 'log_counts'], show=False)

# # # On all genes
# for data in tqdm([cll, cll_b, cll_bp]):
#     # sc.pp.scale(data)
#     # sc.pp.pca(data)
#     # sc.pp.neighbors(data)
#     # sc.tl.diffmap(data)
#     sc.tl.umap(data)

# sc.pl.diffmap(cll, color=['patient_id', 'timepoint', 'log_counts'], show=False)
# sc.pl.diffmap(cll_b, color=['patient_id', 'timepoint', 'log_counts'], show=False)
# sc.pl.diffmap(cll_bp, color=['patient_id', 'timepoint', 'log_counts'], show=False)

# # # On signature genes
# common = pd.read_csv(os.path.join("results", "offtarget_signature.csv"), index_col=0, header=None, squeeze=True)

# cll = cll[:, common.index]
# cll_b = cll_b[:, common.index]
# cll_bp = cll_bp[:, common.index]

# sc.pp.neighbors(cll, use_rep="X")
# sc.pp.neighbors(cll_b, use_rep="X")
# sc.pp.neighbors(cll_bp, use_rep="X")

# sc.tl.diffmap(cll)
# sc.tl.diffmap(cll_b)
# sc.tl.diffmap(cll_bp)

# sc.pl.diffmap(cll, color=['patient_id', 'timepoint', 'log_counts'], show=False)
# sc.pl.diffmap(cll_b, color=['patient_id', 'timepoint', 'log_counts'], show=False)
# sc.pl.diffmap(cll_bp, color=['patient_id', 'timepoint', 'log_counts'], show=False)


# sc.tl.umap(cll)
# sc.tl.umap(cll_b)
# sc.tl.umap(cll_bp)

# sc.pl.umap(cll, color=['patient_id', 'timepoint', 'log_counts'], show=False)
# sc.pl.umap(cll_b, color=['patient_id', 'timepoint', 'log_counts'], show=False)
# sc.pl.umap(cll_bp, color=['patient_id', 'timepoint', 'log_counts'], show=False)


genes = [
    "HMGN2", "CXCR4", "RPL22", "RP11−160E2.6", "H3F3A", "ZFP36L1", "IER2", "NFKBIA",
    "FOS", "LGALS1", "S100A4", "CRIP1", "RPL36A", "RPS28", "RPL17", "ZFP36",
    "DUSP1", "ZFP36L2", "CD27", "CCDC88A", "RPL41", "RPL39", "HIF1A", "CXCL8",
    "CCL3L3", "EEF1A1", "NBEAL1", "CCL3", "AC090498.1", "LAPTM5", "DDX5", "JUN"]

for cell_type in ["CLL", "CD8", "CD4", "Mono"]:
    data = adata[adata.obs['cell_type'] == cell_type, :]
    data.strings_to_categoricals()
    sc.pl.tracksplot(data, var_names=[x for x in genes if x in data.var_names], groupby="timepoint", show=False)
    savefig(plt.gca().figure, os.path.join("results", "scRNA-seq.differential.fig1_diff_genes.tracksplot.{}.svg"
                                           .format(cell_type)))
    plt.close('all')


data = adata[adata.obs['cell_type'] == "CLL", :]

for g in [x for x in genes if x in data.var_names][:4]:
    print(g)
    cell_type = "CLL"
    sc.pl.violin(data, keys=g, groupby="timepoint", use_raw=False, show=False)
    # for child in plt.gca().get_children():
    #     if isinstance(child, matplotlib.collections.PolyCollection):
    #         child.set_rasterized(True)
    savefig(ax.figure, os.path.join("results", "scRNA-seq.differential.fig1_diff_genes.violin.{}.{}.svg"
                                               .format(cell_type, g)))
    plt.close('all')
