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
from collections import OrderedDict

from tqdm import tqdm

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


def save_model(model):
    try:
        import pickle
        file = os.path.join(
            "dca_models",
            prefix + ".dca-zinb_model.pickle")
        with open(file, 'wb') as f:
            pickle.dump(model, f)
    except pickle.PicklingError:
        pass


def chunks(l, n):
    """
    Partition iterable in chunks of size `n`.

    Attributes:

    :param iterable l:
        Iterable (e.g. list or numpy array).

    :param int n:
        Size of chunks to generate.
    """
    n = max(1, n)
    return list(l[i:i + n] for i in range(0, len(l), n))


def add_tracking_to_anndata():
    from inspect import getcallargs
    from functools import wraps

    def tracker(f):
        @wraps(f)
        def wrapper(*args, **kwds):
            adata = args[0]
            attr = "_processing_path"
            if not hasattr(adata, attr):
                setattr(adata, attr, list())
            setattr(adata, attr, getattr(adata, attr) + [f.__name__])
            return f(*args, **kwds)
        return wrapper

    for module in [sc.pp, sc.tl]:
        for function in filter(lambda x: not x.startswith("_"), module.__dict__.keys()):
            i = 0
            max_i = 1
            failing = True
            while failing:
                i += 1
                if i == max_i:
                    continue
                print(module.__name__, function, i)
                try:
                    args = getcallargs(getattr(module, function), [None] * i).keys()
                except TypeError:
                    continue
                failing = False
                i += 1
            if ("adata" in args, "data" in args, "copy" in args):
                setattr(module, function, tracker(getattr(module, function)))


# add_tracking_to_anndata()


mark1 = ['CD3D', 'CD3G', 'CD4', 'CD8A', 'NKG7', 'CD14', "CST3", "CD79A", "TCL1A", "CD19", 'GZMB', 'CD68', 'CD27', 'CD20', 'CD24', "NCR1"]
mark2 = ["IL32", "IFNG", "IFNGR1", "IFNGR1", "IL4R", "IL4", "JUN", "JUNB", "JUND", "JAK1", "JAK2", "GATA1", "JARID2", "KRAS", "MYC"]
mark3 = ["BTK", "LCK", "E2F4", "CXCR4", "ITGA4", "HBA1", "PTPRC"]
red_mark = ["CST3", "TCL1A", "GZMB", "NKG7", "CD3D"]
pbmc_markers = mark1 + mark2 + mark3


prj = Project(os.path.join("metadata", "project_config.scRNA-seq.yaml"))
data_type = "scRNA-seq"

group_variables = [
    "patient_id", "timepoint", "sex", "binet_stage",
    "number_of_prior_treatments", "ttft", "response_at_120", "response", "batch"]

for sample in prj.samples:
    if sample.protocol == "XXGenomics":
        sample.processed_dir = os.path.join(prj.output_dir, "results", "cellranger_count", sample.BSF_name, "outs")
        sample.matrix_dir_raw = os.path.join(sample.processed_dir, "raw_gene_bc_matrices", "GRCh38")
        sample.matrix_dir_filtered = os.path.join(sample.processed_dir, "filtered_gene_bc_matrices", "GRCh38")
        sample.matrix_raw = os.path.join(sample.processed_dir, "raw_gene_bc_matrices_h5.h5")
        sample.matrix_filtered = os.path.join(sample.processed_dir, "filtered_gene_bc_matrices_h5.h5")

# Read up samples and concatenate AnnData objects
adatas = list()
for sample in prj.samples:
    print(sample)
    s = sc.read_10x_h5(sample.matrix_raw, genome="GRCh38")
    for v in group_variables:
        s.obs.loc[:, v] = getattr(sample, v)
    adatas.append(s)
adata = sc.AnnData.concatenate(*adatas, join="inner", batch_key="sample")
adata.var = adata.var.iloc[:, [0]].rename(columns={"gene_ids-0": "gene_ids"})

sc.write("cll-time_course-scRNA-seq.all_samples.h5ad", adata)


# Start preprocessing
prefix = "cll-time_course-scRNA-seq.all_samples.250-50_filter"

adata_ = sc.pp.filter_cells(
    adata, min_counts=250, copy=True)
sc.pp.filter_genes(adata_, min_counts=50, copy=False)
adata_.obs.to_csv(prefix + ".obs.csv")
adata_.obs = pd.DataFrame(index=adata_.obs.index)
sc.write(prefix + ".h5ad", adata_)

# observe NB vs ZINB
adata_tmp = adata_.copy()
adata_tmp.X = np.array(adata_tmp.X.todense())
fig, axis = plt.subplots(1, 1, figsize=(3, 3))
# plot_mean_dropout(adata_tmp, title="Distribution", ax=axis)
savefig(fig, os.path.join("results", "scRNA-seq.original.distribution_fits.svg"))
del adata_tmp

# Denoise with ZINB
model = sc.pp.dca(
    adata,
    mode="denoise", ae_type="zinb-conddisp",
    return_info=True, return_model=True)

adata.obs.to_csv(prefix + ".dca_denoised-zinb.obs.csv")
adata.obs = pd.DataFrame(index=adata.obs.index)
sc.write(prefix + ".dca_denoised-zinb.h5ad", adata)

# Regress out batch
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
# sc.pp.combat(adata, keys=['batch'])

adata.obs.to_csv(prefix + ".dca_denoised-zinb.obs.csv")
adata.obs = pd.DataFrame(index=adata.obs.index)
sc.write(prefix + ".dca_denoised-zinb.h5ad", adata)
sc.pp.scale(adata)

sc.pp.highly_variable_genes(adata, flavour="cell_ranger", n_top_genes=500)
sc.pp.pca(adata, use_highly_variable=True)
sc.pp.neighbors(adata)
sc.tl.diffmap(adata)
sc.tl.umap(adata)

add_cell_type(adata)

adata.obs.to_csv(prefix + ".dca_denoised-zinb.processed.obs.csv")
adata.obs = pd.DataFrame(index=adata.obs.index)
sc.write(prefix + ".dca_denoised-zinb.processed.h5ad", adata)
adata.obs = pd.read_csv(prefix + ".dca_denoised-zinb.processed.obs.csv", index_col=0)
add_cell_type(adata)

sc.pl.pca(adata, color=["cell_type", "patient_id", "timepoint", "batch"], sort_order=False, show=False)
sc.pl.diffmap(adata, color=["cell_type", "patient_id", "timepoint", "batch"], sort_order=False, show=False)
sc.pl.umap(adata, color=["cell_type", "patient_id", "timepoint", "batch"], sort_order=False, show=False)


adata.obs.loc[:, 'log_counts'] = np.log10(adata.obs['n_counts'])

for method in ['pca', 'diffmap', 'umap']:
    print(method)
    axes = getattr(sc.pl, method)(
        adata,
        color=["cell_type", "patient_id", "timepoint", "batch", "log_counts"],
        ncols=5,
        sort_order=False, show=False)
    savefig(
        axes[0].figure,
        os.path.join(
            "results",
            prefix + ".dca_denoised-zinb.{}.svg"
            .format(method)))
plt.close("all")


# Add cell type
# # de novo:
sc.tl.leiden(adata)
sc.pl.pca(adata, colors='leiden')
