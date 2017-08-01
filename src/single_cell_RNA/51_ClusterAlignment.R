require("project.init")
# install.packages("nloptr",repos="https://cloud.r-project.org")
require(Seurat)
require(methods)
project.init2("cll-time_course")


out <- "51_ClusterAlignment/"
dir.create(dirout(out))

# arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Need two arguments: 1 - sample")
}
sample.x <- args[1]

# Hardcoded settings
cluster.resolution <- c(0.8, 1.6)
set.seurat2 <- list(
  PT = list(s1 = "PT_d0", s2 = "PT_d120"),
  PT2 = list(s1 = "PT_d0", s2 = "PT_d280"),
  VZS = list(s1 = "VZS_d0", s2 = "LiveBulk_10x_VZS7_120d"),
  FE = list(s1 = "LiveBulk_10x_FE_FE1_d0_10xTK078", s2 = "LiveBulk_10x_FE7_120d"),
  PBGY = list(s1 = "LiveBulk_10x_PBGY1_0d", s2 = "LiveBulk_10x_PBGY7_150d")
  )

# set up further settings
sample1 <- set.seurat2[[sample.x]]$s1
sample2 <- set.seurat2[[sample.x]]$s2
outS <- paste0(out, sample.x, "/")
dir.create(dirout(outS))

# Process sample1
path1 <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/",sample1,"/outs/","filtered","_gene_bc_matrices/GRCh38/")
obj1 <- Read10X(path1)
colnames(obj1) <- paste0(colnames(obj1), "-1")
obj1 <- CreateSeuratObject(raw.data = obj1)
obj1 <- NormalizeData(object = obj1)
obj1 <- ScaleData(object = obj1)
obj1 <- FindVariableGenes(object = obj1, do.plot = FALSE)

# Process sample2
path2 <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/",sample2,"/outs/","filtered","_gene_bc_matrices/GRCh38/")
obj2 <- Read10X(path2)
colnames(obj2) <- paste0(colnames(obj2), "-2")
obj2 <- CreateSeuratObject(raw.data = obj2)
obj2 <- NormalizeData(object = obj2)
obj2 <- ScaleData(object = obj2)
obj2 <- FindVariableGenes(object = obj2, do.plot = FALSE)

# prep cca
var.genes.union <- union(rownames(head(obj1@hvg.info,2000)), rownames(head(obj2@hvg.info,2000)))

obj1@meta.data[, "sample"] <- gsub("_10xTK078", "", gsub("LiveBulk_10x_", "", sample1))
obj2@meta.data[, "sample"] <- gsub("_10xTK078", "", gsub("LiveBulk_10x_", "", sample2))


# CCA ---------------------------------------------------------------------
pbmc <- RunCCA(object = obj1, object2 = obj2, genes.use = var.genes.union)

pdf(dirout(outS, "CCA_1_9.pdf", width=15, height=15))
DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)
dev.off()

pdf(dirout(outS, "CCA_10_18.pdf", width=15, height=15))
DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 500, dim.use = 10:18, do.balanced = TRUE)
dev.off()

pbmc <- CalcVarExpRatio(object = pbmc, reduction.type = "pca", grouping.var = "sample", dims.use = 1:13)

pbmc.all.save <- pbmc
pbmc <- SubsetData(object = pbmc, subset.name = "var.ratio.pca", accept.low = 0.5)
pbmc.discard <- SubsetData(object = pbmc.all.save, subset.name = "var.ratio.pca", accept.high = 0.5)

pbmc <- AlignSubspace(object = pbmc, reduction.type = "cca", grouping.var = "sample",dims.align = 1:13)

# tsne
pbmc <- RunTSNE(object = pbmc, reduction.use = "cca.aligned", dims.use = 1:13, do.fast = TRUE)
TSNEPlot(object = pbmc, group.by = "sample", do.return = TRUE, pt.size = 0.5)
ggsave(dirout(outS, "sample.pdf"))

# clustering
for(clres in cluster.resolution){
  pbmc <- FindClusters(object = pbmc, resolution=clres, reduction.type = "cca.aligned", dims.use = 1:13, save.SNN = TRUE)
  #   TSNEPlot(object = pbmc, do.return = TRUE, pt.size = 0.5)
  #   ggsave(dirout(outS, "Clusters_",clres,".pdf"))
}

save(pbmc, file=dirout(outS, sample.x, ".RData"))


max.nr.of.cores <- 3

source("src/single_cell_RNA/95_Seurat2.R")