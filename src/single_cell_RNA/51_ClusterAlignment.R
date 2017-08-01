require("project.init")
require(Seurat)
require(methods)
project.init2("cll-time_course")

out <- "51_ClusterAlignment/"
dir.create(dirout(out))


sample1 <- "PT_d0"
path1 <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/",sample1,"/outs/","filtered","_gene_bc_matrices/GRCh38/")
obj1 <- Read10X(path1)
colnames(obj1) <- paste0(colnames(obj1), "-1")
obj1 <- CreateSeuratObject(raw.data = obj1)
obj1 <- NormalizeData(object = obj1)
obj1 <- ScaleData(object = obj1)
obj1 <- FindVariableGenes(object = obj1, do.plot = FALSE)


sample2 <- "PT_d120"
path2 <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/",sample2,"/outs/","filtered","_gene_bc_matrices/GRCh38/")
sample1 <- "PT_d0"
obj2 <- Read10X(path2)
colnames(obj2) <- paste0(colnames(obj2), "-2")
obj2 <- CreateSeuratObject(raw.data = obj2)
obj2 <- NormalizeData(object = obj2)
obj2 <- ScaleData(object = obj2)
obj2 <- FindVariableGenes(object = obj2, do.plot = FALSE)


var.genes.union <- union(rownames(head(obj1@hvg.info,2000)), rownames(head(obj2@hvg.info,2000)))

obj1@meta.data[, "sample"] <- "PT_d0"
obj2@meta.data[, "sample"] <- "PT_d120"


# CCA ---------------------------------------------------------------------
pbmc <- RunCCA(object = obj1, object2 = obj2, genes.use = var.genes.union)

pdf(dirout(out, "CCA_1_9.pdf", width=15, height=15))
DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)
dev.off()

pdf(dirout(out, "CCA_10_18.pdf", width=15, height=15))
DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 500, dim.use = 10:18, do.balanced = TRUE)
dev.off()

pbmc <- CalcVarExpRatio(object = pbmc, reduction.type = "pca", grouping.var = "sample", dims.use = 1:13)

pbmc.all.save <- pbmc
pbmc <- SubsetData(object = pbmc, subset.name = "var.ratio.pca", accept.low = 0.5)
pbmc.discard <- SubsetData(object = pbmc.all.save, subset.name = "var.ratio.pca", accept.high = 0.5)

pbmc <- AlignSubspace(object = pbmc, reduction.type = "cca", grouping.var = "sample",dims.align = 1:13)

pbmc <- RunTSNE(object = pbmc, reduction.use = "cca.aligned", dims.use = 1:13, do.fast = TRUE)
TSNEPlot(object = pbmc, group.by = "sample", do.return = TRUE, pt.size = 0.5)
ggsave(dirout(out, "sample.pdf"))

pbmc <- FindClusters(object = pbmc, reduction.type = "cca.aligned", dims.use = 1:13, save.SNN = TRUE)
TSNEPlot(object = pbmc, do.return = TRUE, pt.size = 0.5)
ggsave(dirout(out, "Clusters.pdf"))


pDat <- data.table(pbmc@dr$tsne@cell.embeddings)
genes <- c("NCR1", "CD14", "CD3D", "NKG7", "CD5", "TCL1A", "CD19", "FCGR3A")
genes %in% rownames(pbmc@data)
for(gene in genes){
  pDat[[gene]] <- pbmc@data[gene,]
  ggplot(pDat, aes_string(x="tSNE_1", y="tSNE_2", color=gene)) + geom_point()
  ggsave(dirout(out, gene, ".pdf"))
}
