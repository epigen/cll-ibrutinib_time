require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("cll-time_course")

R.version

seurat.diff.test <- "negbinom"
seurat.min.pct <- 0.1
seurat.thresh.use <- 0.1

clustering.precision <- seq(0.5, 2.5, 0.4)

out <- "25_Patient_CLL/"
dir.create(dirout(out))


cell = "PBGY"
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Need arguments: 1 - celltype")
} else {
  cell <- args[1]
}

sample.x <- cell

# CREATE FOLDER FOR SPECIFIC CELLTYPE ------------------------------------
outS <- paste0(out, cell, "/")
dir.create(dirout(outS))

if(!file.exists(dirout(outS, cell,".RData"))){
  (load(dirout("11_CellTypes_inclDay30/CLL/CLL.RData")))
  
  pbmcOrig <- pbmc
  
  stopifnot(all(rownames(pbmc@meta.data) == colnames(pbmc@data)))
  
  # SUBSET DATA -------------------------------------------------------------
  cellToKeep.idx <- which(pbmc@meta.data$patient == cell)
  
  stopifnot(!is.na(cellToKeep.idx))
  str(cellToKeep <- rownames(pbmcOrig@meta.data)[cellToKeep.idx])
  pbmc <- SubsetData(pbmcOrig, cells.use = cellToKeep)
  
  pbmc@meta.data <- pbmc@meta.data[,!grepl("ClusterNames", colnames(pbmc@meta.data))]
  pbmc@meta.data <- pbmc@meta.data[,!grepl("res", colnames(pbmc@meta.data))]
  str(pbmc@meta.data)
  
  pDat <- data.table(pbmcOrig@dr$tsne@cell.embeddings)
  pDat$selected <- "no"
  pDat$selected[cellToKeep.idx] <- "yes"
  ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=selected)) + geom_point()
  ggsave(dirout(outS, "SelectedCells.pdf"))
  
  
  # COUNT CELLS -------------------------------------------------------------
  cellCounts <- rbind(
    data.table(sample=pbmc@meta.data$sample, group=cell),
    data.table(sample=pbmcOrig@meta.data$sample, group="all")
  )[,.N, by=c("sample", "group")]
  cellCounts <- dcast.data.table(cellCounts, sample ~ group,value.var="N")
  cellCounts[is.na(get(cell)), c(cell) := 0]
  cellCounts[,fraction := get(cell)/all]
  ggplot(cellCounts, aes(x=sample, y=fraction)) + geom_bar(stat="identity") + coord_flip() + theme_bw(24) + ggtitle(cell)
  ggsave(dirout(outS, "Cellcounts.pdf"), height=8, width=8)
  write.table(cellCounts, dirout(outS, "Cellcounts.tsv"), sep="\t", quote=F, row.names=F)
  
  # PREP DATASET ------------------------------------------------------------
  pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
  pbmc <- FindVariableGenes(pbmc ,fxn.x = expMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  pbmc <- RunPCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5)
  pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = TRUE)
  
  # Clustering
  for(x in clustering.precision){
    pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
    pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
  }
  
  save(pbmc, file=dirout(outS, cell,".RData"))
} else {
  load(file=dirout(outS, cell,".RData"))
  
  if(!.hasSlot(pbmc, "version")){
      pbmc <- UpdateSeuratObject(pbmc)
      update <- TRUE
  }
  
  update <- FALSE
  for(x in clustering.precision){
    if(is.null(pbmc@meta.data[[paste0("ClusterNames_", x)]])){
      update <- TRUE
      pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
      pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
    }
  }
  
  if(update){
    save(pbmc, file=dirout(outS, cell,".RData"))
  }
}
# 
# str(pbmc@pca.rot)
# str(pbmc@meta.data)
# pcDat <- data.table(pbmc@pca.rot, sample=pbmc@meta.data$sample)
# pcDat2 <- melt(pcDat,id.vars="sample")[variable %in% paste0("PC", 1:10)]
# ggplot(pcDat2, aes(x=sample, y=value)) + geom_violin() + facet_grid(variable ~ .) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ggsave(dirout(out, "PC_Distr_",cell,".pdf"), height=29, width=29)
#   
# for(pc in unique(pcDat2$variable)){
#   ggplot(pcDat2[variable == pc], aes(x=sample, y=value)) + geom_violin() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   ggsave(dirout(out, "PC_Distr_",cell,"_", pc,".pdf"), height=15, width=15)
# }

pbmc@meta.data <- pbmc@meta.data[,c("ClusterNames_0.5", "nUMI", "nGene")]

# Cluster analysis
source("src/single_cell_RNA/FUNC_Seurat2.R", echo=TRUE)