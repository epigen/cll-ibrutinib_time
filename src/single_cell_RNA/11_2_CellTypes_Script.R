require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("cll-time_course")
out <- "11_CellTypes/"
dir.create(dirout(out))

# NEEDS TO BE SET TO SUBCLUSTER ---------------------------------------------------------
grouping <- list(
  "allDataBest" = list(
    "Monos" = c(10,5),
    "NurseLikeCells" = c(11),
    "Tcells" = c(2,3,7,8),
    "Bcells" = c(0,1,4,6,9,12)
    ),
  "allDataBest_NoDownSampling_noIGH" = list(
    "Monos" = c(5,10),
    "NurseLikeCells" = c(12),
    "Tcells" = c(2,3,6,9),
    "Bcells" = c(0,1,4,7,8,11,13)
    ),
  "allDataBest_noIGH" = list(
    "Monos" = c(5,10,13),
    "NurseLikeCells" = c(11),
    "Tcells" = c(2,3,7,8),
    "Bcells" = c(0,1,4,6,9,12,14)
    )
)

# ARGUMENTS ---------------------------------------------------------------
sample.x <- "allDataBest"
cell <- "Monos"
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Need two arguments: 1 - sample, 2 - celltype")
} else {
  sample.x <- args[1]
  cell <- args[2]
}
if (!sample.x %in% names(grouping)) {
  stop("Celltypes not defined for dataset ", sample.x)
}
if (!cell %in% names(grouping[[sample.x]])) {
  stop("Celltype ",cell, " not defined for dataset ", sample.x)
}



# LOAD DATA FROM ORIGINAL DATASET -----------------------------------------
groups <- grouping[[sample.x]]

# CREATE DATASET FOR SPECIFIC CELLTYPE ------------------------------------
outS <- paste0(out, sample.x, "_", cell, "/")
dir.create(dirout(outS))

if(!file.exists(dirout(outS, cell,".RData"))){
  outOrig <- paste0("10_Seurat/", sample.x,"/")
  load(dirout(outOrig, sample.x,".RData"))
  pbmcOrig <- pbmc
  
  # SUBSET DATA -------------------------------------------------------------
  str(cellToKeep <- rownames(pbmcOrig@data.info)[which(pbmcOrig@data.info[["ClusterNames_0.5"]] %in% as.character(groups[[cell]]))])
  pbmc <- SubsetData(pbmcOrig, cells.use = cellToKeep)
  
  # COUNT CELLS -------------------------------------------------------------
  cellCounts <- rbind(
    data.table(sample=pbmc@data.info$sample, group=cell),
    data.table(sample=pbmcOrig@data.info$sample, group="all")
  )[,.N, by=c("sample", "group")]
  cellCounts <- dcast.data.table(cellCounts, sample ~ group,value.var="N")
  cellCounts[is.na(get(cell)), c(cell) := 0]
  cellCounts[,fraction := get(cell)/all]
  ggplot(cellCounts, aes(x=sample, y=fraction)) + geom_bar(stat="identity") + coord_flip() + theme_bw(24) + ggtitle(cell)
  ggsave(dirout(outS, "Cellcounts.pdf"), height=8, width=8)
  write.table(cellCounts, dirout(outS, "Cellcounts.tsv"), sep="\t", quote=F, row.names=F)
  
  # PREP DATASET ------------------------------------------------------------
  pbmc <- MeanVarPlot(pbmc ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
  pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
  pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = T)
  
  # Clustering
  for(x in c(seq(0.5,0.9,0.1), 0.95)){
    pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
    pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
  }
  
  save(pbmc, file=dirout(outS, cell,".RData"))
} else {
  load(file=dirout(outS, cell,".RData"))
}

source("src/single_cell_RNA/10_2_Seurat_Script_2.R")
