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

clustering.precision <- seq(0.5, 2.5, 0.2)

out <- "11_CellTypes_inclDay30/"
dir.create(dirout(out))


cell = "Mono"
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Need arguments: 1 - celltype")
} else {
  cell <- args[1]
}
# if (!cell %in% pbmc@data.info$CellType) {
#   stop("Celltypes not defined ", cell)
# }

sample.x <- cell

# CREATE FOLDER FOR SPECIFIC CELLTYPE ------------------------------------
outS <- paste0(out, cell, "/")
dir.create(dirout(outS))

if(!file.exists(dirout(outS, cell,".RData"))){
  load(dirout("10_Seurat_raw/inclDay30_noIGHLK_negbinom/", "inclDay30_noIGHLK.RData"))
	
  pbmcOrig <- pbmc
  
  stopifnot(all(rownames(pbmc@data.info) == colnames(pbmc@data)))
  
  # SUBSET DATA -------------------------------------------------------------
  cellToKeep.idx <- which(pbmc@data.info$CellType == cell)
  
  stopifnot(!is.na(cellToKeep.idx))
  str(cellToKeep <- rownames(pbmcOrig@data.info)[cellToKeep.idx])
  pbmc <- SubsetData(pbmcOrig, cells.use = cellToKeep)
  
  pDat <- pbmcOrig@tsne.rot
  pDat$selected <- "no"
  pDat$selected[cellToKeep.idx] <- "yes"
  ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=selected)) + geom_point()
  ggsave(dirout(outS, "SelectedCells.pdf"))
  
  
  # ASSIGN PATIENT -------------------------------------------------------------  
  if(is.null(pbmc@data.info[["patient"]])){
    pbmc@data.info[["patient"]] <- gsub("_.*", "", gsub("\\d", "", substr(gsub("LiveBulk_10x_", "", pbmc@data.info[["sample"]]), 0,6)))
    
    for(pat in unique(pbmc@data.info[["patient"]])){
      res <- pbmc@data.info[["sample"]]
      res[!grepl(paste0(pat, "\\d?_(d\\d+|\\d+d)"), pbmc@data.info[["sample"]])] <- "IGNORED"
      pbmc@data.info[[paste0("pat_",pat)]] <- res
    }
  }
  
  
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
  pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = F)
  
  # Clustering
  for(x in clustering.precision){
    pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
    pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
  }
  
  # timepoint zero annotation
  if(is.null(pbmc@data.info[["TZero"]])){
    pbmc@data.info[["TZero"]] <- pbmc@data.info$sample
    pbmc@data.info$TZero[!grepl("_0d", pbmc@data.info$TZero) & !grepl("d0", pbmc@data.info$TZero)] <- "IGNORED"
    #     with(pbmc@data.info, table(sample, TZero))
  }
  
  save(pbmc, file=dirout(outS, cell,".RData"))
} else {
  load(file=dirout(outS, cell,".RData"))
  
  update <- FALSE
    
  if(!.hasSlot(pbmc, "version")){
      pbmc <- UpdateSeuratObject(pbmc)
      update <- TRUE
  }

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

# str(pbmc@pca.rot)
# str(pbmc@data.info)
# pcDat <- data.table(pbmc@pca.rot, sample=pbmc@data.info$sample)
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

pbmc@meta.data <- cbind(
  pbmc@meta.data[,grepl("pat_", colnames(pbmc@meta.data)), drop=F],
  pbmc@meta.data[,c("nGene", "nUMI")]
  )

source("src/single_cell_RNA/FUNC_Seurat2.R", echo=TRUE)