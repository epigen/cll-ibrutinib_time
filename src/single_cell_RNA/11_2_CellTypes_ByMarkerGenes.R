require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("cll-time_course")
out <- "11_CellTypes/"
dir.create(dirout(out))

clustering.precision <- seq(0.5, 2.5, 0.2)

# ARGUMENTS ---------------------------------------------------------------
sample.x <- "allDataBest_NoDownSampling_noIGH"
filter = list(
  CD8 = list(genes="CD8A", clusters=c(2,3,6))
  )

cell = "CD8"
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Need two arguments: 1 - filter to use")
} else {
  cell <- args[1]
}
if (!cell %in% names(filter)) {
  stop("Celltypes not defined ", cell)
}


# LOAD DATA FROM ORIGINAL DATASET -----------------------------------------
filter.x <- filter[[cell]]

# CREATE DATASET FOR SPECIFIC CELLTYPE ------------------------------------
outS <- paste0(out, sample.x, "_", cell, "/")
dir.create(dirout(outS))

if(!file.exists(dirout(outS, cell,".RData"))){
  outOrig <- paste0("10_Seurat/", sample.x,"/")
  load(dirout(outOrig, sample.x,".RData"))
  pbmcOrig <- pbmc
  
  stopifnot(all(rownames(pbmc@data.info) == colnames(pbmc@data)))
  
  # SUBSET DATA -------------------------------------------------------------
  cellToKeep.idx <- which(
    pbmcOrig@data.info[["ClusterNames_0.5"]] %in% as.character(filter.x$clusters) & apply(pbmc@data[filter.x$genes,,drop=F]==0,2,sum)==0
    )
  str(cellToKeep <- rownames(pbmcOrig@data.info)[cellToKeep.idx])
  pbmc <- SubsetData(pbmcOrig, cells.use = cellToKeep)
  
  pDat <- pbmcOrig@tsne.rot
  pDat$selected <- "no"
  pDat$selected[cellToKeep.idx] <- "yes"
  ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=selected)) + geom_point()
  ggsave(dirout(outS, "SelectedCells.pdf"))
  
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
  pbmc <- RunTSNE(pbmc, dims.use = 1:11, do.fast = T)
  
  # Clustering
  for(x in clustering.precision){
    pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
    pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
  }
  
  save(pbmc, file=dirout(outS, cell,".RData"))
} else {
  load(file=dirout(outS, cell,".RData"))
  
  update <- FALSE
  for(x in clustering.precision){
    if(is.null(pbmc@data.info[[paste0("ClusterNames_", x)]])){
      update <- TRUE
      pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
      pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
    }
  }
  
  if(update){
    save(pbmc, file=dirout(outS, cell,".RData"))
  }
}

if(is.null(pbmc@data.info[["patient"]])){
  pbmc@data.info[["patient"]] <- gsub("_.*", "", gsub("\\d", "", substr(gsub("LiveBulk_10x_", "", pbmc@data.info[["sample"]]), 0,6)))
  
  for(pat in unique(pbmc@data.info[["patient"]])){
    res <- pbmc@data.info[["sample"]]
    res[!grepl(paste0(pat, "\\d?_(d\\d+|\\d+d)"), pbmc@data.info[["sample"]])] <- "IGNORED"
    pbmc@data.info[[paste0("pat_",pat)]] <- res
  }
}

# READ GENE LISTS ---------------------------------------------------------
geneSetFiles <- list(
  hallmark = "/data/groups/lab_bock/shared/resources/gene_sets/MSigDB/6.0/Human/h.all.v6.0.symbols.gmt",
  immuno = "/data/groups/lab_bock/shared/resources/gene_sets/MSigDB/6.0/Human//c7.all.v6.0.symbols.gmt"
)

balance.barcodes <- TRUE

file.nam <- "hallmark"
source("src/single_cell_RNA/91_Signatures.R", echo=TRUE)

file.nam <- "immuno"
source("src/single_cell_RNA/91_Signatures.R", echo=TRUE)

balance.barcodes <- FALSE

file.nam <- "hallmark"
source("src/single_cell_RNA/91_Signatures.R", echo=TRUE)

file.nam <- "immuno"
source("src/single_cell_RNA/91_Signatures.R", echo=TRUE)

source("src/single_cell_RNA/10_2_Seurat_Script_3.R", echo=FALSE)
source("src/single_cell_RNA/90_fscLVM.R", echo=TRUE)
