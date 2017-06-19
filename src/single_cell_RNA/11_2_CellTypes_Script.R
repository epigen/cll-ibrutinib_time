require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("cll-time_course")

seurat.diff.test <- "tobit"

out <- "11_CellTypes/"
if(seurat.diff.test == "tobit"){
  out <- gsub("/", "_tobit/", out)
}
dir.create(dirout(out))

# NEEDS TO BE SET TO SUBCLUSTER ---------------------------------------------------------
grouping <- list(
  "allDataBest" = list(
    "Monos" = c(10,5),
    "NurseLikeCells" = c(11),
    "Tcells" = c(2,3,7,8),
    "Bcells" = c(0,1,4,6,9,12),
    "clustering" = "ClusterNames_0.5"
    ),
  "allDataBest_NoDownSampling_noIGH" = list(
    "Monos" = c(6,11),
    "NurseLikeCells" = c(16),
    "NKcells" = c(10),
    "TcellsAll" = c(1,2,8,13),
    "Tcells1" = c(14),
    "Tcells2" = c(8,13),
    "Tcells3" = c(1,2),
    "Bcells" = c(0,3,4,5,7,9,15,18),
    "clustering" = "ClusterNames_0.95"
    ),
  "allDataBest_noIGH" = list(
    "Monos" = c(6,12),
    "NurseLikeCells" = c(14),
    "NKcells" = c(10),
    "Tcells1" = c(9),
    "Tcells2" = c(2,7,5),
    "Bcells" = c(0,1,3,4,8,11,13,15),
    "clustering" = "ClusterNames_0.95"
    )
)

clustering.precision <- seq(0.5, 2.5, 0.2)

# ARGUMENTS ---------------------------------------------------------------
sample.x <- "allDataBest_NoDownSampling_noIGH"
cell <- "Bcells"
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
  str(cellToKeep <- rownames(pbmcOrig@data.info)[which(pbmcOrig@data.info[[groups[["clustering"]]]] %in% as.character(groups[[cell]]))])
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
  for(x in clustering.precision){
    if(is.null(pbmc@data.info[[paste0("ClusterNames_", x)]])){
      update <- TRUE
      pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
      pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
    }
  }
  
  if(is.null(pbmc@data.info[["TZero"]])){
    pbmc@data.info[["TZero"]] <- pbmc@data.info$sample
    pbmc@data.info$TZero[!grepl("_0d", pbmc@data.info$TZero) & !grepl("d0", pbmc@data.info$TZero)] <- "IGNORED"
    update <- TRUE
    #     with(pbmc@data.info, table(sample, TZero))
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

# Signature Gene Analysis ---------------------------------------------------------
# geneSetFiles <- list(
#   hallmark = "/data/groups/lab_bock/shared/resources/gene_sets/MSigDB/6.0/Human/h.all.v6.0.symbols.gmt",
#   immuno = "/data/groups/lab_bock/shared/resources/gene_sets/MSigDB/6.0/Human//c7.all.v6.0.symbols.gmt"
# )
# balance.barcodes <- TRUE
# file.nam <- "hallmark"
# source("src/single_cell_RNA/91_Signatures.R", echo=TRUE)
# file.nam <- "immuno"
# source("src/single_cell_RNA/91_Signatures.R", echo=TRUE)
# balance.barcodes <- FALSE
# file.nam <- "hallmark"
# source("src/single_cell_RNA/91_Signatures.R", echo=TRUE)
# file.nam <- "immuno"
# source("src/single_cell_RNA/91_Signatures.R", echo=TRUE)

# analysis of clusters
# seurat.diff.test <- "bimod"
source("src/single_cell_RNA/10_2_Seurat_Script_3.R", echo=TRUE)

# specific over time / t zero analysis
source("src/single_cell_RNA/92_OverTime.R", echo=TRUE)
source("src/single_cell_RNA/93_TimepointZero.R", echo=TRUE)

# fscLVM
# source("src/single_cell_RNA/90_fscLVM.R", echo=TRUE)

print("PIPELINE DONE!")
