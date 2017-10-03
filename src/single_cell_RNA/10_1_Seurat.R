require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("cll-time_course")

# Either analyze all data (if no args given) or jsut the sample given as args
# message("Analyzing all datasets")
f <- c("PT_d0", "PT_d120", "PT_d280", "VZS_d0", 
       "LiveBulk_10x_FE_FE1_d0_10xTK078","LiveBulk_10x_KI_KI1_d0_10xTK078",
       "LiveBulk_10x_PBGY1_0d", "LiveBulk_10x_FE7_120d", "LiveBulk_10x_PBGY7_150d", "LiveBulk_10x_VZS7_120d", 
       "allData", "allData2", "allDataBest", "allDataBest_NoDownSampling",
       "allDataBest_noIGH", "allDataBest_NoDownSampling_noIGH",
       "allDataBest_NoDownSampling_tissueGenes"
       )

args = commandArgs(trailingOnly=TRUE)
cellranger_filtered <- "raw"
f <- 'inclDay30_noIGHLK'
if (length(args) > 0) {
  cellranger_filtered <- args[1]
  if (length(args) > 1) {
    f <- args[2:length(args)]
  }
}



seurat.min.genes <- 200
seurat.diff.test <- "negbinom"

# This tells Seurat whether to use the filtered or raw data matrices
# cellranger_filtered <- "filtered"
# cellranger_filtered <- "raw"

if(!cellranger_filtered %in% c("raw", "filtered")){
  cellranger_filtered <- "filtered"
}

out <- "10_Seurat/"
if(cellranger_filtered == "raw"){
  out <- "10_Seurat_raw/"
}
dir.create(dirout(out))


# some parameters
mito.cutoff <- 0.15
nGene.cutoff <- 3000

# Loop over samples and do analysis
sample.x <- f[1]
for(sample.x in f){
  outS <- paste0(out, sample.x, "_", seurat.diff.test, "/")
  dir.create(dirout(outS))
  message(sample.x)
  tryCatch({
    if(!file.exists(dirout(outS, sample.x,".RData"))){
      data.path <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/",sample.x,"/outs/",cellranger_filtered,"_gene_bc_matrices/GRCh38/")
      if(!dir.exists(data.path)){ # if this fails use the "_mex" directory (for aggregated datasets)
        data.path <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/",sample.x,"/outs/",cellranger_filtered,"_gene_bc_matrices_mex/GRCh38/")
      }
      if(sample.x == "allData2"){ # in this case the same data is loaded by IGH genes are removed
        data.path <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/", "allData","/outs/",cellranger_filtered,"_gene_bc_matrices_mex/GRCh38/")
      }
      if(grepl("\\_noIGH$", sample.x)){ # in this case the same data is loaded by IGH genes are removed
        data.path <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/", gsub("\\_noIGH$", "", sample.x),"/outs/",cellranger_filtered,"_gene_bc_matrices_mex/GRCh38/")
      }
      if(grepl("\\_noIGHLK$", sample.x)){ # in this case the same data is loaded by IGH genes are removed
        data.path <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/", gsub("\\_noIGHLK$", "", sample.x),"/outs/",cellranger_filtered,"_gene_bc_matrices_mex/GRCh38/")
      }
      if(grepl("\\_noRP$", sample.x)){ # in this case the same data is loaded by IGH genes are removed
        data.path <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/", gsub("\\_noRP$", "", sample.x),"/outs/",cellranger_filtered,"_gene_bc_matrices_mex/GRCh38/")
      }
      if(grepl("\\_noRPstrict$", sample.x)){ # in this case the same data is loaded by IGH genes are removed
        data.path <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/", gsub("\\_noRPstrict$", "", sample.x),"/outs/",cellranger_filtered,"_gene_bc_matrices_mex/GRCh38/")
      }
      if(grepl("\\_tissueGenes$", sample.x)){ # in this case the same data is loaded by IGH genes are removed
        data.path <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/", gsub("\\_tissueGenes$", "", sample.x),"/outs/",cellranger_filtered,"_gene_bc_matrices_mex/GRCh38/")
      }
      
      # set up seurat object ----------------------------------------------------
      pbmc.data <- Read10X(data.path)
      
      # Remove IGH genes
      if((sample.x %in% c("allData2")) | grepl("\\_noIGH$", sample.x)){
        allGenes <- rownames(pbmc.data)
        ighGenes <- allGenes[grepl("^IGH", allGenes) & !grepl("^IGHMBP", allGenes)]
        pbmc.data <- pbmc.data[allGenes[!allGenes %in% ighGenes],]
      }
      
      if(grepl("\\_noIGHLK$", sample.x)){
        allGenes <- rownames(pbmc.data)
        ighGenes <- allGenes[(grepl("^IGH", allGenes) & !grepl("^IGHMBP", allGenes)) | grepl("^IGL", allGenes) | grepl("^IGK", allGenes)]
        pbmc.data <- pbmc.data[allGenes[!allGenes %in% ighGenes],]
      }
      
      # Remove IGH and RP*-* genes
      if(grepl("\\_noRP$", sample.x)){
        allGenes <- rownames(pbmc.data)
        rm.genes <- allGenes[grepl("^RP.*?-.*", allGenes)]
        rm.genes <- c(rm.genes,allGenes[grepl("^IGH", allGenes) & !grepl("^IGHMBP", allGenes)])
        pbmc.data <- pbmc.data[allGenes[!allGenes %in% rm.genes],]
      }
      
      # Remove RP*-* genes
      if(grepl("\\_noRPstrict$", sample.x)){
        allGenes <- rownames(pbmc.data)
        rm.genes <- allGenes[grepl("^RP.*?-.*", allGenes) | grepl("^RP[LS]\\d+[AX]?$", allGenes) | grepl("^MT-", allGenes)]
        pbmc.data <- pbmc.data[allGenes[!allGenes %in% rm.genes],]
      }
      
      # Make Seurat object
      pbmc <- new("seurat", raw.data = pbmc.data)
      pbmc <- Setup(pbmc, min.cells = 3, min.genes = seurat.min.genes, do.logNormalize = T, total.expr = 1e4, project = sample.x)
      
      # Get samples in aggregated datasets
      if(grepl("^allData$", sample.x) | grepl("^allData\\_", sample.x)){
        pbmc@data.info[["sample"]] <- fread("metadata/Aggregate_all.csv")$library_id[as.numeric(gsub("^[A-Z]+\\-(\\d+)$", "\\1", colnames(pbmc@data)))]
        pbmc@data.info[["sample"]] <- gsub("_10xTK078", "", gsub("LiveBulk_10x_", "", pbmc@data.info[["sample"]]))
      }
      if(grepl("^allDataBest", sample.x)){
        pbmc@data.info[["sample"]] <- fread("metadata/Aggregate_best.csv")$library_id[as.numeric(gsub("^[A-Z]+\\-(\\d+)$", "\\1", colnames(pbmc@data)))]
        pbmc@data.info[["sample"]] <- gsub("_10xTK078", "", gsub("LiveBulk_10x_", "", pbmc@data.info[["sample"]]))
      }
      
      if(grepl("^inclDay30", sample.x)){
        pbmc@data.info[["sample"]] <- fread("metadata/Aggregate_inclDay30.csv")$library_id[as.numeric(gsub("^[A-Z]+\\-(\\d+)$", "\\1", colnames(pbmc@data)))]
        pbmc@data.info[["sample"]] <- gsub("_10xTK078", "", gsub("LiveBulk_10x_", "", pbmc@data.info[["sample"]]))
      }
      
      if(grepl("\\_noRP$", sample.x)){
        pbmc <- SubsetData(pbmc, cells.use=row.names(subset(pbmc@data.info, sample != "KI_KI1_d0")))
      }
      
      # Mito genes --------------------------------------------------------------
      mito.genes <- grep("^MT-", rownames(pbmc@data), value = T)
      percent.mito <- colSums(expm1(pbmc@data[mito.genes, ]))/colSums(expm1(pbmc@data))
      pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")
      ggplot(pbmc@data.info, aes(y=percent.mito, x=nUMI)) + geom_point() + geom_hline(yintercept=mito.cutoff)
      ggsave(dirout(outS, "MitoPlot.pdf"))
      
      # Number of genes (many genes probably duplets) ---------------------------
      ggplot(pbmc@data.info, aes(x=nUMI, y=nGene)) + geom_point() + geom_hline(yintercept=nGene.cutoff)
      ggsave(dirout(outS, "GenePlot.pdf"))
      
      # Cut genes / mitochondrial -----------------------------------------------
      pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.high = nGene.cutoff) # are you sure you want to use "accept.high" instead of "accept.low" here?
      pbmc <- SubsetData(pbmc, subset.name = "percent.mito", accept.high = mito.cutoff)
      
      # regress out unwanted sources of variation (UMI normalization), also normalize for mitochondrial content
      pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "percent.mito"))
      
      # Identify variable genes from downstream analysis
      pbmc <- MeanVarPlot(pbmc ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
      length(pbmc@var.genes)
      
      if(grepl("\\_tissueGenes$", sample.x)){
        # PCA on Tissue specific genes
        markerGenes <- fread("metadata/Marker_gene_list_for_single_cell_analysis.csv")
        markerGenes <- unique(do.call(c, markerGenes[2:nrow(markerGenes)]))
        markerGenes <- markerGenes[!grepl("^IGH", markerGenes) & markerGenes != ""]
        pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes[pbmc@var.genes %in% markerGenes], do.print = TRUE, pcs.print = 5, genes.print = 5)
      } else {
        # PCA on most variable genes
        pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
      }
      
      # Project other genes? 
      pbmc <- ProjectPCA(pbmc)
      
      # genes associated
      #VizPCA(pbmc, 1:2)
      
      # Jack Straw significant genes
      # pbmc <- JackStraw(pbmc, num.replicate = 100, do.print = FALSE)
      # JackStrawPlot(pbmc, PCs = 1:12)
      
      save(pbmc, file=dirout(outS, sample.x,".RData"))
      
      # t_SNE
      pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = TRUE,check_duplicates = FALSE)
      # ifelse(grepl("\\_tissueGenes$", sample.x), FALSE, TRUE))
      
      # plot clusters
      # write clusters
      # Clustering
      for(x in seq(0.7,0.9,0.2)){
        pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
        pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
      }
      
      save(pbmc, file=dirout(outS, sample.x,".RData"))
    } else {
      load(dirout(outS, sample.x,".RData"))
      
      # Get samples in aggregated datasets - where it wasn't previously done
      if(grepl("^allData$", sample.x) | grepl("^allData\\_", sample.x)){
        pbmc@data.info[["sample"]] <- fread("metadata/Aggregate_all.csv")$library_id[as.numeric(gsub("^[A-Z]+\\-(\\d+)$", "\\1", colnames(pbmc@data)))]
      }
      if(grepl("^allDataBest$", sample.x) | grepl("^allDataBest\\_", sample.x)){
        pbmc@data.info[["sample"]] <- fread("metadata/Aggregate_best.csv")$library_id[as.numeric(gsub("^[A-Z]+\\-(\\d+)$", "\\1", colnames(pbmc@data)))]
      }
      
      if(!is.null(pbmc@data.info[["sample"]])){
        pbmc@data.info[["sample"]] <- gsub("_10xTK078", "", gsub("LiveBulk_10x_", "", pbmc@data.info[["sample"]]))
        save(pbmc, file=dirout(outS, sample.x,".RData"))
      }
    }
    
    max.nr.of.cores <- 3
    
    # #     pbmc@data.info[["sample"]] <- NULL
    source("src/single_cell_RNA/10_2_Seurat_Script_3.R")
    
  }, error = function(e){
    print(e)
  })
}