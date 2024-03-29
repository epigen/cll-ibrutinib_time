require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("cll-time_course")

# Either analyze all data (if no args given) or jsut the sample given as args
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  message("Analyzing all datasets")
  f <- c("PT_d0", "PT_d120", "PT_d280", "VZS_d0", 
         "LiveBulk_10x_FE_FE1_d0_10xTK078","LiveBulk_10x_KI_KI1_d0_10xTK078",
         "LiveBulk_10x_PBGY1_0d", "LiveBulk_10x_FE7_120d", "LiveBulk_10x_PBGY7_150d", "LiveBulk_10x_VZS7_120d", 
         "allData", "allData2", "allDataBest", "allDataBest_NoDownSampling")
  f <- c("PT_d0", "PT_d280", "VZS_d0", 
         "LiveBulk_10x_FE_FE1_d0_10xTK078","LiveBulk_10x_KI_KI1_d0_10xTK078",
         "LiveBulk_10x_PBGY1_0d", "LiveBulk_10x_FE7_120d", "LiveBulk_10x_PBGY7_150d", "LiveBulk_10x_VZS7_120d", 
         "allData", "allData2", "allDataBest", "allDataBest_NoDownSampling")
} else {
  f <- args[1]
}


# This tells Seurat whether to use the filtered or raw data matrices
# cellranger_filtered <- "filtered"
cellranger_filtered <- "raw"

out <- "10_Seurat/"
if(cellranger_filtered == "raw"){
  out <- "10_Seurat_raw/"
}
dir.create(dirout(out))


# some parameters
mito.cutoff <- 0.15
nGene.cutoff <- 3000

# Loop over samples and do analysis
sample.x <- "PT_d120"
for(sample.x in f){
  outS <- paste0(out, sample.x,"/")
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
      
      # set up seurat object ----------------------------------------------------
      pbmc.data <- Read10X(data.path)
      
      # Remove IGH genes
      if(sample.x %in% c("allData2")){
        allGenes <- rownames(pbmc.data)
        ighGenes <- allGenes[grepl("^IGH", allGenes) & !grepl("^IGHMBP", allGenes)]
        pbmc.data <- pbmc.data[allGenes[!allGenes %in% ighGenes],]
      }
      pbmc <- new("seurat", raw.data = pbmc.data)
      pbmc <- Setup(pbmc, min.cells = 3, min.genes = 200, do.logNormalize = T, total.expr = 1e4, project = sample.x)
      
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
      pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.high = 3000)
      pbmc <- SubsetData(pbmc, subset.name = "percent.mito", accept.high = 0.15)
      
      # regress out unwanted sources of variation (UMI normalization), also normalize for mitochondrial content
      pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "percent.mito"))
      
      # Identify variable genes from downstream analysis
      pbmc <- MeanVarPlot(pbmc ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
      length(pbmc@var.genes)
      
      # PCA on most variable genes
      pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
      
      # Project other genes? 
      pbmc <- ProjectPCA(pbmc)
      
      # genes associated
      #VizPCA(pbmc, 1:2)
      
      # Jack Straw significant genes
      # pbmc <- JackStraw(pbmc, num.replicate = 100, do.print = FALSE)
      # JackStrawPlot(pbmc, PCs = 1:12)
      
      # t_SNE
      pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = T)
      
      # plot clusters
      # write clusters
      # Clustering
      for(x in c(seq(0.5,0.9,0.1), 0.95)){
        pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
        pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
      }
      
      save(pbmc, file=dirout(outS, sample.x,".RData"))
    } else {
      load(dirout(outS, sample.x,".RData"))
    }
    
    source("src/single_cell_RNA/10_2_Seurat_Script.R")
    
  }, error = function(e){
    print(e)
  })
}