require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(org.Hs.eg.db)
require(scran)

project.init2("cll-time_course")

seurat.min.genes <- 200
mito.cutoff <- 0.15
nGene.cutoff <- 3000
cluster.precision <- seq(0.5,1.5,0.2)

cellranger_filtered <- "filtered"



out <- "10_Seurat/"
dir.create(dirout(out))


sample.x <- "allDataBest_NoDownSampling_noIGH_RmSample"
outS <- paste0(out, sample.x,"/")
dir.create(dirout(outS))
message(sample.x)
if(!file.exists(dirout(outS, sample.x,".RData"))){
  data.path <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/", "allDataBest_NoDownSampling","/outs/",cellranger_filtered,"_gene_bc_matrices_mex/GRCh38/")
  
  # set up seurat object ----------------------------------------------------
  pbmc.data <- Read10X(data.path)
  
  # Remove IGH genes
  allGenes <- rownames(pbmc.data)
  ighGenes <- allGenes[grepl("^IGH", allGenes) & !grepl("^IGHMBP", allGenes)]
  pbmc.data <- pbmc.data[allGenes[!allGenes %in% ighGenes],]
      
  # Make Seurat object
  pbmc <- new("seurat", raw.data = pbmc.data)
  pbmc <- Setup(pbmc, min.cells = 3, min.genes = seurat.min.genes, do.logNormalize = T, total.expr = 1e4, project = sample.x)
  
  # Get samples in aggregated datasets
  pbmc@data.info[["sample"]] <- fread("metadata/Aggregate_best.csv")$library_id[as.numeric(gsub("^[A-Z]+\\-(\\d+)$", "\\1", colnames(pbmc@data)))]
  pbmc@data.info[["sample"]] <- gsub("_10xTK078", "", gsub("LiveBulk_10x_", "", pbmc@data.info[["sample"]]))
  
  # Mito genes --------------------------------------------------------------
  mito.genes <- grep("^MT-", rownames(pbmc@data), value = T)
  percent.mito <- apply(expm1(pbmc@data[mito.genes, ]), 2, sum)/apply(expm1(pbmc@data), 2, sum)
  pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")
  ggplot(pbmc@data.info, aes(y=percent.mito, x=nUMI)) + geom_point() + geom_hline(yintercept=mito.cutoff)
  ggsave(dirout(outS, "MitoPlot.pdf"))
  
  # Number of genes (many genes probably duplets) ---------------------------
  ggplot(pbmc@data.info, aes(x=nUMI, y=nGene)) + geom_point() + geom_hline(yintercept=nGene.cutoff)
  ggsave(dirout(outS, "GenePlot.pdf"))
  
  # Cut genes / mitochondrial -----------------------------------------------
  pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.high = nGene.cutoff)
  pbmc <- SubsetData(pbmc, subset.name = "percent.mito", accept.high = mito.cutoff)
    
  # regress out unwanted sources of variation (UMI normalization), also normalize for mitochondrial content
  pbmc <- RegressOut(pbmc, latent.vars = c("sample", "nUMI", "percent.mito"))
  
  # Identify variable genes from downstream analysis
  pbmc <- MeanVarPlot(pbmc ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
  length(pbmc@var.genes)
  
  # PCA on most variable genes
  pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
  
  # Project other genes? 
  pbmc <- ProjectPCA(pbmc)
  
  # t_SNE
  pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = TRUE,check_duplicates = FALSE)
  # ifelse(grepl("\\_tissueGenes$", sample.x), FALSE, TRUE))
  
  # plot clusters
  # write clusters
  # Clustering
  for(x in cluster.precision){
    pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
    pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
  }
  
  save(pbmc, file=dirout(outS, sample.x,".RData"))
} else {
  load(dirout(outS, sample.x,".RData"))
}

# pbmc@data.info[["sample"]] <- NULL
source("src/single_cell_RNA/10_2_Seurat_Script_3.R")