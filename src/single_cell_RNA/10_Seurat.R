require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("cll-time_course")

out <- "10_Seurat/"
dir.create(dirout(out))



sample.x <- "PT_d280"
mito.cutoff <- 0.15
nGene.cutoff <- 3000

# f2 <- list.files(paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/"))
# f[!grepl("\\.log$", f)]
f <- c("PT_d0", "PT_d280", "PT_d120", "VZS_d0", "LiveBulk_10x_FE_FE1_d0_10xTK078", "LiveBulk_10x_FE_FE7_d120_10xTK078","LiveBulk_10x_KI_KI1_d0_10xTK078")
f <- c("LiveBulk_10x_PBGY1_0d", "LiveBulk_10x_FE7_120d", "LiveBulk_10x_PBGY7_150d", "LiveBulk_10x_VZS7_120d")

for(sample.x in f){
  outS <- paste0(out, sample.x,"/")
  dir.create(dirout(outS))
  message(sample.x)
  tryCatch({
    if(!file.exists(dirout(outS, sample.x,".RData"))){
      data.path <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/",sample.x,"/outs/filtered_gene_bc_matrices/GRCh38/")
      
      # set up seurat object ----------------------------------------------------
      pbmc.data <- Read10X(data.path)
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

    # PLOT CLUSTERS
    pDat <- data.table(pbmc@tsne.rot)
    for(x in c(seq(0.5,0.9,0.1), 0.95)){
      pDat[[paste0("Cluster_", x)]] <-pbmc@data.info[[paste0("ClusterNames_", x)]]
      ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color=paste0("Cluster_", x))) + geom_point() + ggtitle(sample.x)
      ggsave(dirout(outS, "Cluster_", x, ".pdf"))
    }
    
    write.table(pDat, dirout(outS,"Cluster.tsv"), sep="\t", quote=F, row.names=F)
    
    x <- 0.5
    for(x in c(seq(0.5,0.9,0.1), 0.95)){
      out.cl <- paste0(outS, "Cluster_", x, "/")
      dir.create(dirout(out.cl))
      pbmc@ident <- factor(pbmc@data.info[[paste0("ClusterNames_", x)]])
      names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
      clusters <- as.character(unique(pbmc@ident))
      cl.i <- "0"
      for(cl.i in clusters){
        cluster.markers <- FindMarkers(pbmc,  ident.1 = cl.i, min.pct = 0.25)    
        pdf(dirout(out.cl,"Markers_Cluster",cl.i,".pdf"), height=15, width=15)
        FeaturePlot(pbmc, row.names(cluster.markers)[1:min(nrow(cluster.markers),9)],cols.use = c("grey","blue"))
        dev.off()
        write.table(cluster.markers, dirout(out.cl, "Markers_Cluster",cl.i, ".tsv"), sep="\t", quote=F, row.names=TRUE)
      }
      i1 <- 1
      i2 <- 2
      for(i1 in 1:(length(clusters)-1)){
        for(i2 in (i1+1):length(clusters)){
          cl1 <- clusters[i1]
          cl2 <- clusters[i2]
          cluster.markers <- FindMarkers(pbmc,  ident.1 = cl1, ident.2 = cl2, min.pct = 0.25)
          mm <- row.names(cluster.markers)
          mm <- mm[mm %in% row.names(pbmc@data)]
          if(length(mm) > 0){
            pdf(dirout(out.cl,"Diff_Cluster",cl1,"vs",cl2,".pdf"), height=15, width=15)
            FeaturePlot(object=pbmc,features.plot=row.names(cluster.markers)[1:min(nrow(cluster.markers),9)],cols.use = c("grey","blue"),
                        cells.use=names(pbmc@ident)[pbmc@ident %in% c(cl1, cl2)])
            dev.off()
          }
          write.table(cluster.markers, dirout(out.cl,"Diff_Cluster",cl1,"vs",cl2,".tsv"), sep="\t", quote=F, row.names=TRUE)
        }
      }
    }
      
    # CLUSTER DIFF GENES
    str(pbmc)
    
    # PLOT MARKERS
    markers <- c("CD14","CST3","CD3D","NKG7","CD8A","FCGR3A","NCAM1","CD79A","TRDC","CD1C", "GZMB", "MMP9")
    
    pdf(dirout(outS,"Markers.pdf"), height=15, width=15)
    FeaturePlot(pbmc, markers[markers %in% rownames(pbmc@data)],cols.use = c("grey","blue"))
    dev.off()
    
  }, error = function(e){
    print(e)
  })
}