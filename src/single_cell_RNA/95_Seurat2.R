# THIS SCRIPT REQUIRES:
# pbmc object (seurat)
# outS (output dir)
# seurat.diff.test (differential test)


if(!exists("enrichrDBs")) enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location")
if(!exists("extra.genes.to.plot")) extra.genes.to.plot <- c()
if(!exists("seurat.diff.test")) seurat.diff.test <- "negbinom"
if(!exists("seurat.thresh.use")) seurat.thresh.use <- 0.25
if(!exists("max.nr.of.cores")) max.nr.of.cores <- 3
if(!exists("seurat.min.pct")) seurat.min.pct <- 0.25


print(dim(pbmc@data))

require(pheatmap)
require(ggplot2)
require(doMC)
require(pryr)

# manage memory and number of tasks used
mem_u <- as.numeric(mem_used())/10**6
cores_u = max(1, min(12, floor(180000/mem_u)-1))
if(is.na(mem_u)) cores_u <- 3
cores_u <- min(cores_u, max.nr.of.cores)
message("Memory used: ", mem_u, " cores: ",cores_u)
registerDoMC(cores=cores_u)


# use colnames from pbmc@meta.data as clusterings
# only those with >= 2 groups with > 1 cell
clusterings <- colnames(pbmc@meta.data)[apply(pbmc@meta.data, 2, function(x) sum(table(x[x!="IGNORED"])>1)>1)]
clusterings <- clusterings[!clusterings %in% c("nUMI", "nGene", "orig.ident", "percent.mito", "var.ratio.pca")]
# reorder to do the kmeans at the end
(clusterings <- c(clusterings[!grepl("res", clusterings)], clusterings[grepl("res", clusterings)]))


# PLOT MARKERS 2
message("Plotting Known marker genes")
if(file.exists("metadata/CellMarkers.csv")){
  markers <- fread("metadata/CellMarkers.csv")[Marker != ""]
  outMarkers <- paste0(outS, "Markers/")
  dir.create(dirout(outMarkers))
  for(i in 1:nrow(markers)){
    if(markers[i]$GeneSymbol %in% rownames(pbmc@data)){
      marDat <- data.table(pbmc@dr$tsne@cell.embeddings, Expression=pbmc@data[markers[i]$GeneSymbol,])
      ggplot(marDat, aes(x=tSNE_1, y=tSNE_2, color=Expression)) + geom_point() + 
        scale_color_gradientn(colors=c("grey", "blue")) + theme(legend.position = 'none') +
        ggtitle(paste0(markers[i]$GeneSymbol, "/", markers[i]$Marker, "\n", markers[i]$CellType))
      ggsave(dirout(outMarkers, markers[i]$GeneSymbol,".pdf"), height=7, width=7)
    }
  }
  for(geneSyn in extra.genes.to.plot){
    marDat <- data.table(pbmc@dr$tsne@cell.embeddings, Expression=pbmc@data[geneSyn,])
    ggplot(marDat, aes(x=tSNE_1, y=tSNE_2, color=Expression)) + geom_point() +
      scale_color_gradientn(colors=c("grey", "blue")) + theme(legend.position = 'none') +
      ggtitle(paste0(geneSyn))
    ggsave(dirout(outMarkers, geneSyn,".pdf"), height=7, width=7)
  }
}

# PLOT UMIS ---------------------------------------------------------------
# message("Plotting UMIs")
# umip <- ggplot(data.table(pbmc@dr$tsne@cell.embeddings, UMIs=pbmc@meta.data$nUMI), aes(x=tSNE_1,y=tSNE_2, color=log10(UMIs))) + 
#   scale_color_gradient(low="blue", high="red") +
#   geom_point() + ggtitle(paste(sample.x, "\n",nrow(pbmc@data), "genes\n", ncol(pbmc@data), "cells")) + theme_bw(24)
# ggsave(dirout(outS, "UMI.pdf"),plot=umip)


# PLOT CLUSTERS
message("Plotting Clusters")
pDat <- data.table(pbmc@dr$tsne@cell.embeddings)
cl.x <- "sample"
for(cl.x in clusterings){
  print(cl.x)
  pDat[[cl.x]] <-pbmc@meta.data[[cl.x]]
  labelCoordinates <- pDat[,.(tSNE_1=median(tSNE_1), tSNE_2=median(tSNE_2)),by=cl.x]
  
  ggplot(pDat[get(cl.x) != "IGNORED"], aes_string(x=cl.x)) + geom_bar() + coord_flip()
  ggsave(dirout(outS, "Cluster_counts_", cl.x, ".pdf"))  
  
  ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color=cl.x)) + geom_point(alpha=0.5) + # ggtitle(sample.x) +
    geom_label(data=labelCoordinates, aes_string(x="tSNE_1", y="tSNE_2", label=cl.x), color="black", alpha=0.5)
  ggsave(dirout(outS, "Cluster_tSNE_", cl.x, ".pdf"))
}
write.table(pDat, dirout(outS,"Cluster.tsv"), sep="\t", quote=F, row.names=F)
clusterCounts <- pDat[,-c("tSNE_1", "tSNE_2"), with=TRUE]
clusterCounts <- do.call(rbind, lapply(names(clusterCounts), function(nam) data.table(clusterCounts[[nam]], nam)))
write.table(clusterCounts[, .N, by=c("V1", "nam")], dirout(outS, "ClusterCounts.tsv"), sep="\t", quote=F,row.names=F)


# Markers for each cluster ------------------------------------------------
message("Plotting cluster Markers")
cl.x <- "sample"
# for(cl.x in clusterings){
foreach(cl.x = clusterings) %dopar% { 
  if(!is.null(pbmc@meta.data[[cl.x]])){
    x <- gsub("ClusterNames_","", cl.x)
    out.cl <- paste0(outS, "Cluster_",x, "/")
    dir.create(dirout(out.cl))
    pbmc@ident <- factor(pbmc@meta.data[[cl.x]])
    names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
    clusters <- names(table(pbmc@ident))[table(pbmc@ident)>1]
    clusters <- clusters[clusters != "IGNORED"]
    cl.i <- "PT_d0"
    for(cl.i in clusters){
      # message(cl.i)
      if(!file.exists(dirout(out.cl, "Markers_Cluster",cl.i, ".tsv"))){
        try({
          cluster.markers <- FindMarkers(pbmc,  ident.1 = cl.i, ident.2 = clusters[clusters != cl.i],
                                         test.use=seurat.diff.test, 
                                         min.pct = seurat.min.pct,
                                         thresh.use=seurat.thresh.use)    
          write.table(cluster.markers, dirout(out.cl, "Markers_Cluster",cl.i, ".tsv"), sep="\t", quote=F, row.names=TRUE)
          
          pdf(dirout(out.cl,"Markers_Cluster",cl.i,".pdf"), height=15, width=15)
          FeaturePlot(pbmc, row.names(cluster.markers)[1:min(nrow(cluster.markers),9)],cols.use = c("grey","blue"))
          dev.off()
        }, silent=T)
      }
    }
  }
}


# Differences between clusters --------------------------------------------
message("Plotting cluster differences")
cl.x <- "ClusterNames_0.5"
# for(cl.x in c("sample", paste0("ClusterNames_",c(seq(0.5,0.9,0.1), 0.95)))){
foreach(cl.x = clusterings) %dopar% {  
  if(!is.null(pbmc@meta.data[[cl.x]])){
    x <- gsub("ClusterNames_","", cl.x)
    out.cl <- paste0(outS, "Cluster_",x, "/")
    dir.create(dirout(out.cl))
    pbmc@ident <- factor(pbmc@meta.data[[cl.x]])
    names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
    clusters <- names(table(pbmc@ident))[table(pbmc@ident)>1]
    clusters <- clusters[clusters != "IGNORED"]
    i1 <- 1
    i2 <- 2
    for(i1 in 1:(length(clusters)-1)){
      for(i2 in (i1+1):length(clusters)){
        cl1 <- clusters[i1]
        cl2 <- clusters[i2]
        if(!file.exists(dirout(out.cl,"Diff_Cluster",cl1,"vs",cl2,".tsv")) | !file.exists(dirout(out.cl,"Diff_Cluster",cl2,"vs",cl1,".tsv"))){
            cluster.markers <- FindMarkers(pbmc,  ident.1 = cl1, ident.2 = cl2, 
                                           test.use=seurat.diff.test, 
                                           min.pct = seurat.min.pct,
                                           thresh.use=seurat.thresh.use)
            write.table(cluster.markers, dirout(out.cl,"Diff_Cluster",cl1,"vs",cl2,".tsv"), sep="\t", quote=F, row.names=TRUE)
            try({
              mm <- row.names(cluster.markers)
              mm <- mm[mm %in% row.names(pbmc@data)]
              if(length(mm) > 0){
                pdf(dirout(out.cl,"Diff_Cluster",cl1,"vs",cl2,".pdf"), height=15, width=15)
                FeaturePlot(object=pbmc,features.plot=row.names(cluster.markers)[1:min(nrow(cluster.markers),9)],cols.use = c("grey","blue"),
                            cells.use=names(pbmc@ident)[pbmc@ident %in% c(cl1, cl2)])
                dev.off()
              }
            }, silent=TRUE)
        }
      }
    }
  }
}


# SECOND WAY OF GETTING CLUSTER MARKERS -----------------------------------
# Those markers are specific from the pairwise comparisons
message("Plotting second type of cluster markers")
cl.x <- "ClusterNames_0.5"
cl.x <- "patient"
for(cl.x in clusterings){
  if(!is.null(pbmc@meta.data[[cl.x]])){
    x <- gsub("ClusterNames_","", cl.x)
    if(!file.exists(dirout(outS, "Cluster_HM_",x, ".pdf"))){
      out.cl <- paste0(outS, "Cluster_",x, "/")
      pbmc@ident <- factor(pbmc@meta.data[[cl.x]])
      names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
      clusters <- names(table(pbmc@ident))[table(pbmc@ident)>1]
      clusters <- clusters[clusters != "IGNORED"]
      i1 <- 1
      i2 <- 2
      allClDiff <- data.table()
      for(i1 in 1:(length(clusters)-1)){
        for(i2 in (i1+1):length(clusters)){
          cl1 <- clusters[i1]
          cl2 <- clusters[i2]
          f1 <- dirout(out.cl,"Diff_Cluster",cl1,"vs",cl2,".tsv")
          f2 <- dirout(out.cl,"Diff_Cluster",cl2,"vs",cl1,".tsv")
          if(file.exists(f1) | file.exists(f2)){
            clusterDiff <- if(file.exists(f1)) fread(f1) else fread(f2)
            namesDiff <- if(file.exists(f1)) c(cl1, cl2) else c(cl2,cl1)
            if(ncol(clusterDiff) == 5){
              colnames(clusterDiff) <- c("gene", "pval", "diff", "pct1", "pct2")
              clusterDiff[diff > 0, c("up", "down") := as.list(namesDiff)]
              clusterDiff[diff < 0, c("up", "down") := as.list(rev(namesDiff))]
              allClDiff <- rbind(allClDiff, clusterDiff)
            }
          }
        }
      }
      i <- 1
      allClDiff[,diff2 := abs(pct1-pct2)]
      # plot for each cluster and write table
      for(i in 1:length(clusters)){
        clx <- clusters[i]
        clMarkers2 <- allClDiff[,.(minDiff = min(diff2),N=.N), by=c("up", "gene")][N == length(clusters)-1][up == clx][order(minDiff, decreasing=TRUE)]
        if(nrow(clMarkers2) > 0){
          pdf(dirout(out.cl,"Markers_Cluster",clx,"_version2.pdf"), height=15, width=15)
          FeaturePlot(object=pbmc,features.plot=clMarkers2[1:min(nrow(clMarkers2),9)]$gene,cols.use = c("grey","blue"))
          dev.off()
          write.table(clMarkers2, dirout(out.cl, "Markers_Cluster",clx, "_version2.tsv"), sep="\t", quote=F, row.names=TRUE)
        }
      }
      # Plot for all clusters heatmap
      cllDiffSummary <- allClDiff[,.(minDiff = min(diff2),N=.N), by=c("up", "gene")][N == length(clusters)-1][,rank := rank(-minDiff), by="up"]
      try({
        pdf(dirout(outS, "Cluster_HM_",x, ".pdf"), width=10, height=min(29, nrow(cllDiffSummary[rank < 20])/10))
        DoHeatmap(pbmc, genes.use=cllDiffSummary[rank < 20][order(up)]$gene,order.by.ident=TRUE,slim.col.label=TRUE,remove.key=TRUE)
        dev.off()
      }, silent=TRUE)
    }
  }
}





# NEW ENRICHR -----------------------------------------------------------------
message("EnrichR on markers")
source("src/single_cell_RNA/FUNC_Enrichr.R") #devtools::install_github("definitelysean/enrichR")
cl.x <- "patient"
cl.x <- clusterings[13]
for(cl.x in clusterings){
  if(!is.null(pbmc@data.info[[cl.x]])){
    x <- gsub("ClusterNames_","", cl.x)
    #     if(!file.exists(dirout(outS, "Enrichr_",x,".tsv"))){
    out.cl <- paste0(outS, "Cluster_",x, "/")
    f <- list.files(dirout(out.cl))
    f <- f[grepl("_version2.tsv$", f)]
    genes <- lapply(f, function(fx) fread(dirout(out.cl, fx)))
    names(genes) <- gsub("Markers_", "", gsub("_version2.tsv","",f))
    genes <- genes[sapply(genes, ncol) == 5]
    genes <- lapply(genes, function(fx) fx[V4 > 0.1]$V3)
    genes <- genes[sapply(genes, length) > 4]
    if(length(genes) > 2){
      enrichRes <- data.table()
      for(grp.x in names(genes)){
        ret=try(as.data.table(enrichGeneList.oddsRatio(genes[[grp.x]],databases = enrichrDBs)),silent = FALSE)
        if(!any(grepl("Error",ret)) && nrow(ret) > 0){
          enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
        }
      }
      
      enrichRes <- enrichRes[hitLength > 3]
      write.table(enrichRes, file=dirout(outS, "EnrichOR_",x,".tsv"), sep="\t", quote=F, row.names=F)
      
      enrichRes$category <- gsub("\\_(\\w|\\d){8}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){12}", "", enrichRes$category)
      enrichRes$category <- abbreviate(enrichRes$category, minlength=50)
      enrichRes[, mLog10Q := pmin(-log10(qval),4)]
      
      # order terms by OR
      enrichRes[,term := paste0(category, "_", dbLength)]
      if(length(unique(enrichRes$term)) >= 2){
        distMT <- dist(t(as.matrix(dcast.data.table(enrichRes, grp ~ term, value.var="oddsRatio")[,-"grp",with=F])))
        distMT[is.na(distMT)] <- 0
        hclustObj <- hclust(distMT)
        enrichRes$term <- factor(enrichRes$term, levels=hclustObj$labels[hclustObj$order])
      }
      
      # order groups by similarity (of OR)
      if(length(unique(enrichRes$grp)) >= 2){
        distMT <- dist(t(as.matrix(dcast.data.table(enrichRes, term ~ grp, value.var="oddsRatio")[,-"term",with=F])))
        distMT[is.na(distMT)] <- 0
        hclustObj <- hclust(distMT)
        enrichRes$grp <- factor(enrichRes$grp, levels=hclustObj$labels[hclustObj$order])
      }
      
      # plot
      ggplot(enrichRes[term %in% enrichRes[,.(min(qval)), by="term"][V1 < 0.05]$term], 
             aes(x=grp, y=term, size=log10(oddsRatio), color=mLog10Q)) + 
        geom_point() + scale_color_gradient(low="white", high="red") + theme_bw(12) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("-log10(q) capped at 4")
      ggsave(dirout(outS, "EnrichOR_",x,".pdf"), width=min(29, 6+ length(unique(enrichRes$grp))*0.3), height=min(29, length(unique(enrichRes$category))*0.3 + 4))
    }
  }
}

message("Analysis pipeline completed successfully!")

