print(dim(pbmc@data))

require(pheatmap)
require(ggplot2)
require(doMC)

registerDoMC(cores=9)
enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location")
clusterings <- colnames(pbmc@data.info)
clusterings <- clusterings[!grepl("res", clusterings) & !clusterings %in% c("nUMI", "nGene", "orig.ident", "percent.mito")]

# PLOT MARKERS 2
message("Plotting Known marker genes")
if(file.exists("metadata/CellMarkers.csv")){
  markers <- fread("metadata/CellMarkers.csv")[Marker != ""]
  outMarkers <- paste0(outS, "Markers/")
  dir.create(dirout(outMarkers))
  for(i in 1:nrow(markers)){
    if(markers[i]$GeneSymbol %in% rownames(pbmc@data)){
      marDat <- data.table(pbmc@tsne.rot, Expression=FetchData(object=pbmc,vars.all=markers[i]$GeneSymbol, use.imputed=FALSE)[,1])
      ggplot(marDat, aes(x=tSNE_1, y=tSNE_2, color=Expression)) + geom_point() + 
        scale_color_gradientn(colors=c("grey", "blue")) + theme(legend.position = 'none') +
        ggtitle(paste0(markers[i]$GeneSymbol, "/", markers[i]$Marker, "\n", markers[i]$CellType))
      ggsave(dirout(outMarkers, markers[i]$GeneSymbol,".pdf"), height=7, width=7)
    }
  }
}

# PLOT UMIS ---------------------------------------------------------------
message("Plotting UMIs")
umip <- ggplot(data.table(pbmc@tsne.rot, UMIs=pbmc@data.info$nUMI), aes(x=tSNE_1,y=tSNE_2, color=log10(UMIs))) + 
  scale_color_gradient(low="blue", high="red") +
  geom_point() + ggtitle(paste(sample.x, "\n",nrow(pbmc@data), "genes\n", ncol(pbmc@data), "cells")) + theme_bw(24)
ggsave(dirout(outS, "UMI.pdf"),plot=umip)


# AGGREGATED DATA ANALYSIS ------------------------------------------------
# Plot samples
message("Plotting samples per dataset")
if(!is.null(pbmc@data.info[["sample"]])){
  marDat <- data.table(
    pbmc@tsne.rot,
    sample=pbmc@data.info[["sample"]])
  ggplot(marDat, aes(x=sample)) + geom_bar() + coord_flip()
  ggsave(dirout(outS, "Samples_counts.pdf"), height=7, width=7)
  ggplot(marDat, aes(x=tSNE_1, y=tSNE_2, color=sample)) + geom_point()
  ggsave(dirout(outS, "Samples_tSNE.pdf"), height=7, width=10)
  marDat$sample2 <- gsub("(\\d|\\_)", "", substr(gsub("LiveBulk_10x_", "", marDat$sample), 0,3))
  ggplot(marDat, aes(x=tSNE_1, y=tSNE_2, color=sample2)) + geom_point(alpha=0.5)
  ggsave(dirout(outS, "Samples2_tSNE.pdf"), height=7, width=7)
}

# PLOT CLUSTERS
message("Plotting Clusters")
pDat <- data.table(pbmc@tsne.rot)
for(x in c(seq(0.5,0.9,0.1), 0.95)){
  label.x <- paste0("Cluster_",x)
  pDat[[label.x]] <-pbmc@data.info[[paste0("ClusterNames_", x)]]
  labelCoordinates <- pDat[,.(tSNE_1=median(tSNE_1), tSNE_2=median(tSNE_2)),by=label.x]
  ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color=label.x)) + geom_point() + ggtitle(sample.x) +
    geom_label(data=labelCoordinates, aes_string(x="tSNE_1", y="tSNE_2", label=label.x), color="black", alpha=0.5)
  ggsave(dirout(outS, label.x, ".pdf"))
}
write.table(pDat, dirout(outS,"Cluster.tsv"), sep="\t", quote=F, row.names=F)
clusterCounts <- pDat[,-c("tSNE_1", "tSNE_2"), with=TRUE]
clusterCounts <- do.call(rbind, lapply(names(clusterCounts), function(nam) data.table(clusterCounts[[nam]], nam)))
write.table(clusterCounts[, .N, by=c("V1", "nam")], dirout(outS, "ClusterCounts.tsv"), sep="\t", quote=F,row.names=F)


# Markers for each cluster ------------------------------------------------
message("Plotting cluster Markers")
cl.x <- "ClusterNames_0.95"
for(cl.x in c("sample", paste0("ClusterNames_",c(seq(0.5,0.9,0.1), 0.95)))){
  if(!is.null(pbmc@data.info[[cl.x]])){
    x <- gsub("ClusterNames_","", cl.x)
    out.cl <- paste0(outS, "Cluster_",x, "/")
    dir.create(dirout(out.cl))
    pbmc@ident <- factor(pbmc@data.info[[cl.x]])
    names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
    clusters <- names(table(pbmc@ident))[table(pbmc@ident)>1]
    clusters <- clusters[clusters != "IGNORED"]
    cl.i <- "0"
    for(cl.i in clusters){
      # message(cl.i)
      if(!file.exists(dirout(out.cl, "Markers_Cluster",cl.i, ".tsv"))){
        try({
          cluster.markers <- FindMarkers(pbmc,  ident.1 = cl.i, ident.2 = clusters[clusters != cl.i], min.pct = 0.25)    
          pdf(dirout(out.cl,"Markers_Cluster",cl.i,".pdf"), height=15, width=15)
          FeaturePlot(pbmc, row.names(cluster.markers)[1:min(nrow(cluster.markers),9)],cols.use = c("grey","blue"))
          dev.off()
          write.table(cluster.markers, dirout(out.cl, "Markers_Cluster",cl.i, ".tsv"), sep="\t", quote=F, row.names=TRUE)
        }, silent=T)
      }
    }
  }
}


# Differences between clusters --------------------------------------------
message("Plotting cluster differences")
cl.x <- "ClusterNames_0.5"
# for(cl.x in c("sample", paste0("ClusterNames_",c(seq(0.5,0.9,0.1), 0.95)))){
foreach(cl.x = c("sample", paste0("ClusterNames_",c(seq(0.5,0.9,0.1), 0.95)))) %dopar% {  
  if(!is.null(pbmc@data.info[[cl.x]])){
    x <- gsub("ClusterNames_","", cl.x)
    out.cl <- paste0(outS, "Cluster_",x, "/")
    dir.create(dirout(out.cl))
    pbmc@ident <- factor(pbmc@data.info[[cl.x]])
    names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
    clusters <- names(table(pbmc@ident))[table(pbmc@ident)>1]
    clusters <- clusters[clusters != "IGNORED"]
    i1 <- 1
    i2 <- 2
    for(i1 in 1:(length(clusters)-1)){
      for(i2 in (i1+1):length(clusters)){
        cl1 <- clusters[i1]
        cl2 <- clusters[i2]
        if(!file.exists(dirout(out.cl,"Diff_Cluster",cl1,"vs",cl2,".tsv")) & !file.exists(dirout(out.cl,"Diff_Cluster",cl2,"vs",cl1,".tsv"))){
          try({
            cluster.markers <- FindMarkers(pbmc,  ident.1 = cl1, ident.2 = cl2, min.pct = 0.25)
            mm <- row.names(cluster.markers)
            mm <- mm[mm %in% row.names(pbmc@data)]
            if(length(mm) > 0){
              pdf(dirout(out.cl,"Diff_Cluster",cl1,"vs",cl2,".pdf"), height=15, width=15)
              FeaturePlot(object=pbmc,features.plot=row.names(cluster.markers)[1:min(nrow(cluster.markers),9)],cols.use = c("grey","blue"),
                          cells.use=names(pbmc@ident)[pbmc@ident %in% c(cl1, cl2)])
              dev.off()
              write.table(cluster.markers, dirout(out.cl,"Diff_Cluster",cl1,"vs",cl2,".tsv"), sep="\t", quote=F, row.names=TRUE)
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
cl.x <- "sample"
for(cl.x in c("sample", paste0("ClusterNames_",c(seq(0.5,0.9,0.1), 0.95)))){
  if(!is.null(pbmc@data.info[[cl.x]])){
    x <- gsub("ClusterNames_","", cl.x)
    if(!file.exists(dirout(outS, "Cluster",x,"_HM", ".pdf"))){
      out.cl <- paste0(outS, "Cluster_",x, "/")
      pbmc@ident <- factor(pbmc@data.info[[cl.x]])
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
      pdf(dirout(outS, "Cluster",x,"_HM", ".pdf"), width=10, height=min(29, nrow(cllDiffSummary[rank < 20])/10))
      DoHeatmap(pbmc, genes.use=cllDiffSummary[rank < 20][order(up)]$gene,order.by.ident=TRUE,slim.col.label=TRUE,remove.key=TRUE)
      dev.off()
    }
  }
}





# ENRICHR -----------------------------------------------------------------
message("EnrichR on markers")
library(enrichR) #devtools::install_github("definitelysean/enrichR")
cl.x <- "sample"
cl.x <- "ClusterNames_0.5"
for(cl.x in c("sample", paste0("ClusterNames_",c(seq(0.5,0.9,0.1), 0.95)))){
  if(!is.null(pbmc@data.info[[cl.x]])){
    x <- gsub("ClusterNames_","", cl.x)
    if(!file.exists(dirout(outS, cl.x, "_Enrichr.tsv"))){
      out.cl <- paste0(outS, "Cluster_",x, "/")
      enrichRes <- data.table()
      f <- list.files(dirout(out.cl))
      f <- f[grepl("_version2.tsv$", f)]
      genes <- lapply(f, function(fx) fread(dirout(out.cl, fx)))
      names(genes) <- gsub("Markers_", "", gsub("_version2.tsv","",f))
      genes <- genes[sapply(genes, ncol) == 5]
      genes <- lapply(genes, function(fx) fx[V4 > 0.1]$V3)
      genes <- genes[sapply(genes, length) > 4]
      if(length(genes) > 2){
        for(grp.x in names(genes)){
          ret=try(as.data.table(enrichGeneList(genes[[grp.x]],databases = enrichrDBs)),silent = FALSE)
          if(!any(grepl("Error",ret)) && nrow(ret) > 0){
            enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
          }
        }
        
        enrichRes$n <- sapply(strsplit(enrichRes$genes,","), length)
        enrichRes <- enrichRes[n > 3][qval < 0.05]
        enrichRes$category <- gsub("\\_(\\w|\\d){8}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){12}", "", enrichRes$category)
        enrichRes$category <- make.unique(substr(enrichRes$category, 0, 50))
        write.table(enrichRes, file=dirout(outS, cl.x, "_Enrichr.tsv"), sep="\t", quote=F, row.names=F)
              
        if(nrow(enrichRes) > 2 & length(unique(enrichRes$grp)) > 1){
          pDat <- dcast.data.table(enrichRes, make.names(category) ~ grp, value.var="qval")
          pDatM <- as.matrix(pDat[,-"category", with=F])
          row.names(pDatM) <- pDat$category
          pDatM[is.na(pDatM)] <- 1
          str(pDatM <- pDatM[apply(pDatM <= 10e5,1,sum)>1,])
          str(pDatM <- pDatM[,apply(pDatM <= 10e5,2,sum)>1])
          pDatM <- -log10(pDatM)
          pDatM[pDatM > 4] <- 4
          pDatM[pDatM < 1.3] <- 0
          pdf(dirout(outS, cl.x, "_Enrichr.pdf"),onefile=FALSE, width=min(29, 6+ ncol(pDatM)*0.3), height=min(29, nrow(pDatM)*0.3 + 4))
          pheatmap(pDatM) #, color=gray.colors(12, start=0, end=1), border_color=NA)
          dev.off()
        }
      }
    }
  }
}

message("Analysis pipeline completed successfully!")

