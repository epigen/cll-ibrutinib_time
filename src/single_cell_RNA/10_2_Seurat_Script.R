print(dim(pbmc@data))

# # PLOT MARKERS
# message("Plotting markers")
# markers.list = list(
#   standard = c("CD14","CST3","CD3D","NKG7","CD8A","FCGR3A","NCAM1","CD79A","TRDC","CD1C", "GZMB", "MMP9"),
#   macs = c("ITGAM","SIGLEC1","CD68","MRC1","CD163"),
#   liver1 = c("TFRC","SLC20A1","SLC48A1","FLVCR1a","FLVCR1b","FLVCR1","SLC49A1","HBB","HBA1"),
#   liver2 = c("HBG1","HBG2","HBA2","HBD","HBE1","HBZ","EPOR","TFR2","CD34"),
#   liver3 = c("GYPA","GYPB","ALAS1","DMT1","GATA1","ALAS2","TFR2"),
#   bdbioscience1 = c("CD3G", 
#                     "CD4", 
#                     "CD8A",
#                     "CD19", 
#                     "MS4A1", # CD20
#                     "ITGAX", # CD11c
#                     "IL3RA", # CD123
#                     "NCAM1", # CD56 
#                     "CD34"),
#   bdbioscience2 = c("CD14", 
#                     "CD33", 
#                     "CEACAM8", # CD66b 
#                     "ITGA2B", # CD41
#                     "ITGB3", # CD61
#                     "SELP", # CD62
#                     "GYPA", # CD 235a
#                     "MCAM", # CD 146
#                     "EPCAM" # CD 326
#   )
# )
# marker.name <- "macs"
# # Original attempt
# rownames(pbmc@data)[grepl("CEA",rownames(pbmc@data))]
# for(marker.name in names(markers.list)){
#   markers <- markers.list[[marker.name]]
#   message(marker.name, " missing:")
#   print(markers[!markers %in% rownames(pbmc@data)])
#   pdf(dirout(outS,"Markers_",marker.name,".pdf"), height=15, width=15)
#   FeaturePlot(pbmc, markers[markers %in% rownames(pbmc@data)],cols.use = c("grey","blue"))
#   dev.off()
# }

# PLOT MARKERS 2
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

# Plot datasets
if(sample.x %in% c("allData", "allData2")){
  marDat <- data.table(
    pbmc@tsne.rot,
    sample=fread("metadata/Aggregate_all.csv")$library_id[as.numeric(gsub("^[A-Z]+\\-(\\d+)$", "\\1", colnames(pbmc@data)))])
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
  label.x <- paste0("Cluster_", x)
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

# HEATMAP PLOTS
# x <- 0.5
# pDat <- pbmc@scale.data
# pDat <- pDat[apply(pDat!=0, 1, sum) > ncol(pDat)/2,]
# pDat <- pDat[names(sort(apply(pDat, 1, var), decreasing=TRUE)[1:2000]),]
# pDat[pDat > 4] <- 4
# pDat[pDat < -4] <- 4
# for(x in c(seq(0.5,0.9,0.1), 0.95)){
#   if(!file.exists(dirout(outS,"Clusters_",x,"_HM.jpg"))){
#     try({
#       colAnnot = data.frame(Cluster=pbmc@data.info[[paste0("ClusterNames_", x)]], row.names=colnames(pbmc@data))
#       pDat <- pDat[,colnames(pbmc@data)[order(pbmc@data.info[[paste0("ClusterNames_", x)]])]] # sorted by cluster
#       jpeg(dirout(outS,"Clusters_",x,"_HM.jpg"), height=1500, width=1500)
#       pheatmap(pDat, scale="none",show_rownames=FALSE,show_colnames=FALSE,annotation_col=colAnnot,cluster_cols=FALSE,color=colorRampPalette(c("blue", "black","red"))(50))
#       dev.off()
#       jpeg(dirout(outS,"Clusters_",x,"_HM_column.jpg"), height=1500, width=1500)
#       pheatmap(pDat, scale="none",show_rownames=FALSE,show_colnames=FALSE,annotation_col=colAnnot,cluster_cols=TRUE,color=colorRampPalette(c("blue", "black","red"))(50))
#       dev.off()
#     }, silent=TRUE)
#   }
# }


# Markers for each cluster ------------------------------------------------
message("Plotting cluster Markers")
x <- 0.5
for(x in c(seq(0.5,0.9,0.1), 0.95)){
  out.cl <- paste0(outS, "Cluster_", x, "/")
  dir.create(dirout(out.cl))
  pbmc@ident <- factor(pbmc@data.info[[paste0("ClusterNames_", x)]])
  names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
  clusters <- as.character(unique(pbmc@ident))
  cl.i <- "0"
  for(cl.i in clusters){
    if(!file.exists(dirout(out.cl, "Markers_Cluster",cl.i, ".tsv"))){
      cluster.markers <- FindMarkers(pbmc,  ident.1 = cl.i, min.pct = 0.25)    
      pdf(dirout(out.cl,"Markers_Cluster",cl.i,".pdf"), height=15, width=15)
      FeaturePlot(pbmc, row.names(cluster.markers)[1:min(nrow(cluster.markers),9)],cols.use = c("grey","blue"))
      dev.off()
      write.table(cluster.markers, dirout(out.cl, "Markers_Cluster",cl.i, ".tsv"), sep="\t", quote=F, row.names=TRUE)
    }
  }
}


# Differences between clusters --------------------------------------------
message("Plotting cluster differences")
x <- 0.5
for(x in c(seq(0.5,0.9,0.1), 0.95)){
  out.cl <- paste0(outS, "Cluster_", x, "/")
  dir.create(dirout(out.cl))
  pbmc@ident <- factor(pbmc@data.info[[paste0("ClusterNames_", x)]])
  names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
  clusters <- as.character(unique(pbmc@ident))
  i1 <- 1
  i2 <- 2
  for(i1 in 1:(length(clusters)-1)){
    for(i2 in (i1+1):length(clusters)){
      cl1 <- clusters[i1]
      cl2 <- clusters[i2]
      if(!file.exists(dirout(out.cl,"Diff_Cluster",cl1,"vs",cl2,".tsv"))){
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
      }
    }
  }
}


# SECOND WAY OF GETTING CLUSTER MARKERS -----------------------------------
# Those markers are specific from the pairwise comparisons
message("Plotting second type of cluster markers")
x <- 0.5
for(x in c(seq(0.5,0.9,0.1), 0.95)){
  if(!file.exists(dirout(outS, "Cluster",x,"_HM", ".pdf"))){
    out.cl <- paste0(outS, "Cluster_", x, "/")
    pbmc@ident <- factor(pbmc@data.info[[paste0("ClusterNames_", x)]])
    names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
    clusters <- as.character(unique(pbmc@ident))
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