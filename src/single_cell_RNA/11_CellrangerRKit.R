library("project.init")
project.init2("cll-time_course")

library(cellrangerRkit)
packageVersion("cellrangerRkit")

out <- "11_CellrangerRKit/"
dir.create(dirout(out))

f <- c("PT_d0", "PT_d280", "PT_d120", "VZS_d0", "LiveBulk_10x_FE_FE1_d0_10xTK078", "LiveBulk_10x_FE_FE7_d120_10xTK078","LiveBulk_10x_KI_KI1_d0_10xTK078")
f <- c("PT_d120", "VZS_d0", "LiveBulk_10x_FE_FE1_d0_10xTK078", "LiveBulk_10x_FE_FE7_d120_10xTK078","LiveBulk_10x_KI_KI1_d0_10xTK078")

m.nam <- "PT_d280"
for(m.nam in f){
  message(m.nam)
  outS <- paste0(out, m.nam, "/")
  dir.create(dirout(outS))
  
  pipestance_path <- paste0("~/projects_shared/cll-time_course/results/cellranger_count/", m.nam)
  gbm <- load_cellranger_matrix(pipestance_path)
  analysis_results<-load_cellranger_analysis_results(pipestance_path)
  
  #   str(exprs(gbm))
  #   str(fData(gbm))
  #   str(pData(gbm))
  #   str(analysis_results)
    
  #   sy <- fData(gbm)$symbol
  #   sy[grepl("NCA", sy)]
  
  use_genes <- get_nonzero_genes(gbm)
  gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
  gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
  #   print(dim(gbm_log))
  
  tsne_proj <- analysis_results$tsne
  
  visualize_umi_counts(gbm,tsne_proj[c("TSNE.1","TSNE.2")],limits=c(3,4),marker_size=0.05)
  ggsave(dirout(outS,"UMIs.pdf"), width=7, height=7)  
  
  markers <- c("CD14","CST3","CD3D","NKG7","CD8A","FCGR3A","NCAM1","CD79A","TRDC","CD1C", "GZMB", "MMP9")
  visualize_gene_markers(gbm_log,markers[markers %in% fData(gbm_log)$symbol],tsne_proj[c("TSNE.1","TSNE.2")],limits=c(0,1.5))
  ggsave(dirout(outS,"Markers.jpg"), width=7, height=7)
  
  n_clu <- 2:10
  km_res <- analysis_results$kmeans
  clu_res <- sapply(n_clu, function(x) km_res[[paste(x,"clusters",sep="_")]]$Cluster)
  colnames(clu_res) <- sapply(n_clu, function(x) paste("kmeans",x,sep="."))
  visualize_clusters(clu_res,tsne_proj[c("TSNE.1","TSNE.2")])
  ggsave(dirout(outS,"Clusters.jpg"), width=7, height=7)
  
  x <- 4
  for(x in n_clu){
    out.cl <- paste0(outS, "Clusters_k", x, "/")
    dir.create(dirout(out.cl))
    cluster_results <- analysis_results$kmeans[[paste(x,"clusters",sep="_")]]
    clusters <- unique(cluster_results$Cluster)
    cluster.markers <- prioritize_top_genes(gbm, cluster_results$Cluster, "sseq", min_mean=0.5)
    cl <- 2
    for(cl in clusters){
      cl.diff <- data.table(cluster.markers[[cl]])[significant == TRUE][order(log2fc)]
      if(nrow(cl.diff) > 0){
        visualize_gene_markers(gbm_log,cl.diff[1:min(nrow(cl.diff), 9)]$gene_name,tsne_proj[c("TSNE.1","TSNE.2")],limits=c(0,1.5))
        ggsave(dirout(out.cl,"Markers_Cluster",cl,".pdf"), height=10, width=10)
        write.table(cl.diff, dirout(out.cl, "Markers_Cluster",cl, ".tsv"), sep="\t", quote=F, row.names=TRUE)
      }
    }
    i1 <- 1
    i2 <- 2
    for(i1 in 1:(length(clusters)-1)){
      for(i2 in (i1+1):length(clusters)){
        cl1 <- clusters[i1]
        cl2 <- clusters[i2]
        message("Cluster ", cl1," vs ",cl2)
        
        de_result <- cellrangerRkit:::sseq_differential_expression(t(exprs(gbm)), 
                            cluster_results$Cluster == cl1,cluster_results$Cluster==cl2, 
                            cellrangerRkit:::compute_sseq_params(t(exprs(gbm))), 
                            gene_ids = fData(gbm)$id, gene_symbols = fData(gbm)$symbol)
        de_result$significant <- with(de_result, common_mean >= 0.5 & p_adj < 0.05)
        de_result <- data.table(de_result)
          
        markers <- de_result[significant == TRUE][gene_name %in% fData(gbm_log)$symbol]
        if(nrow(markers) > 0){
          #visualize_gene_markers(gbm_log,markers[1:min(nrow(markers), 9)]$gene_name,tsne_proj[,c("TSNE.1","TSNE.2")],limits=c(0,1.5))
          
          # Recoded to limit to specific cells:
          limits <- c(0,1.5)
          gene_values <- as.matrix(exprs(gbm)[markers[1:min(nrow(markers), 9)]$gene_id,cluster_results$Cluster %in% c(cl1, cl2)])
          if(nrow(markers) > 1){
            gene_values <- t(gene_values)
          }
          
          gene_values[gene_values < limits[1]] <- limits[1]
          gene_values[gene_values > limits[2]] <- limits[2]
          colnames(gene_values) <- markers[1:min(nrow(markers), 9)]$gene_name
          projection <- tsne_proj[cluster_results$Cluster %in% c(cl1, cl2),c("TSNE.1","TSNE.2")]
          projection_names <- colnames(projection)
          colnames(projection) <- c("Component.1", "Component.2")
          proj_gene <- data.frame(cbind(projection, gene_values))
          proj_gene_melt <- melt(proj_gene, id.vars = c("Component.1", "Component.2"))
          p <- ggplot(proj_gene_melt, aes(Component.1, Component.2)) + 
            geom_point(aes(colour = value), size = 0.1) + 
            facet_wrap(~variable) + scale_colour_gradient(low = "grey", high = "red", name = "val") + labs(x = projection_names[1], y = projection_names[2])
          p + theme_bw() + theme(plot.title = element_text(hjust = 0.5),  panel.grid.major = element_blank(), panel.grid.minor = element_blank())          
          ggsave(dirout(out.cl,"Diff_Cluster",cl1,"vs",cl2,".pdf"), height=10, width=10)
          
          write.table(markers, dirout(out.cl,"Diff_Cluster",cl1,"vs",cl2,".tsv"), sep="\t", quote=F, row.names=TRUE)
        }
      }
    }
  }
}