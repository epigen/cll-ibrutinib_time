require("project.init")
require("Seurat")
require("gplots")
require(methods)
require(pheatmap)
require(enrichR) #devtools::install_github("definitelysean/enrichR")


project.init2("cll-time_course")

out <- "13_3_Overtime_Together_tobit/"
dir.create(dirout(out))

sample.x <- "allDataBest_NoDownSampling_noIGH"

inDir <- "11_CellTypes_tobit/"

cells <- list.files(dirout(inDir))
cells <- cells[grepl(paste0(sample.x, "_"), cells)]
cells <- cells[!grepl(".log",cells)]
cells <- gsub(paste0(sample.x, "_"), "", cells)

res <- data.table()

cell <- "Bcells"
for(cell in cells){
  message(cell)  
  datasetName <- paste0(sample.x, "_", cell)
  
  pat.dirs <- list.dirs(dirout(inDir, datasetName))
  pat.dirs <- pat.dirs[grepl("_pat_", pat.dirs)]
  
  eff.sizes <- list()
  
  for(pat.dir in pat.dirs){
    comp.files <- list.files(paste0(pat.dir, "/"))
    comp.files <- comp.files[grepl("Diff_Cluster.+\\.tsv", comp.files)]
    comp.files <- comp.files[!grepl("d280", comp.files)]
    
    for(comp.file in comp.files){
      nams <- strsplit((gsub("Diff_Cluster", "", gsub("\\.tsv", "", comp.file))), "vs")[[1]]
      pat <- gsub("([A-Z]+)\\d?_.+", "\\1", nams[1])
      direction <- ifelse(grepl("d0", gsub("(\\d+)d", "d\\1", nams)[1]), "early_vs_late", "late_vs_early")
      
      
      hits <- fread(paste0(pat.dir, "/",comp.file))[order(V3)]
      colnames(hits) <- c("gene", "pvalue", "logFC", "pct1", "pct2")
      hits[,pct.diff := pct1 - pct2]
      hits$cellType <- cell
      hits$patient <- pat
      
      if(direction == "early_vs_late"){
        hits[,logFC := -logFC]
        hits[,pct.diff := -pct.diff]
      }
      res <- rbind(res, hits[,c("gene", "pvalue", "logFC", "pct.diff", "patient", "cellType")])
    }
  }
}



# PREPARE THE DATA --------------------------------------------------------
res[,qvalue  := p.adjust(pvalue, method="BH")]
write.table(res, file=dirout(out, "SigGenes_overTime.tsv"), quote=F, row.names=F,sep="\t")

res$qvalue2  <- p.adjust(res$pvalue, method="BH")
with(res[pvalue < 0.05], table(patient, cellType))
with(res[qvalue < 0.05], table(patient, cellType))
with(res[qvalue2 < 0.05], table(patient, cellType))


res[,cellPat := paste0(cellType, "_", patient)]
res[,logqval := pmin(5, -log(qvalue))]


res.sig <- res[qvalue < 0.05 & abs(logFC) > 0.3]



# ANALYSIS OVER ALL CELLS -------------------------------------------------
for(cell in c("All", unique(res$cellType))){
  for(filter in c("")){
    res2 <- res.sig
    
    if(cell != "All"){
      res2 <- res.sig[cellType == cell]
    }

    gene.cnt <- res2[,.N, by="gene"]
    gene.cnt <- gene.cnt[order(N, decreasing=TRUE)][!grepl("RP", gene)]
    genes <- gene.cnt[1:min(100, nrow(gene.cnt))]$gene
    res3 <- res2[gene %in% genes]
    
    # Heatmap over all changes
    pDT <- dcast.data.table(res2, gene ~ cellPat, value.var="logFC")
    pMT <- as.matrix(pDT[, -"gene", with=F])
    rownames(pMT) <- pDT$gene
    distx <- dist(pMT[genes,])
    distx[is.na(distx)] <- max(distx,na.rm=T)
    res3$gene <- factor(res3$gene, levels=rownames(pMT[genes,])[hclust(distx)$order])
    
    # Dot plot for top significant genes in all samples
    ggplot(
      res3, 
      aes(y=gene, x=cellPat, color=logFC, size=logqval)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      geom_point() + scale_color_gradient2(low="blue", mid="white", high="red")
    ggsave(dirout(out, cell, "_Genes_Dot",filter,".pdf"), height=length(genes)*0.2+3, width=7)
    
    # Heatmap over all changes
    pDT <- dcast.data.table(res2, gene ~ cellPat, value.var="logFC")
    pMT <- as.matrix(pDT[, -"gene", with=F])
    pCor <- cor(pMT, use="pairwise.complete.obs", method="spearman")
    annot <- data.frame(do.call(rbind, strsplit(colnames(pMT), "_")),row.names=colnames(pMT))
    try({
      pdf(dirout(out, cell, "_CorrHeatmap",filter,".pdf"), onefile=F, height=7, width=8)
      pheatmap(pCor,annotation_row=annot)
      dev.off()
    }, silent=T)

    
    
    # ENRICHR
    enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location")
    enrichRes <- data.table()
    res2[,group := paste0(cellPat, "_", ifelse(logFC > 0, "up", "down"))]
    hitSets <- split(res2$gene, factor(res2$group))
    for(grp.x in names(hitSets)){
      ret=try(as.data.table(enrichGeneList(hitSets[[grp.x]],databases = enrichrDBs)),silent = FALSE)
      if(!any(grepl("Error",ret)) && nrow(ret) > 0){
        enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
      }
    }
    enrichRes$n <- sapply(strsplit(enrichRes$genes,","), length)
    enrichRes <- enrichRes[n > 3]
    write.table(enrichRes[qval < 0.05], file=dirout(out, cell, "_EnrichR",filter, ".tsv"), sep="\t", quote=F, row.names=F)
    
    if(nrow(enrichRes) > 2 & length(unique(enrichRes$grp)) > 1){
      pDat <- dcast.data.table(enrichRes, make.names(category) ~ grp, value.var="qval")
      pDatM <- as.matrix(pDat[,-"category", with=F])
      pDat$category <- gsub("\\_(\\w|\\d){8}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){12}", "", pDat$category)
      pDat$category <- substr(pDat$category, 0, 50)
      row.names(pDatM) <- pDat$category
      pDatM[is.na(pDatM)] <- 1
      str(pDatM <- pDatM[apply(pDatM <= 5e-2,1,sum)>=1,apply(pDatM <= 5e-2,2,sum)>=1, drop=F])
      if(nrow(pDatM) >=2 & ncol(pDatM) >= 2){
        pDatM <- -log10(pDatM)
        pDatM[pDatM > 4] <- 4
        # pDatM[pDatM < 1.3] <- 0
        pdf(dirout(out, cell, "_EnrichR",filter, ".pdf"),onefile=FALSE, width=min(29, 6+ ncol(pDatM)*0.3), height=min(29, nrow(pDatM)*0.3 + 4))
        pheatmap(pDatM) #, color=gray.colors(12, start=0, end=1), border_color=NA)
        dev.off()
      }
    }
  }
}




