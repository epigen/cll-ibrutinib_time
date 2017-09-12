require("project.init")
require("Seurat")
require("gplots")
require(methods)
require(pheatmap)
# require(enrichR) #devtools::install_github("definitelysean/enrichR")

project.init2("cll-time_course")

source("src/single_cell_RNA/FUNC_Enrichr.R") #devtools::install_github("definitelysean/enrichR")


out <- "13_4_Overtime_noRP/"
dir.create(dirout(out))



sample.x <- "allDataBest_NoDownSampling_noRP"

inDir <- "11_CellTypes_noRP_tobit/"

cells <- list.files(dirout(inDir))
cells <- cells[grepl(paste0(sample.x, "_"), cells)]
cells <- cells[!grepl(".log",cells)]
cells <- gsub(paste0(sample.x, "_"), "", cells)

res <- data.table()

cell <- "NKcells"
for(cell in cells){
  message(cell)  
  datasetName <- paste0(sample.x, "_", cell)
  
  pat.dirs <- list.dirs(dirout(inDir, datasetName))
  pat.dirs <- pat.dirs[grepl("_pat_", pat.dirs)]
  
  eff.sizes <- list()
  
  pat.dir <- pat.dirs[1]
  for(pat.dir in pat.dirs){
    comp.files <- list.files(paste0(pat.dir, "/"))
    comp.files <- comp.files[grepl("Diff_Cluster.+\\.tsv", comp.files)]
    comp.files <- comp.files[!grepl("d280", comp.files)]
    
    comp.file <- comp.files[1]
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


res[,qvalue  := p.adjust(pvalue, method="BH")]
write.table(res, file=dirout(out, "SigGenes_overTime.tsv"), quote=F, row.names=F,sep="\t")


# PREPARE THE DATA --------------------------------------------------------
res$qvalue2  <- p.adjust(res$pvalue, method="BH")
with(res[pvalue < 0.05], table(patient, cellType))
with(res[qvalue < 0.05], table(patient, cellType))
with(res[qvalue2 < 0.05], table(patient, cellType))

res[,cellPat := paste0(cellType, "_", patient)]
res[,logqval := pmin(5, -log(qvalue))]

res.sig <- res[qvalue < 0.05 & abs(logFC) > 0.3]


# COMPARE TO PREVIOUS -----------------------------------------------------
res.old <- fread(dirout("13_3_Overtime_Together_tobit/SigGenes_overTime.tsv"))
res.old <- res.old[qvalue < 0.05 & abs(logFC) > 0.3]
pat <- "FE"
cellT <- "CD4"
res.old.comp <- data.table()
for(pat in unique(res.sig$patient)){
  for(cellT in unique(res.sig[patient == pat]$cellType)){
    s1 <- res.old[patient == pat & cellType == cellT]$gene
    s2 <- res.sig[patient == pat & cellType == cellT]$gene
    res.old.comp <- rbind(res.old.comp, data.table(
      patient = pat,
      cellType = cellT,
      x = paste0(cellT, "_", pat),
      old = sum(!s1 %in% s2),
      new = sum(!s2 %in% s1),
      both = sum(s1 %in% s2)
        ))
  }
}
ggplot(melt(res.old.comp, id.vars=c("patient", "cellType", "x")), 
       aes(x=x, y=value, fill=variable)) + geom_bar(stat="identity", position="stack")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "Comparison_to_old.pdf"))




# ANALYSIS OVER ALL CELLS -------------------------------------------------
cell <- "CD8"
for(cell in c("All", unique(res$cellType))){
  print(cell)
  filter <- ""
  for(filter in c("")){
    res2 <- res.sig
    
    if(cell != "All"){
      res2 <- res.sig[cellType == cell]
    }

#     gene.cnt <- res2[,.N, by="gene"]
#     gene.cnt <- gene.cnt[order(N, decreasing=TRUE)][!grepl("RP", gene)]
#     genes <- gene.cnt[1:min(100, nrow(gene.cnt))]$gene
#     res3 <- res2[gene %in% genes]
#     
#     # Heatmap over all changes
#     pDT <- dcast.data.table(res2, gene ~ cellPat, value.var="logFC")
#     pMT <- as.matrix(pDT[, -"gene", with=F])
#     rownames(pMT) <- pDT$gene
#     distx <- dist(pMT[genes,])
#     distx[is.na(distx)] <- max(distx,na.rm=T)
#     res3$gene <- factor(res3$gene, levels=rownames(pMT[genes,])[hclust(distx)$order])
#     
#     # Dot plot for top significant genes in all samples
#     ggplot(
#       res3, 
#       aes(y=gene, x=cellPat, color=logFC, size=logqval)) +
#       theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#       geom_point() + scale_color_gradient2(low="blue", mid="white", high="red")
#     ggsave(dirout(out, cell, "_Genes_Dot",filter,".pdf"), height=length(genes)*0.2+3, width=7)
#     
#     # Heatmap over all changes
#     pDT <- dcast.data.table(res2, gene ~ cellPat, value.var="logFC")
#     pMT <- as.matrix(pDT[, -"gene", with=F])
#     pCor <- cor(pMT, use="pairwise.complete.obs", method="spearman")
#     annot <- data.frame(do.call(rbind, strsplit(colnames(pMT), "_")),row.names=colnames(pMT))
#     try({
#       pdf(dirout(out, cell, "_CorrHeatmap",filter,".pdf"), onefile=F, height=7, width=8)
#       pheatmap(pCor,annotation_row=annot)
#       dev.off()
#     }, silent=T)

    
    
    # ENRICHR
    #     enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location")
    #     enrichRes <- data.table()
    #     res2[,group := paste0(cellPat, "_", ifelse(logFC > 0, "up", "down"))]
    #     hitSets <- split(res2$gene, factor(res2$group))
    #     for(grp.x in names(hitSets)){
    #       ret=try(as.data.table(enrichGeneList(hitSets[[grp.x]],databases = enrichrDBs)),silent = FALSE)
    #       if(!any(grepl("Error",ret)) && nrow(ret) > 0){
    #         enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
    #       }
    #     }
    #     enrichRes$n <- sapply(strsplit(enrichRes$genes,","), length)
    #     enrichRes <- enrichRes[n > 3]
    #     write.table(enrichRes[qval < 0.05], file=dirout(out, cell, "_EnrichR",filter, ".tsv"), sep="\t", quote=F, row.names=F)
    #     
    #     if(nrow(enrichRes) > 2 & length(unique(enrichRes$grp)) > 1){
    #       pDat <- dcast.data.table(enrichRes, make.names(category) ~ grp, value.var="qval")
    #       pDatM <- as.matrix(pDat[,-"category", with=F])
    #       pDat$category <- gsub("\\_(\\w|\\d){8}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){12}", "", pDat$category)
    #       pDat$category <- substr(pDat$category, 0, 50)
    #       row.names(pDatM) <- pDat$category
    #       pDatM[is.na(pDatM)] <- 1
    #       str(pDatM <- pDatM[apply(pDatM <= 5e-2,1,sum)>=1,apply(pDatM <= 5e-2,2,sum)>=1, drop=F])
    #       if(nrow(pDatM) >=2 & ncol(pDatM) >= 2){
    #         pDatM <- -log10(pDatM)
    #         pDatM[pDatM > 4] <- 4
    #         # pDatM[pDatM < 1.3] <- 0
    #         pdf(dirout(out, cell, "_EnrichR",filter, ".pdf"),onefile=FALSE, width=min(29, 6+ ncol(pDatM)*0.3), height=min(29, nrow(pDatM)*0.3 + 4))
    #         pheatmap(pDatM) #, color=gray.colors(12, start=0, end=1), border_color=NA)
    #         dev.off()
    #       }
    #     }
    
    # ENRICHR
    enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location")
    enrichRes <- data.table()
    res2[,group := paste0(cellPat, "_", ifelse(logFC > 0, "up", "down"))]
    hitSets <- split(res2$gene, factor(res2$group))
    
    for(grp.x in names(hitSets)){
      ret=try(as.data.table(enrichGeneList.oddsRatio(hitSets[[grp.x]],databases = enrichrDBs,fdr.cutoff=NULL)),silent = FALSE)
      if(!any(grepl("Error",ret)) && nrow(ret) > 0){
        enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
      }
    }
    if(nrow(enrichRes) > 0){
      enrichRes <- enrichRes[hitLength > 3]
      write.table(enrichRes[qval < 0.05], file=dirout(out, cell, "_EnrichOR_",filter,".tsv"), sep="\t", quote=F, row.names=F)
      
      enrichRes$category <- gsub("\\_(\\w|\\d){8}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){12}", "", enrichRes$category)
      enrichRes$category <- abbreviate(enrichRes$category, minlength=50) # substr(enrichRes$category2, 0, 50)
      enrichRes[, mLog10Q := pmin(-log10(qval),4)]
      enrichRes
      
      # order terms by similarity (of OR)
      enrichRes[,term := paste0(category, "_", dbLength)]
      if(length(unique(enrichRes$term)) >= 2){
        try({
          orMT <- t(as.matrix(dcast.data.table(enrichRes, grp ~ term, value.var="oddsRatio")[,-"grp",with=F]))
          orMT[is.na(orMT)] <- 1
          hclustObj <- hclust(dist(orMT))
          enrichRes$term <- factor(enrichRes$term, levels=hclustObj$labels[hclustObj$order])
        },silent=T)
      }
      
      # order groups by similarity (of OR)
      if(length(unique(enrichRes$grp)) >= 2){
        try({
          orMT <- t(as.matrix(dcast.data.table(enrichRes, term ~ grp, value.var="oddsRatio")[,-"term",with=F]))
          orMT[is.na(orMT)] <- 1
          hclustObj <- hclust(dist(orMT))
          enrichRes$grp <- factor(enrichRes$grp, levels=hclustObj$labels[hclustObj$order])
        }, silent=T)
      }
      
      # plot
      ggplot(enrichRes[term %in% enrichRes[,.(min(qval)), by="term"][V1 < 0.05]$term], 
             aes(x=grp, y=term, size=log10(oddsRatio), color=mLog10Q)) + 
        geom_point() + scale_color_gradient(low="grey", high="red") + theme_bw(12) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("-log10(q) capped at 4")
      ggsave(dirout(out, cell,  "_EnrichOR_", filter,".pdf"), width=min(29, 6+ length(unique(enrichRes$grp))*0.3), height=min(29, length(unique(enrichRes$category))*0.3 + 4))
  
    }
  }
}




atac.signaling.genes <- fread(dirout("cll-time_course.ligand-receptor_repertoire.CLL.gene_level.sig_only.timepoint_mean.clustermap.csv"))$V1
ramilowski.rec_lig <- fread("~/resources_nfortelny/PairsLigRec_Ramilowski_NatComm_2015.txt")
signaling.genes <- unique(c(atac.signaling.genes, ramilowski.rec_lig$Ligand.ApprovedSymbol, ramilowski.rec_lig$Receptor.ApprovedSymbol)) 
signaling.genes <- unique(c(signaling.genes, unique(res[grepl("^CC[LR]\\d+$", gene) | grepl("^CXC[RL]\\d+$", gene) | grepl("^CD\\d+\\w?$", gene)]$gene)))

res.reclig <- res[gene %in% signaling.genes]
res.reclig[, significant := any(qvalue < 0.05) | abs(logFC) > 0.3, by= "gene"]
res.reclig[, significant2 := any(qvalue < 0.05) | abs(logFC) > 1, by= "gene"]


# order groups by similarity (of OR)
if(length(unique(res.reclig$gene)) >= 2){
  try({
    orMT <- t(as.matrix(dcast.data.table(res.reclig, cellPat ~ gene, value.var="logFC")[,-"cellPat",with=F]))
    orMT[is.na(orMT)] <- 1
    hclustObj <- hclust(dist(orMT))
    res.reclig$gene <- factor(res.reclig$gene, levels=hclustObj$labels[hclustObj$order])
  }, silent=T)
}


ggplot(res.reclig[significant == TRUE], aes(x=cellPat, y=gene, color=logFC, size=logqval)) + geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_gradient2(high="red", mid="white", low="blue")
ggsave(dirout(out, "Receptor_Ligand.pdf"), height=29, width=15)

ggplot(res.reclig[significant2 == TRUE], aes(x=cellPat, y=gene, color=logFC, size=logqval)) + geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_gradient2(high="red", mid="white", low="blue")
ggsave(dirout(out, "Receptor_Ligand_Strict.pdf"), height=25, width=15)