require("project.init")
require("Seurat")
require("gplots")
require(methods)
require(pheatmap)
# require(enrichR) #devtools::install_github("definitelysean/enrichR")

project.init2("cll-time_course")

source("src/single_cell_RNA/FUNC_Enrichr.R") #devtools::install_github("definitelysean/enrichR")



list.files(dirout())

sample.x = "nUMI_Cutoff"
inDir <- "11_2_3_CellTypes_UMI_Cutoff/"

list.files(dirout(inDir))

(out <- paste0("13_4_Overtime_",sample.x,"/"))
dir.create(dirout(out))



(cells <- list.files(dirout(inDir)))
(cells <- cells[!grepl(".log",cells) & !grepl("PC_Distr",cells)])

res <- data.table()

cell <- "NK"
for(cell in cells){
  message(cell)  
  
  pat.dirs <- list.dirs(dirout(inDir, cell))
  (pat.dirs <- pat.dirs[grepl("_pat_", pat.dirs)])
  
  eff.sizes <- list()
  
  (pat.dir <- pat.dirs[3])
  print(pat.dirs)
  for(pat.dir in pat.dirs){
    comp.files <- list.files(paste0(pat.dir, "/"))
    comp.files <- comp.files[grepl("Diff_Cluster.+\\.tsv", comp.files)]
    (comp.files <- comp.files[grepl("d0", comp.files)])
    
    comp.file <- comp.files[1]
    for(comp.file in comp.files){
      (nams <- strsplit((gsub("Diff_Cluster", "", gsub("\\.tsv", "", comp.file))), "vs")[[1]])
      (pat <- gsub("([A-Z]+)\\d?_.+", "\\1", nams[1]))
      (direction <- ifelse(grepl("d0", gsub("(\\d+)d", "d\\1", nams)[1]), "early_vs_late", "late_vs_early"))
      
      
      hits <- fread(paste0(pat.dir, "/",comp.file))[order(V3)]
      colnames(hits) <- c("gene", "pvalue", "logFC", "pct1", "pct2")
      hits[,pct.diff := pct1 - pct2]
      hits$cellType <- cell
      hits$patient <- ifelse(direction == "early_vs_late", nams[2], nams[1])
      
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

res.sig[,direction := ifelse(logFC > 0, "up", "down")]
res.sig$direction <- factor(res.sig$direction, levels = c("up","down"))
res.sig[,patient2 := gsub("d30", "d030", patient)]

ggplot(res.sig, aes(x=patient2, fill=direction)) + geom_bar(position="dodge") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ cellType) + ylab("Number of genes") + xlab("")
ggsave(dirout(out, "NumberOfGenesChanging.pdf"), height=7, width=14)

# COMPARE TO PREVIOUS -----------------------------------------------------
# res.old <- fread(dirout("13_3_Overtime_Together_tobit/SigGenes_overTime.tsv"))
# res.old <- res.old[qvalue < 0.05 & abs(logFC) > 0.3]
# pat <- "FE"
# cellT <- "CD4"
# res.old.comp <- data.table()
# for(pat in unique(res.sig$patient)){
#   for(cellT in unique(res.sig[patient == pat]$cellType)){
#     s1 <- res.old[patient == pat & cellType == cellT]$gene
#     s2 <- res.sig[patient == pat & cellType == cellT]$gene
#     res.old.comp <- rbind(res.old.comp, data.table(
#       patient = pat,
#       cellType = cellT,
#       x = paste0(cellT, "_", pat),
#       old = sum(!s1 %in% s2),
#       new = sum(!s2 %in% s1),
#       both = sum(s1 %in% s2)
#         ))
#   }
# }
# ggplot(melt(res.old.comp, id.vars=c("patient", "cellType", "x")), 
#        aes(x=x, y=value, fill=variable)) + geom_bar(stat="identity", position="stack")+ 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ggsave(dirout(out, "Comparison_to_old.pdf"))




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
    enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "Transcription_Factor_PPIs")
    enrichRes <- data.table()
    res2[,group := paste0(cellPat, "_", ifelse(logFC > 0, "up", "down"))]
    hitSets <- split(res2$gene, factor(res2$group))
    
    for(grp.x in names(hitSets)){
      ret=try(as.data.table(enrichGeneList.oddsRatio(hitSets[[grp.x]],databases = enrichrDBs,fdr.cutoff=NULL)),silent = FALSE)
      if(!any(grepl("Error",ret)) && nrow(ret) > 0){
        enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
      }
    }      
	  write.table(enrichRes, file=dirout(out, cell, "_Enrich_",filter,".tsv"), sep="\t", quote=F, row.names=F)

    if(nrow(enrichRes) > 0){
      # plot
      enrClass <- enrichRes$database[1]
      for(enrClass in unique(enrichRes$database)){
        enrichResX <- enrichRes[database == enrClass]
        if(nrow(enrichResX) > 0){
          enrichr.plot(enrichResX)
          ggsave(dirout(out, cell, "_Enrich_", enrClass,".pdf"), width=min(29, 6+ length(unique(enrichResX$grp))*0.3), height=min(29, length(unique(enrichResX$category))*0.3 + 4))
        }
      }
    } else {
      message("No enrichments for ", cell)
    }
  }
}




atac.signaling.genes <- fread(dirout("cll-time_course.ligand-receptor_repertoire.CLL.gene_level.sig_only.timepoint_mean.clustermap.csv"))$V1
ramilowski.rec_lig <- fread("~/resources_nfortelny/PairsLigRec_Ramilowski_NatComm_2015.txt")

outSig <- paste0(out, "Signaling/")
dir.create(dirout(outSig))
# ANALYSIS FOR EACH GENE SEPARATELY --------------------------------------

signaling.genes <- unique(c(atac.signaling.genes, ramilowski.rec_lig$Ligand.ApprovedSymbol, ramilowski.rec_lig$Receptor.ApprovedSymbol)) 
signaling.genes <- unique(c(signaling.genes, unique(res[grepl("NFK",gene) | grepl("^CC[LR]\\d+$", gene) | grepl("^CXC[RL]\\d+$", gene) | grepl("^CD\\d+\\w?$", gene)]$gene)))
res.reclig <- res[gene %in% signaling.genes]
res.reclig[, significant := any(qvalue < 0.05) | abs(logFC) > 0.5, by= "gene"]
res.reclig[, significant2 := any(qvalue < 0.05) | abs(logFC) > 1, by= "gene"]
res.reclig[,patient := gsub("d30", "d030", patient)]
res.reclig[,patient := gsub("d0$", "d000", patient)]
write.table(res.reclig, file=dirout(outSig, "RecLig.tsv"), sep="\t", quote=F,row.names=F)

# order groups by similarity (of OR)
if(length(unique(res.reclig$gene)) >= 2){
  try({
    orMT <- t(as.matrix(dcast.data.table(res.reclig, cellPat ~ gene, value.var="logFC")[,-"cellPat",with=F]))
    orMT[is.na(orMT)] <- 1
    hclustObj <- hclust(dist(orMT))
    res.reclig$gene <- factor(res.reclig$gene, levels=hclustObj$labels[hclustObj$order])
  }, silent=T)
}
ggplot(res.reclig[significant == TRUE], aes(x=patient, y=gene, color=logFC, size=logqval)) + geom_point() + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_gradient2(high="red", mid="white", low="blue") +
  facet_grid(. ~ cellType, scale="free", space="free")
ggsave(dirout(outSig, "Receptor_Ligand_all.pdf"), height=29, width=15)

ggplot(res.reclig[significant2 == TRUE], aes(x=patient, y=gene, color=logFC, size=logqval)) + geom_point() + 
  theme_bw(24) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_gradient2(high="red", mid="white", low="blue") +
  facet_grid(. ~ cellType, scale="free", space="free")
ggsave(dirout(outSig, "Receptor_Ligand_Strict_all.pdf"), height=29, width=15)

for(ct in unique(res.reclig$cellType)){
  ggplot(res.reclig[significant == TRUE & cellType == ct], aes(x=cellPat, y=gene, color=logFC, size=logqval)) + geom_point() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_gradient2(high="red", mid="white", low="blue")
  ggsave(dirout(outSig, "Receptor_Ligand_",ct,".pdf"), height=29, width=6)
  
  ggplot(res.reclig[significant2 == TRUE & cellType == ct], aes(x=cellPat, y=gene, color=logFC, size=logqval)) + geom_point() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_gradient2(high="red", mid="white", low="blue")
  ggsave(dirout(outSig, "Receptor_Ligand_Strict_",ct,".pdf"), height=25, width=6)
}



# Plot expression values of these genes -------------------------------------------

genes <- levels(res.reclig[significant2 == TRUE][order(logFC)]$gene)
(load(dirout("10_Seurat_raw/inclDay30_noIGHLK_negbinom/", "inclDay30_noIGHLK.RData")))
metaDat <- data.table(pbmc@data.info, keep.rownames=T)[nUMI > 1000 & nUMI < 3000]
metaDat[, sample := gsub("d0", "d000", sample)]
metaDat[, sample := gsub("d30", "d030", sample)]
dat2 <- pbmc@data[,metaDat$rn]
dat2 <- cbind(metaDat[,c("sample", "CellType"), with=T], data.table(as.matrix(t(dat2[unique(genes),]))))
dat3 <- melt(dat2, id.vars=c("sample", "CellType"))
dat3 <- dat3[,mean(value), by=c("sample", "CellType", "variable")]
dat3[, V1_norm := scale(V1), by="variable"]
dat3$variable <- factor(dat3$variable, levels=genes)

ggplot(dat3, aes(x=sample, y=variable, fill=V1)) + geom_tile() + facet_grid(. ~ CellType,scales="free", space="free") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + theme(axis.text.y=element_text(size=9)) +
  scale_fill_gradient2(low="blue", mid="white", high="red")
ggsave(dirout(outSig, "Rec_Expr.pdf"), height=29, width=15)

ggplot(dat3, aes(x=sample, y=variable, fill=V1_norm)) + geom_tile() + facet_grid(. ~ CellType,scales="free", space="free") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + theme(axis.text.y=element_text(size=9)) +
  scale_fill_gradient2(low="blue", mid="white", high="red")
ggsave(dirout(outSig, "Rec_Expr_norm.pdf"), height=29, width=15)


# Look at interactions ----------------------------------------------------
dat2 <- pbmc@data[signaling.genes[signaling.genes %in% row.names(pbmc@data)],metaDat$rn]
dat2 <- cbind(metaDat[,c("sample", "CellType"), with=T], data.table(as.matrix(t(dat2))))
dat3 <- melt(dat2, id.vars=c("sample", "CellType"))
dat3 <- dat3[,mean(value), by=c("sample", "CellType", "variable")]
dat3 <- dat3[!is.na(CellType)]
dat3 <- dat3[,V1 := V1-min(V1,na.rm=TRUE), by="variable"]
dat3 <- dat3[,V1 := V1/max(V1,na.rm=TRUE), by="variable"]

recHits <- res.reclig[significant == TRUE]
i <- 1
cellScores <- data.table()
for(i in 1:nrow(recHits)){  
  g <- recHits[i]$gene
  int <- ramilowski.rec_lig[Receptor.ApprovedSymbol == g | Ligand.ApprovedSymbol == g]
  int <- unique(union(int$Receptor.ApprovedSymbol, int$Ligand.ApprovedSymbol))
  int <- int[!int == g]
  int <- int[int %in% rownames(pbmc@data)]
  if(length(int) > 0){
    cellScores <- rbind(cellScores,data.table(
      dat3[variable %in% int & sample == recHits[i]$patient][,mean(V1,na.rm=TRUE),by="CellType"]
      , row=i))
  } else {
    cellScores <- rbind(cellScores,data.table(CellType = "Mono", V1=NA, row=i))
  }
}
cellScores <- dcast.data.table(cellScores, row ~ CellType, value.var="V1")
id.vars <- c("gene", "logFC", "patient", "cellType")
recHits <- cbind(recHits[,id.vars, with=F], cellScores[,-"row",with=F])
recHits <- melt(recHits, id.vars=id.vars)
recHits <- recHits[!is.na(value)]
recHits[, direction := ifelse(logFC > 0, "up", "down")]

recSum <- recHits[,mean(value), by=c("patient", "cellType", "variable", "direction")]
recSum[,pat := gsub("([A-Z]+)_d(\\d+)", "\\1", patient)]
recSum[,timepoint := gsub("([A-Z]+)_d(\\d+)", "\\2", patient)]
recSum[timepoint == "150", timepoint := "120"]

cellTypes <- c("CD4", "CD8", "Mono", "CLL")
recSum <- recSum[variable %in% cellTypes & cellType %in% cellTypes]

ggplot(recSum[direction== "up"], aes(x=cellType, y=variable, fill=V1))+ 
  geom_tile() + facet_grid(pat ~ timepoint) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(outSig, "SignalingNetwork_up.pdf"), height=7, width=7)

ggplot(recSum[direction== "down"], aes(x=cellType, y=variable, fill=V1))+ 
  geom_tile() + facet_grid(pat ~ timepoint) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(outSig, "SignalingNetwork_down.pdf"), height=7, width=7)



# EXPORT MATRIX OF MEANS --------------------------------------------------

metaDat <- data.table(pbmc@data.info, keep.rownames=T)[nUMI > 1000 & nUMI < 3000]
dat2 <- pbmc@data[,metaDat$rn]
dat2 <- cbind(metaDat[,c("sample", "CellType"), with=T], data.table(as.matrix(t(dat2))))
dat3 <- melt(dat2, id.vars=c("sample", "CellType"))
dat3 <- dat3[,mean(value), by=c("sample", "CellType", "variable")]
dat3[, sampleCell := paste0(sample, "_", CellType)]
dat4 <- dcast.data.table(dat3, variable ~ sampleCell, value.var="V1")
dMT <- as.matrix(dat4[,-"variable",with=F])
row.names(dMT) <- dat4$variable
write.csv(dMT, dirout("MeanMatrix.csv"))

# table(ramilowski.rec_lig$Pair.Evidence)
# table(is.na(ramilowski.rec_lig$Pair.Evidence))
# 
# rlDT <- data.table()
# for(pat in unique(res$patient)){
#   resPat <- res[patient == pat][qvalue < 0.05]
#   for(ct1 in unique(resPat$cellType)){
#     for(ct2 in unique(resPat$cellType)){
#       # rec in ct1
#       resCT1 <- resPat[cellType == ct1]
#       recDat <- resCT1[match(ramilowski.rec_lig[Pair.Evidence == "literature supported"]$Receptor.ApprovedSymbol, resCT1$gene)]
#       colnames(recDat) <- paste0("rec_",colnames(recDat))
#       # lig in ct2
#       resCT2 <- resPat[cellType == ct2]
#       ligDat <- resCT2[match(ramilowski.rec_lig[Pair.Evidence == "literature supported"]$Ligand.ApprovedSymbol, resCT2$gene)]
#       colnames(ligDat) <- paste0("lig_",colnames(ligDat))
#       
#       rlDat <- cbind(recDat, ligDat)[!is.na(rec_gene) & !is.na(lig_gene)]
#       rlDT <- rbind(rlDT, rlDat)
#     }
#   }
# }
# 
# (rlDT2 <- rlDT[, c("rec_patient", "rec_cellType", "rec_gene", "rec_logFC", "lig_cellType", "lig_gene", "lig_logFC"), with=F])
# rlDT2[sign(rec_logFC) != sign(lig_logFC)]
# rlDT3 <- rlDT2[sign(rec_logFC) == sign(lig_logFC)]
# ramilowski.rec_lig[Receptor.ApprovedSymbol == "KLRC1" & Ligand.ApprovedSymbol == "HLA-E"]
# rlDT3[order(lig_gene)]
# rlDT3[!grepl("HLA", lig_gene)]
# rlDT3[grepl("HLA", lig_gene)]
# 
# write.table(rlDT3[order(lig_gene)], dirout(out, "Ligand_Interactions.tsv"), sep="\t", quote=F, row.names=F)

message("Completed successfully")