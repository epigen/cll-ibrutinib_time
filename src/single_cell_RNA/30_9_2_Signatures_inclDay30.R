require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(pheatmap)
require(gplots)

project.init2("cll-time_course")
out <- "30_9_Signatures_inclDay30_downLists2/"
dir.create(dirout(out))


sample.x <- "inclDay30_noIGHLK"
load(file=dirout("10_Seurat_raw/", sample.x, "_negbinom/",sample.x,".RData"))

pbmc@data.info$sample <- gsub("d30", "d030", gsub("d0", "d000", pbmc@data.info$sample))
pbmc@data.info$patient <- gsub("([A-Z]+)\\d?_.+", "\\1", pbmc@data.info$sample)

pDat <- data.table(pbmc@data.info, keep.rownames=TRUE)
# pDat[,sample := gsub("d0", "d000", sample)]
# pDat[,sample := gsub("d30", "d030", sample)]
pDat$sample_cell <- paste0(pDat$sample, "_", pDat$CellType)
pDat$sample <- gsub("([0-9]+)d", "d\\1", pDat$sample)
pDat$timepoint <- gsub(".+_(d[0-9]+)", "\\1", pDat$sample)


# counts for cell types
for(ct in unique(pDat$CellType)){
  ggplot(pDat[CellType == ct][,.N, by=c("sample", "CellType")], aes(x=sample, y=N)) + geom_bar(stat='identity') + coord_flip()
  ggsave(dirout(out, "Counts_", ct,".pdf"), width=7, height=7)  
  
  ggplot(pDat[CellType == ct], aes(x=sample, y=nUMI)) + geom_violin() + geom_jitter(width=0.05) +coord_flip()
  ggsave(dirout(out, "CountsJitter_", ct,".pdf"), width=7, height=7)  
}

# A problem is the difference in number of UMIs between samples (in particular the later time points have less)
ggplot(pDat, aes(x=sample, y=nUMI)) + geom_boxplot(outlier.shape=NA) + coord_flip()
ggsave(dirout(out, "UMI_boxplot.pdf"))
ggplot(
  pDat,
  aes_string(x="CellType", y="nUMI", group="sample_cell", fill="patient", alpha="timepoint")) + 
  geom_boxplot(outlier.shape=NA) + coord_flip() + 
  scale_color_manual(values=c("grey", "black")) + scale_alpha_manual(values=c(0.30,0.5,0.85,0.85,1))
ggsave(dirout(out, "UMI_boxplot_byCell.pdf"))

# Read genesets
file <- "/data/groups/lab_bock/shared/resources/gene_sets/MSigDB/6.0/Human/h.all.v6.0.symbols.gmt"
lines <- readLines(file)
genesets <- list()
for(line in lines){
  x <- strsplit(line, "\t")[[1]]
  genesets[[x[1]]] <- x[3:length(x)]
}
load(dirout("20_AggregatedLists/lists_plusBcells.RData"))
genesets <- c(genesets, cll.lists2)

cllDownGenes <- fread(dirout("13_4_Overtime_inclDay30/SigGenes_overTime.tsv"))
cllDownGenes <- cllDownGenes[cellType == "CLL"][qvalue < 0.05 & logFC < 0]
genesets <- c(genesets, split(cllDownGenes$gene, factor(cllDownGenes$patient)))

save(genesets, file=dirout(out, "Genesets.RData"))

file.remove(dirout(out, "Genesets.tsv"))
lapply(names(genesets), function(lnam){
  write(paste(c(lnam, genesets[[lnam]]), collapse="\t"), file=dirout(out, "Genesets.tsv"), append=TRUE, ncolumns = max(sapply(genesets, length)))
  })


# Jaccard of genesets
genesets.jac <- lapply(genesets, function(x) x[x %in% rownames(pbmc@data)])
genesets.jac <- genesets.jac[sapply(genesets.jac, length) > 0]
jac <- matrix(NA, ncol=length(genesets.jac), nrow=length(genesets.jac))
colnames(jac)  <- names(genesets.jac)
rownames(jac)  <- names(genesets.jac)
for(i1 in 1:(length(genesets.jac))){
  for(i2 in i1:length(genesets.jac)){
     if(i1 != i2){
       jac.score <- length(intersect(genesets.jac[[i1]],genesets.jac[[i2]]))/length(union(genesets.jac[[i1]],genesets.jac[[i2]]))
       jac[i1,i2] <- jac.score
       jac[i2,i1] <- jac.score
     }
  }
}
pdf(dirout(out, "1_jaccard.pdf"), height=15, width=15,onefile=F)
pheatmap(jac)
dev.off()

ggplot(data.table(length=sapply(genesets.jac, length), geneset = names(genesets.jac)), aes(y=length, x=geneset)) + geom_bar(stat="identity") + coord_flip()
ggsave(dirout(out, "1_SetLengths.pdf"), height=12, width=10)

if(!file.exists(dirout(out,"Scores.RData"))){  
  # CALCULATE SCORE ---------------------------------------------------------
  data <- pbmc@data
  data@x <- exp(data@x) - 1
  data <- data - apply(data, 1, min)
  data <- data / apply(data, 1, max)
  
  Dat1 <- pDat
  
  n.draws <- 500
  
  (set.nam <- names(genesets)[1])
  for(set.nam in names(genesets)){
    genes <- genesets[[set.nam]]
    genes <- genes[genes %in% rownames(data)]
    
    score.x <- apply(data[genes,], 2, mean)
    
    other.genes <- rownames(data)[!rownames(data) %in% genes]
    # Background distribution
    bg.distr <- matrix(NA, ncol=ncol(data), nrow=n.draws)
    for(i in 1:n.draws){
      bg.distr[i,] <- apply(data[sample(other.genes, length(genes),replace=TRUE),], 2, mean)
    }
    bg.mean <- apply(bg.distr, 2, mean)
    bg.sd <- apply(bg.distr, 2, sd)
    
    Dat1[[set.nam]] <- (score.x - bg.mean)/bg.sd
  }
  qplot(sapply(names(genesets), function(nam) cor(Dat1[[nam]], Dat1$nUMI)), bins=10) + xlab("Correlation")
  ggsave(dirout(out, "CorrUMI_Raw.pdf"))
  
  names(Dat1)
  ggplot(Dat1[patient=="PT" & CellType == "Mono"], aes(x=nUMI, y=HALLMARK_TNFA_SIGNALING_VIA_NFKB)) + geom_point()
  ggsave(dirout(out, "CorrUMI_AggLogCounts_PT_Mono.pdf"))
  ggplot(Dat1, aes(x=nUMI, y=HALLMARK_TNFA_SIGNALING_VIA_NFKB)) + geom_hex() + 
    ggtitle(paste("cor = ", round(cor(Dat1$nUMI, Dat1$HALLMARK_TNFA_SIGNALING_VIA_NFKB),3)))
  ggsave(dirout(out, "CorrUMI_AggLogCounts_hex.pdf"))
  
  
  # REGRESS OUT nUMI --------------------------------------------------------
  # Dat1.cor <- Dat1
  # set.nam <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
  # for(set.nam in names(genesets)){
  #   lm.fit <- lm(data=Dat1.cor, get(set.nam) ~ nUMI)
  #   Dat1.cor[[set.nam]] <- lm.fit$residuals + coef(lm.fit)[1] # add intercept back
  # }
  # qplot(sapply(names(genesets), function(nam) cor(Dat1.cor[[nam]], Dat1.cor$nUMI)), bins=10) + xlab("Correlation")
  # ggsave(dirout(out, "CorrUMI_Regressed.pdf"))
  # 
  # plot(Dat1$HALLMARK_TNFA_SIGNALING_VIA_NFKB, Dat1$nUMI)
  # 
  # plot(Dat1.cor$HALLMARK_TNFA_SIGNALING_VIA_NFKB, Dat1.cor$nUMI)
  
  # DO NOT USE REGRESSED SCORES, then can have < 0
  save(Dat1, file=dirout(out,"Scores.RData"))
  
} else {
  load(file=dirout(out,"Scores.RData"))
}


# ANALYSIS OVER TIME ----------------------------------------------------------------
pDat <- Dat1
pDat2 <- pDat
pDat2 <- pDat2[CellType %in% c("CD8", "CD4", "Mono", "CLL")]# & grepl("PT", sample)]
pDat2[,sample_cell := paste(patient, timepoint, CellType, sep="_")]
out.details <- paste0(out, "Details/")
dir.create(dirout(out.details))


geneset <- "HALLMARK_MYC_TARGETS_V1"
for(geneset in names(genesets)){
  
  # Plot all values as boxplots
  p <- ggplot(
    pDat2,
    aes_string(x="CellType", y=geneset, group="sample_cell", fill="patient", alpha="timepoint")) + 
    scale_color_manual(values=c("grey", "black")) + scale_alpha_manual(values=c(0.35, 0.5, 0.85, 0.85,1))
  ggsave(dirout(out.details, geneset, "_AllData.pdf"), width=7, height=7, plot=p + geom_boxplot(outlier.shape=NA) + coord_flip())
  ggsave(dirout(out.details, geneset, "_AllData2.pdf"), width=7, height=7, plot=p + geom_violin() + coord_flip())
  
  # Only t0 as boxplot
  ggplot(
    pDat2[timepoint == "d000"],
    aes_string(x="CellType", y=geneset, group="sample_cell", fill="patient")) + 
    geom_boxplot(outlier.shape=NA) + coord_flip()
  ggsave(dirout(out.details, geneset, "_t0.pdf"), width=7, height=7)
  
  # just the change over time as heatmap
  #   sDat <- pDat2[,median(get(geneset)), by=c("patient", "timepoint", "CellType")]
  #   sDat <- dcast.data.table(sDat, patient + CellType ~ timepoint, value.var="V1")
  #   sDat[, timechange := d120 - d0]
  #   ggplot(sDat[!is.na(timechange)], aes(x=CellType, y=patient, fill=timechange)) + 
  #     geom_tile() + scale_fill_gradient2(low="blue", high="red", mid="black")
  #   ggsave(dirout(out.details, geneset, "_HM.pdf"), width=5, height=4)
}

venn(list(alpha=genesets$HALLMARK_INTERFERON_ALPHA_RESPONSE, gamma=genesets$HALLMARK_INTERFERON_GAMMA_RESPONSE))


# Make one large heatmap over time
overTime.sig <- data.table()
pat="VZS"
cell="CD4"
geneset <- "Mono_up"
for(pat in unique(pDat2$patient)){
  for(cell in unique(pDat2$CellType)){
    x <- pDat2[patient == pat & CellType == cell]
    tp <- unique(x[timepoint != "d000"]$timepoint)[1]
    for(tp in unique(x[timepoint != "d000"]$timepoint)){
      if(nrow(x[timepoint == "d000"]) > 5 & nrow(x[timepoint == tp]) >5){
        for(geneset in names(genesets)){
          p <- t.test(x[timepoint == "d000"][[geneset]], x[timepoint == tp][[geneset]])$p.value
          ef <- mean(x[timepoint == tp][[geneset]]) - mean(x[timepoint == "d000"][[geneset]])
          overTime.sig <- rbind(overTime.sig, data.table(patient=pat, CellType=cell, pvalue=p, Diff=ef,timepoint=tp, geneset=geneset))
        }
      }
    }
  }
}
# overTime.sig[is.nan(logFC), pvalue := 1]
# overTime.sig[is.nan(logFC), logFC := 0]
# # overTime.sig[logFC == Inf, logFC := max(overTime.sig$logFC[overTime.sig$logFC != Inf])]
# # overTime.sig[logFC == -Inf, logFC := min(overTime.sig$logFC[overTime.sig$logFC != -Inf])]
# overTime.sig[abs(logFC) > 1, logFC := 1 * sign(logFC)]
overTime.sig[,pvalue2 := pmin(5, -1*log10(pvalue))]
save(overTime.sig, file=dirout(out, "OverTime_Tests.RData"))
ggplot(
  overTime.sig, 
  aes(x=paste0(CellType, "_", patient, "_", timepoint), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave(dirout(out, "0_Overview.pdf"),height=15, width=15)

pDat2[,.N, by=c("sample", "CellType")]

# 
# # PLOT INDIVIDUAL EXAMPLES ------------------------------------------------
# 
pat <- "FE"
set.nam <- "CLL_up"
n <- 100
for(pat in unique(pDat$patient)){
  for(set.nam in names(genesets)){
    genes <- genesets[[set.nam]]
    dat <- pbmc@data
    sampleAnnot <- subset(pbmc@data.info, grepl(pat, sample) & CellType %in% c("CD8", "CD4", "Mono", "CLL"))
    cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(with(sampleAnnot, paste0(sample, CellType)))), function(x) sample(x, min(length(x), n))))]
    
    genes <- genes[genes %in% rownames(dat)]
    dat <- dat[genes, cells]
    dat <- dat[apply(dat, 1, max) != 0,,drop=F]
    dat <- dat - apply(dat, 1, min)
    dat <- dat / apply(dat,1, max)
    #     apply(dat, 1, quantile, na.rm=T)
    dat <- dat[, order(with(sampleAnnot[cells,],paste0(CellType, sample))), drop=F]
    
    pdf(dirout(out.details, set.nam, "_", pat, ".pdf"), height=min(29, nrow(dat) * 0.3 + 3), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
    pheatmap(dat, cluster_rows=F, cluster_cols=F, annotation_col=pbmc@data.info[,c("sample", "CellType", "nUMI"), drop=F],show_colnames=FALSE)
    dev.off()
  }
}

# one plot for each celltype
n <- 100
for(ct in c("CD8", "CD4", "Mono", "CLL")){
  message(ct)
  for(set.nam in names(genesets)){
    genes <- genesets[[set.nam]]
    dat <- pbmc@data
    sampleAnnot <- subset(pbmc@data.info, CellType == ct)
    cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(sampleAnnot$sample)), function(x) sample(x, min(length(x), n))))]
    
    genes <- genes[genes %in% rownames(dat)]
    dat <- dat[genes, cells]
    dat <- dat[apply(dat, 1, max) != 0,,drop=F]
    dat <- dat - apply(dat, 1, min)
    dat <- dat / apply(dat,1, max)
    #     apply(dat, 1, quantile, na.rm=T)
    dat <- dat[, order(with(sampleAnnot[cells,],paste0(CellType, sample))), drop=F]
    
    pdf(dirout(out.details, set.nam, "_", ct, ".pdf"), height=min(29, nrow(dat) * 0.3 + 3), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
    pheatmap(dat, cluster_rows=F, cluster_cols=F, annotation_col=pbmc@data.info[,c("sample", "patient", "nUMI"), drop=F],show_colnames=FALSE)
    dev.off()
  }
}


# # ANALYSIS T_ZERO ----------------------------------------------------------------
pDat2 <- Dat1
pDat2 <- pDat2[timepoint %in% c("d000")]
pDat2 <- pDat2[CellType %in% c("CD8", "CD4", "Mono", "CLL")]# & grepl("PT", sample)]

# Make one large heatmap
overTime.sig <- data.table()
pat="FE"
cell="CD4"
geneset <- "Allograft.rejection"
for(cell in unique(pDat2$CellType)){
  x <- pDat2[CellType == cell]
  x <- x[patient %in% x[,.N, by="patient"][N>5]$patient]
  for(pat in unique(x$patient)){
    for(geneset in names(genesets)){
      p <- t.test(x[patient == pat][[geneset]], x[patient != pat][[geneset]])$p.value
      ef <- mean(x[patient == pat][[geneset]])- mean(x[patient != pat][[geneset]])
      overTime.sig <- rbind(overTime.sig, data.table(patient=pat, CellType=cell, pvalue=p, meanDiff=ef, geneset=geneset))
    }
  }
}
overTime.sig[,pvalue2 := pmin(5, -1*log10(pvalue))]
ggplot(
  overTime.sig, 
  aes(x=paste0(CellType, "_", patient), y=geneset, color=meanDiff, size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "0_Overview_TZero.pdf"),height=15, width=15)

# 
# 
# PLOT Gene Sets t0
n <- 100
set.nam <- "P53.pathway"
for(set.nam in names(genesets)){
  genes <- genesets[[set.nam]]
  dat <- pbmc@data
  sampleAnnot <- subset(pbmc@data.info, CellType %in% c("CD8", "CD4", "Mono", "CLL") & (grepl("d0", sample) | grepl("_0d", sample)))
  cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(with(sampleAnnot, paste0(sample, CellType)))), function(x) sample(x, min(length(x), n))))]
  
  genes <- genes[genes %in% rownames(dat)]
  dat <- dat[genes, cells]
  dat <- dat[apply(dat, 1, max) != 0,, drop=F]
  if(nrow(dat) > 0){
    dat <- dat - apply(dat, 1, min)
    dat <- dat / apply(dat,1, max)
    apply(dat, 1, quantile, na.rm=T)
    dat <- dat[, order(with(sampleAnnot[cells,],paste0(CellType, sample))),drop=F]
    
    if(nrow(dat) > 1 & ncol(dat) >1){
      pdf(dirout(out.details, set.nam, "_TZero", ".pdf"), height=min(29, nrow(dat) * 0.3 + 3), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
      pheatmap(dat, cluster_rows=F, cluster_cols=F, annotation_col=pbmc@data.info[,c("sample", "CellType", "nUMI"), drop=F],show_colnames=FALSE)
      dev.off()
    }
  }
}