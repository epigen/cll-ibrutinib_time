require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(pheatmap)
require(gplots)

project.init2("cll-time_course")
out <- "42_1_CNVs_Sigs_nUMI_Cutoff_min50/"
dir.create(dirout(out))


(load(dirout("10_Seurat_raw/inclDay30_noIGHLK_negbinom/inclDay30_noIGHLK.RData")))
pbmc <- UpdateSeuratObject(pbmc)
pbmc <- SubsetData(pbmc, cells.use=rownames(subset(pbmc@meta.data, nUMI > 1000 & nUMI < 3000)))

pbmc@meta.data$sample <- gsub("d30", "d030", gsub("d0", "d000", pbmc@meta.data$sample))
pbmc@meta.data$patient <- gsub("([A-Z]+)\\d?_.+", "\\1", pbmc@meta.data$sample)

pDat <- data.table(pbmc@meta.data, keep.rownames=TRUE)
pDat$sample_cell <- paste0(pDat$sample, "_", pDat$CellType)
pDat$sample <- gsub("([0-9]+)d", "d\\1", pDat$sample)
pDat$timepoint <- gsub(".+_(d[0-9]+)", "\\1", pDat$sample)


# A problem is the difference in number of UMIs between samples (in particular the later time points have less)
ggplot(pDat, aes(x=sample, y=nUMI)) + geom_boxplot(outlier.shape=NA) + coord_flip()
ggsave(dirout(out, "UMI_boxplot.pdf"))
ggplot(
  pDat,
  aes_string(x="cellType", y="nUMI", group="sample_cell", fill="patient", alpha="timepoint")) + 
  geom_boxplot(outlier.shape=NA) + coord_flip() + 
  scale_color_manual(values=c("grey", "black")) + scale_alpha_manual(values=c(0.30,0.85,1))
ggsave(dirout(out, "UMI_boxplot_byCell.pdf"))

# Read genesets
(gene.annot <- fread(paste0(Sys.getenv("RESOURCES"), "/genomes/hg38/10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf")))
genesets <- list()
table(gene.annot$V3)
a2 <- gene.annot[V3 == "gene"]
a2[,gene := gsub('.+\\; gene_name \\"(.+?)\\"\\; .+', "\\1", V9)]
length(unique(a2$gene))
a2
geneOrder <- a2[,c("V1", "V4", "gene"),with=F]
colnames(geneOrder) <- c("chr", "start", "gene")
geneOrder <- geneOrder[chr %in% c(as.character(1:22), "X", "Y")]
geneOrder

geneOrder2 <- geneOrder[gene %in% rownames(pbmc@data)]
slide.stepsize <- 50
chr.x <- "Y"
i2 <- 0
for(chr.x in unique(geneOrder2$chr)){
  geneOrder.chr <- geneOrder2[chr == chr.x]
  nr.sets <- floor((nrow(geneOrder.chr))/slide.stepsize)
  if(nr.sets == 0) nr.sets <- 1
  slide.stepsize2 <- ceiling(nrow(geneOrder.chr)/nr.sets)
  for(i in 1:nr.sets){
    i2 <- i2 + 1
	message("Chromosome ", chr.x, "----------")
	print(quantile(((i-1)*slide.stepsize2+1):min(i*slide.stepsize2, nrow(geneOrder.chr))))
    genesets[[paste0(chr.x, "_", i2)]] <- geneOrder.chr[((i-1)*slide.stepsize2+1):min(i*slide.stepsize2, nrow(geneOrder.chr))]$gene
  }
}
save(genesets, file=dirout(out, "genesets.RData"))

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


if(!file.exists(dirout(out,"Scores.RData"))){  
  # CALCULATE SCORE ---------------------------------------------------------
  data <- pbmc@data
  data@x <- exp(data@x) - 1
  data <- data - apply(data, 1, min)
  data <- data / apply(data, 1, max)
  
  Dat1 <- pDat
  
  # Background distribution
  n.draws <- 500
  bg.distr <- matrix(NA, ncol=ncol(data), nrow=n.draws)
  for(i in 1:n.draws){
    bg.distr[i,] <- apply(data[sample(1:nrow(data), slide.stepsize),], 2, mean)
  }
  bg.mean <- apply(bg.distr, 2, mean)
  bg.sd <- apply(bg.distr, 2, sd)
  
  
  (set.nam <- names(genesets)[1])
  for(set.nam in names(genesets)){
    genes <- genesets[[set.nam]]
    genes <- genes[genes %in% rownames(data)]
    
    score.x <- apply(data[genes,,drop=F], 2, mean)    
    Dat1[[set.nam]] <- (score.x - bg.mean)/bg.sd
  }
  qplot(sapply(names(genesets), function(nam) cor(Dat1[[nam]], Dat1$nUMI)), bins=10) + xlab("Correlation")
  ggsave(dirout(out, "CorrUMI_Raw.pdf"))
  
  names(Dat1)
  
  # DO NOT USE REGRESSED SCORES, then can have < 0
  save(Dat1, file=dirout(out,"Scores.RData"))
  
} else {
  load(file=dirout(out,"Scores.RData"))
}



# ANALYSIS T_ZERO ----------------------------------------------------------------
pDat2 <- Dat1
pDat2 <- pDat2[timepoint %in% c("d0")]
pDat2 <- pDat2[cellType %in% c("CD8", "CD4", "Mono", "CLL")]# & grepl("PT", sample)]

# Make one large heatmap
tzero.sig <- data.table()
pat="FE"
cell="CD4"
geneset <- "Allograft.rejection"
for(cell in unique(pDat2$cellType)){
  x <- pDat2[cellType == cell]
  x <- x[patient %in% x[,.N, by="patient"][N>5]$patient]
  for(pat in unique(x$patient)){
    for(geneset in names(genesets)){
      p <- t.test(x[patient == pat][[geneset]], x[patient != pat][[geneset]])$p.value
      ef <- mean(x[patient == pat][[geneset]])- mean(x[patient != pat][[geneset]])
      tzero.sig <- rbind(tzero.sig, data.table(patient=pat, cellType=cell, pvalue=p, meanDiff=ef, geneset=geneset))
    }
  }
}
tzero.sig[,pvalue2 := pmin(5, -1*log10(pvalue))]
save(tzero.sig, file=dirout(out, "Tzero.sig.RData"))
ggplot(
  tzero.sig[geneset %in% tzero.sig[,min(pvalue), by="geneset"][order(V1)][1:50]$geneset], 
  aes(x=paste0(cellType, "_", patient), y=geneset, color=meanDiff, size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "0_Overview_TZero.pdf"),height=15, width=15)

tzero.sig$chr <- sapply(strsplit(tzero.sig$geneset, "_"), function(split) split[1])
tzero.sig$pos <- sapply(strsplit(tzero.sig$geneset, "_"), function(split) as.numeric(split[2]))
tzero.sig <- tzero.sig[order(chr, pos)]
tzero.sig$geneset <- factor(tzero.sig$geneset, levels=unique(tzero.sig$geneset))

ggplot(tzero.sig[cellType=="CLL"], aes(x=geneset, y=meanDiff, color=patient, alpha = pvalue2)) + 
  geom_point() + theme_bw(9) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "0_TZero_CLL.pdf"), width=29)


ggplot(tzero.sig[cellType=="CLL"], aes(x=geneset, y=meanDiff, color=patient, group=patient)) + 
  geom_line() + theme_bw(9) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "0_TZero_CLL_lines.pdf"), width=29)


# 
# # PLOT Gene Sets t0
# n <- 100
# set.nam <- "P53.pathway"
# for(set.nam in names(genesets)){
#   genes <- genesets[[set.nam]]
#   dat <- pbmc@data
#   sampleAnnot <- subset(pbmc@data.info, cellType %in% c("CD8", "CD4", "Mono", "CLL") & (grepl("d0", sample) | grepl("_0d", sample)))
#   cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(with(sampleAnnot, paste0(sample, cellType)))), function(x) sample(x, min(length(x), n))))]
#   
#   genes <- genes[genes %in% rownames(dat)]
#   dat <- dat[genes, cells]
#   dat <- dat[apply(dat, 1, max) != 0,, drop=F]
#   if(nrow(dat) > 0){
#     dat <- dat - apply(dat, 1, min)
#     dat <- dat / apply(dat,1, max)
#     apply(dat, 1, quantile, na.rm=T)
#     dat <- dat[, order(with(sampleAnnot[cells,],paste0(cellType, sample))),drop=F]
#     
#     if(nrow(dat) > 1 & ncol(dat) >1){
#       pdf(dirout(out.details, set.nam, "_TZero", ".pdf"), height=min(29, nrow(dat) * 0.3 + 3), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
#       pheatmap(dat, cluster_rows=F, cluster_cols=F, annotation_col=pbmc@data.info[,c("sample", "cellType", "nUMI"), drop=F],show_colnames=FALSE)
#       dev.off()
#     }
#   }
# }



# # ANALYSIS OVER TIME ----------------------------------------------------------------
pDat <- Dat1
pDat2 <- pDat
pDat2 <- pDat2[timepoint %in% c("d0", "d120")]
pDat2 <- pDat2[cellType %in% c("CD8", "CD4", "Mono", "CLL")]# & grepl("PT", sample)]

# Make one large heatmap over time
overTime.sig <- data.table()
pat="VZS"
cell="CD4"
geneset <- "Mono_up"
for(pat in unique(pDat2$patient)){
  for(cell in unique(pDat2$cellType)){
    x <- pDat2[patient == pat & cellType == cell]
    if(nrow(x[timepoint == "d0"]) > 5 & nrow(x[timepoint == "d120"]) >5){
      for(geneset in names(genesets)){
        p <- t.test(x[timepoint == "d0"][[geneset]], x[timepoint == "d120"][[geneset]])$p.value
        ef <- mean(x[timepoint == "d120"][[geneset]]) - mean(x[timepoint == "d0"][[geneset]])
        overTime.sig <- rbind(overTime.sig, data.table(patient=pat, cellType=cell, pvalue=p, Diff=ef, geneset=geneset))
      }
    }
  }
}

overTime.sig[,pvalue2 := pmin(5, -1*log10(pvalue))]
overTime.sig[,qvalue := p.adjust(pvalue, method="BH")]
overTime.sig[,qvalue2 := pmin(5, -1*log10(qvalue))]

save(overTime.sig, file=dirout(out, "overTime.sig.RData"))
# overTime.sig[is.nan(logFC), pvalue := 1]
# overTime.sig[is.nan(logFC), logFC := 0]
# # overTime.sig[logFC == Inf, logFC := max(overTime.sig$logFC[overTime.sig$logFC != Inf])]
# # overTime.sig[logFC == -Inf, logFC := min(overTime.sig$logFC[overTime.sig$logFC != -Inf])]
# overTime.sig[abs(logFC) > 1, logFC := 1 * sign(logFC)]
otPlot <- overTime.sig[geneset %in% overTime.sig[,min(pvalue), by="geneset"][order(V1, decreasing=FALSE)][1:50]$geneset]
max22 <- min(abs(min(otPlot$Diff)), max(otPlot$Diff))
otPlot$Diff2 <- otPlot$Diff
otPlot[Diff < -max22, Diff := -max22]
otPlot[Diff > max22, Diff := max22]
ggplot(
  otPlot,
  aes(x=paste0(cellType, "_", patient), y=geneset, color=Diff, size=qvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave(dirout(out, "0_Overview.pdf"),height=15, width=15)


overTime.sig$chr <- sapply(strsplit(overTime.sig$geneset, "_"), function(split) split[1])
overTime.sig$pos <- sapply(strsplit(overTime.sig$geneset, "_"), function(split) as.numeric(split[2]))
overTime.sig <- overTime.sig[order(chr, pos)]
overTime.sig$geneset <- factor(overTime.sig$geneset, levels=unique(overTime.sig$geneset))

ggplot(overTime.sig[cellType=="CLL"], aes(x=geneset, y=Diff, color=patient, alpha = pvalue2)) + 
  geom_point() + theme_bw(9) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "0_OverTime_CLL.pdf"), width=29)


ggplot(overTime.sig[cellType=="CLL"], aes(x=geneset, y=Diff, color=patient, group=patient)) + 
  geom_line() + theme_bw(9) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "0_OverTime_CLL_lines.pdf"), width=29)



# 
# 
# 
# out.details <- paste0(out, "Details/")
# 
# geneset <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"
# for(geneset in names(genesets)){
#   
#   # Plot all values as boxplots
#   ggplot(
#     pDat2,
#     aes_string(x="cellType", y=geneset, group="sample_cell", fill="patient", alpha="timepoint")) + 
#     geom_boxplot(outlier.shape=NA) + coord_flip() + 
#     scale_color_manual(values=c("grey", "black")) + scale_alpha_manual(values=c(0.35,1))
#   ggsave(dirout(out.details, geneset, "_AllData.pdf"), width=7, height=7)
#   
#   # Only t0 as boxplot
#   ggplot(
#     pDat2[timepoint == "d0"],
#     aes_string(x="cellType", y=geneset, group="sample_cell", fill="patient")) + 
#     geom_boxplot(outlier.shape=NA) + coord_flip()
#   ggsave(dirout(out.details, geneset, "_t0.pdf"), width=7, height=7)
#   
#   # just the change over time as heatmap
#   sDat <- pDat2[,median(get(geneset)), by=c("patient", "timepoint", "cellType")]
#   sDat <- dcast.data.table(sDat, patient + cellType ~ timepoint, value.var="V1")
#   sDat[, timechange := d120 - d0]
#   ggplot(sDat[!is.na(timechange)], aes(x=cellType, y=patient, fill=timechange)) + 
#     geom_tile() + scale_fill_gradient2(low="blue", high="red", mid="black")
#   ggsave(dirout(out.details, geneset, "_HM.pdf"), width=5, height=4)
# }
# 
# venn(list(alpha=genesets$HALLMARK_INTERFERON_ALPHA_RESPONSE, gamma=genesets$HALLMARK_INTERFERON_GAMMA_RESPONSE))
# 
# 
# # PLOT INDIVIDUAL EXAMPLES ------------------------------------------------
# 
# pat <- "FE"
# set.nam <- "CLL_up"
# n <- 100
# for(pat in unique(pDat$patient)){
#   for(set.nam in names(genesets)){
#     genes <- genesets[[set.nam]]
#     dat <- pbmc@data
#     sampleAnnot <- subset(pbmc@data.info, grepl(pat, sample) & cellType %in% c("CD8", "CD4", "Mono", "CLL"))
#     cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(with(sampleAnnot, paste0(sample, cellType)))), function(x) sample(x, min(length(x), n))))]
#     
#     genes <- genes[genes %in% rownames(dat)]
#     dat <- dat[genes, cells]
#     dat <- dat[apply(dat, 1, max) != 0,,drop=F]
#     dat <- dat - apply(dat, 1, min)
#     dat <- dat / apply(dat,1, max)
#     #     apply(dat, 1, quantile, na.rm=T)
#     dat <- dat[, order(with(sampleAnnot[cells,],paste0(cellType, sample))), drop=F]
#     
#     pdf(dirout(out.details, set.nam, "_", pat, ".pdf"), height=min(29, nrow(dat) * 0.3 + 3), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
#     pheatmap(dat, cluster_rows=F, cluster_cols=F, annotation_col=pbmc@data.info[,c("sample", "cellType", "nUMI"), drop=F],show_colnames=FALSE)
#     dev.off()
#   }
# }