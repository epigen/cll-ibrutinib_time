require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(pheatmap)
require(gplots)

project.init2("cll-time_course")
out <- "30_5_SignaturesMagic/"
dir.create(dirout(out))


sample.x <- "allDataBest_NoDownSampling_noIGH"

load(file=dirout("10_Seurat/", sample.x, "/",sample.x,".RData"))

pDat <- data.table(pbmc@data.info, keep.rownames=TRUE)
pDat$sample_cell <- paste0(pDat$sample, "_", pDat$cellType)
pDat$patient <- gsub("([A-Z]+)\\d?_.+", "\\1", pDat$sample)
pDat$sample <- gsub("([0-9]+)d", "d\\1", pDat$sample)
pDat$sample <- gsub("d150", "d120", pDat$sample)
pDat$timepoint <- gsub(".+_(d[0-9]+)", "\\1", pDat$sample)


# A problem is the difference in number of UMIs between samples (in particular the later time points have less)
ggplot(pDat, aes(x=sample, y=nUMI)) + geom_boxplot(outlier.shape=NA) + coord_flip()
ggsave(dirout(out, "UMI_boxplot.pdf"))

# Read genesets
file <- "/data/groups/lab_bock/shared/resources/gene_sets/MSigDB/6.0/Human/h.all.v6.0.symbols.gmt"
lines <- readLines(file)
genesets <- list()
for(line in lines){
  x <- strsplit(line, "\t")[[1]]
  genesets[[x[1]]] <- x[3:length(x)]
}
load(dirout("20_AggregatedLists/lists.RData"))
genesets <- c(genesets, cll.lists)


# CALCULATE SCORE ---------------------------------------------------------
magicDat <- fread(dirout("15_Magic/OLD/magic.matrix.csv"))
magicMT <- as.matrix(magicDat[,-"V1", with=F])
rownames(magicMT) <- magicDat$V1
data <- t(magicMT)

stopifnot(all(colnames(data) == rownames(pbmc@data.info$sample)))

# data <- t(scale(t(data)))
data <- data - apply(data, 1, min)
data <- data / apply(data, 1, max)
boxplot(t(data[1:10,]))

# 1 calculate aggregated log counts
Dat1 <- pDat
for(set.nam in names(genesets)){
  genes <- genesets[[set.nam]]
  genes <- genes[genes %in% rownames(data)]
  Dat1[[set.nam]] <- apply(data[genes,], 2, mean)
}
qplot(sapply(names(genesets), function(nam) cor(Dat1[[nam]], Dat1$nUMI)), bins=10) + xlab("Correlation")
ggsave(dirout(out, "CorrUMI_Raw.pdf"))

names(Dat1)
ggplot(Dat1[patient=="PT" & cellType == "Mono"], aes(x=nUMI, y=HALLMARK_TNFA_SIGNALING_VIA_NFKB)) + geom_point()
ggsave(dirout(out, "CorrUMI_AggLogCounts_PT_Mono.pdf"))
ggplot(Dat1, aes(x=nUMI, y=HALLMARK_TNFA_SIGNALING_VIA_NFKB)) + geom_hex() + 
  ggtitle(paste("cor = ", round(cor(Dat1$nUMI, Dat1$HALLMARK_TNFA_SIGNALING_VIA_NFKB),3)))
ggsave(dirout(out, "CorrUMI_AggLogCounts_hex.pdf"))



# REGRESS OUT nUMI --------------------------------------------------------
Dat1.cor <- Dat1
set.nam <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
for(set.nam in names(genesets)){
  lm.fit <- lm(data=Dat1.cor, get(set.nam) ~ nUMI)
  Dat1.cor[[set.nam]] <- lm.fit$residuals + coef(lm.fit)[1] # add intercept back
}
qplot(sapply(names(genesets), function(nam) cor(Dat1.cor[[nam]], Dat1.cor$nUMI)), bins=10) + xlab("Correlation")
ggsave(dirout(out, "CorrUMI_Regressed.pdf"))

plot(Dat1$HALLMARK_TNFA_SIGNALING_VIA_NFKB, Dat1$nUMI)
plot(Dat1.cor$HALLMARK_TNFA_SIGNALING_VIA_NFKB, Dat1.cor$nUMI)



pDat <- Dat1

# ANALYSIS ----------------------------------------------------------------
pDat2 <- pDat
pDat2 <- pDat2[timepoint %in% c("d0", "d120")]
pDat2 <- pDat2[cellType %in% c("CD8", "CD4", "Mono", "CLL")]# & grepl("PT", sample)]

geneset <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"
for(geneset in names(genesets)){
  
  # Plot all values as boxplots
  ggplot(
    pDat2,
    aes_string(x="cellType", y=geneset, group="sample_cell", fill="patient", alpha="timepoint")) + 
    geom_boxplot(outlier.shape=NA) + coord_flip() + 
    scale_color_manual(values=c("grey", "black")) + scale_alpha_manual(values=c(0.35,1))
  ggsave(dirout(out, geneset, "_AllData.pdf"), width=7, height=7)
  
  # Only t0 as boxplot
  ggplot(
    pDat2[timepoint == "d0"],
    aes_string(x="cellType", y=geneset, group="sample_cell", fill="patient")) + 
    geom_boxplot(outlier.shape=NA) + coord_flip()
  ggsave(dirout(out, geneset, "_t0.pdf"), width=7, height=7)
  
  # just the change over time as heatmap
  sDat <- pDat2[,median(get(geneset)), by=c("patient", "timepoint", "cellType")]
  sDat <- dcast.data.table(sDat, patient + cellType ~ timepoint, value.var="V1")
  sDat[, timechange := d120 - d0]
  ggplot(sDat[!is.na(timechange)], aes(x=cellType, y=patient, fill=timechange)) + 
    geom_tile() + scale_fill_gradient2(low="blue", high="red", mid="black")
  ggsave(dirout(out, geneset, "_HM.pdf"), width=5, height=4)
}

venn(list(alpha=genesets$HALLMARK_INTERFERON_ALPHA_RESPONSE, gamma=genesets$HALLMARK_INTERFERON_GAMMA_RESPONSE))


# Make one large heatmap
overTime.sig <- data.table()
pat="PT"
cell="Mono"
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
ggplot(
  overTime.sig, 
  aes(x=paste0(cellType, "_", patient), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw()
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "0_Overview.pdf"),height=15, width=15)




colnames(data) <- rownames(pbmc@data.info)
# PLOT INDIVIDUAL EXAMPLES ------------------------------------------------
pat <- "PT"
set.nam <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
n <- 100
for(pat in unique(pDat$patient)){
  for(set.nam in names(genesets)){
    genes <- genesets[[set.nam]]
    dat <- data
    sampleAnnot <- subset(pbmc@data.info, grepl(pat, sample) & cellType %in% c("CD8", "CD4", "Mono", "CLL"))
    cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(with(sampleAnnot, paste0(sample, cellType)))), function(x) sample(x, min(length(x), n))))]
    
    genes <- genes[genes %in% rownames(dat)]
    dat <- dat[genes, cells]
    dat <- dat[apply(dat, 1, max) != 0,,drop=F]
    dat <- dat - apply(dat, 1, min)
    dat <- dat / apply(dat,1, max)
    apply(dat, 1, quantile, na.rm=T)
    dat <- dat[, order(with(sampleAnnot[cells,],paste0(cellType, sample))),drop=F]
    
    pdf(dirout(out, set.nam, "_", pat, ".pdf"), height=min(29, nrow(dat) * 0.3 + 3), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
    pheatmap(dat, cluster_rows=F, cluster_cols=F, annotation_col=pbmc@data.info[,c("sample", "cellType", "nUMI"), drop=F],show_colnames=FALSE)
    dev.off()
  }
}



# RANDOM MATRIX -----------------------------------------------------------
dat <- pbmc@data
sampleAnnot <- subset(pbmc@data.info, grepl(pat, sample) & cellType %in% c("CD8", "CD4", "Mono", "CLL"))
cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(with(sampleAnnot, paste0(sample, cellType)))), function(x) sample(x, min(length(x), n))))]

dat <- dat[sample(1:nrow(dat), 100), cells]
dat <- dat[apply(dat, 1, max) != 0,,drop=F]
dat <- dat - apply(dat, 1, min)
dat <- dat / apply(dat,1, max)
apply(dat, 1, quantile, na.rm=T)
dat <- dat[, order(with(sampleAnnot[cells,],paste0(cellType, sample))),drop=F]

pdf(dirout(out, "RandomGenes", ".pdf"), height=min(29, nrow(dat) * 0.3 + 1), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
pheatmap(dat, cluster_rows=F, cluster_cols=F, annotation_col=pbmc@data.info[,c("sample", "cellType", "nUMI"), drop=F],show_colnames=FALSE)
dev.off()


# RANDOM MATRIX -----------------------------------------------------------
dat <- data
sampleAnnot <- subset(pbmc@data.info, grepl(pat, sample) & cellType %in% c("CD8", "CD4", "Mono", "CLL"))
cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(with(sampleAnnot, paste0(sample, cellType)))), function(x) sample(x, min(length(x), n))))]

dat <- dat[sample(1:nrow(dat), 100), cells]
dat <- dat[apply(dat, 1, max) != 0,,drop=F]
dat <- dat - apply(dat, 1, min)
dat <- dat / apply(dat,1, max)
apply(dat, 1, quantile, na.rm=T)
dat <- dat[, order(with(sampleAnnot[cells,],paste0(cellType, sample))),drop=F]

pdf(dirout(out, "RandomGenes_Magic", ".pdf"), height=min(29, nrow(dat) * 0.3 + 1), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
pheatmap(dat, cluster_rows=F, cluster_cols=F, annotation_col=pbmc@data.info[,c("sample", "cellType", "nUMI"), drop=F],show_colnames=FALSE)
dev.off()