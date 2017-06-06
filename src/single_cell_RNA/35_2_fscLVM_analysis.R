require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(pheatmap)
require(gplots)


# READ fscLVM files -------------------------------------------------------
str(terms <- make.names(read.csv(dirout(out.fscLVM.x, "Terms.csv"), header=F)$V1))
# Relevance
rel <- fread(dirout(out.fscLVM.x, "Relevance.csv"))
colnames(rel)[1] <- "value"
rel$term <- terms
rel <- rel[order(value, decreasing=TRUE)]
rel$term <- factor(rel$term, levels=rel$term)
# Plot relevance of terms
ggplot(rel[1:15], aes(x=term, y=value)) + geom_point() + coord_flip() + 
  theme_bw(24) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out.fscLVM.x, "Relevance.pdf"))
# association of cells with terms (X)
str(x <- as.matrix(fread(dirout(out.fscLVM.x, "X.csv"))))
colnames(x) <- terms


# Prepare sample annotation ------------------------------------------------------------
pDat <- data.table(pbmc@data.info, keep.rownames=TRUE)
pDat$sample_cell <- paste0(pDat$sample, "_", pDat$cellType)
pDat$patient <- gsub("([A-Z]+)\\d?_.+", "\\1", pDat$sample)
pDat$sample <- gsub("([0-9]+)d", "d\\1", pDat$sample)
pDat$sample <- gsub("d150", "d120", pDat$sample)
pDat$timepoint <- gsub(".+_(d[0-9]+)", "\\1", pDat$sample)


# Merge sample annotation and fscLVM output ------------------------------------------------------------
stopifnot(all(colnames(pbmc@data) == rownames(pbmc@data.info)))
pDat <- data.table(pDat, data.table(x))
genesets <- terms
names(genesets) <- terms
out.analysis <- paste0(out.fscLVM.x, "analysis/")
dir.create(dirout(out.analysis))




# ANALYSIS ----------------------------------------------------------------
pDat2 <- pDat
pDat2 <- pDat2[timepoint %in% c("d0", "d120")]
pDat2 <- pDat2[cellType %in% c("CD8", "CD4", "Mono", "CLL")]# & grepl("PT", sample)]

# geneset <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"
for(geneset in names(genesets)){
  
  # Plot all values as boxplots
  ggplot(
    pDat2,
    aes_string(x="cellType", y=geneset, group="sample_cell", fill="patient", alpha="timepoint")) + 
    geom_boxplot(outlier.shape=NA) + coord_flip() + 
    scale_color_manual(values=c("grey", "black")) + scale_alpha_manual(values=c(0.35,1))
  ggsave(dirout(out.analysis, geneset, "_AllData.pdf"), width=7, height=7)
  
  # Only t0 as boxplot
  ggplot(
    pDat2[timepoint == "d0"],
    aes_string(x="cellType", y=geneset, group="sample_cell", fill="patient")) + 
    geom_boxplot(outlier.shape=NA) + coord_flip()
  ggsave(dirout(out.analysis, geneset, "_t0.pdf"), width=7, height=7)
  
  # just the change over time as heatmap
  sDat <- pDat2[,median(get(geneset)), by=c("patient", "timepoint", "cellType")]
  sDat <- dcast.data.table(sDat, patient + cellType ~ timepoint, value.var="V1")
  sDat[, timechange := d120 - d0]
  ggplot(sDat[!is.na(timechange)], aes(x=cellType, y=patient, fill=timechange)) + 
    geom_tile() + scale_fill_gradient2(low="blue", high="red", mid="black")
  ggsave(dirout(out.analysis, geneset, "_HM.pdf"), width=5, height=4)
}

# venn(list(alpha=genesets$HALLMARK_INTERFERON_ALPHA_RESPONSE, gamma=genesets$HALLMARK_INTERFERON_GAMMA_RESPONSE))


# Make one large heatmap
overTime.sig <- data.table()
pat="FE"
cell="CD4"
geneset <- "Allograft.rejection"
for(pat in unique(pDat2$patient)){
  for(cell in unique(pDat2$cellType)){
    x <- pDat2[patient == pat & cellType == cell]
    if(nrow(x[timepoint == "d0"]) > 5 & nrow(x[timepoint == "d120"]) >5){
      for(geneset in names(genesets)){
        p <- t.test(x[timepoint == "d0"][[geneset]], x[timepoint == "d120"][[geneset]])$p.value
        ef <- median(x[timepoint == "d120"][[geneset]])- median(x[timepoint == "d0"][[geneset]])
        overTime.sig <- rbind(overTime.sig, data.table(patient=pat, cellType=cell, pvalue=p, medianDiff=ef, geneset=geneset))
      }
    }
  }
}
overTime.sig[,pvalue2 := pmin(5, -1*log10(pvalue))]
ggplot(
  overTime.sig, 
  aes(x=paste0(cellType, "_", patient), y=geneset, color=medianDiff, size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out.analysis, "0_Overview.pdf"),height=15, width=15)





# # PLOT INDIVIDUAL EXAMPLES ------------------------------------------------
# 
# pat <- "PT"
# set.nam <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
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
#     dat <- dat[apply(dat, 1, max) != 0,]
#     dat <- dat - apply(dat, 1, min)
#     dat <- dat / apply(dat,1, max)
#     apply(dat, 1, quantile, na.rm=T)
#     dat <- dat[, order(with(sampleAnnot[cells,],paste0(cellType, sample)))]
#     
#     pdf(dirout(out, set.nam, "_", pat, ".pdf"), height=min(29, nrow(dat) * 0.3 + 3), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
#     pheatmap(dat, cluster_rows=F, cluster_cols=F, annotation_col=pbmc@data.info[,c("sample", "cellType", "nUMI"), drop=F],show_colnames=FALSE)
#     dev.off()
#   }
# }