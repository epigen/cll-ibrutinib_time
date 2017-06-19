require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(pheatmap)
require(gplots)
require(enrichR) #devtools::install_github("definitelysean/enrichR")


normalization = ""
# normalization = "normalize.by.sample"
# normalization = "regress.out.nUMI"
normalization = "normalize.by.timepoint.cell"
# stopifnot((normalize.by.sample | regress.out.nUMI) | !(normalize.by.sample | regress.out.nUMI))

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
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out.fscLVM.x, "Relevance.pdf"))
# association of cells with terms (X)
str(x <- as.matrix(fread(dirout(out.fscLVM.x, "X.csv"))))
colnames(x) <- terms

# Genes associated with terms
str(z <- as.matrix(fread(dirout(out.fscLVM.x, "Z.csv"))))
str(genes <- read.csv(dirout(out.fscLVM.x, "Genes_IDs.csv"), header=F)$V1)
colnames(z) <- terms
rownames(z) <- genes
term <- "hidden0"
hitSets <- list()
for(term in as.character(rel$term)){
  hitSets[[term]] <- genes[z[,term] == 1]
}

# ENRICHR
enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location")
enrichRes <- data.table()
for(grp.x in names(hitSets)[1:10]){
  ret=try(as.data.table(enrichGeneList(hitSets[[grp.x]],databases = enrichrDBs)),silent = FALSE)
  if(!any(grepl("Error",ret)) && nrow(ret) > 0){
    enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
  }
}
enrichRes$n <- sapply(strsplit(enrichRes$genes,","), length)
enrichRes <- enrichRes[n > 3]
write.table(enrichRes[qval < 0.05], file=dirout(out.fscLVM.x, "0_EnrichR", ".tsv"), sep="\t", quote=F, row.names=F)
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
    pDatM <- pDatM[order(apply(pDatM, 1, sum), decreasing=TRUE),][1:50,]
    pdf(dirout(out.fscLVM.x, "0_EnrichR", ".pdf"),onefile=FALSE, width=min(29, 6+ ncol(pDatM)*0.3), height=min(29, nrow(pDatM)*0.3 + 4))
    pheatmap(pDatM) #, color=gray.colors(12, start=0, end=1), border_color=NA)
    dev.off()
  }
}
dir.create(dirout(out.fscLVM.x, "Enrichr/"))
for(grp.x in names(hitSets)[1:10]){
  pDat <- enrichRes[grp == grp.x][order(qval)]
  pDat$category <- gsub("\\_(\\w|\\d){8}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){12}", "", pDat$category)
  pDat$category <- factor(pDat$category, levels=rev(pDat$category))
  ggplot(pDat[1:min(20, nrow(pDat))], aes(x=-log10(qval), y=category)) + geom_point() + ggtitle(grp.x)
  ggsave(dirout(out.fscLVM.x, "Enrichr/", grp.x, ".pdf"), height=5, width=13)
}




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
out.analysis <-  paste0(out.fscLVM.x, "analysis", ifelse(normalization != "", "_", ""),normalization,"/")
# if(normalization == "normalize.by.sample") out.analysis <- paste0(out.fscLVM.x, "analysis_normalized/") 
# if(normalization == "regress.out.nUMI") out.analysis <- paste0(out.fscLVM.x, "analysis_regressedOutUMIs/")
dir.create(dirout(out.analysis))

# Normalize if specified
if(normalization == "normalize.by.sample"){
  for(term.x in terms){
    pDat[, (term.x) := scale(get(term.x))[,1], by="sample"]    
  }
}
if(normalization == "regress.out.nUMI"){
  for(term.x in terms){
    pDat[["xxx"]] <- pDat[[term.x]]
    lm <- lm(data=pDat, xxx ~ nUMI)
    pDat[[term.x]] <- lm$residuals + lm$coefficients[1]
    pDat[["xxx"]] <- NULL
  }
}
if(normalization == "normalize.by.timepoint.cell"){
  for(term.x in terms){
    pDat[, (term.x) := scale(get(term.x))[,1], by=c("timepoint", "cellType")]
  }
}

# ANALYSIS OVER TIME ----------------------------------------------------------------
pDat2 <- pDat
pDat2 <- pDat2[timepoint %in% c("d0", "d120")]
pDat2 <- pDat2[cellType %in% c("CD8", "CD4", "Mono", "CLL")]# & grepl("PT", sample)]


# normalize by sample to remove strange effects?
# ggplot(pDat2, aes(x=sample, y=Allograft.rejection)) + geom_violin() + coord_flip()
# pDat2[,Allograft.rejection2 := scale(Allograft.rejection)[,1], by="sample"]
# ggplot(pDat2, aes(x=sample, y=Allograft.rejection)) + geom_violin() + coord_flip()



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
overTime.sig[,cellPat := paste0(cellType, "_", patient)]
ggplot(
  overTime.sig, 
  aes(x=cellPat, y=geneset, color=medianDiff, size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw(24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave(dirout(out.analysis, "0_Overview.pdf"),height=15, width=15)

# Ordered Heatmap
orderDT <- dcast.data.table(overTime.sig, cellPat ~ geneset, value.var="medianDiff")
orderMT <- as.matrix(orderDT[,-"cellPat"])
rownames(orderMT) <- orderDT$cellPat
pheatmap(orderMT)
overTime.sig$cellPat <- factor(overTime.sig$cellPat, levels=rownames(orderMT)[hclust(dist(orderMT))$order])
overTime.sig$geneset <- factor(overTime.sig$geneset, levels=colnames(orderMT)[hclust(dist(t(orderMT)))$order])
ggplot(
  overTime.sig, 
  aes(x=cellPat, y=geneset, color=medianDiff, size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw(24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave(dirout(out.analysis, "0_Overview_Ordered.pdf"),height=15, width=15)



# PLOT Gene Sets over time
pat <- "PT"
n <- 100
set.nam <- "P53.pathway"
for(set.nam in names(hitSets)){
  for(pat in unique(pDat$patient)){
    genes <- hitSets[[set.nam]]
    dat <- pbmc@data
    sampleAnnot <- subset(pbmc@data.info, grepl(pat, sample) & cellType %in% c("CD8", "CD4", "Mono", "CLL"))
    cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(with(sampleAnnot, paste0(sample, cellType)))), function(x) sample(x, min(length(x), n))))]
    
    genes <- genes[genes %in% rownames(dat)]
    dat <- dat[genes, cells]
    dat <- dat[apply(dat, 1, max) != 0,, drop=F]
    if(nrow(dat) > 0){
      dat <- dat - apply(dat, 1, min)
      dat <- dat / apply(dat,1, max)
      apply(dat, 1, quantile, na.rm=T)
      dat <- dat[, order(with(sampleAnnot[cells,],paste0(cellType, sample))),drop=F]
      
      if(nrow(dat) > 1 & ncol(dat) >1){
        pdf(dirout(out.analysis, set.nam, "_", pat, ".pdf"), height=min(29, nrow(dat) * 0.3 + 3), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
        pheatmap(dat, cluster_rows=F, cluster_cols=F, annotation_col=pbmc@data.info[,c("sample", "cellType", "nUMI"), drop=F],show_colnames=FALSE)
        dev.off()
      }
    }
  }
}



# ANALYSIS T_ZERO ----------------------------------------------------------------
pDat2 <- pDat
pDat2 <- pDat2[timepoint %in% c("d0")]
pDat2 <- pDat2[cellType %in% c("CD8", "CD4", "Mono", "CLL")]# & grepl("PT", sample)]

# Make one large heatmap
overTime.sig <- data.table()
pat="FE"
cell="CD4"
geneset <- "Allograft.rejection"
for(cell in unique(pDat2$cellType)){
  x <- pDat2[cellType == cell]
  x <- x[patient %in% x[,.N, by="patient"][N>5]$patient]
  for(pat in unique(x$patient)){
    for(geneset in names(genesets)){
      p <- t.test(x[patient == pat][[geneset]], x[patient != pat][[geneset]])$p.value
      ef <- median(x[patient == pat][[geneset]])- median(x[patient != pat][[geneset]])
      overTime.sig <- rbind(overTime.sig, data.table(patient=pat, cellType=cell, pvalue=p, medianDiff=ef, geneset=geneset))
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
ggsave(dirout(out.analysis, "0_Overview_TZero.pdf"),height=15, width=15)



# PLOT Gene Sets over time
n <- 100
set.nam <- "P53.pathway"
for(set.nam in names(hitSets)){
    genes <- hitSets[[set.nam]]
    dat <- pbmc@data
    sampleAnnot <- subset(pbmc@data.info, cellType %in% c("CD8", "CD4", "Mono", "CLL") & (grepl("d0", sample) | grepl("_0d", sample)))
    cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(with(sampleAnnot, paste0(sample, cellType)))), function(x) sample(x, min(length(x), n))))]
    
    genes <- genes[genes %in% rownames(dat)]
    dat <- dat[genes, cells]
    dat <- dat[apply(dat, 1, max) != 0,, drop=F]
    if(nrow(dat) > 0){
      dat <- dat - apply(dat, 1, min)
      dat <- dat / apply(dat,1, max)
      apply(dat, 1, quantile, na.rm=T)
      dat <- dat[, order(with(sampleAnnot[cells,],paste0(cellType, sample))),drop=F]
      
      if(nrow(dat) > 1 & ncol(dat) >1){
        pdf(dirout(out.analysis, set.nam, "_TZero", ".pdf"), height=min(29, nrow(dat) * 0.3 + 3), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
        pheatmap(dat, cluster_rows=F, cluster_cols=F, annotation_col=pbmc@data.info[,c("sample", "cellType", "nUMI"), drop=F],show_colnames=FALSE)
        dev.off()
      }
    }
}