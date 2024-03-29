require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(gplots)

project.init2("cll-time_course")
out <- "11_CellTypes/"
dir.create(dirout(out))

clustering.precision <- seq(0.5, 2.5, 0.2)

# ARGUMENTS ---------------------------------------------------------------
sample.x <- "allDataBest_NoDownSampling_noIGH"

outS <- paste0(out, "Signatures/")
dir.create(dirout(outS))

# Load original data
outOrig <- paste0("10_Seurat/", sample.x,"/")
load(dirout(outOrig, sample.x,".RData"))
if(is.null(pbmc@data.info[["patient"]])){
  pbmc@data.info[["patient"]] <- gsub("_.*", "", gsub("\\d", "", substr(gsub("LiveBulk_10x_", "", pbmc@data.info[["sample"]]), 0,6)))
  
  for(pat in unique(pbmc@data.info[["patient"]])){
    res <- pbmc@data.info[["sample"]]
    res[!grepl(paste0(pat, "\\d?_(d\\d+|\\d+d)"), pbmc@data.info[["sample"]])] <- "IGNORED"
    pbmc@data.info[[paste0("pat_",pat)]] <- res
  }
}
pbmcOrig <- pbmc
stopifnot(all(rownames(pbmc@data.info) == colnames(pbmc@data)))

# Gene signatures
file <- "/data/groups/lab_bock/shared/resources/gene_sets/MSigDB/6.0/Human//c7.all.v6.0.symbols.gmt"
lines <- readLines(file)
genesets <- list()
for(line in lines){
  x <- strsplit(line, "\t")[[1]]
  genesets[[x[1]]] <- x[3:length(x)]
}

# T-cell data
tcell.barcodes <- data.table(pbmcOrig@data.info, keep.rownames=TRUE)[ClusterNames_0.5 %in% c(6,2,3,9)]$rn
pDat <- pbmcOrig@tsne.rot
pDat$tcell <- rownames(pbmcOrig@data.info) %in% tcell.barcodes
pDat$nUMI <- pbmcOrig@data.info$nUMI
ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=tcell)) + geom_point()
ggsave(dirout(outS, "SigAnalysis_tCells.pdf"))

ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=nUMI)) + geom_point()

dat <- pbmcOrig@data[,tcell.barcodes]
# colSums <- apply(exp(dat), 2, sum)


# SOME TESTS --------------------------------------------------------------
x <- names(genesets)
x[grepl("CD4",x) & grepl("CD8",x)]
tcell.sets <- x[grepl("CD4",x) & grepl("CD8",x)]


# some inital trials
ex <- "GSE6259_CD4_TCELL_VS_CD8_TCELL"
for(ex in c("GSE6259_CD4_TCELL_VS_CD8_TCELL", "GSE8835_CD4_VS_CD8_TCELL", "GSE27786_CD4_VS_CD8_TCELL")){
  genes.cd8 <- genesets[[paste0(ex, "_UP")]]
  genes.cd8 <- genes.cd8[genes.cd8 %in% rownames(dat)]
  cd8.sig <- log(apply(exp(dat[genes.cd8,]),2,sum))
  
  genes.cd4 <- genesets[[paste0(ex, "_DN")]]
  genes.cd4 <- genes.cd4[genes.cd4 %in% rownames(dat)]
  cd4.sig <- log(apply(exp(dat[genes.cd4,]),2,sum))
  
  pdf(dirout(outS, "SigAnalysis_", ex, "_venn.pdf"))
  venn(list(cd8=genes.cd8, cd4=genes.cd4))
  dev.off()
  
  qplot(x=cd8.sig, y=cd4.sig) + geom_hex() + ggtitle(ex)
  ggsave(dirout(outS, "SigAnalysis_", ex, "_hex.pdf"))
  
  qplot(x=cd8.sig, y=cd4.sig) + geom_point() + ggtitle(ex)
  ggsave(dirout(outS, "SigAnalysis_", ex, "_points.pdf"))
  
  #   pDat <- data.table(cd8=cd8.sig, cd4=cd4.sig)
  #   pDat$label <- log10(pbmcOrig@data.info[colnames(dat),]$nUMI)
  #   pDat <- pDat[label > 3]
  #   ggplot(pDat, aes(x=cd8, y=cd4, color=label)) + geom_point()
}

pdf(dirout(outS, "SigAnalysis_venn.pdf"))
venn(list(
  GSE6259_up=genesets$GSE6259_CD4_TCELL_VS_CD8_TCELL_UP, 
  GSE8835_up=genesets$GSE8835_CD4_VS_CD8_TCELL_UP, 
  GSE27786_up=genesets$GSE27786_CD4_VS_CD8_TCELL_UP))
dev.off()


# Get signatures for all gene sets
X=matrix(0, nrow=length(tcell.barcodes), ncol=length(tcell.sets))
colnames(X) <- tcell.sets
rownames(X) <- tcell.barcodes
for(set.x in tcell.sets){
  genes <- genesets[[set.x]]
  genes <- genes[genes %in% rownames(dat)]
  X[,set.x] <- log(apply(exp(dat[genes,]),2,sum))
}


# PC and tSNE on those signatures
pr <- prcomp(X,scale=TRUE)
str(pr$rotation)
str(pr$x)
pDat <- data.table(pr$x, keep.rownames=TRUE)
ggplot(pDat, aes(x=PC1, y=PC2))+geom_hex()
ggsave(dirout(outS, "SigAnalysis_fromSignatures_PCs.pdf"))

tsne <- Rtsne::Rtsne(pr$x[,1:10],dims=2)
str(tsne)
pDat <- data.table(tsne$Y)
pDat$sample <- pbmcOrig@data.info[tcell.barcodes,]$sample 
pDat$CD8 <- dat["CD8A",]
ggplot(pDat, aes(x=V1, y=V2))+geom_hex()
ggplot(pDat, aes(x=V1, y=V2, color=sample))+geom_point() + geom_abline(intercept=0, slope=0.8)
ggsave(dirout(outS, "SigAnalysis_fromSignatures_TSNEs.pdf"))

# What's the difference between the two groups?
pDat[, grp := (V2 + (-0.8*V1)) > 0]
pDat$Var1 <- rownames(pr$x)
ggplot(pDat, aes(x=V1, y=V2, color=grp))+geom_point() + geom_abline(intercept=0, slope=0.8)

# boxDat <- data.table(melt(X))
# boxDat <- merge(boxDat, pDat, by="Var1")
# boxDat[,boxgrp := paste0(Var2,"_", grp)]
# sig = list()
# for(var in unique(boxDat$Var2)){
#   sig[[var]] <- t.test(boxDat[Var2 == var & grp]$value, boxDat[Var2 == var & !grp]$value)$p.value
# }
# boxDat2 <- boxDat[Var2 %in% names(sort(unlist(sig))[1:10])]
# ggplot(boxDat2, aes(x=boxgrp, y=value, group=boxgrp, color=grp)) + geom_boxplot() + coord_flip()

# PC and tSNE on the genes in the signatures
genes <- unique(do.call(c, genesets))
genes <- genes[genes %in% rownames(dat)]
genes <- genes[apply(dat[genes,]==0, 1, sum) > 0]
pr2 <- prcomp(scale(t(pbmcOrig@scale.data[genes,tcell.barcodes])), scale=FALSE, center=FALSE)
biplot(pr2)
pDat <- data.table(pr2$x, keep.rownames=TRUE)
ggplot(pDat, aes(x=PC1, y=PC2))+geom_hex()
ggsave(dirout(outS, "SigAnalysis_fromSigGenes_PCs.pdf"))

tsne2 <- Rtsne::Rtsne(pr2$x[,1:10],dims=2)
str(tsne2)
pDat <- data.table(tsne2$Y)
pDat$sample <- pbmcOrig@data.info[tcell.barcodes,]$sample 
ggplot(pDat, aes(x=V1, y=V2))+geom_hex()
ggplot(pDat, aes(x=V1, y=V2,color=sample))+geom_point()
ggsave(dirout(outS, "SigAnalysis_fromSigGenes_TSNEs.pdf"))


# Show for one example GSE8835
pDat <- pbmcOrig@tsne.rot[tcell.barcodes,]
pDat$score <- X[,"GSE8835_CD4_VS_CD8_TCELL_UP"]
ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=score)) + geom_point()
ggsave(dirout(outS, "SigAnalysis_GSE8835.pdf"))
