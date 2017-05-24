require(data.table)
require(ggplot2)

outX <- paste0(outS, "fscLVM/fscLVM_files_3hidden/")

options(stringsAsFactors=FALSE)


# READ DATA
pDat <- data.table(pbmc@tsne.rot)
pDat$barcode <- colnames(pbmc@data)
pDat$sample <- pbmc@data.info$sample
pDat$patient <- sapply(strsplit(pDat$sample, "_"), function(x) return(gsub("\\d", "", x[1])))
pDat$time <- sapply(strsplit(pDat$sample, "_"), function(x) return(gsub("d150", "d120", gsub("(\\d+)d", "d\\1", x[length(x)]))))

# LIMIT TO TWO TIMEPOINTS
pDat <- pDat[time != "d280"]
pDat[time == "d0", time := "early"]
pDat[time == "d120", time := "late"]
pDat$time <- factor(pDat$time, levels=c("early", "late"))
pats <- unique(pDat$patient)
pDat$patient <- factor(pDat$patient, levels = c("FE", pats[pats != "FE"]))


# READ fscLVM files -------------------------------------------------------
rel <- fread(dirout(outX, "Relevance.csv"))

str(x <- as.matrix(fread(dirout(outX, "X.csv"))))
str(z <- as.matrix(fread(dirout(outX, "Z.csv"))))
str(zChanged <- as.matrix(fread(dirout(outX, "Zchanged.csv"))))
str(w <- as.matrix(fread(dirout(outX, "W.csv"))))

str(terms <- make.names(read.csv(dirout(outX, "Terms.csv"), header=F)$V1))
str(genes <- read.csv(dirout(outX, "Genes_IDs.csv"), header=F)$V1)
str(barcodes <- colnames(pbmc@data))

# label matrices...
colnames(x) <- terms
rownames(x) <- barcodes
x <- x[pDat$barcode,]
colnames(z) <- terms
rownames(z) <- genes
colnames(zChanged) <- terms
rownames(zChanged) <- genes
colnames(w) <- terms
rownames(w) <- genes
colnames(rel)[1] <- "value"
rel$term <- terms

rel <- rel[order(value, decreasing=TRUE)]
rel$term <- factor(rel$term, levels=rel$term)


ggplot(rel[1:15], aes(x=term, y=value)) + geom_point() + coord_flip() + 
  theme_bw(24) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(outX, "Relevance.pdf"))



# OVER TIME ---------------------------------------------------------------
terms2 <- rel[1:15]$term
res <- data.table(terms=terms2)
if(!file.exists(dirout(outX, "Pvalues.csv"))){
  for(term in terms2){
    pDat$score <- x[,term]
    lm.coef <- summary(lm(data=pDat, score ~ patient / time))$coefficients
    
    lm.coef <- data.table(lm.coef, keep.rownames=TRUE)[rn != "(Intercept)"]
    lm.coef$pvalue <- lm.coef$"Pr(>|t|)" * sign(lm.coef$"t value")
    for(rn.x in lm.coef$rn){
      res[terms==term, eval(rn.x) := lm.coef[rn.x == rn]$pvalue]
    }
  }
  
  write.csv(res, file=dirout(outGeneSets, "Pvalues.csv"), row.names=F)
} else {
  res <- fread(dirout(outGeneSets, "Pvalues.csv"))
}





# TIMEPOINT ZERO  -------------------------------------------------------
mat <- as.matrix(res[,!grepl(":", colnames(res)) & colnames(res) != "terms", with=F])
rownames(mat) <- res$terms

# heatmap
logMat <- -log10(abs(mat))
logMat[logMat >= 5] <- 5
logMat <- sign(mat) * logMat

pdf(dirout(outX, "TimepointZero.pdf"), width=10, height=16,onefile=F)
pheatmap(logMat)
dev.off()

# boxplots for individual genesets
dir.create(dirout(outX, "TimepointZero/"))
for(term in rownames(logMat)){
  pDat$score <- x[,term]
  pDat2 <- pDat[time == "early"]
  
  ggplot(pDat2, aes(y=score, x=sample)) + geom_jitter(color="grey", alpha=0.5)+ geom_boxplot(fill="NA")  + coord_flip() + 
    ggtitle(paste(term))
  ggsave(dirout(outX, "TimepointZero/", term, ".pdf"))
}



mat <- as.matrix(res[,!grepl(":", colnames(res)) & colnames(res) != "terms", with=F])
mat <- as.matrix(res[,grepl(":", colnames(res)), with=F])
rownames(mat) <- res$terms

# heatmap
logMat <- -log10(abs(mat))
logMat[logMat >= 5] <- 5
logMat <- sign(mat) * logMat

pdf(dirout(outX, "Overtime.pdf"), width=10, height=16,onefile=F)
pheatmap(logMat)
dev.off()

# boxplots for individual genesets
dir.create(dirout(outX, "OverTime/"))
for(term in rownames(logMat)){
  pDat$score <- x[,term]
  pDat2 <- pDat[patient != "KI"]
  
  ggplot(pDat2, aes(y=score, x=sample)) + geom_jitter(color="grey", alpha=0.5)+ geom_boxplot(fill="NA")  + coord_flip() + 
    ggtitle(paste(term))
  ggsave(dirout(outX, "OverTime/", term, ".pdf"))
}