require(data.table)
require(ggplot2)

# outFSCLVM <- paste0(outS, "fscLVM/fscLVM_files_3hidden/")
stopifnot(dir.exists(dirout(outFSCLVM)))


options(stringsAsFactors=FALSE)

outX <- paste0(outFSCLVM, "results/")
if(exists("balance.barcodes") && balance.barcodes){
  outX <- paste0(outFSCLVM, "results_balanced/")
}
dir.create(dirout(outX))


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


# balance samples if selected
table(pDat$sample)
if(exists("balance.barcodes") && balance.barcodes){
  barcodes.max <- min(table(pDat$sample))
  barcodes.keep <- do.call(c, lapply(split(pDat$barcode, factor(pDat$sample)), function(barcodes) sample(barcodes, barcodes.max)))
  pDat <- pDat[barcode %in% barcodes.keep]
}
table(pDat$sample)


# READ fscLVM files -------------------------------------------------------
rel <- fread(dirout(outFSCLVM, "Relevance.csv"))

str(x <- as.matrix(fread(dirout(outFSCLVM, "X.csv"))))
str(z <- as.matrix(fread(dirout(outFSCLVM, "Z.csv"))))
str(zChanged <- as.matrix(fread(dirout(outFSCLVM, "Zchanged.csv"))))
str(w <- as.matrix(fread(dirout(outFSCLVM, "W.csv"))))

str(terms <- make.names(read.csv(dirout(outFSCLVM, "Terms.csv"), header=F)$V1))
str(genes <- read.csv(dirout(outFSCLVM, "Genes_IDs.csv"), header=F)$V1)
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

# Plot relevance of terms
ggplot(rel[1:15], aes(x=term, y=value)) + geom_point() + coord_flip() + 
  theme_bw(24) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(outX, "Relevance.pdf"))





# CALCULATE SCORES AND FIT MODEL --------------------------------------------------------
genesets <- as.character(rel[1:15]$term)
lm.pvalues <- data.table(geneset=genesets)
lm.effect <- data.table(geneset=genesets)

if(!file.exists(dirout(outX, "Pvalues.csv"))){
  for(set.nam in genesets){    
    pDat$score <- x[,set.nam]
    lm.coef <- summary(lm(data=pDat, score ~ patient / time))$coefficients
    lm.coef <- data.table(lm.coef, keep.rownames=TRUE)[rn != "(Intercept)"]
    
    for(rn.x in lm.coef$rn){
      lm.pvalues[geneset==set.nam, eval(rn.x) := lm.coef[rn.x == rn]$"Pr(>|t|)"]
      lm.effect[geneset==set.nam, eval(rn.x) := lm.coef[rn.x == rn]$Estimate]
    }
  }
  
  write.csv(lm.pvalues, file=dirout(outX, "Pvalues.csv"), row.names=F)
  write.csv(lm.effect, file=dirout(outX, "EffectSize.csv"), row.names=F)
  
} else {
  lm.pvalues <- fread(dirout(outX, "Pvalues.csv"))
  lm.effect <- fread(dirout(outX, "EffectSize.csv"))
}


# this time without multiple testing correction
mat.eff <- as.matrix(lm.effect[,-"geneset",with=F])
mat.pval <- as.matrix(lm.pvalues[,-"geneset", with=F])
rownames(mat.eff) <- lm.effect$geneset
rownames(mat.pval) <- lm.pvalues$geneset
# test that they are the same
stopifnot(nrow(mat.eff) == nrow(mat.pval) & ncol(mat.eff) == ncol(mat.pval))


mat.pval.signed <- sign(mat.eff) * mat.pval
mat.pval.ind <- mat.pval < 0.05
class(mat.pval.ind) <- "integer"
mat.eff.significant <- mat.eff * mat.pval.ind






# TIMELINE ANALYSIS -------------------------------------------------------

# Venn diagrams (use pvalues for this)
mat <- mat.pval.signed[,grepl(":", colnames(mat.pval.signed))]
hitLists <- list()
for(col in colnames(mat)){
  val <- mat[,col]
  hitLists[[paste0(col, "_up")]] <- rownames(mat)[val < 0.05 & val > 0]
  hitLists[[paste0(col, "_down")]] <- rownames(mat)[val > -0.05 & val < 0]
}
names(hitLists) <- gsub(":timelate", "", gsub("patient", "", names(hitLists)))
try({
  pdf(dirout(outX, "Overtime_Venn_up.pdf"))
  venn(hitLists[grepl("_up", names(hitLists))])# | grepl("PT", names(hitLists))])
  dev.off()
}, silent=TRUE)
try({
  pdf(dirout(outX, "Overtime_Venn_down.pdf"))
  venn(hitLists[grepl("_down", names(hitLists))])#  | grepl("PT", names(hitLists))])
  dev.off()
}, silent=TRUE)
ggplot(data.table(list=names(hitLists), size=sapply(hitLists, length)), aes(x=list, y=size))+geom_bar(stat="identity") + coord_flip()
ggsave(dirout(outX, "Overtime_Listsize.pdf"), width=9, height=7)


# Heatmap (plot FCs)
mat <- mat.eff.significant[,grepl(":", colnames(mat.eff.significant))]
mat <- mat[apply(abs(mat), 1, sum) > 0,]

if(nrow(mat) > 50){
  rows <- order(apply(mat, 1, sum), decreasing=T)[1:25]
  rows <- c(rows, order(apply(mat, 1, sum), decreasing=F)[1:25])
  mat <- mat[rows,]
}

if(nrow(mat) > 2){
  mat2 <- mat
  cutval <- min(
    max(abs(mat2[mat2 < 0])), # smallest neg value
    max(mat2[mat2 > 0]) # largest positive value
  )
  mat2[mat2 > cutval] <- cutval
  mat2[mat2 < -cutval] <- -cutval
  
  pdf(dirout(outX, "Overtime.pdf"), width=10, height=min(16, nrow(mat2)*0.25+3),onefile=F)
  pheatmap(mat2)
  dev.off()
  
  # boxplots for individual genesets
  dir.create(dirout(outX, "OverTime/"))
  for(set.nam in rownames(mat)){
    pDat$score <- x[,set.nam]
    pDat2 <- pDat[patient != "KI"]
    
    ggplot(pDat2, aes(y=score, x=sample)) + geom_jitter(color="grey", alpha=0.5)+ geom_boxplot(fill="NA")  + coord_flip() + 
      ggtitle(paste(set.nam))
    ggsave(dirout(outX, "OverTime/", set.nam, ".pdf"))
  }
}




# TIMEPOINT ZERO ANALYSIS -------------------------------------------------------

# Venn diagrams (use pvalues for this)
mat <- mat.pval.signed[,!grepl(":", colnames(mat.pval.signed)) & colnames(mat.pval.signed) != "geneset"]
hitLists <- list()
for(col in colnames(mat)){
  val <- mat[,col]
  hitLists[[paste0(col, "_up")]] <- rownames(mat)[val < 0.05 & val > 0]
  hitLists[[paste0(col, "_down")]] <- rownames(mat)[val > -0.05 & val < 0]
}
names(hitLists) <- gsub(":timelate", "", gsub("patient", "", names(hitLists)))
try({
  pdf(dirout(outX, "TimepointZero_Venn_up.pdf"))
  venn(hitLists[grepl("_up", names(hitLists))])# | grepl("PT", names(hitLists))])
  dev.off()
}, silent=TRUE)
try({
  pdf(dirout(outX, "TimepointZero_Venn_down.pdf"))
  venn(hitLists[grepl("_down", names(hitLists))])#  | grepl("PT", names(hitLists))])
  dev.off()
}, silent=TRUE)
ggplot(data.table(list=names(hitLists), size=sapply(hitLists, length)), aes(x=list, y=size))+geom_bar(stat="identity") + coord_flip()
ggsave(dirout(outX, "TimepointZero_Listsize.pdf"), width=9, height=7)

# Heatmap (plot FCs)
mat <- mat.eff.significant[,!grepl(":", colnames(mat.eff.significant)) & colnames(mat.eff.significant) != "geneset"]
mat <- mat[apply(abs(mat), 1, sum) > 0,]

if(nrow(mat) > 50){
  rows <- order(apply(mat, 1, sum), decreasing=T)[1:25]
  rows <- c(rows, order(apply(mat, 1, sum), decreasing=F)[1:25])
  mat <- mat[rows,]
}

if(nrow(mat) > 2){
  mat2 <- mat
  cutval <- min(
    max(abs(mat2[mat2 < 0])), # smallest neg value
    max(mat2[mat2 > 0]) # largest positive value
  )
  mat2[mat2 > cutval] <- cutval
  mat2[mat2 < -cutval] <- -cutval
  
  pdf(dirout(outX, "TimepointZero.pdf"), width=10, height=min(16, nrow(mat2)*0.25+3),onefile=F)
  pheatmap(mat2)
  dev.off()
  
  # boxplots for individual genesets
  dir.create(dirout(outX, "TimepointZero/"))
  for(set.nam in rownames(mat)){
    pDat$score <- x[,set.nam]
    pDat2 <- pDat[time == "early"]
    
    ggplot(pDat2, aes(y=score, x=sample)) + geom_jitter(color="grey", alpha=0.5)+ geom_boxplot(fill="NA")  + coord_flip() + 
      ggtitle(paste(set.nam))
    ggsave(dirout(outX, "TimepointZero/", set.nam, ".pdf"))
  }
}