
require(pheatmap)
require(gplots)


outGeneSets <- paste0(outS, "Genesets2/", file.nam, "/")
dir.create(dirout(outGeneSets))

file <- geneSetFiles[[file.nam]]


lines <- readLines(file)
genesets <- list()
for(line in lines){
  x <- strsplit(line, "\t")[[1]]
  genesets[[x[1]]] <- x[3:length(x)]
}


# PREPARE SC DATA ---------------------------------------------------------
pDat <- data.table(pbmc@tsne.rot)
pDat$barcode <- colnames(pbmc@data)
pDat$sample <- pbmc@data.info$sample
table(pDat$sample)
pDat$patient <- sapply(strsplit(pDat$sample, "_"), function(x) return(gsub("\\d", "", x[1])))
table(pDat$patient)
pDat$time <- sapply(strsplit(pDat$sample, "_"), function(x) return(gsub("d150", "d120", gsub("(\\d+)d", "d\\1", x[length(x)]))))
table(pDat$time)

dat <- pbmc@data

# LIMIT TO TWO TIMEPOINTS
dat <- dat[,pDat$time != "d280"]
pDat <- pDat[time != "d280"]
pDat[time == "d0", time := "early"]
pDat[time == "d120", time := "late"]
pDat$time <- factor(pDat$time, levels=c("early", "late"))
pats <- unique(pDat$patient)
pDat$patient <- factor(pDat$patient, levels = c("FE", pats[pats != "FE"]))


# Plot number of genes and UMIs ----------------------------------------------------
ggplot(pbmc@data.info, aes(y=nUMI, x=sample)) + geom_jitter(color="grey", alpha=0.5)+ geom_boxplot(fill="NA")  + coord_flip()
ggsave(dirout(outGeneSets, "UMIs.pdf"))
ggplot(pbmc@data.info, aes(y=nGene, x=sample)) + geom_jitter(color="grey", alpha=0.5)+ geom_boxplot(fill="NA")  + coord_flip()
ggsave(dirout(outGeneSets, "GENES.pdf"))



# CALCULATE SCORES AND FIT MODEL --------------------------------------------------------
res <- data.table(geneset=names(genesets))
colSums <- apply(dat, 2, sum)

stopifnot(all(apply(dat, 2, sum)/colSums==1))

if(!file.exists(dirout(outGeneSets, "Pvalues.csv"))){
	for(set.nam in names(genesets)){
	  genes <- genesets[[set.nam]]
	  genes <- genes[genes %in% rownames(dat)]
  
	  pDat$score <- apply(dat[genes,], 2, sum)/colSums  
	  lm.coef <- summary(lm(data=pDat, score ~ patient / time))$coefficients
  
	  lm.coef <- data.table(lm.coef, keep.rownames=TRUE)[rn != "(Intercept)"]
	  lm.coef$pvalue <- lm.coef$"Pr(>|t|)" * sign(lm.coef$"t value")
	  for(rn.x in lm.coef$rn){
	    res[geneset==set.nam, eval(rn.x) := lm.coef[rn.x == rn]$pvalue]
	  }
	}

	write.csv(res, file=dirout(outGeneSets, "Pvalues.csv"), row.names=F)
} else {
	res <- fread(dirout(outGeneSets, "Pvalues.csv"))
}




# TIMELINE ANALYSIS -------------------------------------------------------
mat <- as.matrix(res[,grepl(":", colnames(res)), with=F])
rownames(mat) <- res$geneset

# Venn diagrams
hitLists <- list()
for(col in colnames(mat)){
  val <- mat[,col]
  hitLists[[paste0(col, "_up")]] <- rownames(mat)[val < 0.05 & val > 0]
  hitLists[[paste0(col, "_down")]] <- rownames(mat)[val > -0.05 & val < 0]
}
names(hitLists) <- gsub(":timelate", "", gsub("patient", "", names(hitLists)))
pdf(dirout(outGeneSets, "Overtime_Venn_up.pdf"))
venn(hitLists[grepl("_up", names(hitLists))])# | grepl("PT", names(hitLists))])
dev.off()
pdf(dirout(outGeneSets, "Overtime_Venn_down.pdf"))
venn(hitLists[grepl("_down", names(hitLists))])#  | grepl("PT", names(hitLists))])
dev.off()

# heatmap
logMat <- -log10(abs(mat))
logMat[logMat >= 5] <- 5
logMat <- sign(mat) * logMat

rows <- order(apply(logMat, 1, sum), decreasing=T)[1:25]
rows <- c(rows, order(apply(logMat, 1, sum), decreasing=F)[1:25])

pdf(dirout(outGeneSets, "Overtime.pdf"), width=10, height=16,onefile=F)
pheatmap(logMat[rows,])
dev.off()

# boxplots for individual genesets
dir.create(dirout(outGeneSets, "OverTime/"))
for(set.nam in rownames(logMat)[rows]){
  genes <- genesets[[set.nam]]
  genes <- genes[genes %in% rownames(dat)]
  
  pDat$score <- apply(dat[genes,], 2, sum)/colSums  
  pDat2 <- pDat[patient != "KI"]
  
  ggplot(pDat2, aes(y=score, x=sample)) + geom_jitter(color="grey", alpha=0.5)+ geom_boxplot(fill="NA")  + coord_flip() + 
    ggtitle(paste(set.nam))
  ggsave(dirout(outGeneSets, "OverTime/", set.nam, ".pdf"))
}




# TIMEPOINT ZERO ANALYSIS -------------------------------------------------
mat <- as.matrix(res[,!grepl(":", colnames(res)) & colnames(res) != "geneset", with=F])
rownames(mat) <- res$geneset

# Venn diagrams
hitLists <- list()
for(col in colnames(mat)){
  val <- mat[,col]
  hitLists[[paste0(col, "_up")]] <- rownames(mat)[val < 0.05 & val > 0]
  hitLists[[paste0(col, "_down")]] <- rownames(mat)[val > -0.05 & val < 0]
}
names(hitLists) <- gsub(":timelate", "", gsub("patient", "", names(hitLists)))
pdf(dirout(outGeneSets, "TimepointZero_Venn_up.pdf"))
venn(hitLists[grepl("_up", names(hitLists))])# | grepl("PT", names(hitLists))])
dev.off()
pdf(dirout(outGeneSets, "TimepointZero_Venn_down.pdf"))
venn(hitLists[grepl("_down", names(hitLists))])#  | grepl("PT", names(hitLists))])
dev.off()

# heatmap
logMat <- -log10(abs(mat))
logMat[logMat >= 5] <- 5
logMat <- sign(mat) * logMat

rows <- order(apply(logMat, 1, sum), decreasing=T)[1:25]
rows <- c(rows, order(apply(logMat, 1, sum), decreasing=F)[1:25])

pdf(dirout(outGeneSets, "TimepointZero.pdf"), width=10, height=16,onefile=F)
pheatmap(logMat[rows,])
dev.off()

# boxplots for individual genesets
dir.create(dirout(outGeneSets, "TimepointZero/"))
for(set.nam in rownames(logMat)[rows]){
  genes <- genesets[[set.nam]]
  genes <- genes[genes %in% rownames(dat)]
  
  pDat$score <- apply(dat[genes,], 2, sum)/colSums  
  pDat2 <- pDat[time == "early"]
  
  ggplot(pDat2, aes(y=score, x=sample)) + geom_jitter(color="grey", alpha=0.5)+ geom_boxplot(fill="NA")  + coord_flip() + 
    ggtitle(paste(set.nam))
  ggsave(dirout(outGeneSets, "TimepointZero/", set.nam, ".pdf"))
}