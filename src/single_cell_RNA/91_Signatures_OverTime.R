# CALCULATE SCORES AND FIT MODEL --------------------------------------------------------
lm.pvalues <- data.table(geneset=names(genesets))
lm.effect <- data.table(geneset=names(genesets))
colSums <- apply(dat, 2, sum)

stopifnot(all(apply(dat, 2, sum)/colSums==1))

if(!file.exists(dirout(outGeneSets, "Pvalues.csv"))){
  for(set.nam in names(genesets)){
    genes <- genesets[[set.nam]]
    genes <- genes[genes %in% rownames(dat)]
    
    pDat$score <- apply(dat[genes,], 2, sum)/colSums  
    lm.coef <- summary(lm(data=pDat, score ~ patient / time))$coefficients
    lm.coef <- data.table(lm.coef, keep.rownames=TRUE)[rn != "(Intercept)"]
    
    for(rn.x in lm.coef$rn){
      lm.pvalues[geneset==set.nam, eval(rn.x) := lm.coef[rn.x == rn]$"Pr(>|t|)"]
      lm.effect[geneset==set.nam, eval(rn.x) := lm.coef[rn.x == rn]$Estimate]
    }
  }
  
  write.csv(lm.pvalues, file=dirout(outGeneSets, "Pvalues.csv"), row.names=F)
  write.csv(lm.effect, file=dirout(outGeneSets, "EffectSize.csv"), row.names=F)
  
} else {
  lm.pvalues <- fread(dirout(outGeneSets, "Pvalues.csv"))
  lm.effect <- fread(dirout(outGeneSets, "EffectSize.csv"))
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
  pdf(dirout(outGeneSets, "Overtime_Venn_up.pdf"))
  venn(hitLists[grepl("_up", names(hitLists))])# | grepl("PT", names(hitLists))])
  dev.off()
}, silent=TRUE)
try({
  pdf(dirout(outGeneSets, "Overtime_Venn_down.pdf"))
  venn(hitLists[grepl("_down", names(hitLists))])#  | grepl("PT", names(hitLists))])
  dev.off()
}, silent=TRUE)
ggplot(data.table(list=names(hitLists), size=sapply(hitLists, length)), aes(x=list, y=size))+geom_bar(stat="identity") + coord_flip()
ggsave(dirout(outGeneSets, "Overtime_Listsize.pdf"), width=9, height=7)


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
  
  pdf(dirout(outGeneSets, "Overtime.pdf"), width=10, height=min(16, nrow(mat2)*0.25+3),onefile=F)
  pheatmap(mat2)
  dev.off()
  
  # boxplots for individual genesets
  dir.create(dirout(outGeneSets, "OverTime/"))
  for(set.nam in rownames(mat)){
    genes <- genesets[[set.nam]]
    genes <- genes[genes %in% rownames(dat)]
    
    pDat$score <- apply(dat[genes,], 2, sum)/colSums  
    pDat2 <- pDat[patient != "KI"]
    
    ggplot(pDat2, aes(y=score, x=sample)) + geom_jitter(color="grey", alpha=0.5)+ geom_boxplot(fill="NA")  + coord_flip() + 
      ggtitle(paste(set.nam))
    ggsave(dirout(outGeneSets, "OverTime/", set.nam, ".pdf"))
  }
}

# Heatmap of correlations between samples
try({
  pdf(dirout(outGeneSets, "Overtime_CorrHM.pdf"), width=5, height=5,onefile=F)
  pheatmap(cor(mat.eff[,grepl(":", colnames(mat.eff))], method="pearson"))
  dev.off()
}, silent=TRUE)

