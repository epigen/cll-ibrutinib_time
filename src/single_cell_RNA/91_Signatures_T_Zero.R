
# CALCULATE SCORES AND FIT MODEL --------------------------------------------------------
lm.pvalues <- data.table(geneset=names(genesets))
lm.effect <- data.table(geneset=names(genesets))
# colSums <- apply(dat, 2, sum)
# stopifnot(all(apply(dat, 2, sum)/colSums==1))

patients <- unique(as.character(pDat[time == "early"]$patient))

file1 <- dirout(outGeneSets, "Pvalues_t0.csv")
file2 <- dirout(outGeneSets, "EffectSize_t0.csv")
if(!file.exists(file1)){
  for(set.nam in names(genesets)){
    genes <- genesets[[set.nam]]
    genes <- genes[genes %in% rownames(dat)]
    
    pDat$score <- log(apply(exp(dat[genes,]), 2, sum))
    
    for(pi1 in 1:(length(patients)-1)){
      for(pi2 in (pi1+1):length(patients)){
        pat1 <- patients[pi1]
        pat2 <- patients[pi2]
        lm.coef <- summary(lm(data=pDat[patient %in% c(pat1, pat2) & time == "early"], score ~ patient))$coefficients
        lm.coef <- data.table(lm.coef, keep.rownames=TRUE)[rn != "(Intercept)"]
        
        for(rn.x in lm.coef$rn){
          lm.pvalues[geneset==set.nam, eval(paste0(pat1, "_vs_", pat2)) := lm.coef[1]$"Pr(>|t|)"]
          lm.effect[geneset==set.nam, eval(paste0(pat1, "_vs_", pat2)) := lm.coef[1]$Estimate]
        }
      }
    }
  }
  
  write.csv(lm.pvalues, file=file1, row.names=F)
  write.csv(lm.effect, file=file2, row.names=F)
  
} else {
  lm.pvalues <- fread(file1)
  lm.effect <- fread(file2)
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
# try({
#   pdf(dirout(outGeneSets, "TimepointZero_Venn_up.pdf"))
#   venn(hitLists[grepl("_up", names(hitLists))])# | grepl("PT", names(hitLists))])
#   dev.off()
# }, silent=TRUE)
# try({
#   pdf(dirout(outGeneSets, "TimepointZero_Venn_down.pdf"))
#   venn(hitLists[grepl("_down", names(hitLists))])#  | grepl("PT", names(hitLists))])
#   dev.off()
# }, silent=TRUE)
ggplot(data.table(list=names(hitLists), size=sapply(hitLists, length)), aes(x=list, y=size))+geom_bar(stat="identity") + coord_flip()
ggsave(dirout(outGeneSets, "TimepointZero_Listsize.pdf"), width=9, height=7)

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
    ifelse(sum(mat2 < 0) > 0, max(abs(mat2[mat2 < 0])), NA), # smallest neg value
    ifelse(sum(mat2 > 0) > 0, max(mat2[mat2 > 0]), NA), # largest positive value
    na.rm=T)
  mat2[mat2 > cutval] <- cutval
  mat2[mat2 < -cutval] <- -cutval
  
  pdf(dirout(outGeneSets, "TimepointZero.pdf"), width=10, height=min(16, nrow(mat2)*0.25+3),onefile=F)
  pheatmap(mat2)
  dev.off()
  
  # boxplots for individual genesets
  dir.create(dirout(outGeneSets, "TimepointZero/"))
  for(set.nam in rownames(mat)){
    genes <- genesets[[set.nam]]
    genes <- genes[genes %in% rownames(dat)]
    
    pDat$score <- log(apply(exp(dat[genes,]), 2, sum))
    pDat2 <- pDat[time == "early"]
    
    ggplot(pDat2, aes(y=score, x=sample)) + geom_jitter(color="grey", alpha=0.5)+ geom_boxplot(fill="NA")  + coord_flip() + 
      ggtitle(paste(set.nam))
    ggsave(dirout(outGeneSets, "TimepointZero/", set.nam, ".pdf"))
  }
}
# Heatmap of correlations between samples
try({
  pdf(dirout(outGeneSets, "TimepointZero_CorrHM.pdf"), width=5, height=5,onefile=F)
  pheatmap(cor(mat.eff[,!grepl(":", colnames(mat.eff)) & colnames(mat.eff) != "geneset"], method="pearson"))
  dev.off()
}, silent=TRUE)