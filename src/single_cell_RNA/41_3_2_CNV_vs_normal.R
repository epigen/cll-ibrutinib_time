require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(pheatmap)
require(gplots)

project.init2("cll-time_course")
inDir <- "41_3_CNVs_Sigs_Once_min50/"
out <- "41_3_2_CNV_vs_normal/" 
dir.create(dirout(out))


load(dirout(inDir, "genesets.RData"))
load(file=dirout(inDir,"Scores.RData"))



# ANALYSIS T_ZERO ----------------------------------------------------------------
pDat2 <- Dat1
# pDat2 <- pDat2[timepoint %in% c("d0")]

# Make one large heatmap
tzero.sig <- data.table()
pat="FE"
geneset <- "1_1"
for(pat in unique(pDat2$sample)){
  for(geneset in names(genesets)){
    p <- t.test(pDat2[sample == pat & cellType == "CLL"][[geneset]], pDat2[sample == pat & cellType != "CLL"][[geneset]])$p.value
    ef <- mean(pDat2[sample == pat & cellType == "CLL"][[geneset]])- mean(pDat2[sample == pat & cellType != "CLL"][[geneset]])
    tzero.sig <- rbind(tzero.sig, data.table(sample=pat, pvalue=p, meanDiff=ef, geneset=geneset))
  }
}

tzero.sig[,pvalue2 := pmin(5, -1*log10(pvalue))]
save(tzero.sig, file=dirout(out, "Tzero.sig.RData"))
ggplot(
  tzero.sig[geneset %in% tzero.sig[,min(pvalue), by="geneset"][order(V1)][1:50]$geneset], 
  aes(x=paste0(sample), y=geneset, color=meanDiff, size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "0_Overview_TZero.pdf"),height=15, width=15)
