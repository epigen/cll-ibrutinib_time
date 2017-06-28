require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(pheatmap)
require(gplots)

project.init2("cll-time_course")
inDir <- "41_3_CNVs_Sigs_Once_min50/"
out <- "41_3_3_CNV_PCs/" 
dir.create(dirout(out))


load(dirout(inDir, "genesets.RData"))
load(file=dirout(inDir,"Scores.RData"))



Dat1

colnames(Dat1)[!grepl("\\d+\\_\\d+", colnames(Dat1))]

sx <- "VZS_d0"
for(sx in unique(Dat1$sample)){
  cnv.mat <- as.matrix(Dat1[cellType == "CLL" & sample == sx,grepl("\\d+\\_\\d+", colnames(Dat1)), with=F])
  rownames(cnv.mat) <- Dat1[cellType == "CLL" & sample == sx]$rn
  pr <- prcomp(cnv.mat,scale=TRUE)
  ggplot(data.table(pr$x), aes(PC1, PC2, color=PC3)) + geom_point()
  ggsave(dirout(out,"PCs_", sx, ".pdf"))
  tsne <- Rtsne::Rtsne(X=pr$x)
  ggplot(data.table(tsne$Y), aes(V1, V2)) + geom_point()
  ggsave(dirout(out,"TSNEs_", sx, ".pdf"))
}

