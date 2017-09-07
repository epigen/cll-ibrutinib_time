require("project.init")
# install.packages("nloptr",repos="https://cloud.r-project.org")
require(Seurat)
require(methods)
project.init2("cll-time_course")


out <- "52_2_PseudoTime/"
dir.create(dirout(out))


genes <- c("ITGA4", #CD49D
           "CD200",
           "CD44",
           "ACTB"
           )


(load(dirout("11_CellTypes_tobit/allDataBest_NoDownSampling_noIGH_Bcells/Bcells.RData")))
pbmc.bcells <- pbmc



file.small <- dirout(out, "small.magic.csv")
# if(!file.exists(file.small)){
# LOAD MAGIC DATA ---------------------------------------------------------
magDat <- fread(dirout("15_Magic/","magic.matrix.csv"))
genes.random <- sample(colnames(magDat)[!colnames(magDat) %in% genes], 200)
magicMT <- as.matrix(magDat[,-"V1", with=F])
rownames(magicMT) <- magDat$V1
magicMT <- t(magicMT)

# Match to pbmc data
sample.x <- "allDataBest_NoDownSampling_noIGH"
(load(file=dirout("10_Seurat/", sample.x, "/",sample.x,".RData")))
pbmc.full <- pbmc

stopifnot(all(colnames(data) == rownames(pbmc.full@data.info$sample)))
colnames(magicMT) <- rownames(pbmc.full@data.info)

# Just verify that we got the right sample labels
pDat <- data.table(pbmc.full@data.info, expr = magicMT["CD79A",])
ggplot(pDat, aes(x=paste0(cellType), y=expr)) + geom_violin() + facet_grid(sample ~ .) + ggtitle("CD79A expression")
ggsave(dirout(out, "Magic_SampleMatch_verification.pdf"), height=15, width=15)

write.csv(magicMT[c(genes, genes.random), colnames(magicMT) %in% rownames(pbmc.bcells@data.info)], file=file.small)
# }

magicB <- magicMT[c(genes, genes.random), colnames(magicMT) %in% rownames(pbmc.bcells@data.info)]


# ANALYZE THIS ------------------------------------------------------------

# correlation of ITGA4 and CD200
cor(t(magicB[genes,]))
cor(t(magicB[genes,]), method="spearman")

cMT <- cor(t(magicB))
str(cMT[,genes])
head(rank(cMT[,1])/nrow(cMT))
head(rank(cMT[,2])/nrow(cMT))

ggplot(data.table(t(magicB[genes,])), aes(x=ITGA4, y=CD200, color=ACTB)) + geom_point()
ggplot(data.table(t(magicB[genes,])), aes(x=ITGA4, y=ACTB)) + geom_point()
ggplot(data.table(t(magicB[genes,])), aes(x=CD200, y=ACTB)) + geom_point()

pDat <- data.table(pbmc.full@data.info, data.table(t(magicMT[genes,])))[cellType == "CLL"]
pDat.raw <- data.table(pbmc.full@data.info, data.table(t(as.matrix(pbmc.full@data[genes,]))))[cellType == "CLL"]


sample.x <- "PT_d0"
for(sample.x in unique(pDat$sample)){

  ggplot(pDat[sample == sample.x], aes(x=ITGA4, y=CD200)) + geom_hex() + scale_x_log10() + scale_y_log10() + 
    ggtitle(paste(sample.x, round(cor(pDat[sample == sample.x]$ITGA4, pDat[sample == sample.x]$CD200, method="spearman"), 3)))
  ggsave(dirout(out, sample.x, "_imputed.pdf"), height = 7, width = 7)
  
  ggplot(pDat.raw[sample == sample.x], aes(x=ITGA4, y=CD200)) + geom_hex() + ggtitle(paste(sample.x, "Raw"))
  ggsave(dirout(out, sample.x, "_raw.pdf"), height = 7, width = 7)
}


xDat <- data.table(data.frame(sapply(unique(pDat$sample), function(sample.x) 
  cor(pDat[sample == sample.x]$ITGA4, pDat[sample == sample.x]$CD200, method="spearman"))), keep.rownames=T)
colnames(xDat) <- c("sample", "Correlation")

ggplot(xDat, aes(y=sample,x=Correlation)) + geom_point() + geom_vline(xintercept=0, size=2, color="grey")
ggsave(dirout(out, "Correlations_summary.pdf"))







bcDat <- data.table(pbmc.bcells@data.info, data.table(t(as.matrix(pbmc.bcells@data[genes,]))))
bcDat2 <- data.table(pbmc.bcells@data.info, data.table(t(as.matrix(pbmc.bcells@raw.data[genes,rownames(pbmc.bcells@data.info)]))))


sample.x <- "PT_d0"
for(sample.x in unique(bcDat$sample)){
  ggplot(bcDat[sample == sample.x], aes(x=CD200, y=CD44)) + geom_hex() + ggtitle(paste(sample.x, "Raw"))
  ggsave(dirout(out, sample.x, "_CD200_CD44_raw.pdf"), height = 7, width = 7)
  
  ggplot(bcDat2[sample == sample.x & CD200 != 0 & CD44 != 0], aes(x=CD200, y=CD44)) + geom_point() + ggtitle(paste(sample.x, "Raw"))
  ggsave(dirout(out, sample.x, "_CD200_CD44_raw2.pdf"), height = 7, width = 7)
  
  
  ggplot(bcDat[sample == sample.x], aes(x=ITGA4, y=CD44)) + geom_hex() + ggtitle(paste(sample.x, "Raw"))
  ggsave(dirout(out, sample.x, "_CD40d_CD44_raw.pdf"), height = 7, width = 7)
  
  ggplot(bcDat2[sample == sample.x & ITGA4 != 0 & CD44 != 0], aes(x=ITGA4, y=CD44)) + geom_point() + ggtitle(paste(sample.x, "Raw"))
  ggsave(dirout(out, sample.x, "_CD49d_CD44_raw2.pdf"), height = 7, width = 7)
}