require("Seurat")
require("ggplot2")
require("pheatmap")
require("project.init")
project.init2("cll-time_course")
out <- "14_scRNA_FACS/"
dir.create(dirout(out))


sample.x <- "allDataBest_NoDownSampling_noIGH"
load(file=dirout("10_Seurat/", sample.x, "/",sample.x,".RData"))

# sc Data with cell types
pDat <- data.table(pbmc@data.info, keep.rownames=TRUE)
pDat <- data.table(pDat, pbmc@tsne.rot)
pDat[ClusterNames_0.95 %in% c(6,11), cellType := "Mono"]
pDat[ClusterNames_0.95 %in% c(16), cellType := "NurseLikeCells"]
pDat[ClusterNames_0.95 %in% c(10), cellType := "NK"]
pDat[ClusterNames_0.95 %in% c(14), cellType := "GDT"]
pDat[ClusterNames_0.95 %in% c(0,3,4,12,5,7,9,15,18), cellType := "CLL"]
pDat[ClusterNames_0.5 %in% c(6), cellType := "CD4"]
pDat[ClusterNames_0.5 %in% c(2,3,9) & !ClusterNames_0.95 %in% c(14,19,10), cellType := "CD8"]
# Plot this to check
ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=cellType)) + geom_point()
ggsave(dirout(out, "Celltypes_tSNE.pdf"), width=7, height=7)
# save this in the big file
# pbmc@data.info$cellType <- pDat$cellType

pDat$sample <- gsub("(\\d+)d", "d\\1", pDat$sample)
pDat$patient <- gsub("([A-Z]+)\\d?\\_.+", "\\1", pDat$sample)
pDat$timepoint <- gsub(".+\\_d([0-9]+)", "\\1", pDat$sample)
pDat$timepoint <- gsub("150", "120", pDat$timepoint)
pDat <- pDat[timepoint %in% c("0", "120")]
pDat$sample2 <- paste0(pDat$patient, "_", pDat$timepoint)

scCounts <- pDat[,.N, by=c("cellType", "sample2")][!is.na(cellType)]
scCounts$m <- paste0(scCounts$sample2, "_",scCounts$cellType)
scCounts[,cellCount := sum(N), by="sample2"]
scCounts[,count := N/cellCount]

table(pDat$cellType)


# FACS Data ---------------------------------------------------------------
facsDat <- fread("metadata/facs_quantification.csv")
colnames(facsDat) <- make.names(gsub("(\\+|\\_| |\\-)", "", colnames(facsDat)))
facsDat <- facsDat[timepoint %in% c(0,120, 150)]
facsDat[,sample2:=paste0(patientid, "_", timepoint)]

facsDat2 <- melt(facsDat[,-c("patientid", "timepoint"),with=F], id.vars="sample2")
facsDat2[,m := paste0(sample2, "_", variable)]
facsDat2[,sum(value), by="sample2"]

mDat <- merge(scCounts, facsDat2, by="m")




# PLOTS -------------------------------------------------------------------
ggplot(mDat, aes(x=count, y=value/100, color=cellType, shape=sample2.x)) + 
  geom_point(size=5) + ylab("FACS") + xlab("scRNA") + scale_shape_manual(values=c(97:104,107))
ggsave(dirout(out, "FACS_vs_scRNA.pdf"), height=5, width=6)

ggplot(mDat, aes(x=count, y=value/100, color=cellType, shape=sample2.x)) + 
  geom_point(size=5) + ylab("FACS") + xlab("scRNA") + scale_shape_manual(values=c(97:104,107)) +
  scale_x_log10() + scale_y_log10()
ggsave(dirout(out, "FACS_vs_scRNA_log.pdf"), height=5, width=6)

ggplot(mDat, aes(x=count, y=value, color=cellType)) + geom_point() + ylab("FACS") + xlab("scRNA")
ggsave(dirout(out, "FACS_vs_scRNA_noShapes.pdf"))

