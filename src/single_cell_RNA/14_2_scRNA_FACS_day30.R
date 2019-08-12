require("Seurat")
require("ggplot2")
require("pheatmap")
require("project.init")
project.init2("cll-time_course")
out <- "14_scRNA_FACS_day30/"
dir.create(dirout(out))

(load(dirout("10_Seurat_raw/inclDay30_noIGHLK_negbinom/inclDay30_noIGHLK.RData")))



# sc Data with cell types
pDat <- data.table(pbmc@data.info, keep.rownames=TRUE)
pDat$sample <- gsub("(\\d+)d", "d\\1", pDat$sample)
pDat$patient <- gsub("([A-Z]+)\\d?\\_.+", "\\1", pDat$sample)
pDat$timepoint <- gsub(".+\\_d([0-9]+)", "\\1", pDat$sample)
pDat$timepoint <- gsub("150", "120", pDat$timepoint)
pDat <- pDat[timepoint %in% c("0", "120")]
pDat$sample2 <- paste0(pDat$patient, "_", pDat$timepoint)

pDat$cellType <- pDat$CellType

scCounts <- pDat[,.N, by=c("cellType", "sample2")][!is.na(cellType)]
scCounts$m <- paste0(scCounts$sample2, "_",scCounts$cellType)
scCounts[,cellCount := sum(N), by="sample2"]
scCounts[,count := N/cellCount]

table(pDat$cellType)


# FACS Data ---------------------------------------------------------------
facsDat <- fread("metadata/facs_quantification.csv")
colnames(facsDat) <- make.names(gsub("(\\+|\\_| |\\-)", "", colnames(facsDat)))
facsDat <- facsDat[time %in% c(0,30,120, 150, 280)]
facsDat[,sample2:=paste0(patientid, "_", time)]

facsDat2 <- melt(facsDat[variable=="percentage",c(1,2,8:17)][,-c("patientid", "time"),with=F], id.vars="sample2")
facsDat2 <- facsDat2[!is.na(value)]
facsDat2[,m := paste0(sample2, "_", variable)]
facsDat2[,sum(value), by="sample2"]


table(facsDat2$m)
table(scCounts$m)
facsDat2[,m := gsub("CD14Myeloidcells", "Mono", m)]
facsDat2[,m := gsub("CD4Tcells", "CD4", m)]
facsDat2[,m := gsub("CD8Tcells", "CD8", m)]
facsDat2[,m := gsub("CD56NKcells", "NK", m)]
facsDat2[,m := gsub("CD5CLL", "CLL", m)]

mDat <- merge(scCounts, facsDat2, by="m")

write.table(mDat, dirout(out, "FACS_Data.tsv"), quote=F, row.names=F, sep="\t")


# PLOTS -------------------------------------------------------------------
ggplot(mDat, aes(x=count, y=value/100, color=cellType, shape=sample2.x)) + 
  geom_point(size=5) + ylab("FACS") + xlab("scRNA") + scale_shape_manual(values=c(97:104,107))
ggsave(dirout(out, "FACS_vs_scRNA.pdf"), height=5, width=6)

ggplot(mDat, aes(x=count, y=value/100, color=cellType, shape=sample2.x)) + 
  geom_point(size=5) + ylab("FACS") + xlab("scRNA") + scale_shape_manual(values=c(97:104,107)) +
  scale_x_log10() + scale_y_log10()
ggsave(dirout(out, "FACS_vs_scRNA_log.pdf"), height=5, width=6)

ggplot(mDat, aes(x=count*100, y=value, color=cellType)) + geom_point() + ylab("FACS") + xlab("scRNA") +
  theme_bw(24)
ggsave(dirout(out, "FACS_vs_scRNA_noShapes.pdf"), width=7, height=7)
    
cor(mDat$count, mDat$value, method="spearman")
cor(mDat$count, mDat$value, method="pearson")
