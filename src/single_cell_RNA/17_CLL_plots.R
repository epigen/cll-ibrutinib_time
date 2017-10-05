require("project.init")
# install.packages("nloptr",repos="https://cloud.r-project.org")
require(Seurat)
require(methods)
project.init2("cll-time_course")


out <- "17_CLL_Plots/"
dir.create(dirout(out))


(inDir <- paste0("13_4_Overtime_","inclDay30","/"))
list.files(dirout(inDir))

source("src/single_cell_RNA/FUNC_Enrichr.R")

cll.enr.path <- fread(dirout(inDir, "CLL_EnrichOR_.tsv"))
cll.enr.path.terms <- cll.enr.path[,sum(abs(log(oddsRatio))), by="category"][order(V1,decreasing=T)][1:20]$category
enrichRes <- cll.enr.path[category %in% cll.enr.path.terms]
enrichRes[,term := paste0(category)]
enrichRes$grp2 <- enrichRes[,gsub("_[up|down]", "",grp)]
enrichRes$direction <- enrichRes[,gsub(".+?0_([up|down])", "\\1",grp)]
if(length(unique(enrichRes$term)) >= 2){
  try({
    orMT <- t(as.matrix(dcast.data.table(enrichRes, grp ~ term, value.var="oddsRatio")[,-"grp",with=F]))
    orMT[is.na(orMT)] <- 1
    hclustObj <- hclust(dist(orMT))
    enrichRes$term <- factor(enrichRes$term, levels=hclustObj$labels[hclustObj$order])
  },silent=T)
}
enrichRes[, mLog10Q := pmin(-log10(qval),4)]
ggplot(enrichRes, 
       aes(x=grp, y=term, size=log10(oddsRatio), color=mLog10Q)) + 
  geom_point() + scale_color_gradient(low="grey", high="red") + theme_bw(12) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(. ~ direction, scales="free") + xlab("") + ylab("")
ggsave(dirout(out, "Enrichr_Pathways.pdf"), width=11, height=8)


cll.enr.tfs <- fread(dirout(inDir, "CLL_EnrichOR_2.tsv"))
cll.enr.tfs.terms <- cll.enr.tfs[,sum(abs(log(oddsRatio))), by="category"][order(V1,decreasing=T)][1:20]$category
enrichRes <- cll.enr.tfs[category %in% cll.enr.tfs.terms]
enrichRes[,term := paste0(category)]
enrichRes$grp2 <- enrichRes[,gsub("_[up|down]", "",grp)]
enrichRes$direction <- enrichRes[,gsub(".+?0_([up|down])", "\\1",grp)]
if(length(unique(enrichRes$term)) >= 2){
  try({
    orMT <- t(as.matrix(dcast.data.table(enrichRes, grp ~ term, value.var="oddsRatio")[,-"grp",with=F]))
    orMT[is.na(orMT)] <- 1
    hclustObj <- hclust(dist(orMT))
    enrichRes$term <- factor(enrichRes$term, levels=hclustObj$labels[hclustObj$order])
  },silent=T)
}
enrichRes[, mLog10Q := pmin(-log10(qval),4)]
ggplot(enrichRes, 
       aes(x=grp, y=term, size=log10(oddsRatio), color=mLog10Q)) + 
  geom_point() + scale_color_gradient(low="grey", high="red") + theme_bw(12) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(. ~ direction, scales="free") + xlab("") + ylab("")
ggsave(dirout(out, "Enrichr_TFs.pdf"), width=11, height=8)



(load(dirout("11_CellTypes_inclDay30/CLL/CLL.RData")))

pDat <- data.table(pbmc@data.info)
pDat <- cbind(pDat, data.table(pbmc@tsne.rot))

ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=patient)) + geom_point(alpha=0.5) + theme_bw(24)
ggsave(dirout(out, "Patient.pdf"), height=7, width=8)

pDat$timepoint <- gsub("^[A-Z]+\\_d(\\d+)$", "\\1", pDat$sample)
pDat[timepoint %in% c("0"), timepoint2 := 1]
pDat[timepoint %in% c("30"), timepoint2 := 2]
pDat[timepoint %in% c("120", "150"), timepoint2 := 3]
pDat[timepoint %in% c("280"), timepoint2 := 3.5]
ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=timepoint2)) + geom_point(alpha=0.5) + theme_bw(24) + scale_color_gradient(low="grey", high="red")
ggsave(dirout(out, "Timepoint.pdf"), height=7, width=8)
