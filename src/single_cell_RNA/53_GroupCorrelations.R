require("project.init")
require(Seurat)
require(methods)
project.init2("cll-time_course")


out <- "53_GroupCorrelations/"
dir.create(dirout(out))


if(!file.exists(dirout(out, "Correlations.RData"))){
  message("Calculating Correlations")
  (load(file=dirout('10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.RData')))
  corrs <- data.table()
  for(sample.x in unique(pbmc@data.info$sample)){
    for(cellType.x in unique(pbmc@data.info$cellType)){
      message(sample.x, " ", cellType.x)
      barcodes.use <- data.table(pbmc@data.info, keep.rownames=TRUE)[sample == sample.x & cellType == cellType.x]$rn
      if(length(barcodes.use) >= 2){
        cMat <- cor(as.matrix(pbmc@data[,barcodes.use]))
        corrs <- rbind(corrs, data.table(
          cor = cMat[upper.tri(cMat)], 
          sample=sample.x,
          cellType = cellType.x))
      }
    }
  }
  
  save(corrs, file=dirout(out, "Correlations.RData"))
} else {
  message("loading")
  load(dirout(out, "Correlations.RData"))
}


corrs2 <- corrs[cellType %in% c("CD4","CD8","CLL","Mono")]

ggplot(corrs2, aes(x=cellType, y=cor)) + geom_boxplot(outlier.shape=NA) + facet_grid(sample ~ .)
ggsave(dirout(out, "AllCorrelations.pdf"), height=29, width =15)

corrMeans <- corrs2[,.(meanCor = mean(cor)), by=c("sample", "cellType")]
corrMeans[,patient := substr(sample,0,2)]
corrMeans[,t := gsub(".+?\\_d?(\\d+)d?", "\\1", sample)]
corrMeans[t == "150", t := "120"]
corrMeans <- corrMeans[t != "280" & patient != "KI"]
corrMeans[,t:= paste0("t", t)]
corrMeansWide <- dcast.data.table(corrMeans, cellType + patient ~ t, value.var="meanCor")
corrMeansWide[,diff := t120 - t0]
ggplot(corrMeansWide, aes(x=cellType, y=diff)) + geom_bar(stat="identity") + facet_grid(patient ~ .)
ggsave(dirout(out, "MeanCorrelations.pdf"), width=10, height=10)