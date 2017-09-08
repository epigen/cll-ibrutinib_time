require("project.init")
require(Seurat)
require(methods)
project.init2("cll-time_course")


out <- "53_GroupCorrelations/"
dir.create(dirout(out))

(load(file=dirout('10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.RData')))

table(pbmc@data.info$sample)
sample.x <- "KI_KI1_d0"

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

