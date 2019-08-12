require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("cll-time_course")

out <- "10_Seurat_raw/inclDay30_noIGHLK_negbinom/"
load(dirout(out, "inclDay30_noIGHLK.RData"))


stopifnot(ncol(pbmc@data) == nrow(pbmc@data.info))
stopifnot(rownames(pbmc@data.info) == colnames(pbmc@data))
str(pbmc@data.info)

pbmc@raw.data <- pbmc@raw.data[rownames(pbmc@data), colnames(pbmc@data)]

pbmc@data.info$CellType <- NA
pbmc@data.info$CellType[pbmc@data.info$ClusterNames_0.5 %in% c(0,1,3,4,7)] <- "CLL"
pbmc@data.info$CellType[pbmc@data.info$ClusterNames_0.5 %in% c(2,5)] <- "CD8"
pbmc@data.info$CellType[pbmc@data.info$ClusterNames_0.5 %in% c(6) & pbmc@data["NKG7",] == 0] <- "CD4"
pbmc@data.info$CellType[pbmc@data.info$ClusterNames_0.5 %in% c(9) & pbmc@data["CD3D",] == 0 & pbmc@data["CD3G",] == 0] <- "NK"
pbmc@data.info$CellType[pbmc@data.info$ClusterNames_0.5 %in% c(8)] <- "Mono"
pbmc@data.info$CellType[pbmc@data.info$ClusterNames_0.5 %in% c(10)] <- "NurseLikeCell"


pDat <- data.table(pbmc@tsne.rot)
pDat$CellType <- pbmc@data.info$CellType
ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=CellType)) + geom_point(alpha=0.3) + theme_bw()
ggsave(dirout(out, "CellTypes_tSNE.pdf"), height=7, width=7)

save(pbmc, file=dirout(outS, sample.x,".RData"))