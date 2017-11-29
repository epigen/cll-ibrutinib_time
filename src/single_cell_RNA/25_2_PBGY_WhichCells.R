require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("cll-time_course")

load(dirout("10_Seurat_raw/inclDay30_noIGHLK_negbinom/", "inclDay30_noIGHLK.RData"))

pbmcO <- pbmc

load(dirout("25_Patient_CLL/PBGY/PBGY.RData"))

pbmcPBGY <- pbmc

with(pbmcPBGY@data.info, table(sample, CellType))
with(pbmcO@data.info[colnames(pbmcPBGY@data),], table(sample, CellType))

ggplot(pbmcO@data.info[colnames(pbmcPBGY@data),], aes(x=nUMI, group=sample, color=sample)) + stat_ecdf()
