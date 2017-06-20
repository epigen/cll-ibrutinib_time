require("project.init")
require(Seurat)
require(methods)
project.init2("cll-time_course")

load(dirout("10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.RData"))


# get exponent of log values
str(x <- pbmc@data[1:20, 1:20])
x@x <- exp(x@x)
x

# from the raw data
colSums <- apply(pbmc@raw.data[rownames(pbmc@data),colnames(pbmc@data)[1:20]], 2, sum)
# only the same genes / cells
rr <- pbmc@raw.data[rownames(pbmc@data)[1:20],colnames(pbmc@data)[1:20]]
x
(1/colSums * 1e4)
(1/colSums * 1e4) + 1
log((1/colSums * 1e4) + 1)
pbmc@data[1:20, 1:20]


# ==> SEURAT ADDS +1 to the values!