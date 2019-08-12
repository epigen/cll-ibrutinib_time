require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

seurat.min.genes <- 200

project.init2("cll-time_course")

cellranger_filtered <- "filtered"

sample.x <- "allDataBest_NoDownSampling/"

data.path <- paste0(getOption("PROCESSED.PROJECT"), "results/cellranger_count/",sample.x,"/outs/",cellranger_filtered,"_gene_bc_matrices_mex/GRCh38/")

pbmc.data <- Read10X(data.path)
pbmc <- new("seurat", raw.data = pbmc.data)
pbmc.prenormalize <- pbmc
pbmc <- Setup(pbmc, min.cells = 3, min.genes = seurat.min.genes, do.logNormalize = T, total.expr = 1e4, project = sample.x)
pbmc.normalize <- pbmc

# Mito genes --------------------------------------------------------------
mito.genes <- grep("^MT-", rownames(pbmc@data), value = T)
percent.mito <- colSums(expm1(pbmc@data[mito.genes, ]))/colSums(expm1(pbmc@data))
pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")


# regress out unwanted sources of variation (UMI normalization), also normalize for mitochondrial content
pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "percent.mito"))
pbmc.regressed <- pbmc

str(pbmc.prenormalize@data)
str(pbmc.normalize@data)
str(pbmc.regressed@data)

str(pbmc.prenormalize@scale.data)
str(pbmc.normalize@scale.data)
str(pbmc.regressed@scale.data)

qplot(x=pbmc.normalize@scale.data[,4], y=pbmc.regressed@scale.data[,4]) + geom_hex()

all(pbmc.normalize@scale.data[,4] == pbmc.regressed@scale.data[,4])
all(pbmc.normalize@data[,4] == pbmc.regressed@data[,4])


str(pbmc@raw.data)
str(pbmc@data)
pbmc@data[1:20, 1:20]
pbmc@raw.data[rownames(pbmc@data)[1:20], colnames(pbmc@data)[1:20]]
pbmc@scale.data[rownames(pbmc@data)[1:20], colnames(pbmc@data)[1:20]]

pbmc2 <- Setup(pbmc.prenormalize, min.cells = 3, min.genes = seurat.min.genes, do.logNormalize = F, total.expr = 1e4, project = sample.x)
str(pbmc2@scale.data)

pbmc2@data[rownames(pbmc@data)[1:20], colnames(pbmc@data)[1:20]]
pbmc2@raw.data[rownames(pbmc@data)[1:20], colnames(pbmc@data)[1:20]]

log(1/10000)

qplot(x=pbmc@raw.data[rownames(pbmc@data),colnames(pbmc@data)[1]], y=pbmc@data[,1]) + geom_point() + 
  geom_line(data=data.table(x=1:500, y=log(1:500*1e4/sum(pbmc@raw.data[rownames(pbmc@data),colnames(pbmc@data)[1]]))), aes(x=x, y=y))

lognorm <- pbmc@raw.data[rownames(pbmc@data),colnames(pbmc@data)[1]]
lognorm2 <- log(lognorm*1e4/sum(lognorm))
qplot(x=lognorm2, y=pbmc@data[,1]) + geom_hex() + geom_abline(intercept=0, slope=1)
