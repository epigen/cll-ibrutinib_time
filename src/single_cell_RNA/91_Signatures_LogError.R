require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("cll-time_course")
out <- "11_CellTypes/"
dir.create(dirout(out))


sample.x <- "allDataBest_NoDownSampling_noIGH"
cell <- "Bcells"

outS <- paste0(out, sample.x, "_", cell, "/")

load(file=dirout(outS, cell,".RData"))

balance.barcodes <- FALSE
file.nam <- "hallmark"

require(pheatmap)
require(gplots)

outGeneSets <- paste0(outS, "Genesets2/", file.nam, "/")

file <- geneSetFiles[[file.nam]]


lines <- readLines(file)
genesets <- list()
for(line in lines){
  x <- strsplit(line, "\t")[[1]]
  genesets[[x[1]]] <- x[3:length(x)]
}


# PREPARE SC DATA ---------------------------------------------------------
pDat <- data.table(pbmc@tsne.rot)
pDat$barcode <- colnames(pbmc@data)
pDat$sample <- pbmc@data.info$sample
table(pDat$sample)
pDat$patient <- sapply(strsplit(pDat$sample, "_"), function(x) return(gsub("\\d", "", x[1])))
table(pDat$patient)
pDat$time <- sapply(strsplit(pDat$sample, "_"), function(x) return(gsub("d150", "d120", gsub("(\\d+)d", "d\\1", x[length(x)]))))
table(pDat$time)

dat <- pbmc@data

# LIMIT TO TWO TIMEPOINTS
pDat <- pDat[time != "d280"]
pDat[time == "d0", time := "early"]
pDat[time == "d120", time := "late"]
pDat$time <- factor(pDat$time, levels=c("early", "late"))
pats <- unique(pDat$patient)
pDat$patient <- factor(pDat$patient, levels = c("FE", pats[pats != "FE"]))

samples.small <- pDat[, .N, by = "sample"][N<5]$sample
pDat <- pDat[!sample %in% samples.small]

# balance samples if selected
if(exists("balance.barcodes") && balance.barcodes){
  barcodes.max <- min(table(pDat$sample))
  barcodes.keep <- do.call(c, lapply(split(pDat$barcode, factor(pDat$sample)), function(barcodes) sample(barcodes, barcodes.max)))
  pDat <- pDat[barcode %in% barcodes.keep]
}

# Also limit data to those samples
dat <- dat[,pDat$barcode]




genes <- genesets$HALLMARK_INTERFERON_ALPHA_RESPONSE
genes <- genes[genes %in% rownames(dat)]

colSums <- apply(dat, 2, sum)
pDat$score.corr <- log(apply(exp(dat[genes,]), 2, sum))
pDat$score.err <- apply(dat[genes,], 2, sum)/colSums

ggplot(pDat, aes(x=score.corr, y=score.err)) + geom_point()

ggplot(pDat, aes(x=sample, y=score.err, group=sample, color=time)) + geom_boxplot() + coord_flip()

ggplot(pDat, aes(x=sample, y=score.corr, group=sample, color=time)) + geom_boxplot() + coord_flip()
