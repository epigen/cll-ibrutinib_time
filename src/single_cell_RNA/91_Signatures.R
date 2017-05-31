
require(pheatmap)
require(gplots)

outGeneSets <- paste0(outS, "Genesets/", file.nam, "/")

if(exists("balance.barcodes") && balance.barcodes){
  outGeneSets <- paste0(outS, "Genesets_balanced/", file.nam, "/")
  dir.create(dirout(outS, "Genesets_balanced/"))
} else {
  dir.create(dirout(outS, "Genesets/"))
}
dir.create(dirout(outGeneSets))

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


# Plot number of genes and UMIs ----------------------------------------------------
ggplot(pbmc@data.info, aes(y=nUMI, x=sample)) + geom_jitter(color="grey", alpha=0.5)+ geom_boxplot(fill="NA")  + coord_flip()
ggsave(dirout(outGeneSets, "UMIs.pdf"))
ggplot(pbmc@data.info, aes(y=nGene, x=sample)) + geom_jitter(color="grey", alpha=0.5)+ geom_boxplot(fill="NA")  + coord_flip()
ggsave(dirout(outGeneSets, "GENES.pdf"))



source("src/single_cell_RNA/91_Signatures_OverTime.R")
source("src/single_cell_RNA/91_Signatures_T_Zero.R")