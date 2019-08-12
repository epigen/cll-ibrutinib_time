require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(pheatmap)
require(gplots)

project.init2("cll-time_course")
out <- "16_PlotSomeGenes/"
dir.create(dirout(out))


sample.x <- "inclDay30_noIGHLK"
load(file=dirout("10_Seurat_raw/", sample.x, "_negbinom/",sample.x,".RData"))

pbmc@data.info$sample <- gsub("d30", "d030", gsub("d0", "d000", pbmc@data.info$sample))
pbmc@data.info$patient <- gsub("([A-Z]+)\\d?_.+", "\\1", pbmc@data.info$sample)

pDat <- data.table(pbmc@data.info, keep.rownames=TRUE)

genes <- c("TNF", "NFKBIA", "IL1B", "CXCL8", "CCL3L3", "NFKBIZ", "CCL3", "CXCR4", "CCL4")

g <- "TNF"
for(g in genes){
  pDat[[g]] <- pbmc@data[g,]
  ct <- "Mono"
  for(ct in c("CD8", "CD4", "Mono", "CLL")){
    pDat2 <- pDat[CellType == ct]
    ggplot(pDat2, aes_string(x='sample', y=g)) + geom_violin(aes_string(fill="patient"), color=NA) + theme_bw() + 
      theme(axis.text.x=element_text(angle=90, hjust=1)) + ggtitle(paste(g, ct)) 
    ggsave(dirout(out, ct, "_", g, ".pdf"), width=7, height=7)
  }
}