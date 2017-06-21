require(project.init)
project.init2("cll-time_course")
require("gplots")
require(methods)
require(pheatmap)


project.init2("cll-time_course")

out <- "40_CNVs/"
dir.create(dirout(out))

sample.x <- "allDataBest_NoDownSampling_noIGH"
inDir <- "11_CellTypes_tobit/"
load(file=dirout("10_Seurat/", sample.x, "/",sample.x,".RData"))

# Read genesets
genesets <- list()
(gene.annot <- fread(paste0(Sys.getenv("RESOURCES"), "/genomes/hg38/10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf")))
table(gene.annot$V3)
a2 <- gene.annot[V3 == "gene"]
a2[,gene := gsub('.+\\; gene_name \\"(.+?)\\"\\; .+', "\\1", V9)]
length(unique(a2$gene))
geneOrder <- a2[,c("V1", "V4", "gene"),with=F]
colnames(geneOrder) <- c("chr", "start", "gene")
geneOrder <- geneOrder[chr %in% c(as.character(1:22), "X", "Y")]

geneOrder2 <- geneOrder[gene %in% rownames(pbmc@data)]
slide.stepsize <- 50
chr.x <- 1
i2 <- 0
for(chr.x in unique(geneOrder2$chr)){
  geneOrder.chr <- geneOrder2[chr == chr.x]
  for(i in 1:ceiling((nrow(geneOrder.chr))/slide.stepsize)){
    i2 <- i2 + 1    
    genesets[[paste0(chr.x, "_", i2)]] <- geneOrder.chr[((i-1)*slide.stepsize+1):min(i*slide.stepsize, nrow(geneOrder.chr))]$gene
  }
}



# Read in diff expr data and calculate scores
cells <- list.files(dirout(inDir))
cells <- cells[grepl(paste0(sample.x, "_"), cells)]
cells <- cells[!grepl(".log",cells)]
cells <- gsub(paste0(sample.x, "_"), "", cells)
scorecards <- data.table()
fullList <- data.table()
cell <- "Monos"
for(cell in cells){
  (load(dirout(inDir, sample.x, "_", cell,"/", "overTime/", "DiffExprRes.RData")))
  for(pat in names(res)){
    message(pat, " ", cell)
    for(gs in names(genesets)){
      scorecards <- rbind(scorecards, data.table(
        cellType = cell, 
        patient = pat, 
        pvalue = tryCatch({ t.test(res[[pat]][rn %in% genesets[[gs]]]$avg_diff)$p.value },error=function(e){1}), 
        effectsize = -1 * median(res[[pat]][rn %in% genesets[[gs]]]$avg_diff, na.rm=TRUE),
        geneset = gs
      ))
      try({
        fullList <- rbind(fullList, data.table(
          cellType = cell, 
          patient = pat, 
          geneset = gs,
          res[[pat]][rn %in% genesets[[gs]]]
        ))}, silent=TRUE)
    }
  }
}

scorecards[cellType == "Bcells", cellType := "CLL"]
scorecards[abs(effectsize) > 0.25, effectsize := 0.25 * sign(effectsize)]

# Make one large heatmap
scorecards <- scorecards[cellType %in% c("CD8", "CD4", "Monos", "CLL")]
scorecards[,pvalue2 := pmin(5, -1*log10(pvalue))]
ggplot(
  scorecards, 
  aes(x=paste0(cellType, "_", patient), y=geneset, color=effectsize, size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw(10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "0_Overview.pdf"),height=29, width=8)



# Plot values for all genesets
fullList[,avg_diff := - avg_diff]
gs <- "CD8_down"
fullList[, logp := pmin(-log10(p_val), 5)]
for(gs in names(genesets)){
  ggplot(fullList[geneset == gs], aes(x=paste0(cellType, "_", patient), group=paste0(cellType, "_", patient), y=avg_diff, color=logp)) + 
    geom_hline(yintercept=0,color="grey") + geom_boxplot(outlier.shape=NA, color="lightgrey") + geom_point(alpha = 0.5) + coord_flip() + ggtitle(gs)
  ggsave(dirout(out, "Details_", gs, ".pdf"))
}
