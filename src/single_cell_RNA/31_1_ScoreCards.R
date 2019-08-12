require(project.init)
project.init2("cll-time_course")
require("gplots")
require(methods)
require(pheatmap)


project.init2("cll-time_course")

out <- "31_ScoreCards_tobit/"
dir.create(dirout(out))

sample.x <- "allDataBest_NoDownSampling_noIGH"
inDir <- "11_CellTypes_tobit/"


# Read genesets
file <- "/data/groups/lab_bock/shared/resources/gene_sets/MSigDB/6.0/Human/h.all.v6.0.symbols.gmt"
lines <- readLines(file)
genesets <- list()
for(line in lines){
  x <- strsplit(line, "\t")[[1]]
  genesets[[x[1]]] <- x[3:length(x)]
}
load(dirout("20_AggregatedLists/lists.RData"))
genesets <- c(genesets, cll.lists)

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
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw(24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "0_Overview.pdf"),height=20, width=20)



# Plot values for all genesets
fullList[,avg_diff := - avg_diff]
gs <- "CD8_down"
fullList[, logp := pmin(-log10(p_val), 5)]
for(gs in names(genesets)){
  ggplot(fullList[geneset == gs], aes(x=paste0(cellType, "_", patient), group=paste0(cellType, "_", patient), y=avg_diff, color=logp)) + 
    geom_hline(yintercept=0,color="grey") + geom_boxplot(outlier.shape=NA, color="lightgrey") + geom_point(alpha = 0.5) + coord_flip() + ggtitle(gs)
  ggsave(dirout(out, "Details_", gs, ".pdf"))
}