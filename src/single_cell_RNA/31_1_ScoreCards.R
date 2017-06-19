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


# Read in diff expr data and calculate scores
cells <- list.files(dirout(inDir))
cells <- cells[grepl(paste0(sample.x, "_"), cells)]
cells <- cells[!grepl(".log",cells)]
cells <- gsub(paste0(sample.x, "_"), "", cells)
scorecards <- data.table()
cell <- "Bcells"
for(cell in cells){
  (load(dirout(inDir, sample.x, "_", cell,"/", "overTime/", "DiffExprRes.RData")))
  for(pat in names(res)){
    for(gs in names(genesets)){
      scorecards <- rbind(scorecards, data.table(
        cellType = cell, 
        patient = pat, 
        pvalue = tryCatch({ t.test(res[[pat]][rn %in% genesets[[gs]]]$avg_diff)$p.value },error=function(e){1}), 
        effectsize = median(res[[pat]][rn %in% genesets[[gs]]]$avg_diff, na.rm=TRUE),
        geneset = gs
        ))
    }
  }
}



# Make one large heatmap
scorecards <- scorecards[cellType %in% c("CD8", "CD4", "Monos", "Bcells")]
scorecards[,pvalue2 := pmin(5, -1*log10(pvalue))]
ggplot(
  scorecards, 
  aes(x=paste0(cellType, "_", patient), y=geneset, color=effectsize, size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "0_Overview.pdf"),height=15, width=15)





# PLOT INDIVIDUAL EXAMPLES ------------------------------------------------
# 
# pat <- "PT"
# set.nam <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
# n <- 100
# for(pat in unique(pDat$patient)){
#   for(set.nam in names(genesets)){
#     genes <- genesets[[set.nam]]
#     dat <- pbmc@data
#     sampleAnnot <- subset(pbmc@data.info, grepl(pat, sample) & cellType %in% c("CD8", "CD4", "Mono", "CLL"))
#     cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(with(sampleAnnot, paste0(sample, cellType)))), function(x) sample(x, min(length(x), n))))]
#     
#     genes <- genes[genes %in% rownames(dat)]
#     dat <- dat[genes, cells]
#     dat <- dat[apply(dat, 1, max) != 0,]
#     dat <- dat - apply(dat, 1, min)
#     dat <- dat / apply(dat,1, max)
#     apply(dat, 1, quantile, na.rm=T)
#     dat <- dat[, order(with(sampleAnnot[cells,],paste0(cellType, sample)))]
#     
#     pdf(dirout(out, set.nam, "_", pat, ".pdf"), height=min(29, nrow(dat) * 0.3 + 3), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
#     pheatmap(dat, cluster_rows=F, cluster_cols=F, annotation_col=pbmc@data.info[,c("sample", "cellType", "nUMI"), drop=F],show_colnames=FALSE)
#     dev.off()
#   }
# }