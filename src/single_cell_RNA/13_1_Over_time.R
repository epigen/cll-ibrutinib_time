require("project.init")
require("Seurat")
require(methods)
require(pheatmap)

project.init2("cll-time_course")

out <- "13_Over_time/"
dir.create(dirout(out))

sample.x <- "allDataBest_NoDownSampling_noIGH"

inDir <- "11_CellTypes/"

cells <- list.files(dirout(inDir))
cells <- cells[grepl(paste0(sample.x, "_"), cells)]
cells <- cells[!grepl(".log",cells)]
cells <- gsub(paste0(sample.x, "_"), "", cells)
cells <- cells[cells != "TcellsAll"]
cells <- cells[cells != "CD8"]


outS <- paste0(out, sample.x, "/")
dir.create(dirout(outS))


cell <- "Bcells"
for(cell in cells){
  message(cell)  
  datasetName <- paste0(sample.x, "_", cell)
  
  pat.dirs <- list.dirs(dirout(inDir, datasetName))
  pat.dirs <- pat.dirs[grepl("_pat_", pat.dirs)]
  
  res <- data.table()
  eff.sizes <- list()
  
  for(pat.dir in pat.dirs){
    comp.files <- list.files(paste0(pat.dir, "/"))
    comp.files <- comp.files[grepl("Diff_Cluster.+\\.tsv", comp.files)]
    comp.files <- comp.files[!grepl("d280", comp.files)]
    
    for(comp.file in comp.files){
      nams <- strsplit((gsub("Diff_Cluster", "", gsub("\\.tsv", "", comp.file))), "vs")[[1]]
      hits <- fread(paste0(pat.dir, "/",comp.file))[order(V3)]
      # get effect sizes into list
      eff <- hits$V3
      names(eff) <- hits$V1
      eff.sizes[[paste0(nams[1], "_vs_", nams[2])]] <- eff
      # get significant genes into list      
      hits <- hits[V2 < 0.05 & abs(V3) > 0.3][order(V3)]
      res <- rbind(res, data.table(hits[V3 > 0]$V1, paste0(nams[1], "_vs_", nams[2])))
      res <- rbind(res, data.table(hits[V3 < 0]$V1, paste0(nams[2], "_vs_", nams[1])))
    }
  }
  names(res) <- c("gene", "comparison")
  res$x <- 1
  res$comparison2 <- sapply(strsplit(res$comparison, "_vs_"), function(vec){
    pat <- gsub("([A-Z]+)\\d?_.+", "\\1", vec[1])
    direction <- ifelse(grepl("d0", gsub("(\\d+)d", "d\\1", vec)[1]), "down", "up")
    return(paste0(pat, "_", direction))
  })
  table(res$comparison, res$comparison2)
  
  # how many hits per cell
  ggplot(res, aes(x=comparison2)) + geom_bar() + coord_flip() + ylab("Significant Genes")
  ggsave(dirout(outS, cell, "_Counts.pdf"))
  
  # Heatmaps
  pDT <- dcast.data.table(res, gene ~ comparison2, value.var="x")
  pMT <- as.matrix(pDT[,-"gene", with=F])
  row.names(pMT) <- pDT$gene
  pMT[is.na(pMT)] <- 0
  cnt <- apply(pMT, 1, sum, na.rm=TRUE)
  pDT$cnt <- cnt
  write.table(pDT[order(cnt, decreasing = TRUE)], file=dirout(outS, cell, "_table.tsv"), quote=F, row.names=F, sep="\t")
  
  
  for(i in 2:4){
    if(nrow(pMT[cnt>=i,,drop=F]) > 2){
      pdf(dirout(outS, cell, "_HM_",i,".pdf"), onefile=FALSE,width=min(29, nrow(pMT[cnt>=i,])*0.3 + 3), height=min(29,ncol(pMT[cnt>=i,])*0.3+2))
      pheatmap(t(pMT[cnt>=i,]))
      dev.off()
    }
  }

  
  all.genes <- unique(do.call(c, lapply(eff.sizes, names)))
  try({
    pdf(dirout(outS, cell, "_CorrHM.pdf"), width=5, height=5, onefile=FALSE)
    pheatmap(cor(do.call(cbind, lapply(eff.sizes, function(x) return(x[all.genes]))), use="pairwise.complete.obs"))
    dev.off()
  }, silent=TRUE)
  
  
  # Venn diagrams
  geneLists <- split(res$gene, factor(res$comparison2))
  try({
    pdf(dirout(outS, cell, "_venn_up.pdf"))
    venn(geneLists[grepl("_up", names(geneLists))])
    dev.off()
  }, silent=T)
  
  try({
    pdf(dirout(outS, cell, "_venn_down.pdf"))
    venn(geneLists[grepl("_down", names(geneLists))])
    dev.off()
  }, silent=T)

}