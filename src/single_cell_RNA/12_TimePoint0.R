require("project.init")
require("Seurat")
require(methods)

project.init2("cll-time_course")

out <- "12_TimePoint0/"
dir.create(dirout(out))

sample.x <- "allDataBest_NoDownSampling_noIGH"

cell <- "Bcells"
inDir <- "11_CellTypes/"

f1 <- list.files(dirout(inDir))
f1 <- f1[grepl(paste0(sample.x, "_"), f1)]
f1 <- f1[!grepl(".log",f1)]
f1 <- gsub(paste0(sample.x, "_"), "", f1)

for(cell in f1){
  message(cell)
  outS <- paste0(out, sample.x, "_", cell, "/")
  dir.create(dirout(outS))
  
  datasetName <- paste0(sample.x, "_", cell)
  
  f <- list.files(dirout(inDir, datasetName, "/Cluster_sample/"))
  
  f <- f[grepl("_d?0d?vs.+_d?0d?\\.tsv", f)]
  
  all.patients <- c()
  
  diffGenes <- list()
  for(f.x in f){
    pats <- strsplit(gsub("\\d", "", gsub("\\_[[:upper:]]+", "", gsub("\\.tsv", "", gsub("_d?0d?", "", gsub("Diff_Cluster", "", f.x))))), "vs")[[1]]
    all.patients <- unique(c(all.patients, pats))
    dat <- fread(dirout(inDir, datasetName, "/Cluster_sample/", f.x))
    diffGenes[[paste0(pats[1], "_vs_", pats[2])]] <- dat[V3>0]$V1
    diffGenes[[paste0(pats[2], "_vs_", pats[1])]] <- dat[V3<0]$V1
  }
  
  responder.yes <- c("KI", "KZ", "PBGY", "VZS")
  responder.no <- c("FE", "PT")
  
  responder.genes <- data.table(melt(diffGenes[names(diffGenes) %in% paste(rep(responder.yes, times=length(responder.no)), rep(responder.no, each=length(responder.yes)), sep="_vs_")]))
  responder.genes <- responder.genes[,.N, by="value"][order(N,decreasing=TRUE)]
  write.table(responder.genes, dirout(outS, "Genes_responder.tsv"), sep="\t", quote=F, row.names=F)
  (responder.genes <- responder.genes[N>3]$value)
  
  responderN.genes <- data.table(melt(diffGenes[names(diffGenes) %in% paste(rep(responder.no, times=length(responder.yes)), rep(responder.yes, each=length(responder.no)), sep="_vs_")]))
  responderN.genes <- responderN.genes[,.N, by="value"][order(N,decreasing=TRUE)]
  write.table(responderN.genes, dirout(outS, "Genes_non_responder.tsv"), sep="\t", quote=F, row.names=F)
  (responderN.genes <- responderN.genes[N>3]$value)
  
  stopifnot(length(intersect(responder.genes, responderN.genes)) == 0)
  
  outG <- paste0(outS, "Genes/")
  dir.create(dirout(outG))
  
  (load(dirout(inDir, datasetName, "/", cell, ".RData")))
  
  for(gene.x in c(responder.genes, responderN.genes)){
    # gene.x <- "CD79A"
    
    pDat <- data.table(Expression=pbmc@data[gene.x,], sample=pbmc@data.info[["sample"]])
    pDat <- pDat[grepl("_d?0d?", sample)]
    pDat$sample2 <- gsub("_d?", "", gsub("\\d", "", gsub("\\_[[:upper:]]+", "", pDat$sample)))
    pDat[sample2 %in% responder.yes, responder := TRUE]
    pDat[sample2 %in% responder.no, responder := FALSE]
    
    grps <- c(
      sapply(diffGenes, function(x) gene.x %in% x)[paste(rep(responder.yes, times=length(responder.no)), rep(responder.no, each=length(responder.yes)), sep="_vs_")],
      sapply(diffGenes, function(x) gene.x %in% x)[paste(rep(responder.no, times=length(responder.yes)), rep(responder.yes, each=length(responder.no)), sep="_vs_")]
    )
    
    #     gene.x %in% responder.genes
    #     gene.x %in% responderN.genes
    
    ggplot(pDat, aes(x=sample, y=Expression, fill=responder)) + geom_violin() + geom_boxplot(fill="white", width=.1, outlier.colour=NA, alpha=0.5) + coord_flip() + 
      ggtitle(paste0(
        ifelse(gene.x %in% responder.genes, "responder", "non-responder"), "\n",
        paste(names(grps[!is.na(grps) & grps]),collapse=" ")
        ))
    ggsave(dirout(outG, gene.x, "_Violin.pdf"),height=7, width=7)
    
    ggplot(pDat[,.(pcts=sum(Expression != 0)/length(Expression)), by="sample2"], aes(x=sample2, y=pcts))+geom_bar(stat="identity")
    ggsave(dirout(outG, gene.x, "_pcts.pdf"))
    
    pdf(dirout(outG, gene.x, "_tsne.pdf"))
    FeaturePlot(pbmc,c(gene.x), cols.use=c("grey", "blue"))
    dev.off()
  }
  
  
  
  
  library(enrichR)
  enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location")
  
  for(x in c("responder", "non_responder")){
    try({
      genes <- if(x == "responder") responder.genes else responderN.genes
      enrichRes <- as.data.table(enrichGeneList(genes,databases = enrichrDBs))
      enrichRes$n <- sapply(strsplit(enrichRes$genes,","), length)
      enrichRes <- enrichRes[n > 3][qval < 0.05]
      enrichRes$category <- gsub("\\_(\\w|\\d){8}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){12}", "", enrichRes$category)
      enrichRes$category <- make.unique(substr(enrichRes$category, 0, 50))
      write.table(enrichRes, file=dirout(outS, "Enrichr_",x,".tsv"), sep="\t", quote=F, row.names=F)
      
      if(nrow(enrichRes) > 0){
        pDat <- enrichRes[order(qval, decreasing=TRUE)][1:min(10, nrow(enrichRes))]
        pDat$category <- factor(pDat$category, levels=pDat$category)
        ggplot(pDat, aes(x=category, y=-log10(qval))) + geom_bar(stat="identity") + coord_flip() + theme_bw(12)
        ggsave(dirout(outS, "Enrichr_",x,".pdf"))
      }
    },silent=TRUE)
  }
}


