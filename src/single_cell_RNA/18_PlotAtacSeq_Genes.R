require("project.init")
require("Seurat")
require("gplots")
require(methods)
require(pheatmap)
# require(enrichR) #devtools::install_github("definitelysean/enrichR")

project.init2("cll-time_course")

sample.x = "nUMI_Cutoff"
(inDir <- paste0("13_4_Overtime_",sample.x,"/"))


out <- "18_1_PlotSigGeneDetails/"
dir.create(dirout(out))

# Load hits
hits <- fread(dirout(inDir, "SigGenes_overTime.tsv"))
hits[,patient := gsub("30", "030", patient)]
hits[,cellPat := paste0(patient ,"_", cellType)]
hits[, significant := any(qvalue < 0.05) | abs(logFC) > 0.5, by= "gene"]
hits[,logqval := pmin(5, -log(qvalue))]

# Load atac lists
genesets <- list()
f <- list.files(dirout("ATAC-seq_signatures"),pattern = "top200\\.csv")
fx <- f[1]
for(fx in f){
  genesets[[gsub("\\.csv", "", fx)]] <- fread(dirout("ATAC-seq_signatures/", fx), header=F)$V1
}
genesets[["TFs"]] <- fread(dirout(out,"TFs.tsv"),header=F)$V1


genesets <- genesets#[sapply(genesets,length) <= 200]
sumRes <- data.table()
gs <- genesets[[1]]
for(gs.nam in names(genesets)){
  gs <- genesets[[gs.nam]]
  res2 <- hits[gene %in% gs]
  if(length(unique(res2$gene)) >= 2){
    try({
      orMT <- t(as.matrix(dcast.data.table(res2, cellPat ~ gene, value.var="logFC")[,-"cellPat",with=F]))
      orMT[is.na(orMT)] <- 1
      hclustObj <- hclust(dist(orMT))
      res2$gene <- factor(res2$gene, levels=hclustObj$labels[hclustObj$order])
    }, silent=T)
  }
  ggplot(res2[significant == TRUE], aes(x=patient, y=gene, color=logFC, size=logqval)) + geom_point() + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_text(size=12)) + 
    scale_color_gradient2(high="red", mid="white", low="blue") +
    facet_grid(. ~ cellType, scale="free", space="free") + 
    ggtitle(gs.nam)
  ggsave(dirout(out, "ATAC_", gs.nam, ".pdf"), height=15, width=15)
  sumRes <- rbind(sumRes,data.table(res2[significant == TRUE][,.(meanLFC = mean(logFC), cnt = .N), by=c("patient","cellType")], grp = gs.nam))
}
table(sumRes$grp)
sumRes[,grp := gsub(".gene_level.top200", "", gsub("atac-seq\\.", "", grp))]
table(sumRes$grp)

ggplot(sumRes, aes(x=patient, y=grp, color=meanLFC, size=cnt))+ geom_point() + 
   theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_text(size=12)) + 
   scale_color_gradient2(high="red", mid="white", low="blue") +
   facet_grid(. ~ cellType, scale="free", space="free")
ggsave(dirout(out, "Summary.pdf"), height=15, width=15)
