require("project.init")
# install.packages("nloptr",repos="https://cloud.r-project.org")
require(Seurat)
require(methods)
project.init2("cll-time_course")


out <- "52_PseudoTime/"
dir.create(dirout(out))

(load(dirout("11_CellTypes_tobit/allDataBest_NoDownSampling_noIGH_Bcells/Bcells.RData")))

genes <- c("ITGA4", #CD49D
           "CD200",
           "CD44"
           )

g <- genes[1]

cor(t(as.matrix(pbmc@data[genes,])))

for(g in genes){
  ggplot(data.table(pbmc@tsne.rot, gene = pbmc@data[g,]), aes(x=tSNE_1,y=tSNE_2, color=gene)) + geom_point()+
    scale_color_gradientn(colors=c("grey", "blue")) + theme(legend.position = 'none') + ggtitle(g)
  ggsave(dirout(out, "Gene_", g, ".pdf"))
  
  pDat <- data.table(sample=pbmc@data.info$sample, expr=pbmc@data[g,])
  pDat[expr == 0, expr := NA]
  ggplot(pDat, aes(x=sample, y=expr)) + geom_violin()  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(g)
  ggsave(dirout(out, "Gene_violin_", g, ".pdf"))
}



