require("project.init")
# install.packages("nloptr",repos="https://cloud.r-project.org")
require(Seurat)
require(methods)
project.init2("cll-time_course")


out <- "52_3_Monocle2/"
dir.create(dirout(out))


# install -----------------------------------------------------------------
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("monocle")
require(monocle)

(load(file=dirout("11_CellTypes_tobit/allDataBest_NoDownSampling_noIGH_Bcells/Bcells.RData")))
pbmc <- UpdateSeuratObject(pbmc)
# pbmc.full <- pbmc

# pbmc <- SubsetData(pbmc.full, cells.use=rownames(pbmc.full@meta.data)[pbmc.full@meta.data$sample == "KI_KI1_d0"])

cell.annot <- pbmc@meta.data
cell.annot$barcode <- rownames(pbmc@meta.data)

mcle <- newCellDataSet(
  cellData=pbmc@raw.data[,rownames(pbmc@meta.data)],
  phenoData = AnnotatedDataFrame(cell.annot),
  featureData = AnnotatedDataFrame(data=data.frame(gene_short_name=rownames(pbmc@raw.data), row.names=rownames(pbmc@raw.data))),
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size())

mcle <- detectGenes(mcle, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(mcle), num_cells_expressed >= 10))

mcle <- setOrderingFilter(mcle, intersect(expressed_genes, pbmc@var.genes))
mcle <- estimateSizeFactors(mcle)
mcle <- estimateDispersions(mcle)
plot_ordering_genes(mcle)
ggsave(dirout(out, "Ordered.genes.pdf"), width=7, height=7)
mcle <- reduceDimension(mcle, max_components = 2, method = 'DDRTree')
mcle <- orderCells(mcle)

save(mcle, file=dirout(out, "MonocleDat.RData"))

plot_cell_trajectory(mcle, color_by = "Pseudotime")
ggsave(dirout(out, "Pseudotime.pdf"), width=7, height=7)

plot_cell_trajectory(mcle, color_by = "State")
ggsave(dirout(out, "State.pdf"), width=7, height=7)

plot_cell_trajectory(mcle, color_by = "sample")
ggsave(dirout(out, "Sample.pdf"), width=7, height=7)

plot_cell_trajectory(mcle, color_by = "sample") + facet_grid(sample ~ .)
ggsave(dirout(out, "Sample_grid.pdf"), width=7, height=29)


str(pData(mcle))

table(gsub(".+\\-(\\d+)", "\\1", pData(mcle)$barcode), pData(mcle)$sample)


dim(pData(mcle))
dim(pbmc@meta.data)

stopifnot(all(rownames(pbmc@meta.data) == pData(mcle)$barcode))
pbmc@meta.data$State <- pData(mcle)$State
pbmc@meta.data$Pseudotime <- pData(mcle)$Pseudotime



pDat <- data.table(pbmc@dr$tsne@cell.embeddings, State=pbmc@meta.data$State)
p <- ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color= State)) 
ggsave(dirout(out, "State_Tsne.pdf"), plot=p + geom_point(), height=7, width=7)
ggsave(dirout(out, "State_Tsne_alpha.pdf"), plot=p + geom_point(alpha=0.5), height=7, width=7)


pbmc@ident <- factor(pbmc@meta.data$State)
names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
clusters <- names(table(pbmc@ident))[table(pbmc@ident)>1]
clusters <- clusters[clusters %in% c(1,7,6)]
cl.i <- 1
for(cl.i in clusters){
  # message(cl.i)
  if(!file.exists(dirout(out, "Markers_Cluster",cl.i, ".tsv"))){
    try({
      cluster.markers <- FindMarkers(pbmc,  ident.1 = cl.i, ident.2 = clusters[clusters != cl.i],test.use="tobit", min.pct = 0.25)    
      pdf(dirout(out,"Markers_Cluster",cl.i,".pdf"), height=15, width=15)
      FeaturePlot(pbmc, row.names(cluster.markers)[1:min(nrow(cluster.markers),9)],cols.use = c("grey","blue"))
      dev.off()
      write.table(cluster.markers, dirout(out, "Markers_Cluster",cl.i, ".tsv"), sep="\t", quote=F, row.names=TRUE)
    },silent=TRUE)
  }
}
      
      

state.x <- as.character(pbmc@meta.data$State)
state.x[state.x %in% c("3", "4")] <- "IGNORED"
pbmc@meta.data$StateNoIntermediate <- state.x
state.x[state.x %in% c("2", "5")] <- "IGNORED"
pbmc@meta.data$StateEndpoints <- state.x

(save(pbmc, file=dirout(out, "Bcells_Monocle.RData")))

outS <- paste0(out, "Bcells/")
dir.create(outS)

pbmc@meta.data <- subset(pbmc@meta.data, select=c(State))

source("src/single_cell_RNA/95_Seurat2.R")



# diff.genes <- differentialGeneTest(mcle, fullModelFormulaStr="~sm.ns(Pseudotime)")
# sig.gene.names <- rownames(subset(diff.genes, qval < 0.1))
# plot_pseudotime_heatmap(mcle[sig.gene.names,])


# BEAM_res <- BEAM(mcle, branch_point = 1, cores = 1)
# BEAM_res <- BEAM_res[order(BEAM_res$qval),]
# BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
# plot_genes_branched_heatmap(lung[row.names(subset(BEAM_res, qval < 1e-4)),],
#                             branch_point = 1,
#                             num_clusters = 4,
#                             cores = 1,
#                             use_gene_short_name = T,
#                             show_rownames = T)