require("project.init")
# install.packages("nloptr",repos="https://cloud.r-project.org")
require(Seurat)
require(methods)
project.init2("cll-time_course")


out <- "52_3_Monocle3/"
dir.create(dirout(out))


# install -----------------------------------------------------------------
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("monocle")

# if(!file.exists(dirout(out, "Bcells_Monocle.RData"))){
require(monocle)

(load(file=dirout("11_CellTypes_tobit/allDataBest_NoDownSampling_noIGH_Bcells/Bcells.RData")))
pbmc <- UpdateSeuratObject(pbmc)
# pbmc.full <- pbmc
# 
# pbmc <- SubsetData(pbmc.full, cells.use=rownames(pbmc.full@meta.data)[pbmc.full@meta.data$sample == "KI_KI1_d0"])



# DO MONOCLE ANALYSES -----------------------------------------------------

if(!file.exists(dirout(out, "MonocleDat.RData"))){
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
} else {
  load(dirout(out, "MonocleDat.RData"))
}



# PLOT RESULTS ------------------------------------------------------------

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



# ANNOTATE SEURAT OBJECT --------------------------------------------------

stopifnot(all(rownames(pbmc@meta.data) == pData(mcle)$barcode))
pbmc@meta.data$State <- pData(mcle)$State
pbmc@meta.data$Pseudotime <- pData(mcle)$Pseudotime

pDat <- data.table(pbmc@dr$tsne@cell.embeddings, State=pbmc@meta.data$State)
p <- ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color= State)) 
ggsave(dirout(out, "State_Tsne.pdf"), plot=p + geom_point(), height=7, width=7)
ggsave(dirout(out, "State_Tsne_alpha.pdf"), plot=p + geom_point(alpha=0.5), height=7, width=7)




# FIND DIFFERENTIAL GENES -------------------------------------------------


# SEURAT WAY

# pbmc@ident <- factor(pbmc@meta.data$State)
# names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
# clusters <- names(table(pbmc@ident))[table(pbmc@ident)>1]
# clusters <- clusters[clusters %in% c(1,7,6)]
# cl.i <- 1
# for(cl.i in clusters){
#   # message(cl.i)
#   if(!file.exists(dirout(out, "Markers_Cluster",cl.i, ".tsv"))){
#     try({
#       cluster.markers <- FindMarkers(pbmc,  ident.1 = cl.i, ident.2 = clusters[clusters != cl.i],test.use="tobit", min.pct = 0.25)    
#       pdf(dirout(out,"Markers_Cluster",cl.i,".pdf"), height=15, width=15)
#       FeaturePlot(pbmc, row.names(cluster.markers)[1:min(nrow(cluster.markers),9)],cols.use = c("grey","blue"))
#       dev.off()
#       write.table(cluster.markers, dirout(out, "Markers_Cluster",cl.i, ".tsv"), sep="\t", quote=F, row.names=TRUE)
#     },silent=TRUE)
#   }
# }

# MONOCLE WAY
diff_test_res <- NA
if(!file.exists(dirout(out, "Monocle_State_DiffGenes.tsv"))){
  cds_subset <- mcle[,subset(pData(mcle), State %in% c("1", "6", "7"))$barcode]
  diff_test_res <- differentialGeneTest(cds_subset,
                                        fullModelFormulaStr = "~State + sample",
                                        reducedModelFormulaStr = "~sample",
                                        cores=12)

  diff_test_res[,c("gene_short_name", "pval", "qval", "num_cells_expressed")]
  # plot_genes_jitter(cds_subset,
  #                   grouping = "State", color_by = "sample2", plot_trend = TRUE) +
  #   facet_wrap( ~ feature_label, scales= "free_y")
  write.table(diff_test_res[,c("gene_short_name", "pval", "qval", "num_cells_expressed")], file=dirout(out, "Monocle_State_DiffGenes.tsv"), sep="\t", quote=F, row.names=F)
} else {
  diff_test_res <- fread(dirout(out, "Monocle_State_DiffGenes.tsv"))
}

str(genes1 <- unique(diff_test_res[qval < 0.05]$gene_short_name))
gene.lists <- list(
  All = genes1,
  noRP = genes1[!grepl("^RP.*?-.*", genes1) & !grepl("^RP[LS]\\d+[AX]?$", genes1) & !grepl("^MT-", genes1)]
  )
for(gnam in names(gene.lists)){
  str(genes <- gene.lists[[gnam]])
  dat <- pbmc@data
  sampleAnnot <- subset(pbmc@meta.data, State %in% c(1,6,7), select=c(State, sample))
  n <- 100
  cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(with(sampleAnnot, paste0(sample, State)))), function(x) sample(x, min(length(x), n))))]
  str(genes <- genes[genes %in% rownames(dat)])
  dat <- dat[genes, cells]
  dat <- dat[apply(dat, 1, max) != 0,,drop=F]
  dat <- dat - apply(dat, 1, min)
  dat <- dat / apply(dat,1, max)
  dat <- dat[, order(with(sampleAnnot[cells,],paste0(sample, State))), drop=F]
  
  pdf(dirout(out, "State_Markers_Monocle_",gnam,".pdf"), height=min(29, nrow(dat) * 0.3 + 3), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
  pheatmap::pheatmap(dat, cluster_rows=T, cluster_cols=F, annotation_col=sampleAnnot,show_colnames=FALSE)
  dev.off()
}



# ANNOTATE SEURAT OBJECT FURTHER AND DO ANALYSES

state.x <- as.character(pbmc@meta.data$State)
state.x[state.x %in% c("3", "4")] <- "IGNORED"
pbmc@meta.data$StateNoIntermediate <- state.x
state.x[state.x %in% c("2", "5")] <- "IGNORED"
pbmc@meta.data$StateEndpoints <- state.x

(save(pbmc, file=dirout(out, "Bcells_Monocle.RData")))
# } else {
#   load(dirout(out, "Bcells_Monocle.RData"))
# }

outS <- paste0(out, "Bcells_minpct0.1/")
dir.create(dirout(outS))

pbmc@meta.data <- subset(pbmc@meta.data, select=c(StateEndpoints, State))

seurat.min.pct <- 0.1
source("src/single_cell_RNA/95_Seurat2.R")





# DIFFERENCES BETWEEN STATES / FOR EACH PATIENT SEPARATELY (SEURAT WAY)----------------

pbmc@meta.data$patient <- substr(pbmc@meta.data$sample, 0, 2)
data.table(pbmc@meta.data)[State %in% c(1,6,7),.N, by=c("State", "patient")][order(patient)]

pbmc@ident <- factor(paste0(pbmc@meta.data$patient, "_", pbmc@meta.data$State))
names(pbmc@ident) <- pbmc@cell.names

cellCounts <- data.table(pbmc@meta.data)[State %in% c(1,6,7),.N, by=c("State", "patient")][order(patient)]
cellCounts <- cellCounts[patient != "KI" & N > 10]

if(!file.exists(dirout(out, "StateBySample_DiffGenes.RData"))){
diffGenes <- list()
for(pat in unique(cellCounts$patient)){
    if(nrow(cellCounts[patient == pat]) > 1){
      states <- cellCounts[patient == pat]$State
      message(pat,": ", states)
      for(i1 in 1:(length(states)-1)){
        for(i2 in (i1 + 1):length(states)){
          diffGenes[[paste0(pat, "_", states[i1],"vs",states[i2])]] <- FindMarkers(
            pbmc,  
            ident.1 = paste0(pat, "_", states[i1]), 
            ident.2 = paste0(pat, "_", states[i2]),
            test.use="negbinom") # This is done with min.pct = 0.1 (default)
        }
      }
    }
  }
  
  save(diffGenes, file=dirout(out, "StateBySample_DiffGenes.RData"))
} else{
  load(dirout(out, "StateBySample_DiffGenes.RData"))
}


diffGenes2 <- list()
lnam <- "PT_7vs1"
for(lnam in names(diffGenes)){
  pat <- gsub("^(\\w\\w)_(\\d)vs(\\d)$", "\\1", lnam)
  s1 <- gsub("^(\\w\\w)_(\\d)vs(\\d)$", "\\2", lnam)
  s2 <- gsub("^(\\w\\w)_(\\d)vs(\\d)$", "\\3", lnam)
  diffGenes2[[paste0(pat, "_", s1, "vs",s2)]] <- data.table(diffGenes[[lnam]], keep.rownames=T)[avg_diff > 0 & p_val < 0.0001]$rn
  diffGenes2[[paste0(pat, "_", s2, "vs",s1)]] <- data.table(diffGenes[[lnam]], keep.rownames=T)[avg_diff < 0 & p_val < 0.0001]$rn
}
sapply(diffGenes2, length)


gplots::venn(diffGenes2[grepl("7vs1", names(diffGenes2))])
gplots::venn(diffGenes2[grepl("1vs7", names(diffGenes2))])
gplots::venn(diffGenes2[grepl("6vs7", names(diffGenes2))])
gplots::venn(diffGenes2[grepl("7vs6", names(diffGenes2))])
gplots::venn(diffGenes2[grepl("6vs1", names(diffGenes2))])
gplots::venn(diffGenes2[grepl("1vs6", names(diffGenes2))])

source("src/single_cell_RNA/FUNC_Enrichr.R")




# PLOT GENES --------------------------------------------------------------
pDiffGenes <- do.call(rbind, lapply(names(diffGenes), function(nam) return(data.table(diffGenes[[nam]], name=nam, keep.rownames=TRUE))))
pDiffGenes[,qval := p.adjust(p_val, method="BH")]
pDiffGenes[,comparison := gsub(".+?_(\\dvs\\d)", '\\1', name)]
pDiffGenes[,patient := substr(name, 0,2)]
pDiffGenes[,name2 := paste0(comparison, "_", patient)]

comp.mins <- pDiffGenes[,.(qval_min = pmin(qval)), by=c("rn", "comparison")]
marker.genes <- unique(do.call(c, lapply(unique(comp.mins$comparison), function(comp){
  comp.mins[comparison == comp][order(qval_min, decreasing=TRUE)][1:20]$rn
})))

# order groups by similarity (of OR)
pDat <- pDiffGenes[rn %in% marker.genes]
if(length(unique(pDat$rn)) >= 2){
  try({
    orMT <- t(as.matrix(dcast.data.table(pDat, name2 ~ rn, value.var="avg_diff")[,-"name2",with=F]))
    orMT[is.na(orMT)] <- 1
    hclustObj <- hclust(dist(orMT))
    pDat$rn <- factor(pDat$rn, levels=hclustObj$labels[hclustObj$order])
  }, silent=T)
}
ggplot(pDat, aes(x=name2, y=rn, color=avg_diff, size=pmin(4, -log10(qval)))) + 
  geom_point() + scale_color_gradient2(low="blue", mid="white", high="red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(dirout(out, "StateBySample_Genes.pdf"), height=15, width=7)


# Get enrichments for intersections
ii <- gplots::venn(diffGenes2, intersections=TRUE,show.plot=F)
ii <- attr(ii, "intersections")
StateBySample_Intersections <- ii
save(StateBySample_Intersections, file=dirout(out, "StateBySample_Intersections.RData"))
counts <- sapply(ii, length)
counts[grepl(":", names(counts)) & counts > 10]
counts[counts > 10]
select <- ii[names(counts[counts > 10])]
sapply(select, length)
enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location")
enrichRes <- data.table()
for(grp.x in names(select)){
  ret=try(as.data.table(enrichGeneList.oddsRatio(select[[grp.x]],databases = enrichrDBs)),silent = FALSE)
  if(!any(grepl("Error",ret)) && nrow(ret) > 0){
    enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
  }
}
enrichRes <- enrichRes[hitLength > 3]
write.table(enrichRes[qval < 0.05], file=dirout(out, "StateBySample_Enrichr",".tsv"), sep="\t", quote=F, row.names=F)

enrichr.plot(enrichRes=enrichRes)
ggsave(dirout(out, "StateBySample_Enrichr.pdf"), height=12, width=7)


# Get enrichments for sepearate lists
enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location")
enrichRes <- data.table()
for(grp.x in names(diffGenes2)){
  ret=try(as.data.table(enrichGeneList.oddsRatio(diffGenes2[[grp.x]],databases = enrichrDBs)),silent = FALSE)
  if(!any(grepl("Error",ret)) && nrow(ret) > 0){
    enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
  }
}
enrichRes <- enrichRes[hitLength > 3]
write.table(enrichRes[qval < 0.05], file=dirout(out, "StateBySample_Enrichr_Separate",".tsv"), sep="\t", quote=F, row.names=F)

enrichr.plot(enrichRes=enrichRes)
ggsave(dirout(out, "StateBySample_Enrichr_Separate.pdf"), height=25, width=10)







# BRANCH POINT ANALYSIS PER PATIENT (MONOCLE WAY) -------------------------------------
# mcle.meta.data <- data.table(pData(mcle))
# mcle.meta.data[,patient := substr(sample, 0,2)]
# cellCounts <- mcle.meta.data[State %in% c(1,6,7),.N, by=c("State", "patient")][order(patient)]
# (cellCounts <- cellCounts[patient != "KI" & N > 10])
# 
# if(!file.exists(dirout(out, "StateBySample_DiffGenes.RData"))){
#   diffGenes <- list()
#   for(pat in unique(cellCounts$patient)){
#     if(nrow(cellCounts[patient == pat]) > 1){
#       states <- cellCounts[patient == pat]$State
#       message(pat,": ", states)
#       bp <- 1
#       for(bp in 1:3){
#         mcle.x <- mcle[,mcle.meta.data[patient == pat]$barcode]
#         BEAM_res <- BEAM(mcle.x, branch_point = bp, cores = 1)
#         BEAM_res <- BEAM_res[order(BEAM_res$qval),]
#         BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
#         }
#       }
#     }
#   }
#   
#   save(diffGenes, file=dirout(out, "StateBySample_DiffGenes.RData"))
# } else{
#   load(dirout(out, "StateBySample_DiffGenes.RData"))
# }
