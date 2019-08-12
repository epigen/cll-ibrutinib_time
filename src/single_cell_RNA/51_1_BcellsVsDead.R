require("project.init")
# install.packages("nloptr",repos="https://cloud.r-project.org")
require(Seurat)
require(methods)
project.init2("cll-time_course")


out <- "51_ClusterAlignment/"
dir.create(dirout(out))

list.files(dirout(out, "PT/Cluster_res.0.8/"))

files <- list(
  FE = "FE/Cluster_res.0.8/Diff_Cluster0vs3.tsv",
  PT = "PT/Cluster_res.0.8/Diff_Cluster0vs4.tsv",
  PT2 = "PT2/Cluster_res.0.8/Diff_Cluster0vs4.tsv",
  VZS = "VZS/Cluster_res.0.8/Diff_Cluster0vs1.tsv",
  PBGY = "PBGY/Cluster_res.0.8/Diff_Cluster1vs8.tsv")


dat <- lapply(files, function(f) fread(dirout(out, f)))

lapply(dat,head)

hits <- lapply(dat, function(DT){
  DT2 <- DT[order(V3)]
  return(DT2[!grepl("^MT\\-", V1) & !grepl("RP", V1) & !grepl("IG", V1) & !grepl("HLA", V1) & V2 < 0.05 & abs(V3) > 0.3]) # 
})

hitSets <- list()
for(nam in names(hits)){
  hitSets[[paste0("CLL_", nam)]] <- hits[[nam]][V3 > 0.3]$V1
  hitSets[[paste0("DeadCells_", nam)]] <- hits[[nam]][V3 < -0.3]$V1
}

hitSets

# ENRICHR
source("src/single_cell_RNA/FUNC_Enrichr.R")
enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location")
enrichRes <- data.table()

for(grp.x in names(hitSets)){
  ret=try(as.data.table(enrichGeneList.oddsRatio(hitSets[[grp.x]],databases = enrichrDBs,fdr.cutoff=NULL)),silent = FALSE)
  if(!any(grepl("Error",ret)) && nrow(ret) > 0){
    enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
  }
}

enrichr.plot(enrichRes)
write.table(enrichRes[qval < 0.05], file=dirout(out, "EnrichOR_clean",".tsv"), sep="\t", quote=F, row.names=F)
enrichRes2 <- enrichRes[qval < 0.05]
ggsave(dirout(out, "EnrichOR_clean",".pdf"), width=min(29, 6+ length(unique(enrichRes2$grp))*0.3), height=min(29, length(unique(enrichRes2$category))*0.3 + 4))



enrichRes2 <- data.table()
hitSets2 <- list()
ii <- gplots::venn(hitSets[grepl("DeadCells_", names(hitSets))], intersections=TRUE)
ii <- attr(ii, "intersections")
hitSets2[["Dead"]] <- ii$'DeadCells_FE:DeadCells_PT:DeadCells_PT2:DeadCells_VZS:DeadCells_PBGY'

ii <- gplots::venn(hitSets[grepl("CLL_", names(hitSets))], intersections=TRUE)
ii <- attr(ii, "intersections")
hitSets2[["CLL"]] <- ii$'CLL_FE:CLL_PT:CLL_PT2:CLL_VZS:CLL_PBGY'

for(grp.x in names(hitSets2)){
  ret=try(as.data.table(enrichGeneList.oddsRatio(hitSets2[[grp.x]],databases = enrichrDBs,fdr.cutoff=NULL)),silent = FALSE)
  if(!any(grepl("Error",ret)) && nrow(ret) > 0){
    enrichRes2 <- rbind(enrichRes2, data.table(ret, grp = grp.x))
  }
}

enrichr.plot(enrichRes2)





# FIND CELLS IN OTHER DATASETS --------------------------------------------
cluster.nrs <- list(
  FE = 3,
  PT = 4,
  PT2 = 4,
  VZS = 1,
  PBGY = 8)

dead.cell.codes <- list()
dead.cell.objects <- list()
cl.nr <- "FE"
for(cl.nr in names(cluster.nrs)){
  load(dirout(out,cl.nr, "/", cl.nr,".RData"))
  dead.cell.objects[[cl.nr]] <- pbmc
  xx <- data.table(pbmc@meta.data, keep.rownames=T)[res.0.8 == cluster.nrs[[cl.nr]]]
  for(samp in unique(xx$sample)){
    if(is.null(dead.cell.codes[[samp]])) dead.cell.codes[[samp]] <- c()
    dead.cell.codes[[samp]] <- unique(c(dead.cell.codes[[samp]], gsub("\\-\\d$", '' , xx[sample == samp]$rn)))
  }
}


# list.files(dirout("11_CellTypes_tobit/allDataBest_NoDownSampling_noIGH_Bcells/Bcells.RData"))
load(dirout("10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.RData"))


for(samp in names(dead.cell.codes)){
  message(samp)
  xx <- data.table(pbmc@data.info, keep.rownames=TRUE)[sample == samp]
  xx[,rn := gsub("\\-\\d+$","", rn)]
  print(table(xx[rn %in% dead.cell.codes[[samp]],]$cellType))
  print(paste(round(sum(dead.cell.codes[[samp]] %in% xx$rn)/length(dead.cell.codes[[samp]]),3) * 100, "% of dead cells kept in the analysis"))
}

