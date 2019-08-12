
require("Seurat")
require("pheatmap")

outX <- paste0(outS, "overTime/")
dir.create(dirout(outX))


pbmc@ident <- pbmc@data.info$sample
names(pbmc@ident) <- pbmc@cell.names

table(pbmc@ident)

comp.list <- list(
  PBGY=c("PBGY1_0d", "PBGY7_150d"),
  PT=c("PT_d0", "PT_d120"),
  VZS=c("VZS_d0", "VZS7_120d"),
  FE=c("FE_FE1_d0", "FE7_120d")
  )

target.file <- dirout(outX,"DiffExprRes.RData")
res <- list()
if(!file.exists(target.file)){
  for(l.nam in names(comp.list)){
    message(l.nam)
    if(all(comp.list[[l.nam]] %in% pbmc@ident)){
      try({
        x <- FindMarkers(pbmc,  ident.1 = comp.list[[l.nam]][1], ident.2 = comp.list[[l.nam]][2], test.use=seurat.diff.test, min.pct=0, thresh.use=0)    
        res[[l.nam]] <- data.table(x, keep.rownames=TRUE)
      }, silent=TRUE)
    }
  } 
  save(res, file=target.file)
} else {
  load(target.file)
}

resL <- lapply(res, function(x) setNames(x$avg_diff, nm=x$rn))
str(resL)
genes.union <- Reduce(union, lapply(resL, names))

mat <- do.call(cbind, lapply(resL, function(x) x[genes.union]))
dt <- data.table(mat)

mat.cor <- cor(mat, use="pairwise.complete.obs")

write.csv(mat.cor, dirout(outX, "0_Correlations.csv"))

try({
  pdf(dirout(outX, "HM.pdf"), width=5, height=5, onefile=F)
  pheatmap(mat.cor)
  dev.off()
}, silent=TRUE)

patients <- colnames(dt)
if(length(patients) > 1){
  for(i1 in 1:(length(patients)-1)){
    for(i2 in (i1+1):length(patients)){
      pat1 <- patients[i1]
      pat2 <- patients[i2]
      message(pat1, " vs ", pat2)
      ggplot(dt, aes_string(x=pat1, y=pat2)) + geom_hex()
      ggsave(dirout(outX, pat1, "_vs_", pat2,".pdf"), width=7, height=7)
    }
  }
}