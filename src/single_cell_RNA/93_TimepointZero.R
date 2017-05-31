
require("Seurat")
require("pheatmap")

outX <- paste0(outS, "TimepointZero/")
dir.create(dirout(outX))


pbmc@ident <- pbmc@data.info$sample
names(pbmc@ident) <- pbmc@cell.names

table(pbmc@ident)

t0.samples <- c("FE_FE1_d0", "KI_KI1_d0", "PBGY1_0d", "PT_d0","VZS_d0")

target.file <- dirout(outX,"DiffExprRes.RData")
res <- list()
if(!file.exists(target.file)){
  for(i1 in 1:(length(t0.samples)-1)){
    for(i2 in (i1+1):length(t0.samples)){
      g1 <- t0.samples[i1]
      g2 <- t0.samples[i2]
      message(g1, " vs ",g2)
      if(all(c(g1, g2) %in% pbmc@ident)){
        try({
          x <- FindMarkers(pbmc,  ident.1 = g1, ident.2 = g2, min.pct=0, thresh.use=0)    
          res[[paste0(g1,"_vs_",g2)]] <- data.table(x, keep.rownames=TRUE)
        }, silent=TRUE)
      }
    }
  }
  save(res, file=target.file)
} else {
  load(target.file)
}

resL <- lapply(res, function(x) setNames(x$avg_diff, nm=x$rn))
str(resL)

for(l.nam in names(resL)){
  nn <- strsplit(l.nam, "_vs_")[[1]]
  resL[[paste0(nn[2], "_vs_", nn[1])]] <- -resL[[l.nam]]
}

genes.union <- Reduce(union, lapply(resL, names))

mat <- do.call(cbind, lapply(resL, function(x) x[genes.union]))
dt <- data.table(mat)

mat.cor <- cor(mat, use="pairwise.complete.obs")

write.csv(mat.cor, dirout(outX, "0_Correlations.csv"))


for(samp in t0.samples){
  try({
    pdf(dirout(outX, "0_HM_",samp,".pdf"), width=5, height=5, onefile=F)
    pheatmap(mat.cor[grepl(paste0("^",samp), colnames(mat.cor)),grepl(paste0("^",samp), colnames(mat.cor))])
    dev.off()
  }, silent=TRUE)
}

try({
  pdf(dirout(outX, "0_HM.pdf"), width=8, height=8, onefile=F)
  pheatmap(mat.cor)
  dev.off()
}, silent=TRUE)


outX2 <- paste0(outX, "scatterplots/")
dir.create(dirout(outX2))

patients <- colnames(dt)
for(i1 in 1:(length(patients)-1)){
  for(i2 in (i1+1):length(patients)){
    pat1 <- patients[i1]
    pat2 <- patients[i2]
    message(pat1, " vs ", pat2)
    ggplot(dt, aes_string(x=pat1, y=pat2)) + geom_hex()
    ggsave(dirout(outX2, pat1, "_vs_", pat2,".pdf"), width=7, height=7)
  }
}