require("project.init")
require("Seurat")
require("gplots")
require(methods)
require(pheatmap)
require(foreach)
# require(enrichR) #devtools::install_github("definitelysean/enrichR")

project.init2("cll-time_course")

source("src/single_cell_RNA/FUNC_Enrichr.R") #devtools::install_github("definitelysean/enrichR")


(out <- paste0("17_Rec_Lig/"))
dir.create(dirout(out))




# LOAD DATA ---------------------------------------------------------------

ramilowski.rec_lig <- fread("~/resources_nfortelny/PairsLigRec_Ramilowski_NatComm_2015.txt")
(load(dirout("10_Seurat_raw/inclDay30_noIGHLK_negbinom/", "inclDay30_noIGHLK.RData")))
metaDat <- data.table(pbmc@data.info, keep.rownames=T)[nUMI > 1000 & nUMI < 3000]
all(metaDat$rn %in% colnames(pbmc@data))
ra2 <- ramilowski.rec_lig[Pair.Evidence == "literature supported"][,c("Receptor.ApprovedSymbol", "Ligand.ApprovedSymbol"), with=F]
colnames(ra2) <- c("A", "B")
ra2 <- ra2[A %in% row.names(pbmc@data) & B %in% row.names(pbmc@data)]
sigProteins <- unique(union(ra2$A, ra2$B))


# PREPARE DATA - split by groups ------------------------------------------

dat <- pbmc@data[sigProteins,metaDat$rn]
dat <- dat[apply(dat,1,sum) > 0,]
dat <- dat - apply(dat,1,min)
dat <- dat / apply(dat,1,max)
apply(dat[sample(1:nrow(dat),50),], 1, quantile)
dat <- as.data.frame(t(as.matrix(dat)))
grps <- factor(paste0(metaDat$sample, "_", metaDat$CellType))
datSplit <- split(dat, grps)
lapply(datSplit, function(mat) unique(gsub("[A-Z]+\\-", "", row.names(mat))))



# FUNCTION TO GET INTERACTION STRENGTH FROM MEAN GENE EXPRESSION ----------

getIntStrength <- function(datMeans, ra2){
  intStrength <- data.table()
  
  cellTypes <- unique(datMeans$celltype)
  for(tt in unique(datMeans$time)){
    for(pat in unique(datMeans$patient)){
      for(ct1 in unique(cellTypes)){
        for(ct2 in unique(cellTypes)){
          xx1 <- datMeans[time == tt & patient == pat & celltype == ct1]
          xx2 <- datMeans[time == tt & patient == pat & celltype == ct2]
          if(nrow(xx1) == 0 | nrow(xx2) == 0) next
          val1 <- xx1$value
          names(val1) <- xx1$Var1
          val2 <- xx2$value
          names(val2) <- xx2$Var1        
          intStrength <- rbind(intStrength, data.table(
            str=mean(val1[ra2$A] * val2[ra2$B], na.rm=TRUE), 
            time=tt, patient=pat,
            cell1 = ct1, cell2 = ct2
          ))
        }
      }
    }
  }
  return(intStrength)
}



# RUN THIS FOR BOOTSTRAPED SAMPLES ----------------------------------------

intStrength.file <- dirout(out, "IntStrength.RData")
if(!file.exists(intStrength.file)){
  require(doMC)
  rm(list=c("pbmc"))
  registerDoMC(cores=10)
  
  boots <- 500
  boot.size <- 50
  intStr2 <-  foreach(boot = 1:boots) %dopar% {
    message("Bootstrap:", boot)
    datSplit2 <- lapply(datSplit, function(mat) colMeans(mat[sample(row.names(mat),size=boot.size,replace=TRUE),]))
    datMeans <- do.call(cbind, datSplit2)
    datMeans <- data.table(melt(datMeans))
    datMeans[,celltype := gsub("(.+)_(.+)_(.+)", "\\3", Var2)]
    datMeans[,time := gsub("(.+)_(.+)_(.+)", "\\2", Var2)]
    datMeans[,patient := gsub("(.+)_(.+)_(.+)", "\\1", Var2)]
    return(getIntStrength(datMeans, ra2))
  }
  intStrength <- do.call(rbind, intStr2)
  intStrength[time == "d150", time := "d120"]
  intStrength[time == "d30", time := "d030"]
  save(intStrength, file=dirout(out, "IntStrength.RData"))
} else {
  (load(dirout(out, "IntStrength.RData")))  
}

# Analyzed BOOTSTRAPED SAMPLES ----------------------------------------

cellTypes <- c("CD4", "CD8", "Mono", "CLL")
intStrength <- intStrength[cell1 %in% cellTypes & cell2 %in% cellTypes]

ggplot(intStrength[patient=="PBGY" & cell1 == "Mono" & cell2 == "Mono"], aes(x=time, y=str)) + geom_jitter(width=0.1)
ggplot(intStrength[patient=="PT" & cell1 == "Mono" & cell2 == "Mono"], aes(x=time, y=str)) + geom_jitter(width=0.1)

pat <- "PT"
ct1 <- "CD4"
ct2 <- "CD8"
ggplot(intStrength[patient==pat & cell1 == ct1 & cell2 == ct2], aes(x=time, y=str)) + geom_jitter(width=0.1)
sigTable <- data.table()
for(pat in unique(intStrength$patient)){
  for(ct1 in unique(intStrength$cell1)){
    for(ct2 in unique(intStrength$cell2)){
      message(pat, " ", ct1, " ", ct2)
      xx <- intStrength[patient == pat & cell1 == ct1 & cell2 == ct2]
      if(length(unique(xx$time)) == 1) next
      xx$time <- factor(xx$time, levels=c("d0", unique(xx[time != "d0"]$time)))
      fit <- lm(data=xx,formula=str ~ time)
      sfit <- summary(fit)
      sigTable <- rbind(sigTable, data.table(data.table(sfit$coefficients, keep.rownames=TRUE),
                 fstat_pval = 1 - pf(sfit$fstatistic[1], sfit$fstatistic[2], sfit$fstatistic[3]),
                 patient = pat,
                 cell1 = ct1,
                 cell2 = ct2))
    }
  }
}
sigTable <- sigTable[rn != "(Intercept)"]
sigTable$pvalue <- sigTable[["Pr(>|t|)"]]
sigTable[,qvalue := p.adjust(pvalue, method="BH")]
colnames(sigTable) <- make.names(colnames(sigTable))
sigTable <- sigTable[,-"Pr...t..",with=F]
sigTable <- sigTable[order(t.value,decreasing=TRUE)]
table(sigTable$qvalue < 0.05)
sigTable[,rn := gsub("timed", "", rn)]


##
## T-VALUES
##
i <- 1
for(i in c(1:10, (nrow(sigTable)-10):nrow(sigTable))){
  xx <- sigTable[i]
  ggplot(intStrength[patient==xx$patient & cell1 == xx$cell1 & cell2 == xx$cell2], aes(x=time, y=str)) + geom_jitter(width=0.1) +
    ggtitle(paste(xx$patient, xx$cell1, xx$cell2))
  ggsave(dirout(out, "Example_",i,"_",xx$patient, xx$cell1, xx$cell2, ".pdf"))
}
ggplot(sigTable, aes(x=cell1, y=cell2, fill=t.value)) + 
  geom_tile() + facet_grid(patient ~ rn) + scale_fill_gradient2() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "IntStrength_Tvalues.pdf"), height=7, width=7) 




# Specific plots on full data ---------------------------------------------
datSplit2 <- lapply(datSplit, function(mat) colMeans(mat))
datMeans <- do.call(cbind, datSplit2)
datMeans <- data.table(melt(datMeans))
datMeans[,celltype := gsub("(.+)_(.+)_(.+)", "\\3", Var2)]
datMeans[,time := gsub("(.+)_(.+)_(.+)", "\\2", Var2)]
datMeans[,patient := gsub("(.+)_(.+)_(.+)", "\\1", Var2)]

i <- 1
sRow <- sigTable[i]
timeDat <- foreach(timep = unique(datMeans[patient == sRow$patient & celltype == sRow$cell1]$time)) %do% {
  xx1 <- datMeans[patient == sRow$patient & celltype == sRow$cell1 & time == timep]
  xx2 <- datMeans[patient == sRow$patient & celltype == sRow$cell2 & time == timep]
  val1 <- xx1$value
  names(val1) <- xx1$Var1
  val2 <- xx2$value
  names(val2) <- xx2$Var1        
  data.table(gene1 = ra2$A, gene2=ra2$B, value=val1[ra2$A] * val2[ra2$B], time=timep)
}
timeDat <- do.call(rbind, timeDat)
timeDat[, names := paste0(gene1, "_",gene2)]
pDT <- dcast.data.table(timeDat, names ~ time, value.var="value")
pMT <- as.matrix(pDT[,-"names",with=F])
row.names(pMT) <- pDT$names
pdf(dirout(out, "Example_HM.pdf"),height=25, width=10,onefile=F)
pheatmap(pMT[rowSums(pMT,na.rm=TRUE) != 0,], scale="row",fontsize_row=5)
dev.off()


pdf(dirout(out, "Example_HM_Small.pdf"),height=5, width=5,onefile=F)
pheatmap(pMT[rowSums(pMT,na.rm=TRUE) > max(rowSums(pMT,na.rm=TRUE)) - 0.6,])
dev.off()

(load(dirout("10_Seurat_raw/inclDay30_noIGHLK_negbinom/", "inclDay30_noIGHLK.RData")))

gene1 <- "CD53"
gene2 <- "IL4"

gene1 <- "KLRD1"
gene2 <- "HLA-B"
metaDat1 <- subset(pbmc@data.info, grepl("PBGY", sample) & CellType == sRow$cell1)
metaDat2 <- subset(pbmc@data.info, grepl("PBGY", sample) & CellType == sRow$cell2)

(pDat1 <- data.table(sample=metaDat1$sample, value = pbmc@data[gene1,row.names(metaDat1)]))
(pDat2 <- data.table(sample=metaDat2$sample, value = pbmc@data[gene2,row.names(metaDat2)]))

pDatMean1 <- pDat1[,.(value=mean(value)), by=c("sample")]
pDatMean2 <- pDat2[,.(value=mean(value)), by=c("sample")]

p <- ggplot(pDat1, aes(x=sample, y=value,)) + geom_violin() + geom_point(data=pDatMean1)
ggsave(dirout(out, "Example_Violin",gene1,".pdf"), height=7, width=7, plot=p)
ggsave(dirout(out, "Example_Violin_log",gene1,".pdf"), height=7, width=7, plot=p + scale_y_log10())

p <- ggplot(pDat2, aes(x=sample, y=value,)) + geom_violin() + geom_point(data=pDatMean2)
ggsave(dirout(out, "Example_Violin",gene2,".pdf"), height=7, width=7, plot=p)
ggsave(dirout(out, "Example_Violin_log",gene2,".pdf"), height=7, width=7, plot=p + scale_y_log10())



# AVERAGE OVER BOOTSTRAPS -------------------------------------------------

avDat <- intStrength[cell1 %in% cellTypes & cell2 %in% cellTypes][,.(str=mean(str)), by=c("time", "patient", "cell1", "cell2")]
avDat[,base := str[time == "d0"], by=c("patient", "cell1", "cell2")]
avDat[,strNorm := str - base]
ggplot(avDat, aes(x=cell1, y=cell2, fill=strNorm)) + 
  geom_tile() + facet_grid(patient ~ time) + scale_fill_gradient2() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "IntStrength_MeansNorm.pdf"), height=7, width=7)

ggplot(avDat, aes(x=cell1, y=cell2, fill=str)) + 
  geom_tile() + facet_grid(patient ~ time) + scale_fill_gradient2() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "IntStrength_Means.pdf"), height=7, width=7)






# ON FULL DATA ----------------------------------------------------------------
intStrength.full <- getIntStrength(datMeans, ra2)
intStrength.full[time == "d150", time := "d120"]
intStrength.full[time == "d30", time := "d030"]

avDat <- intStrength.full[cell1 %in% cellTypes & cell2 %in% cellTypes]
avDat[,base := str[time == "d0"], by=c("patient", "cell1", "cell2")]
avDat[,strNorm := str - base]
ggplot(avDat, aes(x=cell1, y=cell2, fill=strNorm)) + 
  geom_tile() + facet_grid(patient ~ time) + scale_fill_gradient2() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "FullData_IntStrength_MeansNorm.pdf"), height=7, width=7)

ggplot(avDat, aes(x=cell1, y=cell2, fill=str)) + 
  geom_tile() + facet_grid(patient ~ time) + scale_fill_gradient2() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "FullData_IntStrength_Means.pdf"), height=7, width=7) 

# compare to Bootstrap averages
ggplot(merge(intStrength.full,intStrength[cell1 %in% cellTypes & cell2 %in% cellTypes][,.(str=mean(str)), by=c("time", "patient", "cell1", "cell2")], by=c("time", "patient", "cell1","cell2")),
       aes(x=str.x, y=str.y)) + geom_point() + xlab("Full data") + ylab("Averaged Bootstraps")
ggsave(dirout(out,"FulLData_Comparison.pdf"))
