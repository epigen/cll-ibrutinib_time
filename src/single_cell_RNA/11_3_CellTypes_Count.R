require("project.init")
project.init2("cll-time_course")
out <- "11_CellTypes/2017_04_26_UnresolvedTcells/"
dir.create(dirout(out))

dirs <-list.dirs(dirout(out), recursive=F)


xxx <- data.table(do.call(rbind, lapply(strsplit(gsub(".*\\/(.+)$", "\\1", dirs), "_"), function(vec) 
  return(c(paste(vec[1:(length(vec)-1)], collapse="_"), vec[length(vec)])))))
  
samples <- unique(xxx$V1)
pDat.scRNA <- data.table()
for(sample.x in samples){
  cells <- xxx[V1 == sample.x]$V2
  
  pDat <- do.call(rbind, lapply(cells, function(cc){
    dat <- fread(dirout(out, sample.x, "_", cc, "/","Cellcounts.tsv"))
    dat$cell <- cc
    colnames(dat)[2] <- "cnt"
    return(dat)
    }))
  
  pDat$sample <- factor(pDat$sample, levels=sort(unique(pDat$sample)))
  
  if(sample.x == "allDataBest_NoDownSampling_noIGH") pDat.scRNA <- pDat
  
  ggplot(pDat, aes(y=sample, x=cell, fill=fraction)) + geom_tile() + theme_bw(24)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(dirout(out, sample.x, "Fractions.pdf"), width=10, height=10)
  
  ggplot(pDat, aes(y=sample, x=cell, fill=log10(cnt+1))) + geom_tile() + theme_bw(24)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(dirout(out, sample.x, "Counts.pdf"), width=10, height=10)
}



# FACS Data ---------------------------------------------------------------
facsDat <- fread("metadata/facs_quantification_17042017.csv")
colnames(facsDat) <- make.names(gsub("(\\+|\\_| )", "", colnames(facsDat)))
facsDat <- facsDat[!apply(is.na(facsDat), 1, sum) > ncol(facsDat)/2]
facsDat <- facsDat[,(!apply(is.na(facsDat), 2, sum) > nrow(facsDat)/2), with=F]

facsDat[, Bcells := CD5CLL + normalBcells]
facsDat[, Monos := CD14Myeloidcells + 0]
facsDat[, Tcells2 := CD3Tcells + CD56NKcells + CD56NKcells + CD4Tcells + NKTcells + Tcells]

facsDat$PatTime <- with(facsDat, paste0(PatientID, "_",time, "d"))

pat <- "PBGY"
for(pat in unique(facsDat$PatientID)){
  pDat <- facsDat[PatientID == pat][,c("PatTime", "Tcells2", "Bcells", "Monos"), with=F]
  pDat$PatTime <- factor(pDat$PatTime, levels=pDat$PatTime)
  pDat <- melt(pDat)
  ggplot(pDat, aes(x=PatTime, color=variable,y=value, group=variable)) + geom_line(size=3) + theme_bw(24)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(dirout(out, "FACS", "_", pat, ".pdf"), width=8, height=8)
  
  pDat.scRNA2 <- pDat.scRNA[,c("sample", "fraction", "cell"), with=F]
  pDat.scRNA2[,sample := gsub("_d(\\d+)", "_\\1d", sample)]
  pDat.scRNA2[,sample := gsub("PBGY\\d", "PBGY", sample)]
  pDat.scRNA2$fraction <- pDat.scRNA2$fraction * 100
  
  colnames(pDat) <- c("sample", "cell", "fraction")
  pDat[cell == "Tcells2", cell := "Tcells"]
  
  ggplot(pDat[cell != "Tcells2"][, Data := "FACS"], aes(x=sample, color=cell,y=fraction, group=paste(cell,Data), linetype=Data)) + 
    geom_line(size=2) + geom_line(data=pDat.scRNA2[sample %in% pDat$sample][cell %in% pDat$cell][,Data := "scRNA"],size=2) +
    theme_bw(24)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(dirout(out, "FACS_vsRNA", "_", pat, ".pdf"), width=8, height=8)
  
}

