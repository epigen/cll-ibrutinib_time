require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(pheatmap)
require(gplots)

project.init2("cll-time_course")
out <- "30_9_4_Clone_Plots/"
dir.create(dirout(out))

list.files(dirout())

(load(dirout("30_9_4_Signatures_nUMI_Cutoff_TFs/","Scores.RData")))
Dat1CLL <- Dat1[CellType == "CLL"]


# Genesets ----------------------------------------------------------------
(load(dirout("30_9_4_Signatures_nUMI_Cutoff_TFs/Genesets.RData")))
genesets.selected <- names(genesets)#[(grepl("d\\d\\d", names(genesets))) | grepl("^Bcells", names(genesets)) | names(genesets) == "Ibrutinib_treatment"]


# CLL patient samples, metadata
pat.dirs <- list.dirs(dirout("25_Patient_CLL"),recursive=F)
names(pat.dirs) <- gsub(dirout("25_Patient_CLL/"), "", pat.dirs)
pat.md <- lapply(pat.dirs, function(d) fread(paste0(d, "/MetaData.tsv")))
pat <- names(pat.md)[1]

pat.md[["VZS"]][,Clone := "CLL1"]
pat.md[["VZS"]][ClusterNames_0.5 == 3, Clone := "CLL2"]
pat.md[["VZS"]][ClusterNames_0.5 == 2, Clone := "CLL3"]
pat.md[["VZS"]][ClusterNames_0.5 == 5, Clone := "CLL4"]

pat.md[["PT"]][,Clone := "CLL1"]
pat.md[["PT"]][ClusterNames_0.5 == 6, Clone := "CLL2"]
pat.md[["PT"]][ClusterNames_0.5 == 5, Clone := "CLL3"]
pat.md[["PT"]][ClusterNames_0.5 == 4, Clone := "CLL4"]
pat.md[["PT"]][ClusterNames_0.5 == 7, Clone := "CLL5"]

pat.md[["FE"]][,Clone := "CLL1"]
pat.md[["FE"]][ClusterNames_0.5 %in% c(5,8), Clone := "CLL2"]
pat.md[["FE"]][ClusterNames_0.5 %in% c(7), Clone := "CLL3"]

pat.md[["PBGY"]][,Clone := "CLL1"]
pat.md[["PBGY"]][ClusterNames_0.5 %in% c(5), Clone := "CLL2"]
pat.md[["PBGY"]][ClusterNames_0.5 %in% c(7), Clone := "CLL3"]

pat.md <- do.call(rbind, pat.md)
pat.md[,timepoint := gsub("\\w+_(.+)", "\\1", sample)]
pat.md[timepoint == "d0",timepoint := "d000"]
pat.md[timepoint == "d30",timepoint := "d030"]

write.table(pat.md, file=dirout(out, "Patient_Clones_Meta.data.tsv"), sep="\t", quote=F, row.names=F)


# Timepoints vs samples ---------------------------------------------------


# ggplot(pat.md[patient == "PT"], aes(x=timepoint, fill=Clone)) + geom_bar(position="dodge")
for(ptx in unique(pat.md$patient)){
  ggplot(pat.md[patient == ptx], aes(x=timepoint, fill=Clone)) + geom_bar(position="fill") + theme_bw(24)
  ggsave(dirout(out, "Counts_", ptx, ".pdf"), height=7, width=7)
}

# ANALYSIS OVER TIME ----------------------------------------------------------------
pDat2 <- merge(Dat1CLL, pat.md[,c("rn", "Clone"), with=F], by="rn")
pDat2[,sample_cell := paste(patient, Clone, timepoint, sep="_")]
out.details <- paste0(out, "Details/")
dir.create(dirout(out.details))


geneset <- "Ibrutinib_treatment"
for(geneset in genesets.selected){
  
  # Plot all values as boxplots
  p <- ggplot(
    pDat2,
    aes_string(x="patient", y=geneset, group="sample_cell", fill="Clone", alpha="timepoint")) + 
    scale_color_manual(values=c("grey", "black")) + scale_alpha_manual(values=c(0.35, 0.5, 0.85, 0.85,1))
  ggsave(dirout(out.details, geneset, "_AllData.pdf"), width=7, height=7, plot=p + geom_boxplot(outlier.shape=NA) + coord_flip())
    
  summaryDat <- pDat2[, .(meanVal = mean(get(geneset)), sdVal = sd(get(geneset)), seVal = sd(get(geneset))/sqrt(.N)), by=c("patient", "Clone", "timepoint")]
  
  p <- ggplot(summaryDat, aes(x= timepoint, y=meanVal, group=paste0(patient, "_", Clone), color=patient)) + ggtitle(geneset)
  ggsave(dirout(out.details, geneset, "_AllData_line.pdf"), width=7, height=7, plot=p + geom_line())
  ggsave(dirout(out.details, geneset, "_AllData_lineSE.pdf"), width=7, height=7, plot=p + geom_line() + geom_errorbar(aes(ymin=meanVal-seVal, ymax=meanVal+seVal), width=0.1))
  ggsave(dirout(out.details, geneset, "_AllData_lineSE_grid.pdf"), width=15, height=15, 
    plot=p + geom_line(size=2) + geom_errorbar(aes(ymin=meanVal-sdVal, ymax=meanVal+sdVal), width=0.3, size=2) + facet_grid(patient ~ Clone))
}


# Make one large heatmap over time
overTime.sig <- data.table()
pat="VZS"
cell="CLL4"
geneset <- "Ibrutinib_treatment"
for(pat in unique(pDat2$patient)){
  for(cell in unique(pDat2$Clone)){
    x <- pDat2[patient == pat & Clone == cell]
    (tp <- unique(x[timepoint != "d000"]$timepoint)[1])
    for(tp in unique(x[timepoint != "d000"]$timepoint)){
      if(nrow(x[timepoint == "d000"]) > 5 & nrow(x[timepoint == tp]) >5){
        for(geneset in genesets.selected){
          p <- t.test(x[timepoint == "d000"][[geneset]], x[timepoint == tp][[geneset]])$p.value
          ef <- mean(x[timepoint == tp][[geneset]]) - mean(x[timepoint == "d000"][[geneset]])
          overTime.sig <- rbind(overTime.sig, data.table(patient=pat, CellType=cell, pvalue=p, Diff=ef,timepoint=tp, geneset=geneset))
        }
      }
    }
  }
}
overTime.sig[,pvalue2 := pmin(5, -1*log10(pvalue))]
save(overTime.sig, file=dirout(out, "OverTime_Tests.RData"))

overTime.sig[,geneset := gsub("d30", "d030", geneset)]

overTime.sig[geneset == "Ibrutinib_treatment"]
ggplot(
  overTime.sig[grepl("Bcells", geneset) | geneset == "Ibrutinib_treatment"], 
  aes(x=paste0(CellType, "_", patient, "_", timepoint), y=geneset, color=log10(Diff), size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(. ~ patient, scale="free")
ggsave(dirout(out, "0_Overview.pdf"),height=8, width=15)


ggplot(
  overTime.sig[!(grepl("Bcells", geneset) | geneset == "Ibrutinib_treatment")], 
  aes(x=paste0(CellType, "_", patient, "_", timepoint), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(. ~ patient, scale="free")
ggsave(dirout(out, "0_Overview_patSigs.pdf"),height=8, width=15)



pDat2[,.N, by=c("sample", "CellType")]

ggplot(
  overTime.sig[grepl("Bcells", geneset) | geneset == "Ibrutinib_treatment"], 
  aes(x=paste0(CellType, "_", timepoint), y=geneset, fill=Diff)) + 
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red") + theme_bw(24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(. ~ patient, scale="free")
ggsave(dirout(out, "0_Overview_Tiles.pdf"),height=4, width=15)

ggplot(
  overTime.sig[!(grepl("Bcells", geneset) | geneset == "Ibrutinib_treatment")], 
  aes(x=paste0(CellType, "_", timepoint), y=geneset, fill=Diff)) + 
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red") + theme_bw(24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(. ~ patient, scale="free")
ggsave(dirout(out, "0_Overview_patSigs_Tiles.pdf"),height=7, width=15)





# Timepoints --------------------------------------------------------------

pat.md[,timepoint := gsub("\\w+_(d\\d+)", "\\1", sample)]
pat.md[timepoint == "d0", timepoint := "d000"]
pat.md[timepoint == "d30", timepoint := "d030"]

pDat <- pat.md[,.N, by=c("patient", "timepoint", "Clone")]
pDat[,sumN := sum(N), by=c("patient", "timepoint")]
pDat[,maxClone := max(N), by=c("patient", "Clone")]
pDat[,fraction := N / sumN * 100]
ggplot(pDat[patient == "PT"], aes(x=timepoint, y=fraction, color=patient, linetype=Clone, group=paste0(Clone, "_", patient))) + 
  geom_line()

ggplot(pDat, aes(x=timepoint, y=fraction, color=patient, linetype=Clone, group=paste0(Clone, "_", patient))) + 
  geom_line(size=2) + ylim(0,100) + ylab("Clone count (%)") + xlab("Timepoint (days)") + theme_bw(24)
ggsave(dirout(out, "CLL_percentage.pdf"), height=7, width=7)









# CNVs --------------------------------------------------------------------


outCNV <- paste0(out, "CNV/")
dir.create(dirout(outCNV))

(load(dirout("42_1_CNVs_Sigs_nUMI_Cutoff_min50/","Scores.RData")))
Dat1Clones <- Dat1
pat.md <- fread(dirout(out, "Patient_Clones_Meta.data.tsv"))
pDat2 <- merge(Dat1Clones, pat.md[,c("rn", "Clone"), with=F], by="rn")
genesets <- colnames(Dat1Clones)[grepl("(X|Y|\\d+)_\\d+", colnames(Dat1Clones))]
# Make one large heatmap over time
overTime.sig <- data.table()
for(pat in unique(pDat2$patient)){
  patDat <- pDat2[patient == pat]
  message(pat)
  print(clones <- unique(patDat$Clone))
  for(cl1 in 1:(length(clones)-1)){
    for(cl2 in (cl1 + 1):length(clones)){
      for(geneset in genesets){
        d1 <- patDat[Clone == clones[cl1]][[geneset]]
        d2 <- patDat[Clone == clones[cl2]][[geneset]]
        if(length(d1) > 5 & length(d2) > 5){
          p <- t.test(d1, d2)$p.value
          ef <- mean(d1) - mean(d2)
          overTime.sig <- rbind(overTime.sig, data.table(patient=pat, Clone1=clones[cl1], Clone2=clones[cl2], pvalue=p, Diff=ef,geneset=geneset))
        }
      }
    }
  }
}

overTime.sig[,qvalue := p.adjust(pvalue, method="BH")]
overTime.sig[,plogq := pmin(5, -1*log10(qvalue))]
overTime.sig[,up_in := ifelse(Diff > 0, Clone1, Clone2)]
overTime.sig[,hits := sum(qvalue < 0.1 & abs(Diff) > 0.5), by=c("geneset", "up_in", "patient")]
overTime.sig[,cloneN := length(unique(c(Clone1, Clone2))), by=c("patient")]

save(overTime.sig, file=dirout(outCNV, "OverTime_Tests.RData"))

with(overTime.sig[cloneN-3 < hits], table(geneset, patient))

(load(dirout("42_1_CNVs_Sigs_nUMI_Cutoff_min50/genesets.RData")))