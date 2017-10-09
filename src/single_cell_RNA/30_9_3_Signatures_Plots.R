require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(pheatmap)
require(gplots)

project.init2("cll-time_course")
out <- "30_9_Plots/"
dir.create(dirout(out))

list.files(dirout())

(load(dirout("30_9_Signatures_inclDay30_downLists2/","Scores.RData")))

Dat1 <- Dat1[CellType == "CLL"]

Dat1


# Genesets ----------------------------------------------------------------
(load(dirout("30_9_Signatures_inclDay30_downLists2/Genesets.RData")))
genesets.selected <- names(genesets)[(grepl("d\\d\\d", names(genesets))) | grepl("^Bcells", names(genesets)) | names(genesets) == "Ibrutinib_treatment"]


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
pat.md[["PT"]][ClusterNames_0.5 == 4, Clone := "CLL5"]

pat.md[["FE"]][,Clone := "CLL1"]
pat.md[["FE"]][ClusterNames_0.5 %in% c(5,8), Clone := "CLL2"]
pat.md[["FE"]][ClusterNames_0.5 %in% c(7), Clone := "CLL3"]

pat.md[["PBGY"]][,Clone := "CLL1"]
pat.md[["PBGY"]][ClusterNames_0.5 %in% c(5), Clone := "CLL2"]
pat.md[["PBGY"]][ClusterNames_0.5 %in% c(7), Clone := "CLL3"]



# ANALYSIS OVER TIME ----------------------------------------------------------------
pDat2 <- merge(Dat1, pat.md[,c("rn", "Clone"), with=F], by="rn")
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
  ggsave(dirout(out.details, geneset, "_AllData2.pdf"), width=7, height=7, plot=p + geom_violin() + coord_flip())
}

# Make one large heatmap over time
overTime.sig <- data.table()
pat="VZS"
cell="CD4"
geneset <- "Mono_up"
for(pat in unique(pDat2$patient)){
  for(cell in unique(pDat2$Clone)){
    x <- pDat2[patient == pat & Clone == cell]
    tp <- unique(x[timepoint != "d000"]$timepoint)[1]
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

ggplot(
  overTime.sig[grepl("Bcells", geneset) | geneset == "Ibrutinib_treatment"], 
  aes(x=paste0(CellType, "_", patient, "_", timepoint), y=geneset, color=Diff, size=pvalue2)) + 
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
