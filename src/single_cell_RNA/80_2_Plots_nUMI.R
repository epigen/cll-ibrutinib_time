require("project.init")
require(Seurat)
require(methods)
require(monocle)

project.init2("cll-time_course")
out <- "80_2_Plots_nUMI_Cutoff/"
dir.create(dirout(out))

source("src/single_cell_RNA/FUNC_Enrichr.R")

celltypes <- c("CD4", "CD8", "Mono", "CLL")



# LOAD DATA ---------------------------------------------------------------

# Load full single cell dataset
if(!"pbmcFull" %in% ls()){
  (load(dirout("10_Seurat_raw/inclDay30_noIGHLK_negbinom/inclDay30_noIGHLK.RData")))
  pbmcFull <- pbmc
  rm(list="pbmc")
  pbmcFull@scale.data <- NULL
  pbmcFull@raw.data <- NULL
}

# Load CLL single cells
if(FALSE){
  (load(dirout("11_2_3_CellTypes_UMI_Cutoff/CLL/CLL.RData")))
  pbmcCLL <- pbmc
  rm(list="pbmc")
  pbmcCLL@scale.data <- NULL
  pbmcCLL@raw.data <- NULL
}



# Load Monocle
# (load(dirout("52_3_Monocle_inclDay30_Endpoints/", "MonocleDat.RData")))

# Load Signatures for whole dataset 
(load(dirout("30_9_4_Signatures_nUMI_Cutoff/","Scores.RData")))
table(Dat1$sample)
sigs.cells <- Dat1

# Load Patient clone data
cloneDat <- fread(dirout("30_9_4_Clone_Plots/Patient_Clones_Meta.data.tsv"))

#############
# CLEAN PATIENT NAMES
#############
cleanPatients <- function(v){
  return(gsub("FE", "P1", gsub("PBGY", "P4", gsub("PT", "P5", gsub("VZS", "P7", v)))))
}
cleanCells <- function(v){
  return(gsub("Mono", "CD14", v))
}


#############
# FIGURE CellTypes + COUNTS  #
#############


# FIGURE 4A: Celltype tsne ------------------------------------------------
pDat <- data.table(Celltype = cleanCells(pbmcFull@data.info$CellType), pbmcFull@tsne.rot)
labelCoordinates <- pDat[,.(tSNE_1=median(tSNE_1), tSNE_2=median(tSNE_2)),by="Celltype"]
ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=Celltype)) + geom_point(alpha=0.2) +
  geom_label(data=labelCoordinates, aes(x=tSNE_1, y=tSNE_2, label=Celltype), color="black", alpha=0.5) +
  theme_bw(24) + guides(color=F) + xlim(-40, 40) +
  xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2")
ggsave(dirout(out, "1_CellTypes_Tsne.pdf"), height=7, width=7)
ggsave(dirout(out, "1_CellTypes_Tsne.jpg"), height=7, width=7)
write.table(dim(pbmcFull@data), file=dirout(out, "1_CellTypes_cellCounts.tsv"))
tabDat <- data.table(pbmcFull@data.info)
tabDat[,Patient := cleanPatients(gsub("(\\w+)_.+", "\\1", sample))]
tabDat[,Timepoint := factor(as.numeric(gsub("(\\w+)_d(.+)", "\\2", sample)))]
write.table(with(tabDat, table(Patient, Timepoint)), file=dirout(out, "1_CellTypes_patient_Time_counts.tsv"), sep="\t", quote=F)


# Cell Counts Figure ------------------------------------------------------
pDat <- data.table(pbmcFull@data.info)[CellType %in% celltypes]
pDat[,timepoint := as.numeric(gsub("\\w+_d(\\d+)", "\\1", sample))]
pDat[,Patient := cleanPatients(gsub("(\\w+)_d\\d+", "\\1", sample))]
pDat[,CellType := cleanCells(CellType)]
ggplot(pDat[nUMI > 1000 & nUMI <3000][,.N, by=c("timepoint", "Patient", "CellType")], 
       aes(x=factor(timepoint), y=N, group=paste0(Patient, "_", CellType), fill=CellType)) + 
  geom_bar(stat="identity") +
  facet_grid( . ~ Patient, scales="free", space="free") + theme_bw(24) +
  scale_fill_discrete(name="Cell type") +
  xlab("Time post ibrutinib (days)") + ylab("Cell count")
ggsave(dirout(out, "1_ExploreCounts_bar.pdf"), width=11, height=7)
ggplot(pDat[,.N, by=c("timepoint", "Patient", "CellType")], aes(x=factor(timepoint), y=N, group=paste0(Patient, "_", CellType), fill=CellType)) + 
  geom_bar(stat="identity") +
  facet_grid( . ~ Patient, scales="free", space="free") + theme_bw(24) +
  scale_fill_discrete(name="Cell type") +
  xlab("Time post ibrutinib (days)") + ylab("Cell count")
ggsave(dirout(out, "1_ExploreCounts_bar_BeforeUMICutoff.pdf"), width=11, height=7)



# PATIENT TSNE --------------------------------------------------------------
pDat <- data.table(Patient = cleanPatients(gsub("(\\w+)_.+", "\\1", pbmcFull@data.info$sample)), pbmcFull@tsne.rot)
ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=Patient)) + geom_point(alpha=0.5) +
  theme_bw(24) + xlim(-40, 40) +
  xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2")
ggsave(dirout(out, "1_Patient_Tsne.pdf"), height=7, width=8)
ggsave(dirout(out, "1_Patient_Tsne.jpg"), height=7, width=8)


# FIGURE S2B UMI Counts --------------------------------------------------------------
pDat <- data.table(pbmcFull@data.info)
pDat[,Patient := cleanPatients(gsub("(\\w+)_.+", "\\1", sample))]
pDat[,Timepoint := as.numeric(gsub("(\\w+)_d(.+)", "\\2", sample))]
ggplot(pDat, aes(x=factor(Timepoint), fill=Patient, y=nUMI)) + geom_boxplot(outlier.colour=NA) +
  facet_grid(. ~ Patient, scales="free", space="free") + ylim(0, 5000) +
  xlab("Time post ibrutinib (days)") + ylab("UMI count") + 
  theme_bw(24) + guides(fill=F) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
ggsave(dirout(out, "1_Umi_Count.pdf"), height=7, width=7)

pDat <- data.table(pbmcFull@data.info)[nUMI > 1000 & nUMI < 3000]
pDat[,Patient := cleanPatients(gsub("(\\w+)_.+", "\\1", sample))]
pDat[,Timepoint := as.numeric(gsub("(\\w+)_d(.+)", "\\2", sample))]
ggplot(pDat, aes(x=factor(Timepoint), fill=Patient, y=nUMI)) + geom_boxplot(outlier.colour=NA) +
  facet_grid(. ~ Patient, scales="free", space="free") +
  xlab("Time post ibrutinib (days)") + ylab("UMI count") + 
  theme_bw(24) + guides(fill=F) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
ggsave(dirout(out, "1_Umi_Count_AfterFiltering.pdf"), height=7, width=7)


#############
# 2. DIFF GENES
#############

# Counts of Diff genes ---------------------------------------------------
res <- fread(dirout("13_4_Overtime_nUMI_Cutoff/", "SigGenes_overTime.tsv"))
res.sig <- res[qvalue < 0.05 & abs(logFC) > 0.3 & cellType %in% celltypes]
res.sig[,direction := ifelse(logFC > 0, "up", "down")]
res.sig$Direction <- factor(res.sig$direction, levels = c("up","down"))
res.sig[,Sample := cleanPatients(gsub("d30", "d030", patient))]
res.sig[,cellType := cleanCells(cellType)]
pDat <- res.sig[, .N, by=c("Sample", "Direction", "cellType")]
pDat[Direction == "down", N := -N]
ggplot(pDat, aes(x=Sample, fill=Direction, y=N)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_bw(24) + ylim(-max(abs(pDat$N)), max(abs(pDat$N))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=12)) +
  facet_grid(. ~ cellType, scales="free", space="free") + 
  ylab("Differentially expressed genes") + xlab("Comparison to day 0") +
  guides(fill = FALSE)
ggsave(dirout(out, "2_DiffGenes_NumberOfGenesChanging2.pdf"), height=7, width=8)

res <- fread(dirout("13_4_Overtime_inclDay30/", "SigGenes_overTime.tsv"))
res.sig <- res[qvalue < 0.05 & abs(logFC) > 0.3 & cellType %in% celltypes]
res.sig[,direction := ifelse(logFC > 0, "up", "down")]
res.sig$Direction <- factor(res.sig$direction, levels = c("up","down"))
res.sig[,Sample := cleanPatients(gsub("d30", "d030", patient))]
res.sig[,cellType := cleanCells(cellType)]
pDat <- res.sig[, .N, by=c("Sample", "Direction", "cellType")]
pDat[Direction == "down", N := -N]
ggplot(pDat, aes(x=Sample, fill=Direction, y=N)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_bw(24) + ylim(-max(abs(pDat$N)), max(abs(pDat$N))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=12)) +
  facet_grid(. ~ cellType, scales="free", space="free") + 
  ylab("Differentially expressed genes") + xlab("Comparison to day 0") +
  guides(fill = FALSE)
ggsave(dirout(out, "2_DiffGenes_NumberOfGenesChanging_Before.pdf"), height=7, width=8)


######
# 3 Enrichments
######

# CLL ENRICHMENTS ---------------------------------------------------------
cell <- "CD8"
for(cell in celltypes){
  cll.enr <- fread(dirout("13_4_Overtime_nUMI_Cutoff/",cell,"_Enrich_.tsv"))
  cll.enr <- cll.enr[!grepl("Mus musculus", category)]
  cll.enr[,Term := abbreviate(gsub("_Homo sapiens.+$", "", category), minlength=50)]
  cll.enr[,grp := gsub("d30", "d030", gsub("CLL_", "", grp))]
  cll.enr$Comparison <- cleanPatients(cll.enr[,gsub("_(up|down)", "",grp)])
  cll.enr$Direction <- cll.enr[,gsub(".+?0_([up|down])", "\\1",grp)]
  cll.enr[,Direction := paste0(toupper(substr(Direction, 0,1)), substr(Direction, 2, 10))]
  cll.enr[, mLog10Q := pmin(-log10(qval),10)]
  qvalCutoff <- 0.05
  
  # Figure D: Paths
  cll.enr.path <- cll.enr[database %in% c("NCI-Nature_2016", "Human_Gene_Atlas", "WikiPathways_2016")]
  cll.enr.path.i <- cll.enr.path[qval < qvalCutoff,.N, by=c("Term", "Direction")]
  cll.enr.path.i <- dcast.data.table(cll.enr.path.i, Term ~ Direction, value.var="N")
  cll.enr.path.i[is.na(Down), Down := 0]
  cll.enr.path.i[is.na(Up), Up := 0]
  cll.enr.path.i[,diff := Up - Down]
  ggplot(cll.enr.path[Term %in% tail(cll.enr.path.i[order(abs(diff))], 20)$Term], 
         aes(x=Comparison, y=Term, size=log10(oddsRatio), color=mLog10Q)) + 
    geom_point() + theme_bw(24) + 
    theme(
      axis.text.x = element_text(size=12, angle = 90, hjust=1, vjust=0.5), 
      axis.text.y=element_text(size=12),
      legend.text=element_text(size=18),
      legend.title=element_text(size=18)
      ) + 
    scale_color_gradient(name=expression(log[10](q-value)), limits=c(-log10(qvalCutoff), 10), low="darkred", high="red", na.value="lightgrey") + 
    scale_size_continuous(name=expression(log[10](Oddsratio))) +
    facet_grid(. ~ Direction) + xlab("Comparison to day 0") + ylab("")
  ggsave(dirout(out, "3_",cell,"_PathEnr.pdf"), width=10, height=7)
  
  # Figure 4D: TFs
  cll.enr.path <- cll.enr[database %in% c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "Transcription_Factor_PPIs")]
  cll.enr.path.i <- cll.enr.path[qval < qvalCutoff,.N, by=c("Term", "Direction")]
  cll.enr.path.i <- dcast.data.table(cll.enr.path.i, Term ~ Direction, value.var="N")
  cll.enr.path.i[is.na(Down), Down := 0]
  cll.enr.path.i[is.na(Up), Up := 0]
  cll.enr.path.i[,diff := Up - Down]
  ggplot(cll.enr.path[Term %in% tail(cll.enr.path.i[order(abs(diff))], 20)$Term], 
         aes(x=Comparison, y=Term, size=log10(oddsRatio), color=mLog10Q)) + 
    geom_point() + theme_bw(24) + 
    theme(
      axis.text.x = element_text(size=12, angle = 90, hjust=1, vjust=0.5), 
      axis.text.y=element_text(size=12),
      legend.text=element_text(size=18),
      legend.title=element_text(size=18)
    ) + 
    scale_color_gradient(name=expression(log[10](q-value)), limits=c(-log10(qvalCutoff), 10), low="darkred", high="red", na.value="lightgrey") + 
    scale_size_continuous(name=expression(log[10](Oddsratio))) +
    facet_grid(. ~ Direction, scales="free") +xlab("Comparison to day 0") + ylab("")
  ggsave(dirout(out, "3_",cell,"_TFEnr.pdf"), width=8, height=7)
}


######
# 4 Signatures
######

(load(file=dirout("30_9_4_Signatures_nUMI_Cutoff2/", "OverTime_Tests.RData")))
overTime.sig[,patient := cleanPatients(patient)]
overTime.sig[,CellType := cleanCells(CellType)]

# MAIN SIGNATURES ---------------------------------------------------------
otsig <- overTime.sig[grepl("Bcells", geneset) | geneset == "Ibrutinib_treatment" | grepl("HALLMARK", geneset)]
otsig[,geneset := gsub("HALLMARK_", "", geneset)]
ggplot(otsig[geneset %in% otsig[, sum(abs(Diff) > 1), by="geneset"][V1 >= 2]$geneset], 
       aes(x=paste0(timepoint, "_", patient), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  ylab("") + xlab("") + theme_bw(24) +  
  scale_color_gradient2(name="Effect size", high="red", mid="white", low="blue") +
  scale_size_continuous(name=expression(log[10](q-value))) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=12), axis.text.y=element_text(size=12)) + 
  facet_grid(. ~ CellType, scales="free", space="free") 
ggsave(dirout(out, "4_CLL_Signatures.pdf"),height=7, width=13)

# FIGURE 4F: Ibrutinib signature ------------------------------------------
sigs.cells[,Timepoint := as.numeric(gsub("d", "", timepoint))]
sigs.cells[,Patient := cleanPatients(patient)]
ggplot(sigs.cells, aes(x=factor(Timepoint), y=Ibrutinib_treatment, fill=Patient)) + geom_boxplot(outlier.colour=NA) +
  facet_grid(.~Patient, scales="free", space="free") + ylim(0,6) + 
  xlab("Time post ibrutinib (days)") + ylab("Ibrutinib signature") +
  guides(fill=FALSE) +
  theme_bw(24) + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
ggsave(dirout(out, "4_Ibrutinib_Signature.pdf"), width=6, height=7)

# B cell signatures -------------------------------------------------------
ggplot(overTime.sig[grepl("d\\d", geneset)], 
       aes(x=paste0(patient, "_", timepoint), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  ylab("") + xlab("") + theme_bw(24) +  
  scale_color_gradient2(name="Effect size", high="red", mid="white", low="blue") +
  scale_size_continuous(name=expression(log[10](q-value))) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=12), axis.text.y=element_text(size=12)) + 
  facet_grid(. ~ CellType, scales="free", space="free") 
ggsave(dirout(out, "4_CLL_Patient_Signatures.pdf"),height=6, width=14)

overTime.sig.new <- overTime.sig
(load(dirout("30_9_Signatures_inclDay30_downLists/","OverTime_Tests.RData")))
overTime.sig[,patient := cleanPatients(patient)]
overTime.sig[,CellType := cleanCells(CellType)]
overTime.sig[, mergeName := paste(patient, CellType, timepoint, geneset, sep="_")]
overTime.sig.new[, mergeName := paste(patient, CellType, timepoint, geneset, sep="_")]
pDat <- merge(overTime.sig, overTime.sig.new, by="mergeName")
ggplot(pDat[grepl("HALLMARK", mergeName)], aes(x=Diff.x, y=Diff.y, color=patient.x, shape=CellType.x)) + 
  geom_abline(intercept=0, slope=1, size=2, color="lightgrey") +
  geom_point() + theme_bw(24) + 
  xlab("Effect size before UMI cutoff") + ylab("Effect size after UMI cutoff") +
  xlim(-6.5,6.5) + ylim(-6.5,6.5) +
  scale_color_discrete(name="Patient") +
  scale_shape_discrete(name="Cell type")
ggsave(dirout(out, "4_Signature_Correlation.pdf"), width=7, height=7)
sink(dirout(out, "4_Signature_correlation.txt"))
with(pDat[grepl("HALLMARK", mergeName)], cor.test(Diff.x, Diff.y, method="spearman"))
with(pDat[grepl("HALLMARK", mergeName)], cor.test(Diff.x, Diff.y, method="pearson"))
sink()

# FIGURE S2C: Random signatures -------------------------------------------
(load(dirout("30_9_4_Random/","OverTime_Tests.RData")))
overTime.sig[,patient := cleanPatients(patient)]
overTime.sig[,CellType := cleanCells(CellType)]
overTime.sig[, size := as.numeric(gsub("\\d*_?set_(\\d+)_\\d+", "\\1", geneset))]
ggplot(overTime.sig, aes(y=Diff, x=as.factor(size))) + 
  geom_boxplot() +
  theme_bw(24) +
  xlab("Size of random geneset") + ylab("Difference over time")
ggsave(dirout(out, "4_Random_Signatures.pdf"), height=7, width=7)
(p <- ggplot(overTime.sig, aes(y=Diff, x=paste0(patient, "_", timepoint))) + 
   geom_boxplot() +
   theme_bw(24) +
   xlab("Comparison to time point 0") + ylab("Difference over time") + 
   theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)))
ggsave(dirout(out, "4_Random_Signatures_Patient.pdf"), height=7, width=7)
ggsave(dirout(out, "4_Random_Signatures_CellType.pdf"), plot = p + facet_grid(.~CellType), height=7, width=14)

(load(dirout("30_9_2_Random/","OverTime_Tests.RData")))
overTime.sig[,patient := cleanPatients(patient)]
overTime.sig[,CellType := cleanCells(CellType)]
(p <- ggplot(overTime.sig, aes(y=Diff, x=paste0(patient, "_", timepoint))) + 
   geom_boxplot() +
   theme_bw(24) +
   xlab("Comparison to time point 0") + ylab("Difference over time") + 
   theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)))
ggsave(dirout(out, "4_Random_Signatures_Patient_before.pdf"), height=7, width=7)


#############
# FIGURE S1 #
#############

# FIGURE S1A: FACS counts ------------------------------------------------
pDat <- fread(dirout("14_scRNA_FACS_day30/FACS_Data.tsv"))
pDat[,cellType := cleanCells(cellType)]
ggplot(pDat, aes(x=count*100, y=value, color=cellType)) + geom_point(size=3) + ylab("FACS") + xlab("scRNA") + theme_bw(24) +
  scale_color_discrete(name="Cell type") + xlab("Single cell RNA percentage") + ylab("FACS percentage")
ggsave(dirout(out, "FigureS1A_FACS.pdf"...=), width=8, height=7)
write.table(data.table(
  Spearman = cor(pDat$count, pDat$value, method="spearman"),
  Pearson = cor(pDat$count, pDat$value, method="pearson")
), file=dirout(out, "FigureS1A_Correlation.tsv"), sep="\t")


# FIGURE S1B-D: SPECIFIC MARKER GENES ---------------------------------------------------
for(gg in c("CD79A", "CD3D", "NKG7", "CD14", "FCGR3A", "CD8A")){
  pDat <- data.table(pbmcFull@tsne.rot, Expression = pbmcFull@data[gg,])
  ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=Expression)) + geom_point(alpha = 0.5) +
    scale_color_gradient(low="grey", high = "blue") +
    theme_bw(24) + guides(color=F) + ggtitle(gg) + xlim(-40, 40) +
    xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2")
  ggsave(dirout(out, "FigureS1_",gg,".pdf"), height=7, width=7)
  ggsave(dirout(out, "FigureS1_",gg,".jpg"), height=7, width=7)
}


#############
# 6 Genes in Monocytes
#############

# FIGURE S4B: Upregulation in Monocytes --------------------------------------------
res <- fread(dirout("13_4_Overtime_nUMI_Cutoff/", "SigGenes_overTime.tsv"))
res[,sum(logFC), by=c("cellType", "gene")][cellType == "Mono"][order(V1, decreasing=TRUE)][1:20]
res[,logqval := pmin(5, -log10(qvalue))]
res[,patient := gsub("d30", "d030", patient)]
ggplot(res[(gene %in% c("TNF", "CCL3", "CCL3L3", "CXCL8") | grepl("NFKBI[Z|A]", gene)) & cellType %in% c("CD8", "CD4", "Mono", "CLL")], 
       aes(x=patient, y=gene, size=logqval, color=logFC)) + geom_point() + 
  theme_bw(24) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) + 
  facet_grid(. ~ cellType) + xlab("") + ylab("") +
  scale_color_gradient2(name=expression(log[10](Foldchange)), high="red", mid="white", low="blue") +
  scale_size_continuous(name=expression(log[10](q-value)))
ggsave(dirout(out, "6_MonoGenes.pdf"), height=7, width=18)



#############
# 7 Changes in Clones
#############

(load(dirout("30_9_4_Clone_Plots/", "OverTime_Tests.RData")))
otsig <- overTime.sig[!grepl("d\\d", geneset)]
otsig <- overTime.sig[grepl("Bcells", geneset) | geneset == "Ibrutinib_treatment" | grepl("HALLMARK", geneset)]
otsig[,geneset := gsub("HALLMARK_", "", geneset)]
otsig[,patient := cleanPatients(patient)]
ggplot(otsig[patient %in% c("P4", "P5", "P7") & geneset %in% otsig[, sum(abs(Diff) > 1), by="geneset"][V1 >= 2]$geneset], 
       aes(x=paste0(timepoint, "_", CellType), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  ylab("") + xlab("") + theme_bw(24) +  
  scale_color_gradient2(name="Effect size", high="red", mid="white", low="blue") +
  scale_size_continuous(name=expression(log[10](q-value))) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=12), axis.text.y=element_text(size=12)) +
  facet_grid(. ~ patient, scales="free", space="free") 
ggsave(dirout(out, "7_Clone_Signatures.pdf"),height=7, width=13)


# PLOT ALL THE TSNEs ------------------------------------------------------
# PT Data
pbmcPTs <- list()
for(pt in c("PT", "FE", "PBGY", "VZS")){
  pt2 <- cleanPatients(pt)
  if(is.null(pbmcPTs[[pt]])){
    (load(dirout("25_Patient_CLL/",pt,"/",pt,".RData")))
    pbmcPTs[[pt]] <- pbmc
  } else {
    pbmc <- pbmcPTs[[pt]]
  }
  pDat <- data.table(pbmc@tsne.rot, Timepoint = as.numeric(gsub("\\w+_d(\\d+)", "\\1", pbmc@data.info$sample)))
  ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=factor(Timepoint))) + geom_point(alpha=0.5) +
    scale_color_manual(name="Time point", values=c("blue4", "deepskyblue2", "seagreen3", "gold4")) +
    xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2") + theme_bw(24)
  ggsave(dirout(out, "7_", pt2, "_Clones_TSNE.pdf"), width=9, height=7)
}

# ggplot(otsig[patient %in% c("P1", "P4", "P7") & geneset %in% otsig[, sum(abs(Diff) > 1), by="geneset"][V1 >= 2]$geneset], 
#        aes(x=paste0(timepoint, "_", CellType), y=geneset, color=Diff, size=pvalue2)) + 
#   geom_point() +
#   ylab("") + xlab("") + theme_bw(24) +  
#   scale_color_gradient2(name="Effect size", high="red", mid="white", low="blue") +
#   scale_size_continuous(name=expression(log[10](q-value))) +
#   theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=12), axis.text.y=element_text(size=12)) + 
#   facet_grid(. ~ patient, scales="free", space="free")
# ggsave(dirout(out, "7_Clone_Signatures_Supp.pdf"),height=7, width=10)
