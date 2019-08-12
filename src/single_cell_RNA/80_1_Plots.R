require("project.init")
require(Seurat)
require(methods)
require(monocle)

project.init2("cll-time_course")
out <- "80_1_Plots/"
dir.create(dirout(out))

source("src/single_cell_RNA/FUNC_Enrichr.R")

celltypes <- c("CD4", "CD8", "Mono", "CLL")



# LOAD DATA ---------------------------------------------------------------

# Load full single cell dataset
(load(dirout("10_Seurat_raw/inclDay30_noIGHLK_negbinom/inclDay30_noIGHLK.RData")))
pbmcFull <- pbmc
rm(list="pbmc")
pbmcFull@scale.data <- NULL
pbmcFull@raw.data <- NULL

# Load CLL single cells
(load(dirout("11_CellTypes_inclDay30/CLL/CLL.RData")))
pbmcCLL <- pbmc
rm(list="pbmc")
pbmcCLL@scale.data <- NULL
pbmcCLL@raw.data <- NULL

# PT Data
(load(dirout("25_Patient_CLL/PT/PT.RData")))
pbmcPT <- pbmc
rm(list.files="pbmc")
pbmcPT@scale.data <- NULL
pbmcPT@raw.data <- NULL

# PBGY Data
(load(dirout("25_Patient_CLL/PBGY/PBGY.RData")))
pbmcPBGY <- pbmc
rm(list.files="pbmc")
pbmcPBGY@scale.data <- NULL
pbmcPBGY@raw.data <- NULL

# Load Monocle
(load(dirout("52_3_Monocle_inclDay30_Endpoints/", "MonocleDat.RData")))

# Load Signatures for whole dataset 
(load(dirout("30_9_Signatures_inclDay30_downLists2/","Scores.RData")))
table(Dat1$sample)
sigs.cells <- Dat1

# Load Patient clone data
cloneDat <- fread(dirout("30_9_Plots/Patient_Clones_Meta.data.tsv"))



#############
# FIGURE 4  #
#############


# FIGURE 4A: Celltype tsne ------------------------------------------------
pDat <- data.table(Celltype = pbmcFull@data.info$CellType, pbmcFull@tsne.rot)
labelCoordinates <- pDat[,.(tSNE_1=median(tSNE_1), tSNE_2=median(tSNE_2)),by="Celltype"]
ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=Celltype)) + geom_point(alpha=0.5) +
  geom_label(data=labelCoordinates, aes(x=tSNE_1, y=tSNE_2, label=Celltype), color="black", alpha=0.5) +
  theme_bw(24) + guides(color=F) + xlim(-40, 40) +
  xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2")
ggsave(dirout(out, "Figure4A_celltype_Tsne.pdf"), height=7, width=7)
ggsave(dirout(out, "Figure4A_celltype_Tsne.jpg"), height=7, width=7)
write.table(dim(pbmcFull@data), file=dirout(out, "Figure4A_cellCounts.tsv"))
tabDat <- data.table(pbmcFull@data.info)
tabDat[,Patient := gsub("(\\w+)_.+", "\\1", sample)]
tabDat[,Timepoint := factor(as.numeric(gsub("(\\w+)_d(.+)", "\\2", sample)))]
write.table(with(tabDat, table(Patient, Timepoint)), file=dirout(out, "Figure4A_patient_Time_counts.tsv"), sep="\t", quote=F)

# FIGURE 4B: Diff genes ---------------------------------------------------
res <- fread(dirout("13_4_Overtime_inclDay30/", "SigGenes_overTime.tsv"))
res.sig <- res[qvalue < 0.05 & abs(logFC) > 0.3]
res.sig[,direction := ifelse(logFC > 0, "up", "down")]
res.sig$Direction <- factor(res.sig$direction, levels = c("up","down"))
res.sig[,Sample := gsub("d30", "d030", patient)]
ggplot(res.sig[cellType %in% celltypes], aes(x=Sample, fill=Direction)) + geom_bar(position="dodge") + 
  theme_bw(24) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=12)) +
  facet_grid(. ~ cellType) + ylab("Number of genes") + xlab("Comparison to day 0")
ggsave(dirout(out, "Figure4B_NumberOfGenesChanging.pdf"), height=7, width=10)

# FIGURE 4C: MONOCLE ----------------------------------------------------
pDat <- data.table(pData(mcle))
pDat[,sample := gsub("d0", "d000", sample)]
pDat[,sample := gsub("d30", "d030", sample)]
pDat[,Timepoint := as.numeric(gsub("\\w+_d", "", sample))]
pDat[,Patient := patient]
ggplot(pDat, aes(y=Pseudotime, x=factor(Timepoint), fill=Patient)) + geom_violin(color=NA) + 
  xlab("Time point (days)") +
  facet_grid(. ~ Patient, scales="free",space="free_x") + theme_bw(24) +
  #theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  guides(fill=FALSE)
ggsave(dirout(out, "Figure4C_Monocle.pdf"), height=7, width=10)

# FIG 4D + 4E: CLL ENRICHMENTS ---------------------------------------------------------
cll.enr <- fread(dirout(("13_4_Overtime_inclDay30/CLL_Enrich_.tsv")))
cll.enr[,Term := gsub("(Homo sapiens|Mus musculus).+$", "\\1", category)]
cll.enr[,grp := gsub("d30", "d030", gsub("CLL_", "", grp))]
cll.enr$Comparison <- cll.enr[,gsub("_(up|down)", "",grp)]
cll.enr$Direction <- cll.enr[,gsub(".+?0_([up|down])", "\\1",grp)]
cll.enr[,Direction := paste0(toupper(substr(Direction, 0,1)), substr(Direction, 2, 10))]
cll.enr[, mLog10Q := pmin(-log10(qval),10)]

# Figure D: Paths
cll.enr.path <- cll.enr[database %in% c("NCI-Nature_2016", "Human_Gene_Atlas", "WikiPathways_2016")]
cll.enr.path.i <- cll.enr.path[qval < 0.05,.N, by=c("Term", "Direction")]
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
  scale_color_gradient(name=expression(log[10](q-value)), limits=c(2, 10), low="pink", high="darkred", na.value="lightgrey") + 
  scale_size_continuous(name=expression(log[10](Oddsratio))) +
  facet_grid(. ~ Direction, scales="free") + xlab("Comparison to day 0") + ylab("")
ggsave(dirout(out, "Figure4D_PathEnr.pdf"), width=10, height=7)

# Figure 4D: TFs
cll.enr.path <- cll.enr[database %in% c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "Transcription_Factor_PPIs")]
cll.enr.path.i <- cll.enr.path[qval < 0.05,.N, by=c("Term", "Direction")]
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
  scale_color_gradient(name=expression(log[10](q-value)), limits=c(2, 10), low="pink", high="darkred", na.value="lightgrey") + 
  scale_size_continuous(name=expression(log[10](Oddsratio))) +
  facet_grid(. ~ Direction, scales="free") +xlab("Comparison to day 0") + ylab("")
ggsave(dirout(out, "Figure4E_TFEnr.pdf"), width=8, height=7)

# FIGURE 4F: Ibrutinib signature ------------------------------------------
sigs.cells[,Timepoint := as.numeric(gsub("d", "", timepoint))]
sigs.cells[,Patient := patient]
ggplot(sigs.cells, aes(x=factor(Timepoint), y=Ibrutinib_treatment, fill=Patient)) + geom_boxplot(outlier.colour=NA) +
  facet_grid(.~Patient) + ylim(0,7) + 
  xlab("Time point (days)") + ylab("Ibrutinib signature") +
  guides(fill=FALSE) +
  theme_bw(24) + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
ggsave(dirout(out, "Figure4F_Ibrutinib_Signature.pdf"), width=7, height=7)



#############
# FIGURE 5  #
#############


# FIGURE 5A: Patient t-SNEs -----------------------------------------------
pDat <-data.table(pbmcPT@tsne.rot, sample = pbmcPT@data.info$sample, Patient = "PT")
pDat <-rbind(pDat, data.table(pbmcPBGY@tsne.rot, sample = pbmcPBGY@data.info$sample, Patient = "PBGY"))
pDat[,Timepoint := gsub("\\w+_d", "Day ", sample)]
pDat$Timepoint <- factor(pDat$Timepoint, levels=unique(pDat$Timepoint[order(as.numeric(gsub("Day ", "", pDat$Timepoint)))]))
ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color= Timepoint)) + geom_point(alpha=0.5) + 
  scale_color_discrete(name="Time point") +
  theme_bw(24) + xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2") +
  facet_grid(. ~ Patient) +
  guides(colour = guide_legend(override.aes=list(alpha=1)))
ggsave(dirout(out, "Figure5A_Pat_tsne.pdf"), width=18, height=7)
ggsave(dirout(out, "Figure5A_Pat_tsne.jpg"), width=18, height=7)

# FIGURE 5B: Signatures in Clones
pDat <- merge(sigs.cells, cloneDat[,c("rn", "Clone"), with=F], by="rn")
pDat2[,sample_cell := paste(patient, Clone, timepoint, sep="_")]

#############
# FIGURE S1 #
#############

# FIGURE S1A: FACS counts ------------------------------------------------
pDat <- fread(dirout("14_scRNA_FACS_day30/FACS_Data.tsv"))
ggplot(pDat, aes(x=count*100, y=value, color=cellType)) + geom_point(size=3) + ylab("FACS") + xlab("scRNA") + theme_bw(24) +
  scale_color_discrete(name="Cell type") + xlab("Single cell RNA percentage") + ylab("FACS percentage")
ggsave(dirout(out, "FigureS1A_FACS.pdf"), width=8, height=7)
write.table(data.table(
  Spearman = cor(pDat$count, pDat$value, method="spearman"),
  Pearson = cor(pDat$count, pDat$value, method="pearson")
), file=dirout(out, "FigureS1A_Correlation.tsv"), sep="\t")


# FIGURE S1B-D: SPECIFIC MARKER GENES ---------------------------------------------------
for(gg in c("CD79A", "CD3D", "NKG7", "CD14", "FCGR3A")){
  pDat <- data.table(pbmcFull@tsne.rot, Expression = pbmcFull@data[gg,])
  ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=Expression)) + geom_point(alpha = 0.5) +
    scale_color_gradient(low="grey", high = "blue") +
    theme_bw(24) + guides(color=F) + ggtitle(gg) + xlim(-40, 40) +
    xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2")
  ggsave(dirout(out, "FigureS1_",gg,".pdf"), height=7, width=7)
  ggsave(dirout(out, "FigureS1_",gg,".jpg"), height=7, width=7)
}



#############
# FIGURE S2 #
#############


# FIGURE S2A: PATIENT TSNE --------------------------------------------------------------
pDat <- data.table(Patient = gsub("(\\w+)_.+", "\\1", pbmcFull@data.info$sample), pbmcFull@tsne.rot)
ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=Patient)) + geom_point(alpha=0.5) +
  theme_bw(24) + xlim(-40, 40) +
  xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2")
ggsave(dirout(out, "FigureS2A_Patient_Tsne.pdf"), height=7, width=8)
ggsave(dirout(out, "FigureS2A_Patient_Tsne.jpg"), height=7, width=8)


# FIGURE S2B UMI Counts --------------------------------------------------------------
pDat <- data.table(pbmcFull@data.info)
pDat[,Patient := gsub("(\\w+)_.+", "\\1", sample)]
pDat[,Timepoint := as.numeric(gsub("(\\w+)_d(.+)", "\\2", sample))]
ggplot(pDat, aes(x=factor(Timepoint), fill=Patient, y=nUMI)) + geom_boxplot(outlier.colour=NA) +
  facet_grid(. ~ Patient) + ylim(0, 5000) +
  xlab("Time point (days)") + ylab("UMI count") + 
  theme_bw(24) + guides(fill=F) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
ggsave(dirout(out, "FigureS2B_Umi_Count.pdf"), height=7, width=7)


# FIGURE S2C: Random signatures -------------------------------------------
(load(dirout("30_9_2_Random/","OverTime_Tests.RData")))
overTime.sig[, size := as.numeric(gsub("\\d*_?set_(\\d+)_\\d+", "\\1", geneset))]
ggplot(overTime.sig, aes(y=Diff, x=as.factor(size))) + geom_violin(fill="lightgrey") + theme_bw(24) +
  xlab("Size of random geneset") + ylab("Difference over time")
ggsave(dirout(out, "FigureS2C_Random_Signatures.pdf"), height=7, width=7)
ggplot(overTime.sig, aes(y=Diff, x=paste0(patient, "_", timepoint))) + geom_violin(fill="lightgrey") + theme_bw(24) +
  xlab("Size of random geneset") + ylab("Difference over time") + facet_grid(. ~ CellType) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
ggsave(dirout(out, "FigureS2C_Random_Signatures_CellType.pdf"), height=7, width=15)


# FIGURE S2C: Monocle Details ---------------------------------------------
pData(mcle)$Timepoint <- factor(as.numeric(gsub("150", "120", gsub("\\w+_d(\\d+)", "\\1", pData(mcle)$sample))))
levels(pData(mcle)$Timepoint)[3] <- "120/150"
plot_cell_trajectory(mcle, color_by = "patient",show_branch_points=FALSE) + ylim(-2,2) +
  facet_grid(Timepoint ~ patient) + guides(color=FALSE) + theme_bw(24) + coord_flip()
ggsave(dirout(out, "FigureS2D_MonocleDetails.pdf"), width=15, height=14)

#############
# FIGURE S3 #
#############


# FIGURE S3A: Clones over time --------------------------------------------
pat.md <- cloneDat
pat.md[,timepoint := as.numeric(gsub("\\w+_d(\\d+)", "\\1", sample))]
pDat <- pat.md[,.N, by=c("patient", "timepoint", "Clone")]
pDat[,sumN := sum(N), by=c("patient", "timepoint")]
pDat[,maxClone := max(N), by=c("patient", "Clone")]
pDat[,fraction := N / sumN * 100]
pDat[, Patient := patient]
ggplot(pDat, aes(x=timepoint, y=fraction, color=Patient, linetype=Clone, group=paste0(Clone, "_", Patient))) + 
  geom_line(size=2) + ylim(0,100) + ylab("Clone count (%)") + xlab("Timepoint (days)") + theme_bw(24) +
  xlim(0,300)
ggsave(dirout(out, "FigureS3A_CLL_Clones.pdf"), height=7, width=10)


#############
# FIGURE S4 #
#############


# FIGURE S4A: Signatures in other cell types --------------------------------------------
(load(file=dirout("30_9_Signatures_inclDay30_downLists/", "OverTime_Tests.RData")))
ggplot(overTime.sig[grepl("HALLMARK", geneset) | geneset == "Ibrutinib_treatment"], 
  aes(x=paste0(patient, "_", timepoint), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  ylab("") + xlab("") + theme_bw(24) +  
  scale_color_gradient2(name="Effect size", high="red", mid="white", low="blue") +
  scale_size_continuous(name=expression(log[10](q-value))) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=12), axis.text.y=element_text(size=12)) + facet_grid(. ~ CellType)
ggsave(dirout(out, "FigureS4A_Signatures_Cells.pdf"),height=14, width=18)

# FIGURE S4B: Upregulation in Monocytes --------------------------------------------
res <- fread(dirout("13_4_Overtime_","inclDay30","/", "SigGenes_overTime.tsv"))
res[,sum(logFC), by=c("cellType", "gene")][cellType == "Mono"][order(V1, decreasing=TRUE)][1:20]
res[,logqval := pmin(5, -log10(qvalue))]
res[,patient := gsub("d30", "d030", patient)]
ggplot(res[(gene %in% c("TNF", "CCL3", "CCL3L3", "CXCL8") | grepl("NFKBI[Z|A]", gene)) & cellType %in% c("CD8", "CD4", "Mono", "CLL")], 
       aes(x=patient, y=gene, color=logqval, size=logFC)) + geom_point() + 
  theme_bw(24) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) + 
  facet_grid(. ~ cellType) + xlab("") + ylab("") +
  scale_color_gradient2(name=expression(log[10](q-value)), high="red", mid="white", low="blue") +
  scale_size_continuous(name=expression(log[10](Foldchange))) +
ggsave(dirout(out, "FigureS4B_MonoGenes.pdf"), height=7, width=18)
