require("project.init")
require(Seurat)
require(methods)
require(monocle)
require(gridExtra)
project.init2("cll-time_course")
out <- "80_3_Plots/"
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
  #pbmcFull@raw.data <- NULL
}

# Load CLL single cells
if(FALSE){
  (load(dirout("11_2_3_CellTypes_UMI_Cutoff/CLL/CLL.RData")))
  pbmcCLL <- pbmc
  rm(list="pbmc")
  pbmcCLL@scale.data <- NULL
  #pbmcCLL@raw.data <- NULL
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
#   theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) + 
  xlab("Time post ibrutinib (days)") + ylab("Cell count")
ggsave(dirout(out, "1_ExploreCounts_bar_BeforeUMICutoff.pdf"), width=11, height=7)


# PATIENT TSNE --------------------------------------------------------------
pDat <- data.table(
  Patient = cleanPatients(gsub("(\\w+)_.+", "\\1", pbmcFull@data.info$sample)), 
  pbmcFull@tsne.rot,
  Timepoint = as.numeric(gsub("[A-Z]+\\_d(\\d+)", "\\1",pbmcFull@data.info$sample))
  )
pDat[Timepoint == 0, color :="blue4"]
pDat[Timepoint == 30, color :="seagreen3"]
pDat[Timepoint == 120, color :="deepskyblue2"]
pDat[Timepoint == 150, color :="deepskyblue2"]
pDat[Timepoint == 280, color :="gold4"]
pDat <- pDat[order(Timepoint)]
pDat$Timepoint <- factor(as.character(pDat$Timepoint), levels=as.character(unique(sort(pDat$Timepoint))))
ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=Patient)) + geom_point(alpha=0.5) +
  theme_bw(24) + xlim(-40, 40) +
  xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2")
ggsave(dirout(out, "1_Patient_Tsne.pdf"), height=7, width=8)
ggsave(dirout(out, "1_Patient_Tsne.jpg"), height=7, width=8)
colDat <- unique(pDat[,c("color", "Timepoint"), with=T])
coll.TIME <- colDat$color
names(coll.TIME) <- as.character(colDat$Timepoint)
ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=Timepoint)) + geom_point(alpha=0.5) +
  theme_bw(24) + xlim(-40, 40) +
  scale_color_manual(values=coll.TIME) +
  xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2")
ggsave(dirout(out, "1_Timepoint_Tsne.pdf"), height=7, width=8)
ggsave(dirout(out, "1_Timepoint_Tsne.jpg"), height=7, width=8)

ptx <- "P1"
gPlots <- list()
for(ptx in unique(pDat$Patient)){
  gPlots[[ptx]] <- ggplot(pDat, aes(x=tSNE_1, y=tSNE_2)) + 
    geom_point(color="grey") + 
    geom_point(data=pDat[Patient == ptx], aes(color=Timepoint), alpha=0.2)  +
    scale_color_manual(values=coll.TIME) +
    theme_bw(12) + xlim(-40,40) + ylim(-40,40) + 
    xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2") + 
    ggtitle(paste("Patient", ptx)) +
    theme(plot.title=element_text(size=24))
}
(p <- grid.arrange(grobs=lapply(gPlots[sort(names(gPlots))], function(p) p + guides(color=FALSE)), ncol=2))
# ggsave(dirout(out, "1_Patient_Tsne2.pdf"), height=7, width=7, plot=p)
ggsave(dirout(out, "1_Patient_Tsne2.jpg"), height=7, width=7, plot=p)
dev.off()

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


##########################
# FIGURE 1_1 MARKER GENES #
##########################

# FIGURE S1A: FACS counts ------------------------------------------------
pDat <- fread(dirout("14_scRNA_FACS_day30/FACS_Data.tsv"))
pDat[,cellType := cleanCells(cellType)]
(facsP <- ggplot(pDat, aes(x=count*100, y=value, color=cellType)) + geom_point(size=3) + ylab("FACS") + xlab("scRNA") +
   scale_color_discrete(name="Cell type") + xlab("Single cell RNA percentage") + ylab("FACS percentage"))
ggsave(dirout(out, "1_1_FACS.pdf"), width=8, height=7, plot= facsP + theme_bw(24))
write.table(data.table(
  Spearman = cor(pDat$count, pDat$value, method="spearman"),
  Pearson = cor(pDat$count, pDat$value, method="pearson")
), file=dirout(out, "1_1_Correlation.tsv"), sep="\t")


# FIGURE S1B-D: SPECIFIC MARKER GENES ---------------------------------------------------
gPlots <- list()
for(gg in c("CD79A", "CD3D", "NKG7", "CD14", "FCGR3A")){
  pDat <- data.table(pbmcFull@tsne.rot, Expression = pbmcFull@data[gg,])
  gPlots[[gg]] <- ggplot(pDat, aes(x=tSNE_1, y=tSNE_2)) + 
    geom_point(alpha = 0.5, aes(color=Expression)) +
    #     geom_point(data=rbind(pDat[Expression == 0][1], pDat[Expression > 0]), aes(color=Expression)) +
    scale_color_gradient(low="grey", high = "blue") +
    theme_bw(12) + guides(color=F) + ggtitle(gg) + xlim(-40, 40) +
    theme(plot.title=element_text(size=24)) +
    xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2")
  ggsave(dirout(out, "1_1_",gg,".pdf"), height=7, width=7)
  ggsave(dirout(out, "1_1_",gg,".jpg"), height=7, width=7)
}
gPlots[["FACS"]] <- facsP + theme_bw(12) + guides(color=F) + ggtitle("  ") + theme(plot.title=element_text(size=24))
(p <- grid.arrange(grobs=gPlots, ncol=3))
# ggsave(dirout(out, "1_Patient_Tsne2.pdf"), height=7, width=7, plot=p)
ggsave(dirout(out, "1_1_Markers.jpg"), height=7, width=10, plot=p)



# AVERAGES OF GENES -------------------------------------------------------
markerGenes <- c("CD79A", "CD3D","CD3G","CD3A","CD3B","CD4", "CD8A","CD19","CD24","CD68","CST3","GZMB","IL32","TCL1A", "NKG7", "CD14", "FCGR3A", "NCR1", "TRDC","CCR7")
markerGenes <- markerGenes[markerGenes %in% row.names(pbmcFull@data)]
pDat <- data.table(pbmcFull@data.info[,c("sample", "CellType")], as.matrix(t(pbmcFull@data[markerGenes,])))
pDat[,sample := gsub("30", "030", cleanPatients(sample))]
pDat[,sample := gsub("d0$", "d000", cleanPatients(sample))]
pDat[,CellType := cleanCells(CellType)]
pDat[CellType == "NurseLikeCell",CellType := "Nurse"]
pDat <- pDat[!is.na(CellType)]
pDat <- data.table(melt(pDat, id.vars=c("sample", "CellType")))
pDat <- pDat[,.(Expression = mean(value)), by=c("sample", "CellType", "variable")]
pDat[,Expression := Expression - min(Expression,na.rm=T), by=c("variable")]
pDat[,Expression := Expression / max(Expression, na.rm=T), by=c("variable")]
orMT <- t(as.matrix(dcast.data.table(pDat, paste0(sample, CellType) ~ variable, value.var="Expression")[,-"sample",with=F]))
orMT[is.na(orMT)] <- 1
hclustObj <- hclust(dist(orMT))
pDat$variable <- factor(pDat$variable, levels=hclustObj$labels[hclustObj$order])
ggplot(pDat, aes(x=sample, y=variable, fill=Expression)) + geom_tile() + facet_grid(. ~ CellType,scales="free", space="free") +
  theme(strip.background=element_rect(fill="grey")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6)) + guides(fill=FALSE) + xlab("") + ylab("") +
  scale_fill_gradient(low="lightgrey", high="blue")
ggsave(dirout(out, "1_1_MarkerMeans.pdf"), height=7, width=7)


# AVERAGES OF GENES (TFs) -------------------------------------------------------
tfs <- fread(dirout("18_1_PlotSigGeneDetails/TFs.tsv"), header=F)$V1
tfs <- tfs[tfs %in% row.names(pbmcFull@data)]
pDat <- data.table(pbmcFull@data.info[,c("sample", "CellType")], as.matrix(t(pbmcFull@data[tfs,])))
pDat[,sample := gsub("30", "030", cleanPatients(sample))]
pDat[,sample := gsub("d0$", "d000", cleanPatients(sample))]
pDat[,patient := gsub("_d\\d+", "", cleanPatients(sample))]
pDat[,CellType := cleanCells(CellType)]
pDat <- pDat[!is.na(CellType)]
pDat <- data.table(melt(pDat, id.vars=c("sample", "CellType")))
pDat <- pDat[,.(Expression = mean(value)), by=c("sample", "CellType", "variable")]
pDat[,Expression := Expression - min(Expression,na.rm=T), by=c("variable")]
pDat[,Expression := Expression / max(Expression, na.rm=T), by=c("variable")]
ggplot(pDat[CellType == "CLL"], aes(x=sample, y=variable, fill=Expression)) + geom_tile() + facet_grid(. ~ patient,scales="free", space="free") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=12)) + guides(fill=FALSE)
ggsave(dirout(out, "1_1_TFs_HM.pdf"), height=5, width=10)


# PLOT SPECIFIC TFs -------------------------------------------------------
pDat <- data.table(pbmcFull@data.info[,c("sample", "CellType")], as.matrix(t(pbmcFull@data[c("NFATC1", "IKZF1"),])))
pDat[,sample := gsub("30", "030", sample)]
pDat <- pDat[CellType == "CLL" & grepl("PT", sample)][,-"CellType",with=F]
(pDat <- melt(pDat,id.vars=c("sample")))
ggplot(pDat, aes(x=sample, y=value)) + geom_jitter(width=0.1, height=0,alpha=0.2) + geom_violin(aes(fill=sample),alpha=0.5) + facet_grid(.~variable)
# ggsave(dirout(out, "1_1_Tf_Expression_Scaled.pdf"), height=5, width=10)
ggsave(dirout(out, "1_1_Tf_Expression_Unscaled.pdf"), height=5, width=10)


metaDatx <- subset(pbmcFull@data.info, CellType == "CLL" & grepl("PT", sample))
pDat <- data.table(metaDatx, t(as.matrix(pbmcFull@raw.data[c("NFATC1", "IKZF1"),row.names(metaDatx)])))
pDat[,sample := gsub("30", "030", sample)]
pDat <- pDat[,c("sample", "NFATC1", "IKZF1"),with=F]
(pDat <- melt(pDat,id.vars=c("sample")))
ggplot(pDat, aes(x=sample, y=value)) + geom_jitter(width=0.2, height=0,alpha=0.2) + geom_violin(aes(fill=sample),alpha=0.5) + facet_grid(.~variable)
ggsave(dirout(out, "1_1_Tf_Expression_Raw.pdf"), height=5, width=10)

normDatOrig <- as.matrix(pbmcFull@data[,row.names(metaDatx)])
rm(list="pbmcFull")
normDat <- normDatOrig
for(col in 1:ncol(normDat)){
  normDat[,col] <- rank(normDat[,col])
}
pDat <- data.table(metaDatx, t(normDat[c("NFATC1", "IKZF1"),]))
pDat[,sample := gsub("30", "030", sample)]
pDat <- pDat[,c("sample", "NFATC1", "IKZF1"),with=F]
(pDat <- melt(pDat,id.vars=c("sample")))
ggplot(pDat, aes(x=sample, y=value)) + geom_jitter(width=0.2, height=0,alpha=0.2) + geom_violin(aes(fill=sample),alpha=0.5) + facet_grid(.~variable)
ggsave(dirout(out, "1_1_Tf_Expression_Ranks.pdf"), height=5, width=10)

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
cell <- "CLL"
for(cell in celltypes){
  cll.enr <- fread(dirout("13_4_Overtime_nUMI_Cutoff/",cell,"_Enrich_.tsv"))
  cll.enr <- cll.enr[!grepl("Mus musculus", category)]
  cll.enr[,Term := abbreviate(gsub("_Homo sapiens.+$", "", category), minlength=50)]
  cll.enr[,grp := gsub("d30", "d030", gsub("CLL_", "", grp))]
  cll.enr$Comparison <- cleanPatients(cll.enr[,gsub("_(up|down)", "",grp)])
  cll.enr$Direction <- cll.enr[,gsub(".+?0_([up|down])", "\\1",grp)]
  cll.enr[,Direction := paste0(toupper(substr(Direction, 0,1)), substr(Direction, 2, 10))]
  cll.enr[, mLog10Q := pmin(-log10(qval),10)]
  cll.enr[,patient := gsub("_d\\d+", "", Comparison)]
  qvalCutoff <- 0.05
  
  if(length(unique(cll.enr$Term)) >= 2){
    try({
      orMT <- t(as.matrix(dcast.data.table(cll.enr, grp ~ Term, value.var="oddsRatio",fun.aggregate=mean)[,-"grp",with=F]))
      orMT[is.na(orMT)] <- 1
      hclustObj <- hclust(dist(orMT))
      cll.enr$Term <- factor(cll.enr$Term, levels=hclustObj$labels[hclustObj$order])
    }, silent=T)
  }
  
  # Figure D: Paths
  cll.enr.path <- cll.enr[database %in% c("NCI-Nature_2016", "Human_Gene_Atlas", "WikiPathways_2016")]
  #cll.enr.path.i <- cll.enr.path[qval < qvalCutoff,.N, by=c("Term", "Direction")]
  cll.enr.path.i <- cll.enr.path[qval < qvalCutoff,.(N=length(unique(patient))), by=c("Term", "Direction")]
  cll.enr.path.i <- dcast.data.table(cll.enr.path.i, Term ~ Direction, value.var="N")
  cll.enr.path.i[is.na(Down), Down := 0]
  cll.enr.path.i[is.na(Up), Up := 0]
  cll.enr.path.i[,diff := Up - Down]
  ggplot(cll.enr.path[Term %in% tail(cll.enr.path.i[order(abs(diff))], 30)$Term], 
         aes(x=Comparison, y=Term, size=log10(oddsRatio), color=mLog10Q)) + 
    geom_point() + theme_bw(24) + 
    theme(
      axis.text.x = element_text(size=12, angle = 90, hjust=1, vjust=0.5), 
      axis.text.y=element_text(size=12),
      legend.text=element_text(size=18),
      legend.title=element_text(size=18)
      ) + 
    scale_color_gradient(name=expression(log[10](q-value)), limits=c(-log10(qvalCutoff), 10), low="darkblue", high="limegreen", na.value="lightgrey") + 
    scale_size_continuous(name=expression(log[10](Oddsratio))) +
    facet_grid(. ~ Direction) + xlab("Comparison to day 0") + ylab("")
  ggsave(dirout(out, "3_",cell,"_PathEnr.pdf"), width=10, height=10)
  
  # Figure 4D: TFs
  cll.enr.path <- cll.enr[database %in% c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "Transcription_Factor_PPIs")]
  #cll.enr.path.i <- cll.enr.path[qval < qvalCutoff,.N, by=c("Term", "Direction")]
  cll.enr.path.i <- cll.enr.path[qval < qvalCutoff,.(N=length(unique(patient))), by=c("Term", "Direction")]
  cll.enr.path.i <- dcast.data.table(cll.enr.path.i, Term ~ Direction, value.var="N")
  cll.enr.path.i[is.na(Down), Down := 0]
  cll.enr.path.i[is.na(Up), Up := 0]
  cll.enr.path.i[,diff := Up - Down]
  ggplot(cll.enr.path[Term %in% tail(cll.enr.path.i[order(abs(diff))], 30)$Term], 
         aes(x=Comparison, y=Term, size=log10(oddsRatio), color=mLog10Q)) + 
    geom_point() + theme_bw(24) + 
    theme(
      axis.text.x = element_text(size=12, angle = 90, hjust=1, vjust=0.5), 
      axis.text.y=element_text(size=12),
      legend.text=element_text(size=18),
      legend.title=element_text(size=18)
    ) + 
    scale_color_gradient(name=expression(log[10](q-value)), limits=c(-log10(qvalCutoff), 10), low="darkblue", high="limegreen", na.value="lightgrey") + 
    scale_size_continuous(name=expression(log[10](Oddsratio))) +
    facet_grid(. ~ Direction, scales="free") +xlab("Comparison to day 0") + ylab("")
  ggsave(dirout(out, "3_",cell,"_TFEnr.pdf"), width=8, height=10)
}


######
# 4 Signatures
######

(load(file=dirout("30_9_4_Signatures_nUMI_Cutoff_TFs/", "OverTime_Tests.RData")))
overTime.sig[,patient := cleanPatients(patient)]
overTime.sig[,CellType := cleanCells(CellType)]
overTime.sig[,qvalue := p.adjust(pvalue, method='BH')]
overTime.sig[,pvalue2 := pmin(5, -log10(qvalue))]

# MAIN SIGNATURES ---------------------------------------------------------
otsig <- overTime.sig[!grepl("_d\\d", geneset) & !grepl("Fereira", geneset)]#[grepl("Bcells", geneset) | geneset == "Ibrutinib_treatment" | grepl("HALLMARK", geneset)]
otsig[,geneset := gsub("HALLMARK_", "", geneset)]
# Selected
ggplot(otsig[geneset %in% otsig[, sum(abs(Diff) > 1), by="geneset"][V1 >= 2]$geneset], 
       aes(x=paste0(timepoint, "_", patient), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  ylab("") + xlab("") + theme_bw(24) +  
  scale_color_gradient2(name="Effect size", high="red", mid="white", low="blue") +
  scale_size_continuous(name=expression(log[10](q-value))) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=12), axis.text.y=element_text(size=12)) + 
  facet_grid(. ~ CellType, scales="free", space="free") 
ggsave(dirout(out, "4_CLL_Signatures2.pdf"),height=7, width=13)
# All
ggplot(otsig, 
       aes(x=paste0(timepoint, "_", patient), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  ylab("") + xlab("") + theme_bw(24) +  
  scale_color_gradient2(name="Effect size", high="red", mid="white", low="blue") +
  scale_size_continuous(name=expression(log[10](q-value))) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=12), axis.text.y=element_text(size=12)) + 
  facet_grid(. ~ CellType, scales="free", space="free") 
ggsave(dirout(out, "4_CLL_Signatures2_All.pdf"),height=15, width=13)


# BCELL SIGNATURE ---------------------------------------------------------
pDat <- sigs.cells
ggplot(pDat[grepl("PT", sample) & CellType == "CLL"], 
       aes(x=Bcells_Human_Gene_Atlas, Bcells_NCI.Nature_2016, color=sample)) + 
  stat_ellipse(alpha=0.5, aes(fill=sample), geom="polygon") +geom_point(alpha=0.5)
ggsave(dirout(out, "4_B_CellSignature.pdf"))


pDat <- sigs.cells[grepl("PT", sample) & CellType == "CLL"][,c("sample", "Bcells_NCI.Nature_2016", "Bcells_WikiPathways_2016", "Bcells_Human_Gene_Atlas")]
pDat[,sample := gsub("PT_d", "day ", sample)]
pDat <- melt(pDat, id.vars="sample")
pDat[,variable := gsub("WikiPathways", "Wiki", gsub("Human ", "", gsub("(\\.|_)", " ", gsub("_20\\d\\d", "", gsub("Bcells_", "", variable)))))]
pDat[,score := scale(value), by="variable"]
(p <- ggplot(pDat,aes(x=score, color=sample,group=sample)) + geom_density(size=1) + facet_grid(variable ~ .) + xlim(-4,4) +
  xlab("Signature") + ylab("Density") + theme_bw(24) + labs(color="P5, CLL"))
ggsave(dirout(out, "4_B_CellSignature2.pdf"),width=8, height=7, plot=p)
ggsave(dirout(out, "4_B_CellSignature2_noGuide.pdf"),width=8, height=7, plot=p + guides(color=FALSE))


# FIGURE 4F: Ibrutinib signature ------------------------------------------
sigs.cells[,Timepoint := as.numeric(gsub("d", "", timepoint))]
sigs.cells[,Patient := cleanPatients(patient)]
ggplot(sigs.cells, aes(x=factor(Timepoint), y=Ibrutinib_treatment, fill=Patient)) + geom_boxplot(outlier.colour=NA) +
  facet_grid(.~Patient, scales="free", space="free") + ylim(0,6) + 
  xlab("Time post ibrutinib (days)") + ylab("Ibrutinib signature") +
  guides(fill=FALSE) +
  theme_bw(24) + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
ggsave(dirout(out, "4_Ibrutinib_Signature.pdf"), width=6, height=7)
# Heatmap
(load(dirout("30_9_4_Signatures_nUMI_Cutoff_TFs/Genesets.RData")))
set.nam <- "Bcells_Human_Gene_Atlas"
set.nam <- "BATF"
set.nam <- "HALLMARK_INFLAMMATORY_RESPONSE"
dat <- pbmcFull@data
sampleAnnot <- subset(pbmcFull@data.info, CellType %in% c("CLL"))
sampleAnnot$sample <- gsub("d30", "d030", gsub("d0", "d000", sampleAnnot$sample))
sampleAnnot$patient <- gsub("(\\w+)_d\\d+", "\\1", sampleAnnot$sample)
sampleAnnot$timepoint <- gsub("(\\w+)_d(\\d+)", "\\2", sampleAnnot$sample)
cells <- rownames(sampleAnnot)[do.call(c, lapply(split(1:nrow(sampleAnnot), factor(sampleAnnot$sample)), function(x) sample(x, min(length(x), 100))))]
dat <- dat[genesets[[set.nam]][genesets[[set.nam]] %in% rownames(dat)], cells]
dat <- dat[apply(dat, 1, max) != 0,,drop=F]
dat <- dat - apply(dat, 1, min)
dat <- dat / apply(dat,1, max)
dat <- dat[, order(with(sampleAnnot[cells,],paste0(patient, timepoint))), drop=F]
pdf(dirout(out, "4_Details2_", set.nam, ".pdf"), height=min(29, nrow(dat) * 0.3 + 3), width=min(ncol(dat)*0.03+3,29), onefile=FALSE)
pheatmap::pheatmap(dat, cluster_rows=TRUE, cluster_cols=F, scale="row",
                   annotation_col=subset(sampleAnnot, select=c(nUMI, timepoint, patient)),
                   show_colnames=FALSE,
                   color=colorRampPalette(c("darkblue", "orange", "yellow"))(8))
dev.off()
# dat <- pbmcFull@data[,rownames(subset(pbmcFull@data.info, CellType == "CLL"))]
# dat <- dat[genesets[[set.nam]][genesets[[set.nam]] %in% rownames(dat)],]
# dat <- dat[rowSums(as.matrix(dat)) > 0,]
# # dat <- dat - apply(dat, 1, min)
# # dat <- dat / apply(dat,1, max)
# sample.ids <- gsub("\\w+\\-(\\d+)", "\\1", colnames(dat))
# mm <- do.call(cbind, lapply(unique(sample.ids), function(s) apply(dat[,sample.ids == s], 1, mean)))
# (colAn <- data.table(subset(pbmcFull@data.info, CellType == "CLL"), keep.rownames=TRUE)[,c("rn", "sample"),with=F])
# colAn <- unique(colAn[,rn := gsub("\\w+\\-(\\d+)", "\\1", rn)])
# colnames(mm) <- colAn$sample[match(unique(sample.ids), colAn$rn)]
# pheatmap::pheatmap(mm, scale="row", cluster_cols=FALSE)
# str(mm)

# patient down signatures -------------------------------------------------------
otsig <- overTime.sig[grepl("d\\d", geneset)]
otsig[,geneset := cleanPatients(geneset)]
ggplot(otsig, 
       aes(x=paste0(patient, "_", timepoint), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  ylab("") + xlab("") + theme_bw(24) +  
  scale_color_gradient2(name="Effect size", high="red", mid="white", low="blue") +
  scale_size_continuous(name=expression(log[10](q-value))) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), axis.text.y=element_text()) + 
  facet_grid(. ~ CellType, scales="free", space="free") 
ggsave(dirout(out, "4_CLL_Patient_Signatures.pdf"),height=6, width=18)

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
overTime.sig[,qvalue := p.adjust(pvalue, method='BH')]
overTime.sig[,pvalue2 := pmin(5, -log10(qvalue))]
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
# 6 Genes in Monocytes
#############

# FIGURE S4B: Upregulation in Monocytes --------------------------------------------
res <- fread(dirout("13_4_Overtime_nUMI_Cutoff/", "SigGenes_overTime.tsv"))
res[,sum(logFC), by=c("cellType", "gene")][cellType == "Mono"][order(V1, decreasing=TRUE)][1:20]
res[,logqval := pmin(5, -log10(qvalue))]
res[,patient := cleanPatients(gsub("d30", "d030", patient))]
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
overTime.sig[,qvalue := p.adjust(pvalue, method='BH')]
overTime.sig[,pvalue2 := pmin(5, -log10(qvalue))]
otsig <- overTime.sig[!grepl("_d\\d", geneset) & !grepl("Fereira", geneset)]#[grepl("Bcells", geneset) | geneset == "Ibrutinib_treatment" | grepl("HALLMARK", geneset)]
# otsig <- overTime.sig[!grepl("d\\d", geneset)]
# otsig <- overTime.sig[grepl("Bcells", geneset) | geneset == "Ibrutinib_treatment" | grepl("HALLMARK", geneset)]
otsig[,geneset := gsub("HALLMARK_", "", geneset)]
otsig[,patient := cleanPatients(patient)]
# Selected
ggplot(otsig[patient %in% c("P4", "P5", "P7") & geneset %in% otsig[, sum(abs(Diff) > 1), by="geneset"][V1 >= 2]$geneset], 
       aes(x=paste0(timepoint, "_", CellType), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  ylab("") + xlab("") + theme_bw(24) +  
  scale_color_gradient2(name="Effect size", high="red", mid="white", low="blue") +
  scale_size_continuous(name=expression(log[10](q-value))) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=12), axis.text.y=element_text(size=12)) +
  facet_grid(. ~ patient, scales="free", space="free") 
ggsave(dirout(out, "7_Clone_Signatures.pdf"),height=7, width=13)
# All
ggplot(otsig[patient %in% c("P4", "P5", "P7")], 
       aes(x=paste0(timepoint, "_", CellType), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  ylab("") + xlab("") + theme_bw(24) +  
  scale_color_gradient2(name="Effect size", high="red", mid="white", low="blue") +
  scale_size_continuous(name=expression(log[10](q-value))) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=12), axis.text.y=element_text(size=12)) +
  facet_grid(. ~ patient, scales="free", space="free") 
ggsave(dirout(out, "7_Clone_Signatures_All.pdf"),height=15, width=13)


# PLOT ALL THE TSNEs ------------------------------------------------------
# PT Data
pbmcPTs <- list()
pList <- list()
for(pt in c("PT", "PBGY")){ #"FE", "PBGY", "VZS")){
  pt2 <- cleanPatients(pt)
  if(is.null(pbmcPTs[[pt]])){
    (load(dirout("25_Patient_CLL_nUMI_Cutoff/",pt,"/",pt,".RData")))
    pbmcPTs[[pt]] <- pbmc
  } else {
    pbmc <- pbmcPTs[[pt]]
  }
  pDat <- data.table(pbmc@dr$tsne@cell.embeddings, Timepoint = as.numeric(gsub("\\w+_d(\\d+)", "\\1", pbmc@meta.data$sample)))
  pList[[pt2]] <- ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=factor(Timepoint))) + geom_point(alpha=0.5) +
    scale_color_manual(name="Time point", values=coll.TIME) +
    xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2") 
  ggsave(dirout(out, "7_", pt2, "_Clones_TSNE.pdf"), width=9, height=7, plot=pList[[pt]] + theme_bw(24))
}
pList <- lapply(pList[sort(names(pList))], function(p) p + guides(color=FALSE) + theme_bw(12))
pList2 <- lapply(names(pList), function(pn) pList[[pn]] + ggtitle(paste("Patient", pn)) + theme(plot.title=element_text(size=24)))
(p <- grid.arrange(grobs=pList2, ncol=1))
ggsave(dirout(out, "7_Clones_All.jpg"), height=7, width=3.5, plot=p)

pbmc <- pbmcPTs[["PT"]]
tfs <- c("NFATC1", "IKZF1")
pDat <- data.table(pbmc@dr$tsne@cell.embeddings, t(as.matrix(pbmc@data[tfs,])), sample=pbmc@meta.data$sample)
pDat[,sample := gsub("30", "030", sample)]
tfsx <- tfs[1]
for(tfsx in tfs){
  ggplot(pDat, aes_string(x="tSNE_1", y="tSNE_2", color=tfsx)) + geom_point() + facet_grid(. ~ sample) + 
    scale_color_gradient(low="lightgrey", high="blue")
  ggsave(dirout(out, "7_",tfsx, ".pdf"), height=7, width=21)
}

for(tfsx in tfs){
  pDat2 <- list()
  for(x in unique(pDat$sample)){
    pDat2[[x]] <- sapply(1:200, function(i) mean(sample(pDat[sample == x][[tfsx]], 30)))
  }
  ggplot(melt(pDat2), aes(x=L1, y=value, fill=L1)) + geom_violin(alpha=0.5) + geom_jitter(width=0.1, height=0, alpha=0.5)
  ggsave(dirout(out, "7_Bootstrap_",tfsx, ".pdf"), height=7, width=7)
}