require("project.init")
require(Seurat)
require(methods)
project.init2("cll-time_course")

out <- "50_CellCycle/"
dir.create(dirout(out))

load(dirout("10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.RData"))

pbmc2 <- UpdateSeuratObject(pbmc)

cc.genes <- readLines(con = "~/resources_nfortelny/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:98]
pbmc2 <- CellCycleScoring(object = pbmc2, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)

data.table(pbmc2@meta.data)

stopifnot(all(rownames(pbmc2@dr$tsne@cell.embeddings) == rownames(pbmc2@meta.data)))

pDat <- data.table(pbmc2@dr$tsne@cell.embeddings, pbmc2@meta.data)
write.table(pDat, file=dirout(out, "Cells.txt"), sep="\t",row.names=F, quote=F)

ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color="Phase")) + geom_point(alpha=0.3) + ggtitle("Cell cycle") +
  scale_color_manual(values=c("red", "green", "blue"))
ggsave(dirout(out, "Cell_Cycle_RGB.pdf"), width=7, height=7)

ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color="Phase")) + geom_point(alpha=0.5) + ggtitle("Cell cycle")
ggsave(dirout(out, "Cell_Cycle.pdf"), width=7, height=7)

ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color="S.Score")) + geom_point() + ggtitle("S Score")
ggsave(dirout(out, "S_Score.pdf"), width=7, height=7)

ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color="G2M.Score")) + geom_point() + ggtitle("G2M.Score")
ggsave(dirout(out, "G2M.Score.pdf"), width=7, height=7)


ggplot(pDat, aes(x=S.Score, y=G2M.Score, color=Phase)) + geom_point() 
ggsave(dirout(out, "ScoreDetails.pdf"))

ggplot(pDat, aes(x=cellType, fill=Phase)) + 
  geom_bar(position="dodge") + 
  #scale_fill_manual(values=c("red", "blue", "yellow")) + 
  facet_grid(sample ~ .)
ggsave(dirout(out, "Proportions_Patients.pdf"), height=15, width=10)

pDat2 <- pDat[,.(cnt = .N), by=c("sample", "Phase", "cellType")]
pDat2[,sampleCellCount  := sum(cnt), by=c("sample", "cellType")]
pDat2[,prop := cnt/sampleCellCount]

ggplot(pDat2[cnt > 50], aes(x=cellType, y=prop, fill=Phase)) + 
  geom_bar(stat="identity", position="dodge") + 
  #scale_fill_manual(values=c("red", "green", "blue")) + 
  facet_grid(sample ~ .)
ggsave(dirout(out, "Proportions_Patients2.pdf"), height=15, width=10)

# ggplot(pDat2, aes(x=sample, y=cellType, size=sampleCellCount)) + 
#   geom_point(aes(alpha=prop), data=pDat2[Phase == "G1"], color="red") + 
#   geom_point(aes(alpha=prop), data=pDat2[Phase == "G2M"], color="blue") + 
#   geom_point(aes(alpha=prop), data=pDat2[Phase == "S"], color="yellow") + 
#   scale_alpha_continuous()
# ggsave(dirout(out, "Proportions_PatientsPoints.pdf"), height=15, width=10)

ggplot(pDat2, aes(y=sample, x=cellType, size=sampleCellCount, alpha=prop, color=Phase)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_alpha_continuous(range=c(0,1)) + scale_color_manual(values=c("red", "green", "blue"))
ggsave(dirout(out, "Proportions_PatientsPoints.pdf"), height=10, width=10)


ggplot(pDat2, aes(x="", y=prop, fill=Phase)) + 
  geom_bar(stat="identity") + coord_polar("y", start=0) +
  facet_grid(sample ~ cellType)
ggsave(dirout(out, "Proportions_pie.pdf"), height=15, width=15)


write.table(pDat2, file=dirout(out, "Numbers.txt"), sep="\t",row.names=F, quote=F)


pDat2 <- fread(dirout(out, "Numbers.txt"))

pDat2[,patient := substr(pDat2$sample, 0,2)]

pDatPT <- pDat2[sample == "PT_d0"]
pDatPT[,patient := "PT2"]
pDat3 <- rbind(pDat2, pDatPT)
pDat3[sample == "PT_d280", patient := "PT2"]
pDat3 <- pDat3[patient != "KI"]

pDat3[grepl("_d?0d?$", sample), timepoint := "early"]
pDat3[is.na(timepoint), timepoint := "late"]
table(pDat3$sample, pDat3$timepoint)
table(pDat3$patient, pDat3$timepoint)

pDat4 <- dcast.data.table(pDat3[,c("patient", "Phase", "cellType","timepoint", "prop"), with=F], patient + Phase + cellType ~ timepoint, value.var="prop")
pDat4[is.na(early), early := 0]
pDat4[is.na(late), late := 0]

pDat4[,sumEarly := sum(early), by=c("patient", "Phase", "cellType")]
pDat4[,sumLate := sum(late), by=c("patient","Phase", "cellType")]

pDat4[sumLate == 0, early := 0]
pDat4[sumEarly == 0, late := 0]


pDat4[,diff := late - early]

ggplot(pDat4, aes(x=cellType, y=diff, fill=Phase)) + geom_bar(stat="identity",position = "dodge") + facet_grid(patient~.)
ggsave(dirout(out, "ProportionChanges.pdf"), height=15, width=15)


# The Other way -----------------------------------------------------------
# not working
# require(org.Hs.eg.db)
# require(scran)
# 
# ccMarkers <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
# anno <- select(org.Hs.eg.db, keys=rownames(pbmc@data), keytype="SYMBOL", column="ENSEMBL")
# ensembl <- anno$ENSEMBL[match(rownames(pbmc@data), anno$SYMBOL)]
# assignments <- cyclone(as.matrix(pbmc@data), ccMarkers, gene.names=ensembl, BPPARAM=MulticoreParam(workers=1), verbose=TRUE)
# 
