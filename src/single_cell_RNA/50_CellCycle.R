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

ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color="Phase")) + geom_point(alpha=0.3) + ggtitle("Cell cycle") +
  scale_color_manual(values=c("red", "green", "blue"))
ggsave(dirout(out, "Cell_Cycle_RGB.pdf"), width=7, height=7)

ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color="Phase")) + geom_point(alpha=0.5) + ggtitle("Cell cycle")
ggsave(dirout(out, "Cell_Cycle.pdf"), width=7, height=7)

ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color="S.Score")) + geom_point() + ggtitle("S Score")
ggsave(dirout(out, "S_Score.pdf"), width=7, height=7)

ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color="G2M.Score")) + geom_point() + ggtitle("G2M.Score")
ggsave(dirout(out, "G2M.Score.pdf"), width=7, height=7)
