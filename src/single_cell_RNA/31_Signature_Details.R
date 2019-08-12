require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(pheatmap)
require(gplots)

project.init2("cll-time_course")
out <- "31_Signature_Details/"
dir.create(dirout(out))


(load(file=dirout("30_9_Signatures_inclDay30/","Scores.RData")))
(load(file=dirout("30_9_Signatures_inclDay30/", "OverTime_Tests.RData")))

Dat1[CellType == "CLL"]

ggplot(
  overTime.sig[grepl("HALLMARK", geneset) | geneset == "Ibrutinib_treatment"], 
  aes(x=paste0(patient, "_", timepoint), y=geneset, color=Diff, size=pvalue2)) + 
  geom_point() +
  scale_color_gradient2(low="blue", mid="white", high="red") + theme_bw() +
  ylab("") + xlab("") + theme_bw(24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(. ~ CellType)
ggsave(dirout(out, "AllCells_Hallmark.pdf"),height=20, width=20)





