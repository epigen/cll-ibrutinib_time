require("project.init")
# install.packages("nloptr",repos="https://cloud.r-project.org")
require(Seurat)
require(methods)
project.init2("cll-time_course")


out <- "17_Mono_Plots/"
dir.create(dirout(out))


(inDir <- paste0("13_4_Overtime_","inclDay30","/"))
res <- fread(dirout(inDir, "SigGenes_overTime.tsv"))

res[,sum(logFC), by=c("cellType", "gene")][cellType == "Mono"][order(V1, decreasing=TRUE)][1:20]
res[,logqval := pmin(5, -log10(qvalue))]

res[,patient := gsub("d30", "d030", patient)]

ggplot(res[(gene %in% c("TNF", "CCL3", "CCL3L3", "CXCL8") | grepl("NFKBI[Z|A]", gene)) & cellType %in% c("CD8", "CD4", "Mono", "CLL")], 
       aes(x=patient, y=gene, color=logFC, size=logqval)) + geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_gradient2(high="red", mid="white", low="blue") +
  facet_grid(. ~ cellType) + xlab("") + ylab("")
ggsave(dirout(out, "Mono_NFKB_Genes.pdf"), height=7, width=10)
