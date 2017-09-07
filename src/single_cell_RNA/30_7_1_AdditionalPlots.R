require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(pheatmap)
require(gplots)

project.init2("cll-time_course")
inDir <- "30_7_Signatures_RowNormalized_Sampled_2/"
out <- paste0(inDir, "AdditionalAnalyses/")
dir.create(dirout(out))

(load(dirout(inDir,"Scores.RData")))


ggplot(Dat1[patient == "PT" & cellType == "CLL"], aes(x=HALLMARK_ALLOGRAFT_REJECTION, y=HALLMARK_MYC_TARGETS_V1, color=timepoint)) + 
  geom_point(alpha=0.5)
ggsave(dirout(out, "CLL_Gradient.pdf"))


ggplot(Dat1[patient == "PT" & cellType == "Mono"], aes(x=HALLMARK_TNFA_SIGNALING_VIA_NFKB, y=timepoint, color=timepoint)) + 
  geom_jitter(alpha=1)
ggsave(dirout(out, "Monos_time.pdf"))


ggplot(Dat1[patient == "PT" & cellType == "Mono"], aes(x=HALLMARK_TNFA_SIGNALING_VIA_NFKB, y=HALLMARK_MYC_TARGETS_V1, color=timepoint)) + 
  geom_point(alpha=0.5)
ggsave(dirout(out, "Monos_Gradient.pdf"))
