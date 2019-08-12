library(project.init)
project.init2("cll-time_course")

out <- "04_SequencingSaturation/"
dir.create(dirout(out))


# ANALYZE METRICS ---------------------------------------------------------
mDat <- data.table()

# Prozessed samples
path2 <- "~/projects_shared/cll-time_course/results/cellranger_count/"
f2 <- list.files(path2)
f2 <- f2[grepl("d\\d?\\d?0$", f2) | grepl("LiveBulk.+?078$", f2) |  grepl("\\d?\\d?0d", f2)]
f2 <- f2[!grepl("\\.log", f2)]
f.x <- f2[1]
for(f.x in f2){
  print(f.x)
  mDat <- rbind(mDat, data.table(t(fread(paste0(path2, f.x,"/outs/metrics_summary.csv"))), sample=gsub(m.nam, "", f.x), keep.rownames=T))
}

mDat[,value := as.numeric(gsub("\\%", "", gsub("\\,", "", V1)))]
mDat[,metric := gsub(" ", "", rn)]
sort(table(mDat$metric))
mDat[metric == "NumberofCells", metric := "EstimatedNumberofCells"]
pDat <- dcast.data.table(mDat, sample ~ metric, value.var = "value")
# should we resequence
pDat$sample2 <- gsub("LiveBulk_10x_", "", pDat$sample)
pDat$sample2 <- gsub("_10xTK078", "", pDat$sample2)
pDat$sample2 <- gsub("FE_FE", "FE", pDat$sample2)

table(pDat$sample2)

ggplot(pDat, 
       aes(x=MedianUMICountsperCell, fill=EstimatedNumberofCells, y=SequencingSaturation, label=sample2)) + 
  geom_label(color="white") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlim(0,3000) +
  theme_bw(24)#+ ylim(0,8000)
ggsave(dirout(out, "x.pdf"), width=12, height=7)
