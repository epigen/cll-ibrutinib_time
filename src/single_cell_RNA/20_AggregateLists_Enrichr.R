require("project.init")
require("Seurat")
require("gplots")
require(methods)
require(pheatmap)
require(gdata)
require(enrichR) #devtools::install_github("definitelysean/enrichR")

out <- "20_AggregatedLists/"
dir.create(dirout(out))

project.init2("cll-time_course")
atacRes <- list()

cell <- "Bcell"
for(cell in c("Bcell", "CD4", "CD8", "CLL","Mono", "NK")){
  x <- fread(paste0(getOption("PROCESSED.PROJECT"), "results/", "cll-time_course_peaks.coverage.joint_qnorm.pca_fix.power.diff_timepoint.limma.",cell,"_diff.gene_signature.weighted.csv"))
  atacRes[[paste0(cell, "_up")]] <- x[V2 > 0.5]$V1
  atacRes[[paste0(cell, "_down")]] <- x[V2 < -0.5]$V1
}


if(!file.exists(dirout(out, "CLL_Signature.xls")){
  system(paste0("wget http://genome.cshlp.org/content/suppl/2013/11/21/gr.152132.112.DC1/Supplemental_File2_diffExpGenes.xls -O ", dirout(out, "CLL_Signature.xls")))
}
if(!file.exists(dirout(out, "Proliferation_Signature.xlsx"))){
  system(paste0("wget http://www.impactjournals.com/oncotarget/index.php?journal=oncotarget&page=article&op=downloadSuppFile&path%5B%5D=16961&path%5B%5D=24097 -O ", dirout(out, "Proliferation_Signature.xlsx")))
}
if(!file.exists(dirout(out, "CLL_Signature2.xlsx"))){
  system(paste0("wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-014-0125-z/MediaObjects/13073_2014_125_MOESM2_ESM.xlsx -O ", dirout(out, "CLL_Signature2.xlsx")))
}


# Fereira et al, Genome Res 2013
x2 <- data.table(read.xls(dirout(out, "CLL_Signature.xls"), sheet=2))
atacRes[["Fereira_normal"]] <- x2[md.tumor < md.normal]$genename
atacRes[["Fereira_tumor"]] <- x2[md.tumor > md.normal]$genename
x2 <- data.table(read.xls(dirout(out, "CLL_Signature.xls"), sheet=3))
atacRes[["Fereira_C1"]] <- x2[md.hcC1 > md.hcC2]$genename
atacRes[["Fereira_C2"]] <- x2[md.hcC1 < md.hcC2]$genename

# Ramaker et al., Oncotarget 2017
x2 <- data.table(read.xls(dirout(out, "Proliferation_Signature.xlsx"), sheet=))
atacRes[["Ramaker_Proliferation"]] <- x2[-1,][[1]]

# two groups
# x2 <- data.table(read.xls(dirout(out, "CLL_Signature2.xlsx"), sheet=3))

# Ibrutinib study
x2 <- fread(dirout(out,"ibrutinib_treatment_expression.timepoint_name.csv"))
atacRes[["Ibrutinib_treatment"]] <- x2[padj < 0.05 & log2FoldChange < 0]$gene


cll.lists <- atacRes
save(cll.lists, file=dirout(out, "lists.RData"))

# ENRICHR
enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location")
enrichRes <- data.table()
hitSets <- atacRes
for(grp.x in names(hitSets)){
  ret=try(as.data.table(enrichGeneList(hitSets[[grp.x]],databases = enrichrDBs)),silent = FALSE)
  if(!any(grepl("Error",ret)) && nrow(ret) > 0){
    enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
  }
}
enrichRes$n <- sapply(strsplit(enrichRes$genes,","), length)
enrichRes <- enrichRes[n > 3]
write.table(enrichRes[qval < 0.05], file=dirout(out, cell, "_EnrichR",filter, ".tsv"), sep="\t", quote=F, row.names=F)

if(nrow(enrichRes) > 2 & length(unique(enrichRes$grp)) > 1){
  pDat <- dcast.data.table(enrichRes, make.names(category) ~ grp, value.var="qval")
  pDatM <- as.matrix(pDat[,-"category", with=F])
  pDat$category <- gsub("\\_(\\w|\\d){8}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){12}", "", pDat$category)
  pDat$category <- substr(pDat$category, 0, 50)
  row.names(pDatM) <- pDat$category
  pDatM[is.na(pDatM)] <- 1
  str(pDatM <- pDatM[apply(pDatM <= 5e-2,1,sum)>=1,apply(pDatM <= 5e-2,2,sum)>=1, drop=F])
  if(nrow(pDatM) >=2 & ncol(pDatM) >= 2){
    pDatM <- -log10(pDatM)
    pDatM[pDatM > 4] <- 4
    # pDatM[pDatM < 1.3] <- 0
    pdf(dirout(out, "EnrichR.pdf"),onefile=FALSE, width=min(29, 6+ ncol(pDatM)*0.3), height=min(29, nrow(pDatM)*0.3 + 4))
    pheatmap(pDatM) #, color=gray.colors(12, start=0, end=1), border_color=NA)
    dev.off()
  }
}