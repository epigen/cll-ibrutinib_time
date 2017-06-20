
require(data.table)
require(pheatmap)

out.cnvs <- paste0(outS, "CNVs/")
dir.create(dirout(out.cnvs))

(gene.annot <- fread(paste0(Sys.getenv("RESOURCES"), "/genomes/hg38/10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf")))
table(gene.annot$V3)
a2 <- gene.annot[V3 == "gene"]
a2[,gene := gsub('.+\\; gene_name \\"(.+?)\\"\\; .+', "\\1", V9)]
length(unique(a2$gene))

geneOrder <- a2[,c("V1", "V4", "gene"),with=F]
colnames(geneOrder) <- c("chr", "start", "gene")

table(geneOrder$chr)
geneOrder <- geneOrder[chr %in% c(as.character(1:22), "X", "Y")]

geneOrder[order(chr, start, decreasing=FALSE)]
plot(geneOrder$start)


table(rownames(pbmc@data) %in% geneOrder$gene)
mat <- pbmc@scale.data
genes.used <- intersect(geneOrder$gene, rownames(mat))

mat <- mat[genes.used,]
geneOrder2 <- geneOrder[gene %in% genes.used]

slide.stepsize <- 10
# slide.windowsize <- 100
# slide.diff <- slide.windowsize - slide.stepsize

# create a matrix with regions * cells
cnv.scores <- matrix(0, ncol=ncol(mat), nrow=geneOrder2[,.(x = ceiling((.N)/slide.stepsize)), by="chr"][, sum(x)])
colnames(cnv.scores) <- colnames(mat)

# scale row of matrix
mat2 <- t(scale(t(mat)))
mat2[mat2 > 3] <- 3
mat2[mat2 < -3] <- -3

# sliding average
chr.x <- 1
i2 <- 0
rowNams <- c()
for(chr.x in unique(geneOrder$chr)){
  geneOrder.chr <- geneOrder2[chr == chr.x]
  for(i in 1:ceiling((nrow(geneOrder.chr))/slide.stepsize)){
    i2 <- i2 + 1
    cnv.scores[i2,] <- apply(mat2[geneOrder.chr[((i-1)*slide.stepsize+1):min(i*slide.stepsize, nrow(geneOrder.chr))]$gene,], 2, mean, na.rm=TRUE)
    rowNams <- c(rowNams, paste0(chr.x, "_", i2))
  }
}
rownames(cnv.scores) <- rowNams
cnv.scores <- cnv.scores[,order(as.numeric(gsub("[A-Z]+\\-(\\d+)","\\1",colnames(cnv.scores))))]


save(cnv.scores, file=dirout(out.cnvs, "scores.RData"))

pdf(dirout(out.cnvs, "hits.pdf"), height=10, width=10, onefile=F)
# str(prows <-  rownames(cnv.scores)[order(abs(apply(cnv.scores, 1, sum)), decreasing=TRUE)[1:200]])
# str(pcols <- order(abs(apply(cnv.scores, 2, sum)), decreasing=TRUE)[1:200])
pheatmap(cnv.scores[,colnames(cnv.scores) %in% sample(colnames(cnv.scores), 500)], 
         cluster_rows=FALSE, cluster_cols=FALSE,
         annotation_col=subset(pbmc@data.info, select=sample),
         annotation_row=data.frame(chromosome=gsub("(.+)\\_.+", "\\1", rownames(cnv.scores)),row.names=rownames(cnv.scores)),
         show_rownames=F,
         show_colnames=F
         )
dev.off()
         
tail(rownames(cnv.scores))
