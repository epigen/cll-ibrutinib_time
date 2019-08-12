require("project.init")
require(Seurat)
require(methods)
project.init2("cll-time_course")


# sbatch --job-name="Seurat fscLVM" --ntasks=32 --mem=180000 --partition=longq --time=5-00:00:00 \
# --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/35_fscLVM.R" \
# --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/35_fscLVM.log


out.fscLVM <- "35_fscLVM_2/"
dir.create(dirout(out.fscLVM))


# LISTS -------------------------------------------------------------------
(load(dirout("20_AggregatedLists/lists.RData")))

file <- "/data/groups/lab_bock/shared/resources/gene_sets/MSigDB/6.0/Human/h.all.v6.0.symbols.gmt"
lines <- readLines(file)
for(ll.nam in names(cll.lists)){
  lines <- c(lines, paste(c(ll.nam, "xx", cll.lists[[ll.nam]]), collapse="\t"))
}

fileConn<-file(dirout(out.fscLVM, "Lists.gmt"))
writeLines(lines, fileConn)
close(fileConn)


# INPUT MATRIX ------------------------------------------------------------
load(dirout("10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.RData"))

# exp.csv <- pbmc@raw.data[row.names(pbmc@data), colnames(pbmc@data)]
exp.csv <- pbmc@data
str(exp.csv)
table(gsub("\\w+\\-(\\d)","\\1", colnames(pbmc@data)), pbmc@data.info$sample) # this is just to see if the rows are ordered like data.info
exp.csvT <- t(as.matrix(exp.csv))
exp.csvT[1:10,1:10]
row.names(exp.csvT) <- pbmc@data.info$sample
exp.csvT[1:10,1:10]
table(apply(exp.csvT, 2, max) > 0)

if(!file.exists(dirout(out.fscLVM, "matrix.csv.gz"))){
  write.csv(exp.csvT, file=dirout(out.fscLVM, "matrix.csv"))
  system(paste('sed -e s/\\"\\",// -i', dirout(out.fscLVM, "matrix.csv"))) # removes the "" from the first column of the matrix
  system(paste("rm", dirout(out.fscLVM, "matrix.csv.gz"))) # removes the old file otherwise the next step won't work
  system(paste("gzip", dirout(out.fscLVM, "matrix.csv")))  # zip
}



nhidden <- 3
for(nhidden in c(3)){
  out.fscLVM.x <- paste0(out.fscLVM, "fscLVM_files_", nhidden, "hidden/")
	dir.create(dirout(out.fscLVM.x))
  if(!file.exists(dirout(out.fscLVM.x, "X.csv"))){
    system(paste("python src/single_cell_RNA/90_fscLVM.py",
  	             dirout(out.fscLVM, "matrix.csv.gz"),
  	             dirout(out.fscLVM.x),
                 dirout(out.fscLVM, "Lists.gmt"),
  	             nhidden
  	  ))
  }
  
  # ANALYZE
  source("src/single_cell_RNA/35_2_fscLVM_analysis.R", echo=TRUE)
}