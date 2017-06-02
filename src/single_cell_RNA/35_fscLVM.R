require("project.init")
require(Seurat)
require(methods)
project.init2("cll-time_course")


# sbatch --job-name="Seurat fscLVM" --ntasks=32 --mem=180000 --partition=longq --time=5-00:00:00 \
# --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/35_fscLVM.R" \
# --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/35_fscLVM.log


out.fscLVM <- "35_fscLVM/"
dir.create(dirout(out.fscLVM))
     
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
for(nhidden in 1:5){
	outDir <- dirout(out.fscLVM, "fscLVM_files_", nhidden, "hidden/")
	dir.create(outDir)

  if(!dir.exists(outDir)){
  	system(paste("python src/single_cell_RNA/90_fscLVM.py",
  	             dirout(out.fscLVM, "matrix.csv.gz"),
  	             outDir,
  	             "$RESOURCES/gene_sets/MSigDB/6.0/Human/h.all.v6.0.symbols.gmt",
  	             nhidden
  	  ))
  }
  
  # ANALYZE
  # 	outFSCLVM <- paste0(out.fscLVM, "fscLVM_files_", nhidden, "hidden/")
  # 	balance.barcodes <- TRUE
  # 	source("src/single_cell_RNA/90_fscLVM_test.R", echo=TRUE)
  # 	balance.barcodes <- FALSE	
  # 	source("src/single_cell_RNA/90_fscLVM_test.R", echo=TRUE)
}