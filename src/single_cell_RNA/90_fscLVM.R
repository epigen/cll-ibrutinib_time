require("project.init")
require(Seurat)
require(methods)
# 
# project.init2("cll-time_course")
# out <- "11_CellTypes/"
# dir.create(dirout(out))
# 
# sample.x <- "allDataBest_NoDownSampling_noIGH"
# cell <- "Monos"
# 
# outS <- paste0(out, sample.x, "_", cell, "/")



out.fscLVM <- paste0(outS, "fscLVM/")
dir.create(dirout(out.fscLVM))
                     
# load(file=dirout(outS, cell,".RData"))


# exp.csv <- pbmc@raw.data[row.names(pbmc@data), colnames(pbmc@data)]
exp.csv <- pbmc@scale.data
str(exp.csv)
table(gsub("\\w+\\-(\\d)","\\1", colnames(pbmc@data)), pbmc@data.info$sample)
exp.csvT <- t(as.matrix(exp.csv))
exp.csvT[1:10,1:10]
row.names(exp.csvT) <- pbmc@data.info$sample
exp.csvT[1:10,1:10]

if(!file.exists(dirout(out.fscLVM, "matrix.csv.gz"))){
  write.csv(exp.csvT, file=dirout(out.fscLVM, "matrix.csv"))
  system(paste('sed -e s/\\"\\",// -i', dirout(out.fscLVM, "matrix.csv")))
  system(paste("rm", dirout(out.fscLVM, "matrix.csv.gz")))
  system(paste("gzip", dirout(out.fscLVM, "matrix.csv")))
}



nhidden <- 3
for(nhidden in 1:5){
	outDir <- dirout(out.fscLVM, "fscLVM_files_", nhidden, "hidden/")
	dir.create(outDir)

	system(paste("python src/single_cell_RNA/90_fscLVM.py",
	             dirout(out.fscLVM, "matrix.csv.gz"),
	             outDir,
	             "$RESOURCES/gene_sets/MSigDB/6.0/Human/h.all.v6.0.symbols.gmt",
	             nhidden
	  ))
}