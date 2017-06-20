require("project.init")
require(Seurat)
require(methods)
project.init2("cll-time_course")

out <- "15_Magic/"
dir.create(dirout(out))


load(dirout("10_Seurat/allDataBest_NoDownSampling_noIGH/allDataBest_NoDownSampling_noIGH.RData"))

exp.csv <- pbmc@data
str(exp.csv)
table(gsub("\\w+\\-(\\d)","\\1", colnames(pbmc@data)), pbmc@data.info$sample) # this is just to see if the rows are ordered like data.info
exp.csvT <- t(as.matrix(exp.csv))
exp.csvT[1:10,1:10]
row.names(exp.csvT) <- pbmc@data.info$sample
exp.csvT[1:10,1:10]
table(apply(exp.csvT, 2, max) > 0)

if(!file.exists(dirout(out, "matrix.csv.gz"))){
  write.csv(exp.csvT, file=dirout(out, "matrix.csv"))
  system(paste('sed -e s/\\"\\",// -i', dirout(out, "matrix.csv"))) # removes the "" from the first column of the matrix
  system(paste("rm", dirout(out, "matrix.csv.gz"))) # removes the old file otherwise the next step won't work
  system(paste("gzip", dirout(out, "matrix.csv")))  # zip
}

if(!file.exists(dirout(out, "magic.matrix.csv"))){
  system("module load python/3.5.0; python3 src/single_cell_RNA/15_Magic.py")
}