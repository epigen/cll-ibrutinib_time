require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("cll-time_course")

R.version

# out <- "25_Patient_CLL_nUMI_Cutoff/"
out <- "25_Patient_CLL/"
dir.create(dirout(out))

pat <- "PT"
(load(dirout(out, pat, "/", pat, ".RData")))
x <- pbmc@data.info
x$Clone <- "CLL1"
x[x$ClusterNames_0.5 == 3,]$Clone <- "CLL2"
x[x$ClusterNames_0.5 == 2,]$Clone <- "CLL3"
x[x$ClusterNames_0.5 == 5,]$Clone <- "CLL4"
pbmc@data.info <- x
save(pbmc, file=dirout(out, pat, "/", pat, ".RData"))

pat.md[["PT"]][,Clone := "CLL1"]
pat.md[["PT"]][ClusterNames_0.5 == 6, Clone := "CLL2"]
pat.md[["PT"]][ClusterNames_0.5 == 5, Clone := "CLL3"]
pat.md[["PT"]][ClusterNames_0.5 == 4, Clone := "CLL4"]
pat.md[["PT"]][ClusterNames_0.5 == 7, Clone := "CLL5"]
  
pat <- "VZS"
(load(dirout(out, pat, "/", pat, ".RData")))
x <- pbmc@data.info
x$Clone <- "CLL1"
x[x$ClusterNames_0.5 == 3,]$Clone <- "CLL2"
x[x$ClusterNames_0.5 == 2,]$Clone <- "CLL3"
x[x$ClusterNames_0.5 == 5,]$Clone <- "CLL4"
pbmc@data.info <- x
save(pbmc, file=dirout(out, pat, "/", pat, ".RData"))

pat <- "FE"
(load(dirout(out, pat, "/", pat, ".RData")))
x <- pbmc@data.info
x$Clone <- "CLL1"
x[x$ClusterNames_0.5 %in% c(5,8),]$Clone <- "CLL2"
x[x$ClusterNames_0.5 == 7,]$Clone <- "CLL3"
pbmc@data.info <- x
save(pbmc, file=dirout(out, pat, "/", pat, ".RData"))

pat <- "PBGY"
(load(dirout(out, pat, "/", pat, ".RData")))
x <- pbmc@data.info
x$Clone <- "CLL1"
x[x$ClusterNames_0.5 == 5,]$Clone <- "CLL2"
x[x$ClusterNames_0.5 == 7,]$Clone <- "CLL3"
pbmc@data.info <- x
save(pbmc, file=dirout(out, pat, "/", pat, ".RData"))

system("bash src/single_cell_RNA/25_1_CLL_Patients_run.sh")