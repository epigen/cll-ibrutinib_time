library("project.init")
project.init2("cll-time_course")

devtools::install_url("https://github.com/satijalab/seurat/releases/download/v1.4.0/Seurat_1.4.0.12.tgz", binary = TRUE)
library(Seurat)

