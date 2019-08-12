setOption(option="ROUT.DIR", value="/scratch/lab_bock/shared/projects/cll-time_course/results/single_cell_RNA/")
if(!dir.exists(getOption("ROUT.DIR"))){
  dir.create(getOption("ROUT.DIR"))
}
