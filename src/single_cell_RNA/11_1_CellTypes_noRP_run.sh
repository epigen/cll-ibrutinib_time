module unload R/3.2.3
module load R/3.3.2

cells=("Monos" "NurseLikeCells" "Bcells" "NKcells" "Tcells1" "CD4" "CD8")

for cell in ${cells[@]}
do
  echo $cell
  sbatch --job-name="Seurat $cell" --ntasks=12 --mem=180000 --partition=longq --time=1-22:00:00 \
      --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_noRP.R $cell" \
      --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes_noRP_tobit/allDataBest_NoDownSampling_noRP_${cell}.log
done