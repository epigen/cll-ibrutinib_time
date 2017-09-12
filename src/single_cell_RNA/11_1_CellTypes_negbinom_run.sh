module unload R/3.2.3
module load R/3.3.2

cells=("Monos" "NurseLikeCells" "Bcells" "NKcells" "Tcells1" "CD4" "CD8")

for cell in ${cells[@]}
do
  echo $cell
  sbatch --job-name="Seurat noIGHKL $cell" --ntasks=5 --mem=180000 --partition=longq --time=4-00:00:00 \
      --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_noIGHLK.R $cell" \
      --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes_noIGHKL/11_CellTypes_noIGHKL_${cell}.log
done

for cell in ${cells[@]}
do
  echo $cell
  sbatch --job-name="Seurat noRPstrict $cell" --ntasks=5 --mem=180000 --partition=longq --time=4-00:00:00 \
      --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_noRPstrict.R $cell" \
      --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes_noRPstrict/11_CellTypes_noRPstrict_${cell}.log
done

for cell in ${cells[@]}
do
  echo $cell
  sbatch --job-name="Seurat negbinom $cell" --ntasks=12 --mem=180000 --partition=longq --time=4-00:00:00 \
      --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_allDataBest_negbinom.R $cell" \
      --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes_negbinom/11_CellTypes_negbinom_${cell}.log
done