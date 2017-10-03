module unload R/3.2.3
module load R/3.3.2

cells=("CLL" "NurseLikeCell" "CD8" "CD4" "NK" "Mono")
cells=("CD8" "CLL")

for cell in ${cells[@]}
do
  echo $cell
  sbatch --job-name="Seurat noIGHKL $cell" --cpus-per-task=5 --mem=180000 --partition=longq --time=7-00:00:00 \
      --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_2_CellTypes_inclDay30_3.R $cell" \
      --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes_inclDay30/${cell}.log
done



## OLD RUNS:
for cell in ${cells[@]}
do
  echo $cell
  sbatch --job-name="Seurat noIGHKL $cell" --ntasks=5 --mem=180000 --partition=longq --time=4-00:00:00 \
      --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_noIGHLK.R $cell" \
      --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes_noIGHKL_negbinom/11_CellTypes_noIGHKL_${cell}.log
done

for cell in ${cells[@]}
do
  echo $cell
  sbatch --job-name="Seurat noRPstrict $cell" --ntasks=5 --mem=180000 --partition=longq --time=4-00:00:00 \
      --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_noRPstrict.R $cell" \
      --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes_noRPstrict_negbinom/11_CellTypes_noRPstrict_${cell}.log
done

for cell in ${cells[@]}
do
  echo $cell
  sbatch --job-name="Seurat negbinom $cell" --ntasks=5 --mem=180000 --partition=longq --time=4-00:00:00 \
      --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_allDataBest_negbinom.R $cell" \
      --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes_negbinom_negbinom/11_CellTypes_negbinom_${cell}.log
done