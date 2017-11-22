module unload R/3.2.3
module unload R/3.3.2
module unload gcc/4.8.2
module unload gcc/6.1.0
module load gcc/7.1.0
module load R/3.4.0

cells=("CLL" "NurseLikeCell" "CD8" "CD4" "NK" "Mono")

for cell in ${cells[@]}
do
  echo $cell
  sbatch --job-name="Seurat noIGHKL $cell" --cpus-per-task=5 --mem=180000 --partition=longq --time=7-00:00:00 \
      --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_3_CellTypes_UMI_Cutoff.R $cell" \
      --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/11_2_3_CellTypes_UMI_Cutoff/${cell}.log
done


for cell in ${cells[@]}
do
  echo $cell
  sbatch --job-name="Seurat noIGHKL $cell" --cpus-per-task=5 --mem=180000 --partition=longq --time=7-00:00:00 \
      --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_2_CellTypes_inclDay30_4.R $cell" \
      --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes_inclDay30/${cell}.log
done