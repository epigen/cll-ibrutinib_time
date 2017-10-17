module unload R/3.2.3
module unload R/3.3.2
module unload gcc/4.8.2
module load gcc/7.1.0
module load R/3.4.0

cells=("PT" "PBGY" "FE" "VZS")

for cell in ${cells[@]}
do
  echo $cell
  sbatch --job-name="Seurat patient CLL $cell" --cpus-per-task=5 --mem=180000 --partition=longq --time=7-00:00:00 \
      --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/25_2_CLL_PatientGroups.R $cell" \
      --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/25_Patient_CLL_nUMI_Cutoff/${cell}.log
done