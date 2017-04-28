module unload R/3.2.3
module load R/3.3.2

array=("PT_d0" "PT_d120" "PT_d280" "VZS_d0" \
  "LiveBulk_10x_FE_FE1_d0_10xTK078" "LiveBulk_10x_KI_KI1_d0_10xTK078" \
  "LiveBulk_10x_PBGY1_0d" "LiveBulk_10x_FE7_120d" "LiveBulk_10x_PBGY7_150d" "LiveBulk_10x_VZS7_120d" \
  "allData" "allData2" "allDataBest" "allDataBest_NoDownSampling" \
  "allDataBest_noIGH" "allDataBest_NoDownSampling_noIGH" \
  "allDataBest_NoDownSampling_tissueGenes"
)

array=("allDataBest_NoDownSampling_tissueGenes")

for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="Seurat $file" --ntasks=12 --mem=180000 --partition=longq --time=3-00:00:00 \
        --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/10_1_Seurat.R filtered ${file}" \
        --output="${PROCESSED}/cll-time_course/results/single_cell_RNA/10_Seurat/${file}.log"
done


# RAW ANALYSIS, not used because it fails often and doesn't make a large difference
# for file in ${array[@]}
# do
#     echo $file
#     sbatch --job-name="Seurat $file raw" --ntasks=1 --mem=180000 --partition=longq --time=3-00:00:00 \
#         --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/10_1_Seurat.R raw ${file}" \
#         --output="${PROCESSED}/cll-time_course/results/single_cell_RNA/10_Seurat_raw/${file}.log"
# done