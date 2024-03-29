module unload R/3.2.3
module load R/3.3.2

array=("PT_d0" "PT_d120" "PT_d280" "VZS_d0" \
  "LiveBulk_10x_FE_FE1_d0_10xTK078" "LiveBulk_10x_KI_KI1_d0_10xTK078" \
  "LiveBulk_10x_PBGY1_0d" "LiveBulk_10x_FE7_120d" "LiveBulk_10x_PBGY7_150d" "LiveBulk_10x_VZS7_120d" \
  "allData" "allData2" "allDataBest" "allDataBest_NoDownSampling" \
  "allDataBest_noIGH" "allDataBest_NoDownSampling_noIGH" \
  "allDataBest_NoDownSampling_tissueGenes" \
  "allDataBest_NoDownSampling_noRP" \
  "allDataBest_NoDownSampling_noIGHLK" "allDataBest_NoDownSampling_noRPstrict" \
  "inclDay30_noIGHLK"
)

array=("inclDay30_noIGHLK")

for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="Seurat $file" --ntasks=1 --mem=100000 --partition=develop --time=05:00:00 \
        --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/10_1_Seurat.R raw ${file}" \
        --output="${PROCESSED}/cll-time_course/results/single_cell_RNA/10_Seurat/${file}2.log"
done


## Run the dataset with the sample effect removed
echo $file
sbatch --job-name="Seurat Removed Sample Effect" --ntasks=12 --mem=180000 --partition=longq --time=3-00:00:00 \
    --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/10_2_Seurat_RmSampleEffect.R" \
    --output="${PROCESSED}/cll-time_course/results/single_cell_RNA/10_Seurat/RmSampleEffect.log"

echo $file
sbatch --job-name="Seurat Removed Sample Effect And CellCyle" --ntasks=12 --mem=180000 --partition=longq --time=3-00:00:00 \
    --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/10_2_Seurat_RmSampleEffect2.R" \
    --output="${PROCESSED}/cll-time_course/results/single_cell_RNA/10_Seurat/RmSampleEffect_CellCycle.log"

echo $file
sbatch --job-name="Seurat Removed CellCycle" --ntasks=12 --mem=180000 --partition=longq --time=3-00:00:00 \
    --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/10_2_Seurat_NoCellCycle.R" \
    --output="${PROCESSED}/cll-time_course/results/single_cell_RNA/10_Seurat/RmCellCycle.log"

# RAW ANALYSIS, not used because it fails often and doesn't make a large difference
# for file in ${array[@]}
# do
#     echo $file
#     sbatch --job-name="Seurat $file raw" --ntasks=1 --mem=180000 --partition=longq --time=3-00:00:00 \
#         --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/10_1_Seurat.R raw ${file}" \
#         --output="${PROCESSED}/cll-time_course/results/single_cell_RNA/10_Seurat_raw/${file}.log"
# done