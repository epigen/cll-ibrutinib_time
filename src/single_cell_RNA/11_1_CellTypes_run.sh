module unload R/3.2.3
module load R/3.3.2

# datasets=("allDataBest" "allDataBest_noIGH" "allDataBest_NoDownSampling_noIGH")
datasets=("allDataBest_NoDownSampling_noIGH")
cells=("Monos" "NurseLikeCells" "Bcells" "NKcells" "Tcells1" "Tcells2" "Tcells3")

for dataset in ${datasets[@]}
do   
  for cell in ${cells[@]}
  do   
    echo $dataset
    echo $cell
    outFile="${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes/${dataset}_${cell}.log"
    echo $outFile
    sbatch --job-name="Seurat $dataset $cell" --ntasks=32 --mem=180000 --partition=longq --time=3-00:00:00 \
        --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_Script.R ${dataset} ${cell}" \
        --output=$outFile
  done
done

sbatch --job-name="Seurat CD8" --ntasks=32 --mem=180000 --partition=longq --time=3-00:00:00 \
    --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_ByMarkerGenes.R CD8" \
    --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes/allDataBest_NoDownSampling_noIGH_CD8.log


  


######
######
######

module unload R/3.2.3
module load R/3.3.2

# datasets=("allDataBest" "allDataBest_noIGH" "allDataBest_NoDownSampling_noIGH")
datasets=("allDataBest_NoDownSampling_noIGH")
cells=("Monos" "NurseLikeCells" "Bcells" "NKcells" "Tcells1" "Tcells2" "Tcells3")

for dataset in ${datasets[@]}
do   
  for cell in ${cells[@]}
  do   
    echo $dataset
    echo $cell
    outFile="${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes/${dataset}_${cell}.log"
    echo $outFile
    sbatch --job-name="Seurat $dataset $cell" --ntasks=1 --mem=80000 --partition=shortq --time=10:00:00 \
        --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_Script.R ${dataset} ${cell}" \
        --output=$outFile
  done
done

sbatch --job-name="Seurat CD8" --ntasks=1 --mem=80000 --partition=develop --time=05:00:00 \
    --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_ByMarkerGenes.R CD8" \
    --output=${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes/allDataBest_NoDownSampling_noIGH_CD8.log