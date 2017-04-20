module unload R/3.2.3
module load R/3.3.2

datasets=("allDataBest" "allDataBest_noIGH" "allDataBest_NoDownSampling_noIGH")
cells=("Monos" "NurseLikeCells" "Bcells" "Tcells")
datasets=("allDataBest")


for dataset in ${datasets[@]}
do   
  for cell in ${cells[@]}
  do   
    echo $dataset
    echo $cell
    outFile="${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes/${dataset}_${cell}.log"
    echo $outFile
    sbatch --job-name="Seurat $dataset $celltype" --ntasks=1 --mem=180000 --partition=longq --time=3-00:00:00 \
        --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_Script.R ${dataset} ${cell}" \
        --output=$outFile
  done
done



#dataset="allDataBest"
#cell="NurseLikeCells"
#outFile="${PROCESSED}/cll-time_course/results/single_cell_RNA/11_CellTypes/${dataset}_${cell}.log"
#echo $outFile
#sbatch --job-name="Seurat $dataset $celltype" --ntasks=1 --mem=180000 --partition=longq --time=3-00:00:00 \
#    --wrap="Rscript ${CODEBASE}/cll-time_course/src/single_cell_RNA/11_2_CellTypes_Script.R ${dataset} ${cell}" \
#    --output=$outFile