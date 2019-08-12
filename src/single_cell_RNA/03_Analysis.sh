# SEURAT
# DOesn't work for some stupid reason
sbatch --job-name="Seurat analysis" --ntasks=1 --cpus-per-task=1 --mem=80000 --partition=longq \
    --wrap="module unload R/3.2.3; module load R/3.3.2; module list; Rscript /home/nfortelny/code/cll-time_course/src/single_cell_RNA/10_Seurat.R" \
    --output="${PROCESSED}/cll-time_course/results/single_cell_RNA/10_Seurat.log"

# Cellranger
sbatch --job-name="Cellranger analysis" --ntasks=1 --cpus-per-task=1 --mem=80000 --partition=longq \
    --wrap="module list; Rscript /home/nfortelny/code/cll-time_course/src/single_cell_RNA/11_CellrangerRKit.R" \
    --output="${PROCESSED}/cll-time_course/results/single_cell_RNA/11_Cellranger.log"