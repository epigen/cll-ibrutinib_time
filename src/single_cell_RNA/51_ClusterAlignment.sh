module unload gcc/4.8.2
module load gcc/7.1.0
module unload R/3.2.3
module load R/3.4.0

pats=("PT" "PT2" "FE" "VZS" "PBGY")

for pat in ${pats[@]}
do
  sbatch --job-name="51_ClusterAlignment.R ${pat}" --ntasks=3 --mem=60000 --partition=mediumq --time=1-12:00:00 \
  --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/51_ClusterAlignment.R ${pat}" \
  --output="$CODEBASE/cll-time_course/src/single_cell_RNA/51_ClusterAlignment_${pat}.log"
done