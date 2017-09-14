################################
######## CLL over TIME #########
################################
module unload R/3.2.3
module load R/3.3.2


# Overtime together
array=("noIGHLK" "noRPstrict" "allDataBest")
array=("noIGHLK" "noRPstrict")
for file in ${array[@]}
do            
  echo $file
  sbatch --job-name="cll-time_course 13_4_OverTime $file" --ntasks=1 --mem=30000 --partition=shortq --time=08:00:00 \
      --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/13_4_OverTime_noRP.R ${file}" \
      --output="$CODEBASE/cll-time_course/src/single_cell_RNA/13_4_OverTime_${file}.log"
done

# MAGIC
sbatch --job-name="15_Magic.R" --ntasks=32 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/15_Magic.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/15_Magic.log"

# FSCLVM
sbatch --job-name="fscLVM" --ntasks=32 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/90_fscLVM.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/90_fscLVM.log"

sbatch --job-name="35_1_fscLVM.R" --ntasks=32 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/35_1_fscLVM.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/35_1_fscLVM.log"

sbatch --job-name="35_1_fscLVM_2.R" --ntasks=32 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/35_1_fscLVM_2.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/35_1_fscLVM_2.log"

# SIGNATURES
sbatch --job-name="30_1_Signatures_Overview.R" --ntasks=12 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/30_1_Signatures_Overview.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/30_1_Signatures_Overview.log"

sbatch --job-name="30_3_Signatures_Regressed.R" --ntasks=12 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/30_3_Signatures_Regressed.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/30_3_Signatures_Regressed.log"

sbatch --job-name="30_4_Signatures_RowNormalized.R" --ntasks=12 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/30_4_Signatures_RowNormalized.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/30_4_Signatures_RowNormalized.log"

sbatch --job-name="30_5_Signatures_Magic.R" --ntasks=12 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/30_5_Signatures_Magic.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/30_5_Signatures_Magic.log"

sbatch --job-name="30_6_Signatures_RowNormalized_Sampled.R" --ntasks=12 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/30_6_Signatures_RowNormalized_Sampled.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/30_6_Signatures_RowNormalized_Sampled.log"

sbatch --job-name="30_7_Signatures_RowNormalized_Sampled_exclGenes.R" --ntasks=12 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/30_7_Signatures_RowNormalized_Sampled_exclGenes.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/30_7_Signatures_RowNormalized_Sampled_exclGenes_2.log"

# Signatures for CNVS
sbatch --job-name="41_CNV_Signatures.R" --ntasks=12 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/41_CNV_Signatures.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/41_CNV_Signatures.log"

sbatch --job-name="41_2_CNV_Signatures_RandomOnce.R" --ntasks=12 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/41_2_CNV_Signatures_RandomOnce.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/41_2_CNV_Signatures_RandomOnce.log"

sbatch --job-name="41_3_CNV_Sigs_Once_min50.R" --ntasks=12 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/41_3_CNV_Sigs_Once_min50.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/41_3_CNV_Sigs_Once_min50.log"

sbatch --job-name="41_4_CNV_Sigs_Once_min20.R" --ntasks=12 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/41_4_CNV_Sigs_Once_min20.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/41_4_CNV_Sigs_Once_min20.log"

sbatch --job-name="41_5_CNV_Sigs_Once_min100.R" --ntasks=12 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/41_5_CNV_Sigs_Once_min100.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/41_5_CNV_Sigs_Once_min100.log"

sbatch --job-name="41_9_Random.R" --ntasks=12 --mem=180000 --partition=longq \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/41_9_Random.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/41_9_Random.log"


# SEURAT2
module unload gcc/4.8.2; module load gcc/7.1.0; module unload R/3.2.3; module load R/3.4.0
sbatch --job-name="51_ClusterAlignment.R" --ntasks=3 --mem=30000 --partition=mediumq --time=3:00:00 \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/51_ClusterAlignment.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/51_ClusterAlignment.log"

sbatch --job-name="52_2_Pseudtime_Magic.R" --ntasks=1 --mem=180000 --partition=longq --time=3:00:00 \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/52_2_Pseudtime_Magic.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/52_2_Pseudtime_Magic.log"

module unload gcc/4.8.2; module load gcc/7.1.0; module unload R/3.2.3; module load R/3.4.0
sbatch --job-name="52_3_Monocle.R" --ntasks=12 --mem=180000 --partition=longq --time=3-00:00:00 \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/52_3_Monocle.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/52_3_Monocle.log"

module unload R/3.2.3; module load R/3.3.2
sbatch --job-name="53_GroupCorrelations.R" --ntasks=1 --mem=70000 --partition=develop --time=08:00:00 \
    --wrap="Rscript $CODEBASE/cll-time_course/src/single_cell_RNA/53_GroupCorrelations.R" \
    --output="$CODEBASE/cll-time_course/src/single_cell_RNA/53_GroupCorrelations.log"