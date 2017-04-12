#!/bin/bash

#######################################################################################
### ASK BSF TO RUN MKFASTQ FOR THE FIRST SAMPLES AND THEN PUT IN THE RIGHT PATHS ######
#######################################################################################



###
### MKFASTQ
###
# sbatch --job-name="10X mkfastq" --ntasks=4 --cpus-per-task=10 --mem=80000 --partition=longq \
#     --wrap="cellranger mkfastq --run=/scratch/lab_bsf/projects/BSF_0284_HH5C7BBXX_runfolder/161230_ST-J00104_0085_AHH5C7BBXX/ --csv=/home/nfortelny/code/10x_first_runs/metadata/samplesheet_cellranger.csv --output-dir=~/projects_shared/cll-time_course/cellranger/" \
#     --output="/home/nfortelny/projects_shared/cll-time_course/cellranger_mkfastq.log"


###
### ORIGINAL CLL TIME SERIES
###
basePW=/home/nfortelny/projects_shared/cll-time_course/
mkdir ${basePW}/results/cellranger_count/
cd ${basePW}/results/cellranger_count/

ll ${basePW}/cellranger/HH5C7BBXX/

array=("PT_d0" "PT_d120" "PT_d280" "PT_d280_15bp" "PT_d280_25bp" "PT_d280_50bp" "VZS_d0" "VZS_d120")

for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --ntasks=4 --mem=100000 --partition=longq \
        --wrap="cellranger count --id=${file} --transcriptome=/home/nfortelny/resources_nfortelny/10X_Genomics/refdata-cellranger-GRCh38-1.2.0/ --fastqs=${basePW}/data/cellranger/HH5C7BBXX/${file}/ --cells=5000 --nopreflight" \
        --output="${file}.log"
done


###
### THE NEW SAMPLES FROM CLL TIME SERIES
###
inPath=/data/groups/lab_bsf/sequences/BSF_0289_HH5MJBBXX_l6/fastq_path/HH5MJBBXX
outPath=/home/nfortelny/projects_shared/cll-time_course/results/cellranger_count
cd $outPath
array=("LiveBulk_10x_FE_FE1_d0_10xTK078" "LiveBulk_10x_FE_FE7_d120_10xTK078" "LiveBulk_10x_KI_KI1_d0_10xTK078")
for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --ntasks=3 --mem=50000 --partition=longq \
        --wrap="cellranger count --id=$file --transcriptome=/home/nfortelny/resources_nfortelny/10X_Genomics/refdata-cellranger-GRCh38-1.2.0/ --fastqs=${inPath}/${file} --cells=5000 --nopreflight" \
            --output="$file.log"
done


###
### THE THIRD BATCH OF SAMPLES
###
inPath=/data/groups/lab_bsf/sequences/BSF_0309_HHNFKBBXX_l6_10x/fastq_path/HHNFKBBXX/
outPath=/home/nfortelny/projects_shared/cll-time_course/results/cellranger_count/
cd $outPath
array=("LiveBulk_10x_FE7_120d" "LiveBulk_10x_PBGY1_0d" "LiveBulk_10x_PBGY7_150d" "LiveBulk_10x_VZS7_120d")
for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --ntasks=9 --mem=180000 --partition=longq \
        --wrap="cellranger count --id=$file --transcriptome=/home/nfortelny/resources_nfortelny/10X_Genomics/refdata-cellranger-GRCh38-1.2.0/ --fastqs=${inPath}/${file} --cells=5000 --localcores=9 --nopreflight" \
            --output="$file.log"
done