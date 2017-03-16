#!/bin/bash

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
    sbatch --job-name="10X $file" --ntasks=4 --cpus-per-task=10 --mem=100000 --partition=longq \
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
    sbatch --job-name="10X $file" --ntasks=3 --cpus-per-task=3 --mem=50000 --partition=longq \
        --wrap="cellranger count --id=$file --transcriptome=/home/nfortelny/resources_nfortelny/10X_Genomics/refdata-cellranger-GRCh38-1.2.0/ --fastqs=${inPath}/${file} --cells=5000 --nopreflight" \
            --output="$file.log"
done