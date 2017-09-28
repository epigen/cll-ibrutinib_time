#!/bin/bash

basePW=/home/nfortelny/projects_shared/cll-time_course/results/
mkdir ${basePW}/cellranger_count/
cd ${basePW}/cellranger_count/

sbatch --job-name="10x aggr PT VZS" --ntasks=2 --cpus-per-task=5 --mem=30000 --partition=longq \
    --wrap="cellranger aggr --id=PT_vs_VZS --csv=$CODEBASE/cll-time_course/metadata/Aggregate_csv_d0_PT_VZS.csv" \
    --output="PT_vs_VZS.log"

sbatch --job-name="10x aggr PT timepoints" --ntasks=2 --cpus-per-task=5 --mem=30000 --partition=longq \
    --wrap="cellranger aggr --id=PT_timepoints --csv=$CODEBASE/cll-time_course/metadata/Aggregate_csv_PT_d0to280.csv" \
    --output="PT_timepoints.log"

sbatch --job-name="10x aggr Timepoint_0" --ntasks=2 --cpus-per-task=5 --mem=30000 --partition=longq \
    --wrap="cellranger aggr --id=Timepoint_0 --csv=$CODEBASE/cll-time_course/metadata/Aggregate_csv_timepoint_0.csv" \
    --output="Timepoint_0.log"

sbatch --job-name="10x aggr cll_time_series" --ntasks=9 --mem=180000 --partition=longq \
    --wrap="cellranger aggr --id=allData --csv=$CODEBASE/cll-time_course/metadata/Aggregate_all.csv" \
    --output="allData.log"

sbatch --job-name="10x aggr cll_time_series" --ntasks=9 --mem=180000 --partition=longq \
    --wrap="cellranger aggr --id=allDataBest --csv=$CODEBASE/cll-time_course/metadata/Aggregate_best.csv" \
    --output="allDataBest.log"

sbatch --job-name="10x aggr cll_time_series" --ntasks=9 --mem=180000 --partition=longq \
    --wrap="cellranger aggr --id=allDataBest_NoDownSampling --normalize=none --csv=$CODEBASE/cll-time_course/metadata/Aggregate_best.csv" \
    --output="allDataBest_NoDownSampling.log"

sbatch --job-name="10x aggr Day 30 dev" --ntasks=3 --mem=180000 --partition=develop \
    --wrap="cellranger aggr --id=inclDay30_2 --normalize=none --csv=$CODEBASE/cll-time_course/metadata/Aggregate_inclDay30.csv" \
    --output="inclDay30_2.log"