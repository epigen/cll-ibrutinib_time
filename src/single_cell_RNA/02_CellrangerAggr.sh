#!/bin/bash

basePW=/home/nfortelny/projects_shared/cll-time_course/results/
mkdir ${basePW}/cellranger_count/
cd ${basePW}/cellranger_count/

sbatch --job-name="10x aggr PT VZS" --ntasks=2 --cpus-per-task=5 --mem=30000 --partition=longq \
    --wrap="cellranger aggr --id=PT_vs_VZS --csv=/home/nfortelny/code/10x_first_runs/metadata/aggregation_csv_d0_PT_VZS.csv" \
    --output="PT_vs_VZS.log"

sbatch --job-name="10x aggr PT timepoints" --ntasks=2 --cpus-per-task=5 --mem=30000 --partition=longq \
    --wrap="cellranger aggr --id=PT_timepoints --csv=/home/nfortelny/code/10x_first_runs/metadata/aggregation_csv_PT_d0to280.csv" \
    --output="PT_timepoints.log"

sbatch --job-name="10x aggr Timepoint_0" --ntasks=2 --cpus-per-task=5 --mem=30000 --partition=longq \
    --wrap="cellranger aggr --id=Timepoint_0 --csv=/home/nfortelny/code/10x_first_runs/metadata/aggregation_csv_timepoint_0.csv" \
    --output="Timepoint_0.log"