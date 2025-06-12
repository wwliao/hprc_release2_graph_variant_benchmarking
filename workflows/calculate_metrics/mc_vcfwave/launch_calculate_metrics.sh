#!/bin/bash
#SBATCH --job-name=launch_calculate_metrics
#SBATCH --output=launch_calculate_metrics-%j.out
#SBATCH --partition=pi_hall
#SBATCH --constraint=nogpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=1-00:00:00

module purge
module load Nextflow/24.04.4

SAMPLE_SHEET="samplesheet.csv"
GRAPH="mc_vcfwave"
OUTDIR="results"

nextflow run main.nf \
    -ansi-log true \
    -profile mccleary \
    --sample_sheet ${SAMPLE_SHEET} \
    --graph ${GRAPH} \
    --outdir ${OUTDIR}
