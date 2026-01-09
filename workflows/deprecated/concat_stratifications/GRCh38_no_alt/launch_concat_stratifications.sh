#!/bin/bash
#SBATCH --job-name=launch_concat_stratifications
#SBATCH --output=launch_concat_stratifications-%j.out
#SBATCH --partition=pi_hall
#SBATCH --constraint=nogpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=7-00:00:00

module purge
module load Nextflow/24.04.4

SAMPLE_SHEET="samplesheet.csv"
PANMASK_A="/vast/palmer/pi/hall/wl474/projects/hprc_r2/variant_benchmarking/stratifications/panmask/hg38.pm151a-v2.easy.bed.gz"
PANMASK_B="/vast/palmer/pi/hall/wl474/projects/hprc_r2/variant_benchmarking/stratifications/panmask/hg38.pm151b-v2.easy.bed.gz"
GIAB_EASY="/vast/palmer/pi/hall/wl474/projects/hprc_r2/variant_benchmarking/stratifications/giab_v3.6/GRCh38_notinalldifficultregions.bed.gz"
GIAB_LOWMAPSD="/vast/palmer/pi/hall/wl474/projects/hprc_r2/variant_benchmarking/stratifications/giab_v3.6/GRCh38_alllowmapandsegdupregions.bed.gz"
GIAB_REPEAT="/vast/palmer/pi/hall/wl474/projects/hprc_r2/variant_benchmarking/stratifications/giab_v3.6/GRCh38_AllTandemRepeats_ge101bp_slop5.bed.gz"
GIAB_OTHERDIFF="/vast/palmer/pi/hall/wl474/projects/hprc_r2/variant_benchmarking/stratifications/giab_v3.6/GRCh38_allOtherDifficultregions.bed.gz"
OUTDIR="results"

nextflow run main.nf \
    -ansi-log true \
    -profile mccleary \
    --sample_sheet ${SAMPLE_SHEET} \
    --panmask_a ${PANMASK_A} \
    --panmask_b ${PANMASK_B} \
    --giab_easy ${GIAB_EASY} \
    --giab_lowmapsd ${GIAB_LOWMAPSD} \
    --giab_repeat ${GIAB_REPEAT} \
    --giab_otherdiff ${GIAB_OTHERDIFF} \
    --outdir ${OUTDIR}
