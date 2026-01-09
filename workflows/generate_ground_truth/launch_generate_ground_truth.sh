#!/bin/bash
#SBATCH --job-name=launch_generate_ground_truth
#SBATCH --output=launch_generate_ground_truth-%j.out
#SBATCH --partition=pi_hall
#SBATCH --constraint=nogpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=28-00:00:00

module purge
module load Nextflow/24.04.4

THRESHOLD=0.96
SAMPLE_SHEET="samplesheet.csv"
OUTDIR="results"
REF_FASTA="/gpfs/gibbs/pi/ycgh/wl474/resources/reference_genomes/GRCh38_no_alt/GRCh38_no_alt.fa"
REF_FAI="${REF_FASTA}.fai"
REGIONS_BED="GRCh38_no_alt.autosomes.bed"
RESOLVE_SV_TO_SEQ_SCRIPT="resolve_sv_to_sequence.py"
ADD_SOURCES_SCRIPT="add_sources_from_basepair_metrics.py"
CONSOLIDATE_CALLERS_SCRIPT="consolidate_callers.py"

nextflow run main.nf \
    -ansi-log true \
    -profile mccleary \
    --threshold ${THRESHOLD} \
    --sample_sheet ${SAMPLE_SHEET} \
    --outdir ${OUTDIR} \
    --ref_fasta ${REF_FASTA} \
    --ref_fai ${REF_FAI} \
    --regions_bed ${REGIONS_BED} \
    --resolve_sv_to_seq_script ${RESOLVE_SV_TO_SEQ_SCRIPT} \
    --add_sources_script ${ADD_SOURCES_SCRIPT} \
    --consolidate_callers_script ${CONSOLIDATE_CALLERS_SCRIPT}
