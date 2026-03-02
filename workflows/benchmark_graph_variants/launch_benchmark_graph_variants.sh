#!/bin/bash
#SBATCH --job-name=launch_benchmark_graph_variants
#SBATCH --output=launch_benchmark_graph_variants-%j.out
#SBATCH --partition=pi_hall
#SBATCH --constraint=nogpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=28-00:00:00

module purge
module load Nextflow/24.04.4

SAMPLE_SHEET="samplesheet.csv"
OUTDIR="results"
QUERY_VCF="/gpfs/gibbs/pi/ycgh/wl474/projects/hprc_r2/variant_benchmarking/graph_variants/traversal_decomposition/genotyping-pipelines/prepare-vcf-MC/hprc-v2.0-mc-grch38_filtered_ids_biallelic.svtyped.vcf.gz"
QUERY_TBI="${QUERY_VCF}.tbi"
REF_FASTA="/gpfs/gibbs/pi/ycgh/wl474/resources/reference_genomes/GRCh38_no_alt/GRCh38_no_alt.fa"
REGIONS_TSV="/gpfs/gibbs/pi/ycgh/wl474/projects/hprc_r2/variant_benchmarking/stratifications/giab_grch38_stratifications.v3.6.tsv"
MIN_GAP=1000

nextflow run main.nf \
    -ansi-log true \
    -profile mccleary \
    --sample_sheet ${SAMPLE_SHEET} \
    --query_vcf ${QUERY_VCF} \
    --query_tbi ${QUERY_TBI} \
    --ref_fasta ${REF_FASTA} \
    --regions_tsv ${REGIONS_TSV} \
    --min_gap ${MIN_GAP} \
    --outdir ${OUTDIR}
