#!/bin/bash
#SBATCH --job-name=launch_preprocess_small_variants
#SBATCH --output=launch_preprocess_small_variants-%j.out
#SBATCH --partition=pi_hall
#SBATCH --constraint=nogpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=1-00:00:00

module purge
module load Nextflow/24.04.4

SAMPLE_SHEET="samplesheet.csv"
VCF="/vast/palmer/pi/hall/wl474/projects/hprc_r2/variant_benchmarking/graph_vcfs/mc/hprc-v2.0-mc-grch38.wave.vcf.gz"
REF_FASTA="/gpfs/gibbs/pi/ycgh/wl474/resources/reference_genomes/GRCh38_no_alt/GRCh38_no_alt.fa"
REF_FAI="/gpfs/gibbs/pi/ycgh/wl474/resources/reference_genomes/GRCh38_no_alt/GRCh38_no_alt.fa.fai"
REF_SDF="/gpfs/gibbs/pi/ycgh/wl474/resources/reference_genomes/GRCh38_no_alt/SDF"
OUTDIR="results"

nextflow run main.nf \
    -ansi-log true \
    -profile mccleary \
    --sample_sheet ${SAMPLE_SHEET} \
    --vcf ${VCF} \
    --ref_fasta ${REF_FASTA} \
    --ref_fai ${REF_FAI} \
    --ref_sdf ${REF_SDF} \
    --outdir ${OUTDIR}
