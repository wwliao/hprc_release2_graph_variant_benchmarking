#!/usr/bin/env nextflow

process BCFTOOLS_NORM_VIEW_ANNOTATE {
    tag "${sample}"

    input:
    tuple val(sample), path(region_bed), path(region_tbi)
    path vcf
    path ref_fasta
    path ref_fai

    output:
    path "${sample}.${ref_fasta.baseName}.mc_vcfwave.stratified.vcf.gz*"

    script:
    """
    bcftools view -a -I -s ${sample} -Ou ${vcf} \
        | bcftools norm -m -any -f ${ref_fasta} -Ou \
        | bcftools view -e 'GT="ref" | GT~"\\." | STRLEN(REF)>1000 | STRLEN(ALT)>1000 | ABS(STRLEN(ALT)-STRLEN(REF)) >= 50' \
        | bcftools annotate \
            -a ${region_bed} \
            -c CHROM,FROM,TO,INFO/REGIONS \
            -h <(echo '##INFO=<ID=REGIONS,Number=.,Type=String,Description="Stratification regions overlapping this variant">') \
            -l REGIONS:unique \
            --write-index=tbi \
            -Oz -o ${sample}.${ref_fasta.baseName}.mc_vcfwave.stratified.vcf.gz
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.region_bed), file(row.region_tbi)) }
        .set { samples_ch }

    BCFTOOLS_NORM_VIEW_ANNOTATE(samples_ch, file(params.vcf), file(params.ref_fasta), file(params.ref_fai))
}
