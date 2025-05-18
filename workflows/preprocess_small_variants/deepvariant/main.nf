#!/usr/bin/env nextflow

process BCFTOOLS_NORM_VIEW_ANNOTATE {
    tag "${sample}"

    input:
    tuple val(sample), path(region_bed), path(region_tbi), path(deepvariant_vcf)
    path ref_fasta
    path ref_fai

    output:
    path "${deepvariant_vcf.baseName.tokenize('.')[0..2].join('.')}.stratified.vcf.gz*"

    script:
    """
    bcftools norm -m -any -f ${ref_fasta} -Ou ${deepvariant_vcf} \
        | bcftools view -f 'PASS' -e 'GT="ref" | GT~"\\." | STRLEN(REF)>1000 | STRLEN(ALT)>1000 | ABS(STRLEN(ALT)-STRLEN(REF)) >= 50' \
        | bcftools annotate \
            -a ${region_bed} \
            -c CHROM,FROM,TO,INFO/REGIONS \
            -h <(echo '##INFO=<ID=REGIONS,Number=.,Type=String,Description="Stratification regions overlapping this variant">') \
            -l REGIONS:unique \
            --write-index=tbi \
            -Oz -o ${deepvariant_vcf.baseName.tokenize('.')[0..2].join('.')}.stratified.vcf.gz
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.region_bed), file(row.region_tbi), file(row.deepvariant_vcf)) }
        .set { samples_ch }

    BCFTOOLS_NORM_VIEW_ANNOTATE(samples_ch, file(params.ref_fasta), file(params.ref_fai))
}
