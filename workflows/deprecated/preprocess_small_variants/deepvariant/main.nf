#!/usr/bin/env nextflow

process BCFTOOLS_NORM {
    tag "${sample}"

    input:
    tuple val(sample), path(deepvariant_vcf)
    path ref_fasta
    path ref_fai

    output:
    tuple val(sample), path("${deepvariant_vcf.baseName.tokenize('.')[0..2].join('.')}.normalized.vcf.gz"), path("${deepvariant_vcf.baseName.tokenize('.')[0..2].join('.')}.normalized.vcf.gz.tbi")

    script:
    """
    bcftools norm -m -any -f ${ref_fasta} -Ou ${deepvariant_vcf} \
        | bcftools view -f 'PASS' -i 'GT~"1" & STRLEN(REF)<=1000 & STRLEN(ALT)<=1000' -Ou \
        | bcftools sort --write-index=tbi -Oz -o ${deepvariant_vcf.baseName.tokenize('.')[0..2].join('.')}.normalized.vcf.gz
    """
}

process RTG_VCFDECOMPOSE {
    tag "${sample}"

    input:
    tuple val(sample), path(normalized_vcf), path(normalized_tbi)
    path ref_sdf

    output:
    tuple val(sample), path("${normalized_vcf.baseName.tokenize('.')[0..2].join('.')}.decomposed.vcf.gz")

    script:
    """
    rtg vcfdecompose --break-indels --break-mnps -t ${ref_sdf} -i ${normalized_vcf} -o ${normalized_vcf.baseName.tokenize('.')[0..2].join('.')}.decomposed.vcf.gz 
    """
}

process BCFTOOLS_ANNOTATE {
    tag "${sample}"

    input:
    tuple val(sample), path(region_bed), path(region_tbi), path(decomposed_vcf)

    output:
    path "${decomposed_vcf.baseName.tokenize('.')[0..2].join('.')}.stratified.vcf.gz*"

    script:
    """
    bcftools sort -Ou ${decomposed_vcf} \
        | bcftools norm -d exact -Ou \
        | bcftools view -i 'ABS(STRLEN(ALT)-STRLEN(REF)) < 50' -Ou \
        | bcftools annotate \
            -a ${region_bed} \
            -c CHROM,FROM,TO,INFO/REGIONS \
            -h <(echo '##INFO=<ID=REGIONS,Number=.,Type=String,Description="Stratification regions overlapping this variant">') \
            -l REGIONS:unique \
            --write-index=tbi \
            -Oz -o ${decomposed_vcf.baseName.tokenize('.')[0..2].join('.')}.stratified.vcf.gz
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.region_bed), file(row.region_tbi), file(row.deepvariant_vcf)) }
        .set { samples_ch }

    BCFTOOLS_NORM(samples_ch.map { sample, region_bed, region_tbi, deepvariant_vcf -> tuple(sample, deepvariant_vcf) }, file(params.ref_fasta), file(params.ref_fai))
    RTG_VCFDECOMPOSE(BCFTOOLS_NORM.out, file(params.ref_sdf))

    // Combine samples_ch with RTG_VCFDECOMPOSE output
    annotate_input = samples_ch
        .join(RTG_VCFDECOMPOSE.out)
        .map { sample, region_bed, region_tbi, deepvariant_vcf, decomposed_vcf -> tuple(sample, region_bed, region_tbi, decomposed_vcf) } 

    BCFTOOLS_ANNOTATE(annotate_input)
}
