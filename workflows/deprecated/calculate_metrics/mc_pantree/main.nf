#!/usr/bin/env nextflow

process CALCULATE_METRICS {
    tag "${sample}"

    input:
    tuple val(sample), path(truth_vcf), path(query_vcf)
    val graph

    output:
    path "${sample}.${graph}.metrics.txt"

    script:
    """
    calculate_metrics.py \
        -s ${sample} \
        -g ${graph} \
        -c Dipcall \
        -o ${sample}.${graph}.metrics.txt \
        ${truth_vcf} \
        ${query_vcf}
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.truth_vcf), file(row.query_vcf)) }
        .set { samples_ch }

    CALCULATE_METRICS(samples_ch, params.graph)
}
