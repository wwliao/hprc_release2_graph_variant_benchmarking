#!/usr/bin/env nextflow

process RTG_VCFEVAL {
    tag "${sample}"

    input:
    tuple val(sample), path(truth_vcf), path(truth_tbi), path(query_vcf), path(query_tbi)
    path ref_sdf
    path region_bed

    output:
    path "${sample}"

    script:
    """
    rtg RTG_MEM=${task.memory.toGiga()}G vcfeval \
        --threads ${task.cpus} \
        --baseline ${truth_vcf} \
        --calls ${query_vcf} \
        --template ${ref_sdf} \
        --bed-regions ${region_bed} \
        --output ${sample} \
        --output-mode annotate \
        --all-records \
        --ref-overlap \
        --no-roc
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.truth_vcf), file(row.truth_tbi), file(row.query_vcf), file(row.query_tbi)) }
        .set { samples_ch }

    RTG_VCFEVAL(samples_ch, file(params.ref_sdf), file(params.region_bed))
}
