#!/usr/bin/env nextflow

process BCFTOOLS_FILTER {
    tag "${sample}"

    input:
    tuple val(sample), path(truth_vcf), path(truth_tbi)

    output:
    tuple val(sample), path("${truth_vcf.baseName.tokenize('.')[0..3].join('.')}.filtered.vcf.gz"), path("${truth_vcf.baseName.tokenize('.')[0..3].join('.')}.filtered.vcf.gz.tbi")

    script:
    """
    bcftools view -i 'INFO/NCALLERS>=2' ${truth_vcf} \
        | sed 's/SVTYPE=CPX/SVTYPE=INS/g' \
        | bgzip > ${truth_vcf.baseName.tokenize('.')[0..3].join('.')}.filtered.vcf.gz

    tabix -p vcf ${truth_vcf.baseName.tokenize('.')[0..3].join('.')}.filtered.vcf.gz
    """

}

process AARDVARK_COMPARE {
    tag "${sample}"

    input:
    tuple val(sample), path(dipcall_bed), path(truth_vcf), path(truth_tbi)
    path query_vcf
    path query_tbi
    path ref_fasta
    path regions_tsv
    val min_gap

    output:
    path "${sample}"

    script:
    """
    awk '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2])\$/' ${dipcall_bed} > ${dipcall_bed.baseName}.autosomes.bed
    aardvark compare \
        --threads ${task.cpus} \
        --reference ${ref_fasta} \
        --truth-vcf ${truth_vcf} \
        --query-vcf ${query_vcf} \
        --compare-label ${sample} \
        --truth-sample ${sample} \
        --query-sample ${sample} \
        --regions ${dipcall_bed.baseName}.autosomes.bed \
        --stratification ${regions_tsv} \
        --output-dir ${sample} \
        --output-debug ${sample} \
        --min-variant-gap ${min_gap}
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.dipcall_bed), file(row.truth_vcf), file(row.truth_tbi)) }
        .set { samples_ch }

    BCFTOOLS_FILTER(samples_ch.map { sample, dipcall_bed, truth_vcf, truth_tbi -> tuple(sample, truth_vcf, truth_tbi) })

    // Combine samples_ch with BCFTOOLS_FILTER output
    filtered_samples_ch = samples_ch
        .join(BCFTOOLS_FILTER.out)
        .map { sample, dipcall_bed, truth_vcf, truth_tbi, filtered_truth_vcf, filtered_truth_tbi -> tuple(sample, dipcall_bed, filtered_truth_vcf, filtered_truth_tbi) } 

    AARDVARK_COMPARE(filtered_samples_ch, file(params.query_vcf), file(params.query_tbi), file(params.ref_fasta), file(params.regions_tsv), params.min_gap)
}
