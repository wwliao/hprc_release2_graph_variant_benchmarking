#!/usr/bin/env nextflow

process RESOLVE_SV_TO_SEQ {
    tag "${sample}: ${caller}"

    input:
    tuple val(sample), val(caller), path(query_vcf)
    path ref_fasta
    path ref_fai
    path resolve_sv_to_seq_script

    output:
    tuple val(sample), val(caller), path("${query_vcf.baseName.tokenize('.')[0..2].join('.')}.resolved.vcf.gz"), path("${query_vcf.baseName.tokenize('.')[0..2].join('.')}.resolved.vcf.gz.tbi")

    script:
    """
    bcftools view -f 'PASS,.' -Ou ${query_vcf} \
        | bcftools norm -m -any -Oz -o ${query_vcf.baseName.tokenize('.')[0..2].join('.')}.biallelic.vcf.gz

    python3 ${resolve_sv_to_seq_script} \
        ${query_vcf.baseName.tokenize('.')[0..2].join('.')}.biallelic.vcf.gz \
        ${ref_fasta} ${query_vcf.baseName.tokenize('.')[0..2].join('.')}.resolved.vcf

    bcftools annotate -x INFO -Ou ${query_vcf.baseName.tokenize('.')[0..2].join('.')}.resolved.vcf \
        | bcftools sort --write-index=tbi -Oz -o ${query_vcf.baseName.tokenize('.')[0..2].join('.')}.resolved.vcf.gz
    """
}

process AARDVARK_MERGE {
    tag "${sample}: ${caller1}, ${caller2}, and ${caller3}"

    input:
    tuple val(sample), val(caller1), path(query_vcf1), path(query_tbi1), val(caller2), path(query_vcf2), path(query_tbi2), val(caller3), path(query_vcf3), path(query_tbi3)
    path ref_fasta
    path regions_bed

    output:
    tuple val(sample), path("${sample}.${regions_bed.baseName}.${caller1}_${caller2}_${caller3}.vcf.gz"), path("${sample}.${regions_bed.baseName}.${caller1}_${caller2}_${caller3}.vcf.gz.tbi")

    script:
    """
    aardvark merge \
        --threads ${task.cpus} \
        --reference ${ref_fasta} \
        --regions ${regions_bed} \
        --input-vcf ${query_vcf1} \
        --vcf-tag ${caller1} \
        --input-vcf ${query_vcf2} \
        --vcf-tag ${caller2} \
        --input-vcf ${query_vcf3} \
        --vcf-tag ${caller3} \
        --vcf-sample ${sample} \
        --min-variant-gap 1000 \
        --merge-strategy all \
        --conflict-select 0 \
        --output-vcfs .

    mv passing.vcf.gz \
        ${sample}.${regions_bed.baseName}.${caller1}_${caller2}_${caller3}.vcf.gz
    mv passing.vcf.gz.tbi \
        ${sample}.${regions_bed.baseName}.${caller1}_${caller2}_${caller3}.vcf.gz.tbi
    """
}

process AARDVARK_COMPARE {
    tag "${sample}: ${caller}"

    input:
    tuple val(sample), path(truth_vcf), path(truth_tbi), val(caller), path(query_vcf), path(query_tbi), val(min_gap)
    path ref_fasta
    path regions_bed

    output:
    tuple val(sample), path(truth_vcf), val(caller), path("${sample}.${regions_bed.baseName}.${caller}.compared.vcf.gz"), path("${sample}.${regions_bed.baseName}.${caller}.region_summary.tsv.gz")

    script:
    """
    aardvark compare \
        --threads ${task.cpus} \
        --reference ${ref_fasta} \
        --regions ${regions_bed} \
        --truth-vcf ${truth_vcf} \
        --query-vcf ${query_vcf} \
        --truth-sample ${sample} \
        --query-sample ${sample} \
        --min-variant-gap ${min_gap} \
        --output-dir . \
        --output-debug .

    mv truth.vcf.gz \
        ${sample}.${regions_bed.baseName}.${caller}.compared.vcf.gz
    mv region_summary.tsv.gz \
        ${sample}.${regions_bed.baseName}.${caller}.region_summary.tsv.gz
    """
}

process ADD_SOURCES_FROM_BASEPAIR_METRICS {
    tag "${sample}: ${caller}"
    
    input:
    tuple val(sample), path(truth_vcf), val(caller), path(compared_vcf), path(region_summary_tsv)
    path add_sources_script
    val threshold

    output:
    tuple val(sample), path("${compared_vcf.baseName.tokenize('.')[0..3].join('.')}.vcf.gz"), path("${compared_vcf.baseName.tokenize('.')[0..3].join('.')}.vcf.gz.tbi")

    script:
    def actual_threshold = caller == 'DeepVariant' ? 1.0 : threshold
    """
    python3 ${add_sources_script} \
        --input ${compared_vcf} \
        --truth ${truth_vcf} \
        --output ${compared_vcf.baseName.tokenize('.')[0..3].join('.')}.vcf.gz \
        --region-summary ${region_summary_tsv} \
        --sample ${sample} \
        --caller ${caller} \
        --threshold ${actual_threshold}
    """
}

process CONSOLIDATE_CALLERS {
    tag "${sample}"

    input:
    tuple val(sample), path(vcfs), path(tbis)
    path regions_bed
    path consolidate_callers_script

    output:
    tuple val(sample), path("${sample}.${regions_bed.baseName}.merged.vcf.gz"), path("${sample}.${regions_bed.baseName}.merged.vcf.gz.tbi")

    script:
    """
    python3 ${consolidate_callers_script} \
        --contigs-bed ${regions_bed} \
        --output ${sample}.${regions_bed.baseName}.merged.vcf.gz \
        ${vcfs}
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.caller, file(row.query_vcf), file(row.query_tbi)) }
        .set { samples_ch }

    RESOLVE_SV_TO_SEQ(samples_ch.map { sample, caller, query_vcf, query_tbi -> tuple(sample, caller, query_vcf) }, file(params.ref_fasta), file(params.ref_fai), file(params.resolve_sv_to_seq_script))

    def primary_callers = ['dipcall', 'PAV', 'longcallD'] as Set
    def primary_callers_order = ['dipcall', 'PAV', 'longcallD']  // Preserve order

    def branched_ch = RESOLVE_SV_TO_SEQ.out
        .branch { sample, caller, vcf, tbi ->
            primary: primary_callers.contains(caller)
            other: true
        }

    primary_grouped_ch = branched_ch.primary
        .groupTuple(by: 0)
        .filter { sample, callers, vcfs, tbis -> 
            callers.size() == primary_callers.size()  // Ensure all primary callers present
        }
        .map { sample, callers, vcfs, tbis ->
            def caller_map = [:]
            callers.eachWithIndex { caller, idx ->
                caller_map[caller] = [vcfs[idx], tbis[idx]]
            }
            
            def result = [sample]
            primary_callers_order.each { caller ->
                result.add(caller)
                result.add(caller_map[caller][0])
                result.add(caller_map[caller][1])
            }
            
            return tuple(*result)
        }

    AARDVARK_MERGE(primary_grouped_ch, file(params.ref_fasta), file(params.regions_bed))

    other_callers_ch = branched_ch.other

    // Combine AARDVARK_MERGE.out with each other caller for the same sample
    AARDVARK_MERGE.out
        .combine(other_callers_ch, by: 0)  // Combine by sample (index 0)
        .map { sample, truth_vcf, truth_tbi, caller, query_vcf, query_tbi ->
            // Set min_gap based on caller
            def min_gap = (caller == 'DeepVariant') ? 50 : 1000
            tuple(sample, truth_vcf, truth_tbi, caller, query_vcf, query_tbi, min_gap) }
        .set { other_paired_ch }

    AARDVARK_COMPARE(other_paired_ch, file(params.ref_fasta), file(params.regions_bed))

    ADD_SOURCES_FROM_BASEPAIR_METRICS(AARDVARK_COMPARE.out, file(params.add_sources_script), params.threshold)

    AARDVARK_MERGE.out
        .concat(ADD_SOURCES_FROM_BASEPAIR_METRICS.out)
        .groupTuple(by: 0)
        .set { grouped_vcfs_ch }

    CONSOLIDATE_CALLERS(grouped_vcfs_ch, file(params.regions_bed), file(params.consolidate_callers_script))
}
