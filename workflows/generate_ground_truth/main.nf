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
    def base_name = query_vcf.baseName.tokenize('.')[0..2].join('.')
    def setgt_cmd = caller == "longcallD"
        ? """
        bcftools +setGT ${base_name}.resolved.vcf.gz --write-index=tbi -Oz -o ${base_name}.resolved.tmp.vcf.gz -- -t a -n u
        mv ${base_name}.resolved.tmp.vcf.gz ${base_name}.resolved.vcf.gz
        mv ${base_name}.resolved.tmp.vcf.gz.tbi ${base_name}.resolved.vcf.gz.tbi
        """
        : ""
    """
    bcftools view -f 'PASS,.' -Ou ${query_vcf} \
        | bcftools norm -m -any -Oz -o ${base_name}.biallelic.vcf.gz

    python3 ${resolve_sv_to_seq_script} \
        ${base_name}.biallelic.vcf.gz \
        ${ref_fasta} ${base_name}.resolved.vcf

    bcftools annotate -x INFO -Ou ${base_name}.resolved.vcf \
        | bcftools sort --write-index=tbi -Oz -o ${base_name}.resolved.vcf.gz

    ${setgt_cmd}
    """
}

process AARDVARK_MERGE {
    tag "${sample}: ${caller1}, ${caller2}, and ${caller3}"

    input:
    tuple val(sample), path(dipcall_bed), val(caller1), path(query_vcf1), path(query_tbi1), val(caller2), path(query_vcf2), path(query_tbi2), val(caller3), path(query_vcf3), path(query_tbi3)
    path ref_fasta

    output:
    tuple val(sample), path(dipcall_bed), path("${sample}.${ref_fasta.baseName}.autosomes.${caller1}_${caller2}_${caller3}.vcf.gz"), path("${sample}.${ref_fasta.baseName}.autosomes.${caller1}_${caller2}_${caller3}.vcf.gz.tbi")

    script:
    """
    awk '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2])\$/' ${dipcall_bed} > ${dipcall_bed.baseName}.autosomes.bed
    aardvark merge \
        --threads ${task.cpus} \
        --reference ${ref_fasta} \
        --regions ${dipcall_bed.baseName}.autosomes.bed \
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
        ${sample}.${ref_fasta.baseName}.autosomes.${caller1}_${caller2}_${caller3}.vcf.gz
    mv passing.vcf.gz.tbi \
        ${sample}.${ref_fasta.baseName}.autosomes.${caller1}_${caller2}_${caller3}.vcf.gz.tbi
    """
}

process SMALL_VARIANT_FILTER {
    tag "${sample}: ${small_caller}"

    input:
    tuple val(sample), val(small_caller), path(small_resolved_vcf), path(small_resolved_tbi)

    output:
    tuple val(sample), val(small_caller), path("${small_resolved_vcf.baseName.tokenize('.')[0..3].join('.')}.lt50.vcf.gz"), path("${small_resolved_vcf.baseName.tokenize('.')[0..3].join('.')}.lt50.vcf.gz.tbi")

    script:
    """
    bcftools view -i 'abs(strlen(REF)-strlen(ALT)) < 50' -Ou ${small_resolved_vcf} \
        | bcftools sort --write-index=tbi -Oz -o ${small_resolved_vcf.baseName.tokenize('.')[0..3].join('.')}.lt50.vcf.gz
    """
}

process BCFTOOLS_CONCAT {
    tag "${sample}: ${small_caller} and ${sv_caller}"

    input:
    tuple val(sample), val(small_caller), path(small_resolved_vcf), path(small_resolved_tbi), val(sv_caller), path(sv_resolved_vcf), path(sv_resolved_tbi)

    output:
    tuple val(sample), val(sv_caller), path("${sv_resolved_vcf.baseName.tokenize('.')[0..3].join('.')}.joint.vcf.gz"), path("${sv_resolved_vcf.baseName.tokenize('.')[0..3].join('.')}.joint.vcf.gz.tbi")

    script:
    """
    bcftools view -i 'abs(strlen(REF)-strlen(ALT)) >= 50' -Ou ${sv_resolved_vcf} \
        | bcftools sort --write-index=tbi -Oz -o ${sv_resolved_vcf.baseName.tokenize('.')[0..3].join('.')}.ge50.vcf.gz

    bcftools concat -a --write-index=tbi -Oz -o ${sv_resolved_vcf.baseName.tokenize('.')[0..3].join('.')}.joint.vcf.gz \
        ${small_resolved_vcf} \
        ${sv_resolved_vcf.baseName.tokenize('.')[0..3].join('.')}.ge50.vcf.gz
    """
}

process AARDVARK_COMPARE {
    tag "${sample}: ${caller}"

    input:
    tuple val(sample), path(dipcall_bed), path(truth_vcf), path(truth_tbi), val(caller), path(query_vcf), path(query_tbi), val(min_gap)
    path ref_fasta

    output:
    tuple val(sample), path(truth_vcf), val(caller), path("${sample}.${ref_fasta.baseName}.autosomes.${caller}.compared.vcf.gz"), path("${sample}.${ref_fasta.baseName}.autosomes.${caller}.region_summary.tsv.gz")

    script:
    """
    awk '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2])\$/' ${dipcall_bed} > ${dipcall_bed.baseName}.autosomes.bed
    aardvark compare \
        --threads ${task.cpus} \
        --reference ${ref_fasta} \
        --regions ${dipcall_bed.baseName}.autosomes.bed \
        --truth-vcf ${truth_vcf} \
        --query-vcf ${query_vcf} \
        --truth-sample ${sample} \
        --query-sample ${sample} \
        --min-variant-gap ${min_gap} \
        --output-dir . \
        --output-debug .

    mv truth.vcf.gz \
        ${sample}.${ref_fasta.baseName}.autosomes.${caller}.compared.vcf.gz
    mv region_summary.tsv.gz \
        ${sample}.${ref_fasta.baseName}.autosomes.${caller}.region_summary.tsv.gz
    """
}

process ADD_CALLERS_FROM_BASEPAIR_METRICS {
    tag "${sample}: ${caller}"
    
    input:
    tuple val(sample), path(truth_vcf), val(caller), path(compared_vcf), path(region_summary_tsv)
    path add_callers_script
    val primary_threshold
    val other_threshold

    output:
    tuple val(sample), path("${compared_vcf.baseName.tokenize('.')[0..3].join('.')}.vcf.gz"), path("${compared_vcf.baseName.tokenize('.')[0..3].join('.')}.vcf.gz.tbi")

    script:
    // Assign thresholds based on caller
    def primary_callers = ['dipcall', 'PAV', 'longcallD']
    def actual_threshold = caller == 'DeepVariant' ? 1.0 :
                          (primary_callers.contains(caller) ? primary_threshold : other_threshold)
    """
    python3 ${add_callers_script} \
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
    def samples_ch = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.caller, file(row.dipcall_bed), file(row.query_vcf), file(row.query_tbi)) }

    RESOLVE_SV_TO_SEQ(samples_ch.map { sample, caller, dipcall_bed, query_vcf, query_tbi -> tuple(sample, caller, query_vcf) }, file(params.ref_fasta), file(params.ref_fai), file(params.resolve_sv_to_seq_script))

    def primary_callers = ['dipcall', 'PAV', 'longcallD'] as Set
    def primary_callers_order = ['dipcall', 'PAV', 'longcallD']  // Preserve order

    def branched_ch = RESOLVE_SV_TO_SEQ.out
        .branch { sample, caller, vcf, tbi ->
            primary: primary_callers.contains(caller)
            small: caller == 'DeepVariant'
            other: true
        }

    def primary_grouped_ch = samples_ch
        .join(branched_ch.primary, by: [0, 1])
        .map {sample, caller, dipcall_bed, query_vcf, query_tbi, vcf, tbi ->
            tuple(sample, dipcall_bed, caller, vcf, tbi)
        }
        .groupTuple(by: [0, 1])
        .filter { sample, dipcall_bed, callers, vcfs, tbis ->
            callers.size() == primary_callers.size()  // Ensure all primary callers present
        }
        .map { sample, dipcall_bed, callers, vcfs, tbis ->
            def caller_map = [:]
            callers.eachWithIndex { caller, idx ->
                caller_map[caller] = [vcfs[idx], tbis[idx]]
            }

            def result = [sample, dipcall_bed]
            primary_callers_order.each { caller ->
                result.add(caller)
                result.add(caller_map[caller][0])
                result.add(caller_map[caller][1])
            }

            return tuple(*result)
        }

    AARDVARK_MERGE(primary_grouped_ch, file(params.ref_fasta))

    SMALL_VARIANT_FILTER(branched_ch.small)

    def small_caller_ch = SMALL_VARIANT_FILTER.out
    def other_callers_ch = branched_ch.other

    def concat_ch = small_caller_ch
        .combine(other_callers_ch, by: 0) // Combine by sample (index 0)

    BCFTOOLS_CONCAT(concat_ch)

    // Combine AARDVARK_MERGE.out with each caller for the same sample
    def paired_ch = AARDVARK_MERGE.out
        .combine(branched_ch.primary.concat(small_caller_ch, BCFTOOLS_CONCAT.out), by: 0)  // Combine by sample (index 0)
        .map { sample, dipcall_bed, truth_vcf, truth_tbi, caller, query_vcf, query_tbi ->
        // Set min_gap based on caller
        def min_gap = (caller == 'DeepVariant') ? 50 : 1000
        tuple(sample, dipcall_bed, truth_vcf, truth_tbi, caller, query_vcf, query_tbi, min_gap) }

    AARDVARK_COMPARE(paired_ch, file(params.ref_fasta))

    ADD_CALLERS_FROM_BASEPAIR_METRICS(AARDVARK_COMPARE.out, file(params.add_callers_script), params.primary_threshold, params.other_threshold)

    def grouped_vcfs_ch = AARDVARK_MERGE.out
        .map { sample, dipcall_bed, vcf, tbi -> tuple(sample, vcf, tbi) }
        .concat(ADD_CALLERS_FROM_BASEPAIR_METRICS.out)
        .groupTuple(by: 0)

    CONSOLIDATE_CALLERS(grouped_vcfs_ch, file(params.regions_bed), file(params.consolidate_callers_script))
}
