#!/usr/bin/env nextflow

process AWK_SORT_BGZIP_TABIX {
    tag "${sample}"

    input:
    tuple val(sample), path(dipcall)
    path panmask_a
    path panmask_b
    path giab_easy
    path giab_lowmapsd
    path giab_repeat
    path giab_otherdiff

    output:
    path "${sample}.strat_regions.bed.gz"
    path "${sample}.strat_regions.bed.gz.tbi"

    script:
    """
    awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "Dipcall"}' ${dipcall} > tmp.strat_regions.bed
    bgzip -dc ${panmask_a} | awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "PM_151a"}' >> tmp.strat_regions.bed
    bgzip -dc ${panmask_b} | awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "PM_151b"}' >> tmp.strat_regions.bed
    bgzip -dc ${giab_easy} | awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "GIAB_Easy"}' >> tmp.strat_regions.bed
    bgzip -dc ${giab_lowmapsd} | awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "GIAB_LowMapSD"}' >> tmp.strat_regions.bed
    bgzip -dc ${giab_repeat} | awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "GIAB_Repeat"}' >> tmp.strat_regions.bed
    bgzip -dc  ${giab_otherdiff} | awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "GIAB_OtherDiff"}' >> tmp.strat_regions.bed

    sort -k1,1 -k2,2n -k3,3n tmp.strat_regions.bed | bgzip -c > ${sample}.strat_regions.bed.gz && tabix -p bed ${sample}.strat_regions.bed.gz
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.dipcall_bed)) }
        .set { samples_ch }

    AWK_SORT_BGZIP_TABIX(samples_ch, file(params.panmask_a), file(params.panmask_b), file(params.giab_easy), file(params.giab_lowmapsd), file(params.giab_repeat), file(params.giab_otherdiff))
}
