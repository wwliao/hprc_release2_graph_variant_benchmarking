# Graph Variant Benchmarking for the HPRC Release 2

This repository contains Nextflow workflows for graph variant benchmarking across 230 of 231 individuals from the Human Pangenome Reference Consortium Release 2 (HPRC R2). HG00272 was excluded due to a likely large-scale misassembly on chromosome X and was therefore not included in the HPRC R2 pangenome graph.

The workflows in this repository support:

- Construction of merged callsets using multiple variant callers 
- Derivation of per-sample joint ground truth callsets and benchmarking of variants derived from the pangenome graph (hereafter referred to as *graph variants*)

## Data Source

Merged callsets used in this repository are provided in the [variant calling repository](https://github.com/wwliao/hprc_release2_variant_calling), with download paths listed in [merged_callsets.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/merged_callsets.index.csv). An additional index file, [variant_callsets.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/variant_callsets.index.csv), provides download paths for related resources such as dipcall BED files.

## Aardvark for Variant Comparison

Variant comparison is challenging because different methods often produce substantially different representations of the same underlying sequence changes. This issue is particularly pronounced for graph variants, whose representations can deviate significantly from those produced by traditional callers. To minimize bias introduced by representation differences, we use [Aardvark](https://github.com/PacificBiosciences/aardvark) for variant comparison.

Aardvark operates at the haplotype sequence level by reconstructing local haplotypes through the joint consideration of small variants and structural variants (SVs). This enables direct comparison of the underlying sequences rather than relying on direct matching of VCF records, making the comparison robust to representation differences.

Formally, let $R$ denote the reference sequence in a region, $T$ the truth haplotype sequence, and $Q$ the query haplotype sequence. Let $d(A, B)$ denote the edit distance between sequences $A$ and $B$. Under this formulation, the relationships between true positives (TP), false negatives (FN), and false positives (FP) can be expressed as:

$$ TP + FN = d(R, T) $$
$$ TP + FP = d(R, Q) $$
$$ FN + FP = d(T, Q) $$

Solving for $TP$ gives:

$$ TP = \frac{d(R, T) + d(R, Q) - d(T, Q)}{2} $$

Basepair-level recall and precision are defined as:

$$ \text{Recall} = \frac{TP}{d(R, T)} $$
$$ \text{Precision} = \frac{TP}{d(R, Q)} $$

In addition to basepair-level metrics, Aardvark derives variant-level metrics by mapping haplotype sequence agreement back to individual variants. Specifically, it identifies the maximal set of alleles that can be incorporated into the reference sequence such that the resulting truth and query haplotypes are identical at the basepair level. Variants contributing to this consistent haplotype reconstruction are labeled as TP, while the remaining variants are assigned as FN or FP depending on their source. Importantly, variant representations do not need to match, as long as the resulting haplotype sequences are identical.

In this work, we apply Aardvark in two contexts. First, during callset merging, we use haplotype sequence comparisons to identify consistent variants across callers and construct a merged callset for each sample. Second, during benchmarking of graph variants, we derive the joint ground truth callset from the merged callset for each sample by requiring support from at least two callers, and compare query and truth callsets through haplotype sequence reconstruction. Specifically, basepair-level metrics define agreement during merging, whereas benchmarking is summarized using variant-level metrics derived from these haplotype sequence comparisons.

## Merged Callset Construction

Existing merged callsets typically separate small variants and SVs, for example:

- Small variant sets derived from DeepVariant
- SV sets derived from PAV, supported by additional SV caller

However, separating variants in this way introduces several issues. First, variant classification depends on representation: a single large deletion may be represented either as one SV or as multiple small variants. Second, representations may differ between callsets: a large deletion in one callset may appear as multiple small variants in another. These inconsistencies complicate downstream comparison and benchmarking.

To address this, we construct per-sample merged callsets that include both small variants and SVs, enabling consistent comparison at the haplotype sequence level using Aardvark.

All steps (including merged callset construction and graph variant benchmarking) are restricted to confident genomic regions defined by the dipcall BED file for each sample. The corresponding files are listed in the Data Source section.

### Preprocessing and haplotype consistency

We first convert VCFs from all 14 callers into sequence-resolved form. Because the merged callset integrates VCFs from assembly-based callers (dipcall and PAV) and a HiFi-based caller (longcallD), we unphase longcallD genotypes prior to merging. Although longcallD reports phased genotypes, the phasing is valid only within local phase blocks and may not be consistent across genomic regions. Since Aardvark clusters variants based on genomic distance rather than phase-block boundaries, retaining these phased genotypes could introduce incorrect haplotype structure across clusters.

This unphasing step does not affect `aardvark merge`, which evaluates haplotype sequences within each cluster. However, it is important when the resulting backbone callset is later treated as the truth callset for `aardvark compare`, where phased genotypes could otherwise impose incorrect haplotype constraints on query callsets.

### Backbone callset construction

We used 14 variant callers in total. Among them, three callers generate joint callsets containing both small variants and SVs: two assembly-based callers (dipcall and PAV) and one HiFi-based caller (longcallD). These three callers form the backbone of the merged callsets, while all 14 callers are incorporated later as supporting evidence.

In addition to producing joint callsets, these backbone callers provide base-level accurate variant representations. In contrast, many SV callers (e.g. SVIM-asm) merge similar heterozygous SVs into homozygous events, reducing base-level accuracy and complicating sequence-based comparisons.

We merged the sequence-resolved VCFs from the three backbone callers using `aardvark merge`. Aardvark groups nearby variants into clusters based on genomic distance, reconstructs haplotype sequences within each cluster, and compares haplotypes across callers.

We apply a priority order (from high to low): dipcall, PAV, and longcallD. This priority is used to select the reported variant representation when multiple callers produce identical haplotype sequences.

For each cluster:

- If at least two callers produce identical haplotype sequences, the representation from the highest-priority caller among them is selected
- If all three callers differ, variants from dipcall are selected
- If only one caller contributes variants, those variants are retained

This strategy produces a backbone callset that is close to a union of calls while preferentially selecting representations supported by consistent haplotype sequences.

### Caller support evaluation

To incorporate evidence from all callers, we treat each of the 14 sequence-resolved callsets as a query callset and compare it to the backbone callset using `aardvark compare`, which computes basepair-level recall and precision for each cluster.

A caller is considered to support a cluster if both recall and precision meet a specified threshold. Thresholds are chosen based on caller type:

- Joint callers (dipcall, PAV, longcallD): 0.9
- Small variant caller (DeepVariant): 1.0
- SV callers (cuteSV-asm, SVIM-asm, cuteSV, DeBreak, DELLY, pbsv, sawfish, Sniffles2, SVDSS, SVIM): 0.5

For joint callers, we require basepair-level recall and precision to be at least 0.9. These callers generate joint callsets containing both small variants and SVs, and generally provide base-level accurate variant representations. However, because `aardvark compare` is used here to assess caller support rather than exact agreement, we allow limited mismatches instead of requiring perfect haplotype sequence matches.

For DeepVariant, we require exact matches (recall = precision = 1.0) because it is a small-variant-only caller and these variants are typically represented precisely. To reduce the impact of missing SVs in the query callset, we adjust the clustering distance from 1,000 bp to 50 bp when running `aardvark compare`, minimizing spurious mismatches between small variants and nearby SVs.

For SV callers, we use a threshold of 0.5. Before comparison, we combine SVs (≥ 50 bp) with small variants (< 50 bp) from DeepVariant to reduce mismatches caused by missing small variants. In addition, many SV callers merge similar heterozygous SVs into homozygous events, which further reduces base-level concordance. A threshold of 0.5 therefore provides a more appropriate criterion for support.

### Final callset construction

Finally, we consolidate caller support information across all outputs from `aardvark merge` and `aardvark compare`. For each variant, the `INFO/SOURCES` and `INFO/CALLERS` values are merged into a single `INFO/CALLERS` field listing supporting callers, and an `INFO/NCALLERS` field indicating the number of supporting callers. We also assign a unique variant ID and annotate `INFO/SVTYPE` and `INFO/SVLEN` for variants ≥ 50 bp.

## Graph Variant Benchmarking

For each sample, we derive a joint ground truth callset from the merged callset by selecting variants supported by at least two callers (i.e., `INFO/NCALLERS ≥ 2`). This filtering step retains variants with consistent support across methods and serves as the truth callset for benchmarking.

Prior to benchmarking, graph variants are preprocessed using a modified version of the [prepare-vcf-MC workflow](https://github.com/eblerjana/genotyping-pipelines/tree/main/prepare-vcf-MC) to decompose bubbles based on allele traversals, producing a representation suitable for comparison.

Graph variants are then compared against the joint ground truth callset using `aardvark compare`. Comparisons are restricted to confident genomic regions defined by the dipcall BED file for each sample.

In addition, GIAB stratification v3.6 BED files are applied to evaluate performance across different genomic contexts, enabling more detailed assessment of variant calling performance in regions such as segmental duplication and tandem repeat regions.
