# Graph Variant Benchmarking for the HPRC Release 2

This repository contains Nextflow workflows for graph variant benchmarking across 230 of 231 individuals from the Human Pangenome Reference Consortium (HPRC) Release 2. HG00272 was excluded due to a likely large-scale misassembly on chromosome X and was therefore not included in the HPRC Release 2 pangenome graph.

The workflows in this repository support:

- Per-sample ground truth generation using multiple variant callers
- Benchmarking variants derived from the pangenome graph, referred to as *graph variants* below

## Why benchmark with Aardvark?

Benchmarking variants is challenging because variant representation differs substantially across methods and tools. This problem is even more severe for graph variants, whose representations can be very different from those produced by traditional linear-reference-based callers.

To minimize bias introduced by representation differences, we use Aardvark, a recently developed benchmarking tool with the following advantages:

- Haplotype-based comparison: Aardvark reconstructs haplotype sequences by *jointly* considering small variants and structural variants (SVs), instead of comparing variants record by record.
- Basepair-level metrics: In addition to variant-level metrics, Aardvark provides *basepair-level recall and precision*, which are particularly informative for benchmarking SVs.

Together, these features allow for a more robust and representation-agnostic comparison.

## Why create joint ground truth VCFs?

Existing ground truth datasets typically separate small variants and SVs, for example:

- Small variant truth sets derived from DeepVariant
- SV truth sets derived from PAV, supported by at least one additional SV caller

However, separating ground truths in this way introduces several issues:

- Variant classification depends on representation: Whether a variant is labeled as a small variant or an SV often depends on how it is represented (e.g., a single large deletion vs. many small variants).
- Inconsistent representations between truth and query: A large deletion in the truth set may appear as multiple small variants in a query callset, or vice versa, complicating fair benchmarking.

As a result, we generate per-sample joint ground truth VCFs that contain both small variants and SVs, enabling consistent benchmarking despite differences in variant representation using Aardvark.

## How are the joint ground truth VCFs created?

We used 14 variant callers in total. However, only three callers produce joint callsets that include both small variants and SVs: two assembly-based callers (dipcall and PAV) and one HiFi-based caller (longcallD). These three callers are therefore used to construct the *backbone* of the joint ground truth VCFs, while the remaining 11 callers are incorporated later as additional supporting evidence. In addition to generating joint callsets, these three callers provide base-level accurate variant representations. In contrast, many SV callers (e.g. SVIM-asm) merge similar heterozygous SVs into homozygous events, which can reduce base-level accuracy and complicate sequence-based comparisons.

We first generate sequence-resolved VCFs from the three backbone callers and merge them using Aardvark. Aardvark groups nearby variants into smaller clusters (sub-regions) based on genomic distance, reconstructs haplotype sequences for each cluster, and compares haplotype sequences across callers. We run `aardvark merge` with the following *priority order* (from high to low): dipcall, PAV, and longcallD. This priority is used only to select the reported variant representation when multiple callers agree at the haplotype level.

For each cluster, if two out of three callers produce identical haplotype sequences, the variant representation from the higher-priority caller among the two is reported (e.g. if PAV and longcallD agree, PAV is used). If all three callers produce different haplotype sequences, variants from dipcall are reported. If only one caller has variants in a cluster, variants from that caller are reported. Using this strategy, the resulting backbone VCF is close to a union callset, while preferentially selecting representations supported by the majority haplotype sequence.

We treat the VCF from each of the remaining 11 callers as a query callset and compare it to the backbone VCF using `aardvark compare`, which computes basepair-level recall and precision for each cluster. For a given cluster, we require both recall and precision to be greater than or equal to a specified *threshold* in order to claim that the caller supports the corresponding backbone sub-region.

Unlike `aardvark merge`, we do not require an exact haplotype sequence match. This is because most of these callers (with the exception of DeepVariant, which calls only small variants) report SVs only. When such callsets are compared against the backbone, which contains both small variants and SVs, missing small variants can lead to haplotype sequence mismatches. In addition, as noted above, many SV callers merge similar heterozygous SVs into homozygous events, which further reduces base-level concordance. Using basepair-level recall and precision thresholds therefore provides a more appropriate criterion for determining caller support.

For the small-variant-only caller, DeepVariant, we require basepair-level recall and precision to be 1.0 (i.e. an exact match), given the small size and precise representation of these variants. To reduce the impact of missing SVs in the query callset, we adjust the genomic distance used for clustering nearby variants, reducing it from 1,000 bp to 50 bp when running `aardvark compare`. This tighter clustering minimizes the chance that small variants are grouped with nearby SVs that are absent from the query callset, which would otherwise lead to spurious haplotype mismatches.

Finally, we consolidate caller support information from the multiple VCFs produced by `aardvark merge` and `aardvark compare`. For each variant, the `INFO/SOURCES` values from all input VCFs are combined into a single `INFO/CALLERS` field, which lists the callers supporting the variant, and an `INFO/NCALLERS` field, which indicates the number of supporting callers. We also assign a unique ID to each variant and annotate `INFO/SVTYPE` and `INFO/SVLEN` for variants ≥ 50 bp. For benchmarking, variants supported by at least two callers are selected as the truth callset.

To tune the *threshold* for balancing basepair-level recall and precision, we used HG002 callsets generated by the same 14 callers under comparable conditions. Specifically, we used assemblies constructed with equivalent pipelines and similar read depth, as well as HiFi data with read depth comparable to the HPRC Release 2 samples, to generate a joint callset for HG002. We then applied the same caller-support criterion (≥ 2 callers) to the HG002 joint callset and compared it against the GIAB T2T-HG002 Q100 v1.1 ground truth. **Thresholds were selected based on the resulting *precision–recall curves*: we chose 0.96 when using GRCh38 as the reference and 0.92 when using T2T-CHM13 v2 as the reference.**
