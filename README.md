# Graph Variant Benchmarking for the HPRC Release 2

This repository contains Nextflow workflows for graph variant benchmarking across 230 of 231 individuals from the Human Pangenome Reference Consortium (HPRC) Release 2, with HG00272 excluded due to a likely large-scale misassembly on chromosome X.

The workflows in this repository support:

- Per-sample ground truth generation using multiple variant callers
- Benchmarking variants derived from the pangenome graph, referred to as graph variants below

## Why benchmark with Aardvark?

Benchmarking variants is challenging because variant representation differs substantially across methods and tools. This problem is even more severe for graph variants, whose representations can be very different from those produced by traditional linear-reference-based callers.

To minimize bias introduced by representation differences, we use Aardvark, a recently developed benchmarking tool with the following advantages:

- Haplotype-based comparison: Aardvark reconstructs haplotype sequences by jointly considering small variants and structural variants (SVs), instead of comparing variants record by record.
- Basepair-level metrics: In addition to variant-level metrics, Aardvark provides basepair-level recall and precision, which are particularly informative for benchmarking SVs.

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

Our goal is to generate per-sample joint ground truth VCFs that contain both small variants and SVs.

We used 14 variant callers in total. However, only three callers produce joint call sets that include both small variants and SVs: two assembly-based callers (dipcall and PAV) and one HiFi-based caller (longcallD). These three callers are therefore used to construct the backbone of the joint ground truth VCFs, while the remaining 11 callers are incorporated later as additional supporting evidence. In addition to generating joint call sets, these three callers provide base-level accurate variant representations. In contrast, many SV callers (e.g. SVIM-asm) merge similar heterozygous SVs into homozygous events, which can reduce base-level accuracy and complicate sequence-based comparisons.

We first generate sequence-resolved VCFs from the three backbone callers and merge them using Aardvark. Aardvark groups nearby variants into smaller clusters (sub-regions) based on genomic distance, reconstructs haplotype sequences for each cluster, and compares haplotype sequences across callers. We run `aardvark merge` with the following priority order (from high to low): dipcall, PAV, and longcallD. This priority is used only to select the reported variant representation when multiple callers agree at the haplotype level.

For each cluster, if two out of three callers produce identical haplotype sequences, the variant representation from the higher-priority caller among the two is reported (e.g. if PAV and longcallD agree, PAV is used). If all three callers produce different haplotype sequences, variants from dipcall are reported. If only one caller has variants in a cluster, variants from that caller are reported. Using this strategy, the resulting backbone VCF is close to a union callset, while preferentially selecting representations supported by the majority haplotype sequence.

