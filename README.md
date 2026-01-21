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

## Why create custom ground truth VCFs?

Existing ground truth datasets typically separate small variants and SVs, for example:

- Small variant truth sets derived from DeepVariant
- SV truth sets derived from PAV, supported by at least one additional SV caller

However, separating ground truths in this way introduces several issues:

- Variant classification depends on representation: Whether a variant is labeled as a small variant or an SV often depends on how it is represented (e.g., a single large deletion vs. many small variants).
- Inconsistent representations between truth and query: A large deletion in the truth set may appear as multiple small variants in a query callset, or vice versa, complicating fair benchmarking.

To address these problems, we perform joint benchmarking without separating variants by type or size. As a result, we generate custom per-sample ground truth VCFs that contain both small variants and SVs, enabling benchmarking that is robust to differences in variant representation using Aardvark.

