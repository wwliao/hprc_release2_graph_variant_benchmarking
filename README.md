# Graph Variant Benchmarking for the HPRC Release 2

This repository contains Nextflow workflows for graph variant benchmarking across 230 of 231 individuals from the Human Pangenome Reference Consortium (HPRC) Release 2.

One individual (HG00272) was excluded from the pangenome graph because of a likely large-scale misassembly on chromosome X.

The workflows in this repository support:

- Generating per-sample ground truth variant sets using multiple variant callers
- Benchmarking variants derived from the pangenome graph (referred to as graph variants below)

## Ground Truth Generation

To evaluate the quality of graph variants, a ground truth variant set is required for each sample. However, no validated ground truth exists for these 230 samples.

To address this limitation, we adopted a multi-caller consensus approach:

- Variants were called using 14 variant callers per sample (4 assembly-based and 10 HiFi-based callers)
- Variants from all callers were merged into a single VCF per sample
- Variants supported by at least two callers were treated as truth calls

