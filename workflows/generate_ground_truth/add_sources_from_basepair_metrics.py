#!/usr/bin/env python3
import argparse
import csv
import gzip
import sys

import pysam


def load_basepair_metrics(region_summary_path):
    """
    Load BASEPAIR recall and precision per region from a region summary file.
    Regions with no truth or query basepairs are treated as uninformative and
    assigned None values.
    """
    metrics_by_region = {}

    with gzip.open(region_summary_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("comparison") != "BASEPAIR":
                continue

            rid = int(row["region_id"])
            truth_tp = int(row["truth_tp"])
            truth_total = int(row["truth_total"])
            query_tp = int(row["query_tp"])
            query_total = int(row["query_total"])

            recall = None if truth_total == 0 else (truth_tp / truth_total)
            precision = None if query_total == 0 else (query_tp / query_total)

            metrics_by_region[rid] = {"recall": recall, "precision": precision}

    return metrics_by_region


def rec_key_biallelic(rec):
    """
    Key a biallelic record by (CHROM, POS, REF, ALT).
    """
    if rec.alts is None or len(rec.alts) != 1:
        raise ValueError(f"Expected biallelic record but got alts={rec.alts}")
    return (rec.contig, rec.pos, rec.ref, rec.alts[0])


def build_input_index(input_vcf_path):
    """
    Build an in-memory dict: (CHROM, POS, REF, ALT) -> pysam.VariantRecord
    """
    vcf = pysam.VariantFile(input_vcf_path)
    idx = {}
    for rec in vcf:
        idx[rec_key_biallelic(rec)] = rec
    vcf.close()
    return idx

parser = argparse.ArgumentParser(
    description=(
        "Drive output by a truth VCF so all truth records are present. If a "
        "truth record also exists in the input VCF, read FORMAT/RI from the "
        "input record, look up BASEPAIR metrics for that region ID in "
        "region_summary.tsv.gz, and set INFO/SOURCES to caller when both "
        "recall and precision ≥ the given threshold. If the truth record is "
        "missing in the input VCF, write basic record info only and do not "
        "set INFO/SOURCES."
    )
)
parser.add_argument(
    "-i", "--input", required=True,
    help="Input VCF (may miss records)."
)
parser.add_argument(
    "-t", "--truth", required=True,
    help="Truth VCF (defines the full record set)."
)
parser.add_argument(
    "-o", "--output", required=True,
    help="Output VCF. If the filename ends with .vcf.gz, a bgzip-compressed VCF and a tabix index (.tbi) will be created."
)
parser.add_argument(
    "--region-summary", required=True,
    help=(
        "region_summary.tsv.gz containing BASEPAIR recall and precision "
        "per region."
    )
)
parser.add_argument(
    "--sample", required=True,
    help="Sample name used to read FORMAT fields."
)
parser.add_argument(
    "--caller", required=True,
    help="Caller name to include in INFO/SOURCES when metrics pass the threshold."
)
parser.add_argument(
    "--threshold", required=True, type=float,
    help=(
        "Minimum value required for both BASEPAIR recall and precision "
        "(>= threshold)."
    )
)
args = parser.parse_args()

# Load BASEPAIR recall/precision by region_id
metrics_by_region = load_basepair_metrics(args.region_summary)

# Validate sample exists in input VCF (RI comes from input)
input_vcf = pysam.VariantFile(args.input)
header = input_vcf.header.copy()
if args.sample not in input_vcf.header.samples:
    sys.exit(f"Error: sample {args.sample!r} not found in input VCF")
input_vcf.close()

# Ensure INFO/SOURCES exists in the output header
if "SOURCES" not in header.info:
    header.add_line(
        '##INFO=<ID=SOURCES,Number=.,Type=String,'
        'Description="List of tools or technologies that called the same record">'
    )

# Validate sample exists in truth VCF
truth_vcf = pysam.VariantFile(args.truth)
if args.sample not in truth_vcf.header.samples:
    sys.exit(f"Error: sample {args.sample!r} not found in truth VCF")

# Index input records for fast lookup
input_index = build_input_index(args.input)

out = pysam.VariantFile(args.output, "w", header=header)

for truth_rec in truth_vcf:
    k = rec_key_biallelic(truth_rec)
    input_rec = input_index.get(k)

    # Create a new record with core fields carried over
    new = out.new_record(
        contig=truth_rec.contig,
        start=truth_rec.start,
        stop=truth_rec.stop,
        alleles=truth_rec.alleles,
    )


    if input_rec is not None:
        # Only set INFO/SOURCES when input record exists and metrics pass
        sample_data = input_rec.samples[args.sample]

        # Region ID from FORMAT/RI for the given sample
        region_id = None
        if "RI" in sample_data:
            try:
                region_id = int(sample_data["RI"])
            except Exception:
                region_id = None

        metrics = metrics_by_region.get(region_id) if region_id is not None else None

        pass_metrics = False
        if metrics is not None:
            pass_metrics = (
                metrics["recall"] is not None
                and metrics["precision"] is not None
                and metrics["recall"] >= args.threshold
                and metrics["precision"] >= args.threshold
            )

        if pass_metrics:
            new.info["SOURCES"] = args.caller


    # Preserve genotype and phasing for the target sample
    new.samples[args.sample]["GT"] = truth_rec.samples[args.sample].get("GT")
    new.samples[args.sample].phased = truth_rec.samples[args.sample].phased

    # If missing in input, do nothing: no INFO/SOURCES set
    out.write(new)

out.close()
truth_vcf.close()

# Create tabix index if output is bgzip-compressed
if args.output.endswith(".gz"):
    pysam.tabix_index(args.output, preset="vcf", force=True)
