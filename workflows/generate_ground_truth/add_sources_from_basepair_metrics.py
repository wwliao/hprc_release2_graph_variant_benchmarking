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


parser = argparse.ArgumentParser(
    description=(
        "Add the INFO/SOURCES field in an input VCF using BASEPAIR recall and "
        "precision from a region_summary.tsv.gz file. For each variant, the "
        "region ID is read from FORMAT/RI and the corresponding BASEPAIR "
        "metrics are evaluated. If both recall and precision ≥ the given "
        "threshold, SOURCES is set to include caller."    
    )
)
parser.add_argument(
    "-i", "--input", required=True,
    help="Input VCF."
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
    "-t", "--threshold", required=True, type=float,
    help=(
        "Minimum value required for both BASEPAIR recall and precision "
        "(>= threshold)."
    )
)
args = parser.parse_args()

# Load BASEPAIR recall/precision by region_id
metrics_by_region = load_basepair_metrics(args.region_summary)

vcf = pysam.VariantFile(args.input)
header = vcf.header

# Ensure INFO/SOURCES exists in the output header
if "SOURCES" not in header.info:
    header.add_line(
        '##INFO=<ID=SOURCES,Number=.,Type=String,'
        'Description="List of tools or technologies that called the same record">'
    )

# Validate the sample exists before processing
if args.sample not in vcf.header.samples:
    sys.exit(f"Error: sample {args.sample!r} not found in VCF")

out = pysam.VariantFile(args.output, "w", header=header)

# Process input VCF records
for rec in vcf:
    # Region ID from FORMAT/RI for the given sample
    sample_data = rec.samples[args.sample]
    region_id = int(sample_data["RI"])

    # Create a new record with core fields carried over
    new = out.new_record(
        contig=rec.contig,
        start=rec.start,
        stop=rec.stop,
        alleles=rec.alleles,
    )

    # Preserve genotype and phasing for the target sample
    new.samples[args.sample]["GT"] = rec.samples[args.sample].get("GT")
    new.samples[args.sample].phased = rec.samples[args.sample].phased

    # Lookup metrics for this region and decide SOURCES
    metrics = metrics_by_region.get(region_id)

    pass_metrics = (
        metrics["recall"] is not None
        and metrics["precision"] is not None
        and metrics["recall"] >= args.threshold
        and metrics["precision"] >= args.threshold
    )

    if pass_metrics:
        new.info["SOURCES"] = f"{args.caller}"

    out.write(new)

out.close()
vcf.close()

# Create tabix index if output is bgzip-compressed
if args.output.endswith(".gz"):
    pysam.tabix_index(args.output, preset="vcf", force=True)
