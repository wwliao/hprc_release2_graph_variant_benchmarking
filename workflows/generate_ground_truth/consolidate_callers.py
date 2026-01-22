#!/usr/bin/env python3
import argparse
import pysam


def read_contigs_from_bed(bedfile):
    contigs = set()
    with open(bedfile) as infile:
        for line in infile:
            chrom = line.strip().split("\t")[0]
            contigs.add(chrom)
    return contigs


def clean_header_ordered(orig_header, allowed_contigs):
    DROP_INFO = {"SOURCES", "MR"}
    KEEP_FILTER = {"PASS"}
    KEEP_FORMAT = {"GT"}

    TYPE_ORDER = ("FILTER", "INFO", "FORMAT", "CONTIG")

    new_header = pysam.VariantHeader()

    # Add fileformat first
    for rec in orig_header.records:
        if rec.type == "GENERIC" and rec.key == "fileformat":
            new_header.add_record(rec)
            break

    # Add records in order
    for wanted_type in TYPE_ORDER:
        for rec in orig_header.records:
            if rec.type != wanted_type:
                continue

            if wanted_type == "FILTER":
                if rec.get("ID") not in KEEP_FILTER:
                    continue

            elif wanted_type == "INFO":
                if rec.get("ID") in DROP_INFO:
                    continue

            elif wanted_type == "FORMAT":
                if rec.get("ID") not in KEEP_FORMAT:
                    continue

            elif wanted_type == "CONTIG":
                if rec.get("ID") not in allowed_contigs:
                    continue

            new_header.add_record(rec)

        # Insert required INFO definitions inside INFO block
        if wanted_type == "INFO":
            if "SVTYPE" not in new_header.info:
                new_header.add_line(
                    '##INFO=<ID=SVTYPE,Number=1,Type=String,'
                    'Description="Type of structural variant">'
                )
            if "SVLEN" not in new_header.info:
                new_header.add_line(
                    '##INFO=<ID=SVLEN,Number=A,Type=Integer,'
                    'Description="Length of structural variant">'
                )
            if "CALLERS" not in new_header.info:
                new_header.add_line(
                    '##INFO=<ID=CALLERS,Number=.,Type=String,'
                    'Description="List of variant callers supporting this record">'
                )
            if "NCALLERS" not in new_header.info:
                new_header.add_line(
                    '##INFO=<ID=NCALLERS,Number=1,Type=Integer,'
                    'Description="Number of variant callers supporting this record">'
                )

    # Add samples
    for sample in orig_header.samples:
        new_header.add_sample(sample)

    return new_header


def parse_sources(val):
    if val is None:
        return []
    if isinstance(val, (tuple, list)):
        out = []
        for x in val:
            out.extend(str(x).split(","))
        return out
    return str(val).split(",")


def record_key_biallelic(rec):
    alt = rec.alts[0] if rec.alts else None
    return (rec.contig, rec.pos, rec.ref, alt)


def infer_svtype_svlen(ref, alt):
    """
    Types:
      SNP: len(ref)==len(alt)==1
      MNP: len(ref)==len(alt)>1 and <50
      CPX: len(ref)==len(alt)>=50     (treated as SV; SVLEN=len(ref))
      INS/DEL: len(ref)!=len(alt)     (SVLEN=len(alt)-len(ref))
    """
    lr, la = len(ref), len(alt)

    if lr == la:
        if lr == 1:
            return "SNP", None
        if lr < 50:
            return "MNP", None
        return "CPX", lr

    d = la - lr
    return ("INS" if d > 0 else "DEL"), d


def make_record_id(chrom, pos_1based, ref, alt):
    svtype, svlen = infer_svtype_svlen(ref, alt)

    if svtype in {"SNP", "MNP"}:
        rid = f"{chrom}-{pos_1based}-{svtype}-{ref}{alt}"
        return rid, svtype, svlen

    if svtype == "CPX":
        rid = f"{chrom}-{pos_1based}-CPX-{len(ref)}"
        return rid, svtype, svlen

    # INS / DEL
    rid = f"{chrom}-{pos_1based}-{svtype}-{abs(svlen)}"
    return rid, svtype, svlen


def uniquify_id(rid, seen):
    """
    Ensure ID uniqueness by appending _2, _3, ... if needed.
    """
    n = seen.get(rid, 0) + 1
    seen[rid] = n
    return rid if n == 1 else f"{rid}_{n}"


parser = argparse.ArgumentParser(
    description=(
        "Consolidate caller support information from multiple record-aligned, "
        "biallelic VCFs produced by `aardvark merge` and `aardvark compare`. "
        "For each variant, INFO/SOURCES values from all input VCFs are "
        "combined into INFO/CALLERS with INFO/NCALLERS indicating the number "
        "of supporting callers. The script also assigns a unique ID to each "
        "variant and annotates INFO/SVTYPE and INFO/SVLEN for variants ≥50 bp."
    )
)
parser.add_argument("-o", "--output", required=True, help="Output VCF. If the filename ends with .vcf.gz, a bgzip-compressed VCF and a tabix index (.tbi) will be created.")
parser.add_argument("--contigs-bed", required=True, help="BED file whose first column defines the contigs to keep in the output header.")
parser.add_argument("vcfs", nargs="+", help="Input VCFs (record-aligned, biallelic, same variants in same order).")
args = parser.parse_args()

ins = [pysam.VariantFile(p) for p in args.vcfs]
allowed_contigs = read_contigs_from_bed(args.contigs_bed)
header = clean_header_ordered(ins[0].header, allowed_contigs)
out = pysam.VariantFile(args.output, "w", header=header)

seen_id_counts = {}

iters = [f.fetch() for f in ins]
for recs in zip(*iters):
    r0 = recs[0]

    # Safety check: records aligned across inputs
    key0 = record_key_biallelic(r0)
    for i, r in enumerate(recs[1:], start=2):
        if record_key_biallelic(r) != key0:
            raise RuntimeError(
                f"Record mismatch: input#1 {key0} vs input#{i} {record_key_biallelic(r)}"
            )

    # Consolidate INFO/SOURCES across inputs
    callers = set()
    for r in recs:
        for x in parse_sources(r.info.get("SOURCES")):
            x = x.strip()
            if x:
                callers.add(x)

    # Build ID + SV annotations
    rid, svtype, svlen = make_record_id(r0.contig, r0.pos, r0.ref, r0.alts[0])
    rid = uniquify_id(rid, seen_id_counts)

    # Create a NEW record under the OUTPUT header (avoids "unknown INFO" errors)
    new = out.new_record(
        contig=r0.contig,
        start=r0.start,
        stop=r0.stop,
        id=rid,
        alleles=r0.alleles,
        qual=None,
        filter="PASS"
    )

    # Add INFO/SVTYPE and INFO/SVLEN
    if svlen is not None:
        if (svtype in {"INS", "DEL"} and abs(svlen) >= 50) or (svtype == "CPX" and svlen >= 50):
            new.info["SVTYPE"] = svtype
            new.info["SVLEN"] = int(svlen)

    # Add INFO/CALLERS and INFO/NCALLERS
    new.info["CALLERS"] = tuple(sorted(callers))
    new.info["NCALLERS"] = len(callers)

    # Keep only FORMAT/GT and change 1/1 -> 1|1 only
    for sample in header.samples:
        gt = r0.samples[sample].get("GT")
        phased = r0.samples[sample].phased

        new.samples[sample]["GT"] = gt
        if gt == (1, 1) and not phased:
            new.samples[sample].phased = True
        else:
            new.samples[sample].phased = phased

    out.write(new)

out.close()
for f in ins:
    f.close()

# Create tabix index if output is bgzip-compressed
if args.output.endswith(".gz"):
    pysam.tabix_index(args.output, preset="vcf", force=True)
