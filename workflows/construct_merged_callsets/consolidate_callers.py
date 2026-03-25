#!/usr/bin/env python3
import argparse
import pysam


def read_contigs_from_bed(bedfile):
    contigs = set()
    with open(bedfile) as infile:
        for line in infile:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            chrom = line.split("\t")[0]
            contigs.add(chrom)
    return contigs


def clean_header_keep_first(orig_header, allowed_contigs):
    """
    Keep the first VCF's header (so INFO/SOURCES and INFO/MR are preserved),
    optionally restrict CONTIG lines by allowed_contigs, and ensure 
    INFO/CALLERS, INFO/NCALLERS, INFO/SVTYPE, and INFO/SVLEN exist in the 
    output header.
    """
    new_header = pysam.VariantHeader()

    KEEP_FILTER = {"PASS"}
    KEEP_INFO = {"SOURCES", "MR"}
    KEEP_FORMAT = {"GT"}

    TYPE_ORDER = ("FILTER", "INFO", "FORMAT", "CONTIG")

    new_header = pysam.VariantHeader()

    # Add records in order
    for wanted_type in TYPE_ORDER:
        for rec in orig_header.records:
            if rec.type != wanted_type:
                continue

            if wanted_type == "FILTER":
                if rec.get("ID") not in KEEP_FILTER:
                    continue

            elif wanted_type == "INFO":
                if rec.get("ID") not in KEEP_INFO:
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
            new_header.add_line(
                '##INFO=<ID=SVTYPE,Number=1,Type=String,'
                'Description="Type of structural variant">'
            )
            new_header.add_line(
                '##INFO=<ID=SVLEN,Number=A,Type=Integer,'
                'Description="Length of structural variant">'
            )
            new_header.add_line(
                '##INFO=<ID=CALLERS,Number=.,Type=String,'
                'Description="List of variant callers supporting this record">'
            )
            new_header.add_line(
                '##INFO=<ID=NCALLERS,Number=1,Type=Integer,'
                'Description="Number of variant callers supporting this record">'
            )

    # Copy samples (order preserved)
    for sample in orig_header.samples:
        new_header.add_sample(sample)

    return new_header


def parse_info_list(val):
    if val is None:
        return []
    return list(val)


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
    n = seen.get(rid, 0) + 1
    seen[rid] = n
    return rid if n == 1 else f"{rid}_{n}"


parser = argparse.ArgumentParser(
    description=(
        "INFO/SOURCES and INFO/MR from the first input VCF are preserved. "
        "INFO/CALLERS in the output is constructed by combining INFO/SOURCES "
        "from the first VCF and INFO/CALLERS from the remaining VCFs; "
        "INFO/NCALLERS records the number of supporting callers. "
        "All inputs must be record-aligned, biallelic, and in the same order."
    )
)
parser.add_argument(
    "-o",
    "--output",
    required=True,
    help="Output VCF. If the filename ends with .vcf.gz, a bgzip-compressed VCF and a tabix index (.tbi) will be created.",
)
parser.add_argument(
    "--contigs-bed",
    required=True,
    help="BED file whose first column defines contigs to keep in the output header.",
)
parser.add_argument(
    "vcfs",
    nargs="+",
    help="Input VCFs. First VCF is the backbone; remaining VCFs provide INFO/CALLERS support.",
)
args = parser.parse_args()

if len(args.vcfs) < 2:
    raise SystemExit("ERROR: Provide at least 2 input VCFs (base VCF + >=1 caller VCF).")

ins = [pysam.VariantFile(p) for p in args.vcfs]

allowed_contigs = read_contigs_from_bed(args.contigs_bed) if args.contigs_bed else None
header = clean_header_keep_first(ins[0].header, allowed_contigs)

out = pysam.VariantFile(args.output, "w", header=header)

seen_id_counts = {}

iters = [f.fetch() for f in ins]
for recs in zip(*iters):
    r0 = recs[0]

    # Safety check: record-aligned across inputs
    key0 = record_key_biallelic(r0)
    for i, r in enumerate(recs[1:], start=2):
        if record_key_biallelic(r) != key0:
            raise RuntimeError(
                f"Record mismatch: input#1 {key0} vs input#{i} {record_key_biallelic(r)}"
            )

    # Consolidate supporting tools: SOURCES from first VCF, CALLERS from others
    callers = set()

    for r, tag in [(recs[0], "SOURCES")] + [(r, "CALLERS") for r in recs[1:]]:
        for x in parse_info_list(r.info.get(tag)):
            x = x.strip()
            if x:
                callers.add(x)

    # Build ID + (optional) SV annotations
    rid, svtype, svlen = make_record_id(r0.contig, r0.pos, r0.ref, r0.alts[0])
    rid = uniquify_id(rid, seen_id_counts)

    # Create a NEW record under the OUTPUT header
    new = out.new_record(
        contig=r0.contig,
        start=r0.start,
        stop=r0.stop,
        id=rid,
        alleles=r0.alleles,
        qual=r0.qual,
        filter=r0.filter.keys() if len(r0.filter.keys()) > 0 else "PASS",
    )

    # Copy INFO from the first VCF (preserve SOURCES, MR, etc.)
    for k in r0.info.keys():
        new.info[k] = r0.info[k]

    # Add/overwrite INFO/SVTYPE and INFO/SVLEN (SV >= 50 bp)
    if svlen is not None:
        if (svtype in {"INS", "DEL"} and abs(svlen) >= 50) or (svtype == "CPX" and svlen >= 50):
            new.info["SVTYPE"] = svtype
            new.info["SVLEN"] = int(svlen)

    # Apply DeepVariant filtering rule
    if "DeepVariant" in callers:
        allowed_with_dv = {"DeepVariant", "dipcall", "PAV", "longcallD"}
        callers = callers & allowed_with_dv

    # Add/overwrite INFO/CALLERS and INFO/NCALLERS
    new.info["CALLERS"] = tuple(sorted(callers))
    new.info["NCALLERS"] = int(len(callers))

    # Keep only FORMAT/GT and change 1/1 -> 1|1 only (preserve original phasing otherwise)
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
