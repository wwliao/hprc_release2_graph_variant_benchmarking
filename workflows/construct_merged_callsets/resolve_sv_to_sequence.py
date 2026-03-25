#!/usr/bin/env python3
import argparse
from cyvcf2 import VCF, Writer
import pysam


def left_anchor(ref, chrom, pos):
    """
    Return (new_pos, anchor_base).

    The anchor base is the single reference base immediately to the left of POS.
    POS is shifted left by one base (1-based coordinates).
    """
    new_pos = pos - 1
    anchor = ref.fetch(chrom, new_pos - 1, new_pos)
    return new_pos, anchor


def resolve_dup(rec, ref):
    """
    Resolve DUP to a sequence-based representation.

    - If REF == "N", left-anchor the variant and set ALT = anchor + duplicated sequence
    - Otherwise, set ALT to the duplicated sequence
    """
    chrom = rec.CHROM
    pos = rec.POS
    end = rec.INFO.get("END")

    dup_seq = ref.fetch(chrom, pos - 1, end)

    if rec.REF == "N":
        new_pos, anchor = left_anchor(ref, chrom, pos)
        rec.set_pos(new_pos - 1)  # cyvcf2 expects 0-based positions
        rec.REF = anchor
        rec.ALT = [anchor + dup_seq]
    else:
        rec.ALT = [dup_seq]


def resolve_dup_tandem(rec, ref):
    """
    Resolve DUP:TANDEM using copy number (CN).

    - If REF == "N":
      REF is set to one copy of the duplicated unit
      ALT is set to CN copies of the duplicated unit
    - Otherwise:
      ALT is set to REF followed by (CN - 1) copies of the duplicated unit
    """
    chrom = rec.CHROM
    pos = rec.POS
    end = rec.INFO.get("END")

    cn = int(rec.format("CN")[0][0])

    if rec.REF == "N":
        unit = ref.fetch(chrom, pos - 1, end)
        rec.REF = unit
        rec.ALT = [unit * cn]
    else:
        # Append repeated sequence after the existing REF
        unit = ref.fetch(chrom, pos, end)
        rec.ALT = [rec.REF + unit * (cn - 1)]


def resolve_del(rec, ref):
    """
    Resolve sequence-based deletions.

    - Only modifies records where ALT == "N"
    - Symbolic deletions (<DEL>) are skipped
    """
    if rec.ALT[0] == "<DEL>":
        # Skip symbolic deletions (e.g., DELLY, sawfish, Sniffles2)
        return

    if rec.ALT[0] == "N":
        chrom = rec.CHROM
        pos = rec.POS
        new_pos, anchor = left_anchor(ref, chrom, pos)
        rec.set_pos(new_pos - 1)  # cyvcf2 expects 0-based positions
        # REF contains the deleted sequence
        rec.REF = anchor + rec.REF
        rec.ALT = [anchor]


def resolve_ins(rec, ref):
    """
    Resolve sequence-based insertions.

    - Only modifies records where REF == "N"
    - Symbolic insertions (<INS>) are skipped
    """
    if rec.ALT[0] == "<INS>":
        # Skip symbolic insertions (e.g., DeBreak, DELLY, Sniffles2, SVIM)
        return

    if rec.REF == "N":
        chrom = rec.CHROM
        pos = rec.POS
        new_pos, anchor = left_anchor(ref, chrom, pos)
        rec.set_pos(new_pos - 1)  # cyvcf2 expects 0-based positions
        ins_seq = rec.ALT[0]
        rec.REF = anchor
        rec.ALT = [anchor + ins_seq]


parser = argparse.ArgumentParser(
    description=(
        "Keep only PASS variants, remove BND/TRA/INV, and convert remaining "
        "records to sequence-resolved representations when possible."
    )
)
parser.add_argument("in_vcf", help="Input VCF")
parser.add_argument("ref_fa", help="Reference FASTA (indexed with .fai)")
parser.add_argument("out_vcf", help="Output VCF")
args = parser.parse_args()

ref = pysam.FastaFile(args.ref_fa)
vcf = VCF(args.in_vcf)
w = Writer(args.out_vcf, vcf)

for rec in vcf:
    if rec.FILTER is not None:
        # Skip records that are not PASS
        continue

    svt = rec.INFO.get("SVTYPE")
    if svt in {"BND", "TRA", "INV"}:
        # Skip unsupported SV types
        continue

    if svt == "DUP":
        resolve_dup(rec, ref)
        w.write_record(rec)

    elif svt == "DUP:TANDEM":
        resolve_dup_tandem(rec, ref)
        w.write_record(rec)

    elif svt == "DEL":
        resolve_del(rec, ref)
        # Write only non-symbolic deletions
        if rec.ALT[0] != "<DEL>":
            w.write_record(rec)

    elif svt == "INS":
        resolve_ins(rec, ref)
        # Write only non-symbolic insertions
        if rec.ALT[0] != "<INS>":
            w.write_record(rec)

    else:
        # SNVs or already sequence-resolved records
        w.write_record(rec)

w.close()
vcf.close()
ref.close()
