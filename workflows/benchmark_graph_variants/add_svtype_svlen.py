#!/usr/bin/env python3
import argparse
from cyvcf2 import VCF, Writer

def trim_with_anchor(ref, alt):
    # Trim identical suffixes
    while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    # Trim identical prefixes
    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
    return ref, alt

def classify_sv(trimmed_ref, trimmed_alt):
    svlen = len(trimmed_alt) - len(trimmed_ref)
    is_sv = abs(svlen) >= 50 or max(len(trimmed_ref), len(trimmed_alt)) >= 50
    if not is_sv:
        return None, svlen
    return ("DEL" if svlen < 0 else "INS"), svlen

def main(in_vcf, out_vcf):
    invcf = VCF(in_vcf)

    if not invcf.contains("SVTYPE"):
        invcf.add_info_to_header({
            "ID": "SVTYPE",
            "Description": "Type of structural variant",
            "Type": "String",
            "Number": "1"
        })
    if not invcf.contains("SVLEN"):
        invcf.add_info_to_header({
            "ID": "SVLEN",
            "Description": "Length of structural variant",
            "Type": "Integer",
            "Number": "A"
        })

    outvcf = Writer(out_vcf, invcf)

    for rec in invcf:
        ref = rec.REF
        alts = rec.ALT or []
        if not alts:
            outvcf.write_record(rec)
            continue

        alt = alts[0]  # biallelic assumption

        # Trimming for measurement (doesn't change REF/ALT in file)
        trimmed_ref, trimmed_alt = trim_with_anchor(ref, alt)

        svtype, svlen = classify_sv(trimmed_ref, trimmed_alt)

        # Set SVTYPE only if it qualifies as SV
        if svtype is not None:
            rec.INFO["SVTYPE"] = svtype
            rec.INFO["SVLEN"] = svlen
        
        outvcf.write_record(rec)

    outvcf.close()
    invcf.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add INFO/SVTYPE and INFO/SVLEN with trimming.")
    parser.add_argument("in_vcf", help="Input biallelic, sequence-resolved VCF")
    parser.add_argument("-o", "--out_vcf", required=True, help="Output VCF")
    args = parser.parse_args()
    main(args.in_vcf, args.out_vcf)
