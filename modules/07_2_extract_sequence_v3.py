"""Compatibility sequence-extraction module for the corrected candidate-gene pipeline.

This module extracts the SNP-centered reference sequence from a reference FASTA
and emits REF/ALT allele sequences.  It is deliberately strict about coordinate
and REF validation so downstream BLAST results are not built from fabricated
sequences.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pipeline_common import clean_dna, first_value, read_table, reverse_complement, write_table

CHR_COLS = ["Chromosome", "Chr", "chromosome", "CHROM", "Local_Chromosome", "Ensembl_Chromosome"]
POS_COLS = ["Position", "POS", "position", "BP", "bp", "Local_Position", "Ensembl_Position"]
REF_COLS = ["REF", "Ref", "Reference", "reference_allele", "Allele1"]
ALT_COLS = ["ALT", "Alt", "Alternate", "alternate_allele", "Allele2"]


def read_fasta(path: str) -> dict[str, str]:
    records: dict[str, str] = {}
    current = None
    chunks: list[str] = []
    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if current is not None:
                    records[current] = clean_dna("".join(chunks))
                current = line[1:].split()[0]
                chunks = []
            elif current is not None:
                chunks.append(line)
    if current is not None:
        records[current] = clean_dna("".join(chunks))
    return records


def get_contig(records: dict[str, str], chrom: object) -> str | None:
    value = str(chrom)
    candidates = [value]
    if value.lower().startswith("chr"):
        candidates.append(value[3:])
    else:
        candidates.append(f"chr{value}")
    for candidate in candidates:
        if candidate in records:
            return candidate
    for name in records:
        if name.split("|")[0] == value:
            return name
    return None


def run(input_path: str, fasta_path: str, output_path: str, flank: int = 250) -> pd.DataFrame:
    if flank < 0:
        raise ValueError("flank must be >= 0")

    table = read_table(input_path).copy()
    if table.empty:
        raise ValueError("Input SNP table is empty")
    if "SNP_ID" not in table.columns:
        raise ValueError("Input table must contain SNP_ID")

    records = read_fasta(fasta_path)
    rows = []

    for _, row in table.iterrows():
        sid = str(row["SNP_ID"])
        ref = clean_dna(first_value(row, REF_COLS, ""))
        alt = clean_dna(first_value(row, ALT_COLS, ""))
        chrom = first_value(row, CHR_COLS, "")
        pos_value = first_value(row, POS_COLS, None)
        result = {
            "SNP_ID": sid,
            "REF": ref,
            "ALT": alt,
            "Sequence_Source": "Reference_FASTA",
            "REF_Sequence": "",
            "ALT_Sequence": "",
            "REF_ReverseComplement": "",
            "ALT_ReverseComplement": "",
            "Variant_Offset": pd.NA,
            "sequence_status": "ERROR",
            "sequence_error": "",
        }

        try:
            if not ref or not alt:
                raise ValueError("Missing REF or ALT allele")
            if ref == alt:
                raise ValueError("REF and ALT alleles are identical")
            if pd.isna(pos_value):
                raise ValueError("Missing SNP position")

            position = int(float(pos_value))
            if position < 1:
                raise ValueError("SNP position must be 1-based and >= 1")

            contig = get_contig(records, chrom)
            if contig is None:
                raise ValueError(f"Reference contig not found for chromosome {chrom}")

            sequence = records[contig]
            start = max(1, position - flank)
            end = min(len(sequence), position + flank)
            ref_sequence = sequence[start - 1 : end]
            offset = position - start

            observed = ref_sequence[offset : offset + len(ref)]
            if observed != ref:
                raise ValueError(f"REF mismatch: FASTA has {observed}, expected {ref}")

            if offset + len(ref) > len(ref_sequence):
                raise ValueError("REF allele extends beyond extracted sequence")

            alt_sequence = ref_sequence[:offset] + alt + ref_sequence[offset + len(ref) :]
            result.update(
                {
                    "REF_Sequence": ref_sequence,
                    "ALT_Sequence": alt_sequence,
                    "REF_ReverseComplement": reverse_complement(ref_sequence),
                    "ALT_ReverseComplement": reverse_complement(alt_sequence),
                    "Variant_Offset": offset,
                    "sequence_status": "OK",
                }
            )
        except Exception as exc:
            result["sequence_error"] = str(exc)

        rows.append(result)

    out = pd.DataFrame(rows)
    write_table(out, output_path)
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract SNP-centered REF/ALT sequences from a reference FASTA")
    parser.add_argument("input")
    parser.add_argument("fasta")
    parser.add_argument("output")
    parser.add_argument("--flank", type=int, default=250)
    args = parser.parse_args()
    run(args.input, args.fasta, args.output, args.flank)
