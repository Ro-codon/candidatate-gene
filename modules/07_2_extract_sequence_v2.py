"""Sequence extraction utilities for candidate gene pipelines.

This module extracts nucleotide or protein sequences for genomic intervals,
transcripts, genes, or variant windows from common reference file formats.

Supported inputs
----------------
- FASTA references
- GFF3/GTF annotations
- BED-like interval specifications
- Direct genomic coordinates

The implementation is intentionally dependency-light and defensive so it can be
used inside larger annotation workflows without forcing a specific framework.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union
import gzip
import os
import re

DNA_COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


@dataclass(frozen=True)
class SequenceRecord:
    """Simple FASTA record container."""

    identifier: str
    sequence: str
    description: str = ""


@dataclass(frozen=True)
class GenomicInterval:
    """1-based closed genomic interval."""

    chrom: str
    start: int
    end: int
    strand: str = "+"
    label: str = ""

    def normalized(self) -> "GenomicInterval":
        if self.start <= self.end:
            return self
        return GenomicInterval(
            chrom=self.chrom,
            start=self.end,
            end=self.start,
            strand=self.strand,
            label=self.label,
        )

    @property
    def length(self) -> int:
        n = self.normalized()
        return n.end - n.start + 1


class SequenceExtractionError(RuntimeError):
    """Raised when a sequence extraction operation fails."""


def _open_text(path: Union[str, os.PathLike]):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def read_fasta(path: Union[str, os.PathLike]) -> Dict[str, SequenceRecord]:
    """Read a FASTA file into memory.

    Parameters
    ----------
    path:
        Path to a FASTA or FASTA.GZ file.

    Returns
    -------
    dict
        Mapping from sequence identifier to SequenceRecord.
    """
    records: Dict[str, SequenceRecord] = {}
    current_id: Optional[str] = None
    current_desc = ""
    chunks: List[str] = []

    with _open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = SequenceRecord(
                        identifier=current_id,
                        sequence="".join(chunks),
                        description=current_desc,
                    )
                header = line[1:].strip()
                if not header:
                    raise SequenceExtractionError("Encountered empty FASTA header.")
                parts = header.split(None, 1)
                current_id = parts[0]
                current_desc = parts[1] if len(parts) > 1 else ""
                chunks = []
            else:
                if current_id is None:
                    raise SequenceExtractionError("FASTA sequence data found before any header.")
                chunks.append(line)

    if current_id is not None:
        records[current_id] = SequenceRecord(
            identifier=current_id,
            sequence="".join(chunks),
            description=current_desc,
        )

    if not records:
        raise SequenceExtractionError(f"No FASTA records found in {path!s}")

    return records


_GFF_ATTR_SPLIT = re.compile(r"[;\n\r]+")
_GFF_KV_SPLIT = re.compile(r"[=\s]")


def parse_gff_attributes(attr_field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    if not attr_field or attr_field == ".":
        return attrs
    for part in _GFF_ATTR_SPLIT.split(attr_field.strip()):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(None, 1)
        else:
            key, value = part, ""
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def read_gff_intervals(
    path: Union[str, os.PathLike],
    feature_types: Sequence[str] = ("gene", "mRNA", "transcript", "CDS", "exon"),
) -> List[GenomicInterval]:
    """Parse intervals from GFF3/GTF annotations."""
    wanted = {ft.lower() for ft in feature_types}
    intervals: List[GenomicInterval] = []

    with _open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, feature_type, start, end, _score, strand, _phase, attrs = fields
            if wanted and feature_type.lower() not in wanted:
                continue
            parsed = parse_gff_attributes(attrs)
            label = parsed.get("ID") or parsed.get("Name") or parsed.get("gene_id") or parsed.get("transcript_id") or feature_type
            try:
                intervals.append(
                    GenomicInterval(
                        chrom=chrom,
                        start=int(start),
                        end=int(end),
                        strand=strand if strand in {"+", "-"} else "+",
                        label=label,
                    )
                )
            except ValueError as exc:
                raise SequenceExtractionError(f"Invalid coordinates in GFF line: {line}") from exc
    return intervals


def reverse_complement(seq: str) -> str:
    return seq.translate(DNA_COMPLEMENT)[::-1]


def sanitize_sequence(seq: str) -> str:
    return re.sub(r"\s+", "", seq).upper()


def extract_interval_sequence(reference: Dict[str, SequenceRecord], interval: GenomicInterval) -> str:
    """Extract 1-based closed interval sequence from an in-memory reference."""
    interval = interval.normalized()
    if interval.chrom not in reference:
        raise SequenceExtractionError(f"Chromosome/contig '{interval.chrom}' not found in reference.")
    seq = sanitize_sequence(reference[interval.chrom].sequence)
    if interval.start < 1 or interval.end < 1:
        raise SequenceExtractionError("Genomic coordinates must be positive and 1-based.")
    if interval.end > len(seq):
        raise SequenceExtractionError(
            f"Interval {interval.chrom}:{interval.start}-{interval.end} exceeds contig length {len(seq)}."
        )
    subseq = seq[interval.start - 1 : interval.end]
    if interval.strand == "-":
        return reverse_complement(subseq)
    return subseq


def extract_fasta_by_intervals(
    fasta_path: Union[str, os.PathLike],
    intervals: Iterable[GenomicInterval],
) -> Dict[str, str]:
    """Return extracted sequences keyed by interval label or coordinate string."""
    reference = read_fasta(fasta_path)
    out: Dict[str, str] = {}
    for iv in intervals:
        key = iv.label or f"{iv.chrom}:{iv.start}-{iv.end}({iv.strand})"
        out[key] = extract_interval_sequence(reference, iv)
    return out


def extract_windows_around_positions(
    fasta_path: Union[str, os.PathLike],
    chrom: str,
    positions: Iterable[int],
    flank: int = 250,
    strand: str = "+",
    prefix: str = "pos",
) -> Dict[str, str]:
    reference = read_fasta(fasta_path)
    if chrom not in reference:
        raise SequenceExtractionError(f"Chromosome/contig '{chrom}' not found in reference.")
    seq_len = len(sanitize_sequence(reference[chrom].sequence))
    out: Dict[str, str] = {}
    for pos in positions:
        start = max(1, int(pos) - int(flank))
        end = min(seq_len, int(pos) + int(flank))
        iv = GenomicInterval(chrom=chrom, start=start, end=end, strand=strand, label=f"{prefix}_{pos}")
        out[iv.label] = extract_interval_sequence(reference, iv)
    return out


def write_fasta(records: Dict[str, str], path: Union[str, os.PathLike], line_width: int = 60) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wt", encoding="utf-8") as handle:
        for header, seq in records.items():
            handle.write(f">{header}\n")
            clean = sanitize_sequence(seq)
            for i in range(0, len(clean), line_width):
                handle.write(clean[i : i + line_width] + "\n")


def intervals_from_bed(
    path: Union[str, os.PathLike],
    one_based: bool = False,
    default_strand: str = "+",
) -> List[GenomicInterval]:
    intervals: List[GenomicInterval] = []
    with _open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(fields[1]) + (0 if one_based else 1)
            end = int(fields[2])
            name = fields[3] if len(fields) >= 4 else ""
            strand = fields[5] if len(fields) >= 6 and fields[5] in {"+", "-"} else default_strand
            intervals.append(GenomicInterval(chrom=chrom, start=start, end=end, strand=strand, label=name))
    return intervals


class SequenceExtractor:
    """Convenience wrapper for repeated extraction tasks."""

    def __init__(self, fasta_path: Union[str, os.PathLike]):
        self.fasta_path = str(fasta_path)
        self.reference = read_fasta(self.fasta_path)

    def extract(self, interval: GenomicInterval) -> str:
        return extract_interval_sequence(self.reference, interval)

    def extract_many(self, intervals: Iterable[GenomicInterval]) -> Dict[str, str]:
        out: Dict[str, str] = {}
        for iv in intervals:
            key = iv.label or f"{iv.chrom}:{iv.start}-{iv.end}({iv.strand})"
            out[key] = self.extract(iv)
        return out

    def extract_gene_regions(
        self,
        gff_path: Union[str, os.PathLike],
        feature_types: Sequence[str] = ("gene",),
    ) -> Dict[str, str]:
        intervals = read_gff_intervals(gff_path, feature_types=feature_types)
        return self.extract_many(intervals)


__all__ = [
    "SequenceRecord",
    "GenomicInterval",
    "SequenceExtractor",
    "SequenceExtractionError",
    "read_fasta",
    "read_gff_intervals",
    "intervals_from_bed",
    "extract_interval_sequence",
    "extract_fasta_by_intervals",
    "extract_windows_around_positions",
    "write_fasta",
    "reverse_complement",
    "sanitize_sequence",
    "parse_gff_attributes",
]
