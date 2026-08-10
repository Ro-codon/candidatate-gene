"""Module 12: end-to-end orchestrator for the corrected candidate-gene pipeline."""
from __future__ import annotations

import argparse
import importlib
import logging
import sys
from pathlib import Path

import pandas as pd

# modules/12_run_pipeline.py lives inside the modules directory.
MODULES_DIR = Path(__file__).resolve().parent
ROOT = MODULES_DIR.parent

if str(MODULES_DIR) not in sys.path:
    sys.path.insert(0, str(MODULES_DIR))

from pipeline_common import read_table, write_table, setup_logging

LOGGER = logging.getLogger("candidate_gene_pipeline")

MODULE_NAMES = [
    "01_prepare_snps",
    "02_generate_allele_sequences",
    "03_local_blast",
    "04_coordinate_validation",
    "05_merge_coordinates",
    "06_local_neighborhood",
    "07_api_neighborhood",
    "08_functional_annotation",
    "09_pubmed_literature",
    "10_score_candidates",
    "11_generate_master_workbook",
]


def resolve_path(value: str) -> str:
    """Resolve a user-supplied path relative to the project root."""
    if not value:
        return ""
    path = Path(value)
    if not path.is_absolute():
        path = ROOT / path
    return str(path)


def merge_unique(base: pd.DataFrame, path) -> pd.DataFrame:
    """Merge one SNP_ID keyed table without overwriting existing columns."""
    path = Path(path)
    if not path.exists():
        LOGGER.warning("Merge file does not exist: %s", path)
        return base

    data = read_table(path)
    if data.empty or "SNP_ID" not in data.columns:
        return base

    data = data.drop_duplicates(subset=["SNP_ID"], keep="first")
    return base.merge(
        data,
        on="SNP_ID",
        how="left",
        suffixes=("", "_extra"),
    )


def load_modules():
    """Load numbered modules in pipeline order."""
    loaded = []
    for name in MODULE_NAMES:
        LOGGER.info("Loading %s", name)
        loaded.append(importlib.import_module(name))
    return tuple(loaded)


def run(args):
    out = Path(args.out)
    if not out.is_absolute():
        out = ROOT / out
    out.mkdir(parents=True, exist_ok=True)

    input_path = Path(resolve_path(args.input))
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    (
        m1, m2, m3, m4, m5, m6, m7,
        m8, m9, m10, m11,
    ) = load_modules()

    p1 = out / "01_prepared.csv"
    p2 = out / "02_alleles.csv"
    fasta = out / "02_alleles.fasta"
    p3 = out / "03_local_blast.csv"
    p4 = out / "04_coordinate_validation.csv"
    p5 = out / "05_coordinates_merged.csv"
    p6 = out / "06_local_neighborhood.csv"
    p7 = out / "07_api_neighborhood.csv"
    p8 = out / "08_annotation.csv"
    p9 = out / "09_pubmed.csv"
    p10 = out / "10_scored.csv"
    final = out / "candidate_gene_master.xlsx"

    LOGGER.info("MODULE 01: prepare SNPs")
    m1.run(str(input_path), str(p1))

    LOGGER.info("MODULE 02: generate REF/ALT sequences")
    m2.run(
        str(p1),
        str(fasta),
        str(p2),
        args.allow_midpoint,
        resolve_path(args.reference),
        args.flank,
    )

    LOGGER.info("MODULE 03: local BLAST")
    m3.run(
        str(fasta),
        str(p3),
        resolve_path(args.ncbi_db),
        resolve_path(args.ensembl_db),
        args.online_ncbi,
    )

    LOGGER.info("MODULE 04: coordinate/API validation")
    m4.run(
        str(p3),
        str(p4),
        args.species,
        args.min_identity,
        args.min_coverage,
    )

    LOGGER.info("MODULE 05: merge coordinates with original SNP rows")
    m5.run(str(p1), str(p4), str(p5))

    LOGGER.info("MODULE 06/07: local and Ensembl API neighborhood")
    if args.gff:
        gff = Path(resolve_path(args.gff))
        if not gff.exists():
            raise FileNotFoundError(f"GFF/GTF file not found: {gff}")
        m6.run(str(p4), str(gff), str(p6), args.window)
        m7.run(str(p4), str(p6), str(p7), args.species, args.window)
    else:
        LOGGER.warning("No --gff supplied; neighborhood tables will be empty.")
        pd.DataFrame(columns=["SNP_ID"]).to_csv(p6, index=False)
        pd.DataFrame(columns=["SNP_ID"]).to_csv(p7, index=False)

    LOGGER.info("MODULE 08: functional annotation")
    annotation_input = out / "annotation_input.csv"
    annotation_data = read_table(p4)
    annotation_data = merge_unique(annotation_data, p6)
    annotation_data = merge_unique(annotation_data, p7)
    write_table(annotation_data, str(annotation_input))
    m8.run(
        str(annotation_input),
        str(p8),
        args.species,
        args.ncbi_email,
    )

    LOGGER.info("MODULE 09: PubMed literature mining")
    m9.run(
        str(p8),
        str(p9),
        args.ncbi_email,
        args.max_articles,
    )

    LOGGER.info("MODULE 10: candidate scoring")
    scoring_input = out / "scoring_input.csv"
    scored_data = read_table(p5)
    scored_data = merge_unique(scored_data, p2)
    scored_data = merge_unique(scored_data, p7)
    scored_data = merge_unique(scored_data, p8)
    scored_data = merge_unique(scored_data, p9)
    write_table(scored_data, str(scoring_input))
    m10.run(str(scoring_input), str(p10))

    LOGGER.info("MODULE 11: master Excel workbook")
    m11.run(
        str(p10),
        str(final),
        blast=str(p3),
        validation=str(p4),
        local_neighborhood=str(p6),
        api_neighborhood=str(p7),
        annotation=str(p8),
        literature=str(p9),
        scoring=str(p10),
    )

    LOGGER.info("Pipeline completed: %s", final)
    return final


def build_parser():
    parser = argparse.ArgumentParser(
        description="End-to-end SNP-to-candidate-gene annotation pipeline."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--reference", default="")
    parser.add_argument("--flank", type=int, default=250)
    parser.add_argument("--ncbi-db", default="")
    parser.add_argument("--ensembl-db", default="")
    parser.add_argument(
        "--online-ncbi",
        action="store_true",
        help="Enable optional online NCBI BLAST cross-check.",
    )
    parser.add_argument("--gff", default="")
    parser.add_argument("--out", default="results")
    parser.add_argument("--species", default="triticum_aestivum")
    parser.add_argument("--window", type=int, default=100000)
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=70.0)
    parser.add_argument("--allow-midpoint", action="store_true")
    parser.add_argument("--ncbi-email", default="")
    parser.add_argument("--max-articles", type=int, default=20)
    return parser


def main():
    setup_logging()
    args = build_parser().parse_args()
    return run(args)


if __name__ == "__main__":
    main()
