"""Module 12: end-to-end orchestrator.

Usage:
python modules/12_run_pipeline.py --input filtered_snps.csv --reference reference.fa \
  --ncbi-db db/ncbi/wheat_dna_db --ensembl-db db/ensembl/wheat_dna_db \
  --gff annotation.gff3 --out results
"""
from __future__ import annotations
import argparse, importlib, sys
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"modules"))

from pipeline_common import read_table, write_table, setup_logging


def run(args):
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    m1=importlib.import_module("01_prepare_snps"); m2=importlib.import_module("02_generate_allele_sequences"); m3=importlib.import_module("03_local_blast"); m4=importlib.import_module("04_coordinate_validation"); m5=importlib.import_module("05_merge_coordinates"); m6=importlib.import_module("06_local_neighborhood"); m7=importlib.import_module("07_api_neighborhood"); m8=importlib.import_module("08_functional_annotation"); m9=importlib.import_module("09_pubmed_literature"); m10=importlib.import_module("10_score_candidates"); m11=importlib.import_module("11_generate_master_workbook")
    p1=out/"01_prepared.csv"; p2=out/"02_alleles.csv"; fasta=out/"02_alleles.fasta"; p3=out/"03_local_blast.csv"; p4=out/"04_coordinate_validation.csv"; p5=out/"05_coordinates_merged.csv"; p6=out/"06_local_neighborhood.csv"; p7=out/"07_api_neighborhood.csv"; p8=out/"08_annotation.csv"; p9=out/"09_pubmed.csv"; p10=out/"10_scored.csv"; final=out/"candidate_gene_master.xlsx"
    m1.run(args.input,str(p1)); m2.run(str(p1),str(fasta),str(p2),args.allow_midpoint); m3.run(str(fasta),str(p3),args.ncbi_db,args.ensembl_db); m4.run(str(p3),str(p4),args.species,args.min_identity,args.min_coverage); m5.run(str(p1),str(p4),str(p5));
    if args.gff: m6.run(str(p4),args.gff,str(p6),args.window); m7.run(str(p4),str(p6),str(p7),args.species,args.window)
    else:
        pd.DataFrame(columns=["SNP_ID"]).to_csv(p6,index=False); pd.DataFrame(columns=["SNP_ID"]).to_csv(p7,index=False)
    # Build annotation input from coordinate + API neighborhood + local neighborhood.
    ann=read_table(p4); 
    for p in [p6,p7]:
        d=read_table(p)
        if not d.empty: ann=ann.merge(d.drop_duplicates("SNP_ID"),on="SNP_ID",how="left")
    write_table(ann,str(out/"annotation_input.csv")); m8.run(str(out/"annotation_input.csv"),str(p8),args.species,args.ncbi_email)
    litin=read_table(p8); m9.run(str(p8),str(p9),args.ncbi_email,args.max_articles)
    scored=read_table(p5)
    for p in [p7,p8,p9]:
        d=read_table(p)
        if not d.empty: scored=scored.merge(d.drop_duplicates("SNP_ID"),on="SNP_ID",how="left",suffixes=("","_extra"))
    write_table(scored,str(out/"scoring_input.csv")); m10.run(str(out/"scoring_input.csv"),str(p10))
    m11.run(str(p10),str(final),blast=str(p3),validation=str(p4),local_neighborhood=str(p6),api_neighborhood=str(p7),annotation=str(p8),literature=str(p9),scoring=str(p10))
    print(final)

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--input",required=True); ap.add_argument("--reference",default=""); ap.add_argument("--ncbi-db",default=""); ap.add_argument("--ensembl-db",default=""); ap.add_argument("--gff",default=""); ap.add_argument("--out",default="results"); ap.add_argument("--species",default="triticum_aestivum"); ap.add_argument("--window",type=int,default=100000); ap.add_argument("--min-identity",type=float,default=90); ap.add_argument("--min-coverage",type=float,default=70); ap.add_argument("--allow-midpoint",action="store_true"); ap.add_argument("--ncbi-email",default=""); ap.add_argument("--max-articles",type=int,default=20); args=ap.parse_args(); setup_logging(); run(args)
