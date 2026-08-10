"""Module 12: end-to-end orchestrator for the corrected candidate-gene pipeline."""
from __future__ import annotations
import argparse, importlib, sys
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"modules"))
from pipeline_common import read_table, write_table, setup_logging

def merge_unique(base,path):
    d=read_table(path)
    if d.empty or "SNP_ID" not in d:return base
    return base.merge(d.drop_duplicates("SNP_ID"),on="SNP_ID",how="left",suffixes=("","_extra"))

def run(args):
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    mods=[importlib.import_module(f"{i:02d}_{name}") for i,name in [(1,"prepare_snps"),(2,"generate_allele_sequences"),(3,"local_blast"),(4,"coordinate_validation"),(5,"merge_coordinates"),(6,"local_neighborhood"),(7,"api_neighborhood"),(8,"functional_annotation"),(9,"pubmed_literature"),(10,"score_candidates"),(11,"generate_master_workbook")]]
    m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11=mods
    p1=out/"01_prepared.csv"; p2=out/"02_alleles.csv"; fasta=out/"02_alleles.fasta"; p3=out/"03_local_blast.csv"; p4=out/"04_coordinate_validation.csv"; p5=out/"05_coordinates_merged.csv"; p6=out/"06_local_neighborhood.csv"; p7=out/"07_api_neighborhood.csv"; p8=out/"08_annotation.csv"; p9=out/"09_pubmed.csv"; p10=out/"10_scored.csv"; final=out/"candidate_gene_master.xlsx"
    m1.run(args.input,str(p1)); m2.run(str(p1),str(fasta),str(p2),args.allow_midpoint,args.reference,args.flank); m3.run(str(fasta),str(p3),args.ncbi_db,args.ensembl_db,args.online_ncbi); m4.run(str(p3),str(p4),args.species,args.min_identity,args.min_coverage); m5.run(str(p1),str(p4),str(p5))
    if args.gff: m6.run(str(p4),args.gff,str(p6),args.window); m7.run(str(p4),str(p6),str(p7),args.species,args.window)
    else: pd.DataFrame(columns=["SNP_ID"]).to_csv(p6,index=False); pd.DataFrame(columns=["SNP_ID"]).to_csv(p7,index=False)
    ann=read_table(p4)
    for p in [p6,p7]: ann=merge_unique(ann,p)
    write_table(ann,str(out/"annotation_input.csv")); m8.run(str(out/"annotation_input.csv"),str(p8),args.species,args.ncbi_email); m9.run(str(p8),str(p9),args.ncbi_email,args.max_articles)
    scored=read_table(p5)
    for p in [p2,p7,p8,p9]: scored=merge_unique(scored,p)
    write_table(scored,str(out/"scoring_input.csv")); m10.run(str(out/"scoring_input.csv"),str(p10)); m11.run(str(p10),str(final),blast=str(p3),validation=str(p4),local_neighborhood=str(p6),api_neighborhood=str(p7),annotation=str(p8),literature=str(p9),scoring=str(p10)); return final

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--input",required=True); ap.add_argument("--reference",default=""); ap.add_argument("--flank",type=int,default=250); ap.add_argument("--ncbi-db",default=""); ap.add_argument("--ensembl-db",default=""); ap.add_argument("--online-ncbi",action="store_true"); ap.add_argument("--gff",default=""); ap.add_argument("--out",default="results"); ap.add_argument("--species",default="triticum_aestivum"); ap.add_argument("--window",type=int,default=100000); ap.add_argument("--min-identity",type=float,default=90); ap.add_argument("--min-coverage",type=float,default=70); ap.add_argument("--allow-midpoint",action="store_true"); ap.add_argument("--ncbi-email",default=""); ap.add_argument("--max-articles",type=int,default=20); args=ap.parse_args(); setup_logging(); print(run(args))
