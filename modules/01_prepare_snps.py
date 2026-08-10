"""Module 01: normalize and validate the previously filtered SNP table."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from pipeline_common import first_value, normalize_snp_id, read_table, write_table

ID_COLS=["SNP_ID","SNP","rsID","rs_id","Marker","Marker_ID","SNPID","ID"]
TRAIT_COLS=["Trait","trait","Phenotype","phenotype"]
P_COLS=["GWAS_P","GWAS_PVALUE","P_value","P-value","p_value","PValue","pvalue"]
REF_COLS=["REF","Ref","Reference","reference_allele","Allele1"]
ALT_COLS=["ALT","Alt","Alternate","alternate_allele","Allele2"]


def run(input_path: str, output_path: str) -> pd.DataFrame:
    df=read_table(input_path).copy()
    if df.empty: raise ValueError("Input SNP table is empty")
    ids=[]
    for i,row in df.iterrows(): ids.append(normalize_snp_id(first_value(row,ID_COLS),i+1))
    df.insert(0,"SNP_ID",ids) if "SNP_ID" not in df.columns else df.__setitem__("SNP_ID",ids)
    df["Trait"]=[first_value(r,TRAIT_COLS) for _,r in df.iterrows()] if "Trait" not in df.columns else df["Trait"]
    df["GWAS_P"]=[first_value(r,P_COLS) for _,r in df.iterrows()] if "GWAS_P" not in df.columns else df["GWAS_P"]
    df["REF"]=[str(first_value(r,REF_COLS,"")).upper() for _,r in df.iterrows()] if "REF" not in df.columns else df["REF"].astype(str).str.upper()
    df["ALT"]=[str(first_value(r,ALT_COLS,"")).upper() for _,r in df.iterrows()] if "ALT" not in df.columns else df["ALT"].astype(str).str.upper()
    df["GWAS_P"] = pd.to_numeric(df["GWAS_P"],errors="coerce")
    df["input_row"]=range(1,len(df)+1)
    df["alleles_valid"]=(df["REF"].str.fullmatch(r"[ACGTN]+",na=False) & df["ALT"].str.fullmatch(r"[ACGTN]+",na=False) & (df["REF"]!="") & (df["ALT"]!=""))
    df["is_biallelic_snp"]=(df["REF"].str.len()==1)&(df["ALT"].str.len()==1)&(df["REF"]!=df["ALT"])
    write_table(df,output_path); return df

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("output"); a=ap.parse_args(); run(a.input,a.output)
