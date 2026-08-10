"""Module 07: Ensembl Plants REST overlap validation of local neighborhood."""
from __future__ import annotations
import argparse
import pandas as pd
from pipeline_common import request_json, read_table, write_table

def run(coords_path, local_path, output_path, species="triticum_aestivum", window=100000):
    c=read_table(coords_path); l=read_table(local_path); rows=[]
    local=l.groupby("SNP_ID") if not l.empty and "SNP_ID" in l else {}
    for _,r in c.drop_duplicates("SNP_ID").iterrows():
        sid=r.SNP_ID; chrom=str(r.get("Ensembl_Chromosome") or r.get("Local_Chromosome") or ""); pos=pd.to_numeric(r.get("Ensembl_Position") if pd.notna(r.get("Ensembl_Position")) else r.get("Local_Position"),errors="coerce")
        base={"SNP_ID":sid,"api_validated":False,"API_Gene_ID":"","API_Gene_Name":"","API_Gene_Start":pd.NA,"API_Gene_End":pd.NA,"API_Gene_Strand":"","API_Distance_bp":pd.NA,"api_validation_status":"NO_COORDINATE","api_error":""}
        if not chrom or pd.isna(pos): rows.append(base); continue
        region=f"{chrom}:{max(1,int(pos)-window)}-{int(pos)+window}"; data,status,err=request_json(f"https://rest.ensembl.org/overlap/region/{species}/{region}",params={"feature":"gene"},headers={"Content-Type":"application/json","Accept":"application/json"})
        if not isinstance(data,list): base["api_validation_status"]="API_UNAVAILABLE"; base["api_error"]=err or f"HTTP {status}"; rows.append(base); continue
        genes=[x for x in data if x.get("feature_type")=="gene" or x.get("object_type")=="Gene"]; genes.sort(key=lambda x: min(abs(int(pos)-int(x.get("start",pos))),abs(int(pos)-int(x.get("end",pos)))))
        if genes:
            x=genes[0]; st=int(x.get("start",0)); en=int(x.get("end",0)); dist=0 if st<=int(pos)<=en else min(abs(int(pos)-st),abs(int(pos)-en)); base.update({"api_validated":True,"API_Gene_ID":x.get("id",""),"API_Gene_Name":x.get("external_name","") or x.get("display_name",""),"API_Gene_Start":st,"API_Gene_End":en,"API_Gene_Strand":x.get("strand",""),"API_Distance_bp":dist,"api_validation_status":"HIT"})
        else: base["api_validation_status"]="NO_API_GENE"
        rows.append(base)
    out=pd.DataFrame(rows); write_table(out,output_path); return out
if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("coordinates"); ap.add_argument("local"); ap.add_argument("output"); ap.add_argument("--species",default="triticum_aestivum"); ap.add_argument("--window",type=int,default=100000); a=ap.parse_args(); run(a.coordinates,a.local,a.output,a.species,a.window)
