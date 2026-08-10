"""Module 04: resolve best local hit and cross-check Ensembl Plants REST.

Ensembl is used as an independent validation source. API failures are recorded,
never converted into a false validation.
"""
from __future__ import annotations
import argparse
import pandas as pd
from pipeline_common import request_json, read_table, write_table

BASE="https://rest.ensembl.org"

def parse_subject(subject):
    import re
    s=str(subject)
    m=re.search(r"(?:chr|chromosome)[^0-9A-Za-z]*(\w+)[^0-9]*(\d+)",s,re.I)
    if m: return m.group(1),int(m.group(2))
    m=re.search(r"[:|_]([0-9]+)$",s)
    return "",int(m.group(1)) if m else None

def ensembl_lookup(species, stable_id):
    url=f"{BASE}/lookup/id/{stable_id}"
    return request_json(url,params={"expand":1},headers={"Content-Type":"application/json","Accept":"application/json"})

def ensembl_overlap(species, chrom, pos, window=1):
    if chrom in {"", "nan", "None"} or pos is None: return None,None,"missing_coordinate"
    region=f"{chrom}:{max(1,pos-window)}-{pos+window}"
    url=f"{BASE}/overlap/region/{species}/{region}"
    return request_json(url,params={"feature":"gene"},headers={"Content-Type":"application/json","Accept":"application/json"})

def run(blast_path, output_path, species="triticum_aestivum", min_identity=90, min_coverage=70):
    b=read_table(blast_path); rows=[]
    if b.empty:
        write_table(pd.DataFrame(columns=["SNP_ID","Best_Allele","Local_Database","Local_Subject","Local_Chromosome","Local_Position","Local_Identity","Local_Coverage","Ensembl_Chromosome","Ensembl_Position","Ensembl_ID","Ensembl_Gene","api_validated","validation_status","api_error"]),output_path); return pd.DataFrame()
    b["Identity"]=pd.to_numeric(b["Identity"],errors="coerce"); b["Coverage"]=pd.to_numeric(b["Coverage"],errors="coerce")
    for sid,g in b.groupby("SNP_ID",dropna=False):
        hits=g[(g.Identity>=min_identity)&(g.Coverage>=min_coverage)].sort_values(["Bitscore","Identity","Coverage"],ascending=False)
        if hits.empty: rows.append({"SNP_ID":sid,"validation_status":"NO_QUALIFYING_LOCAL_HIT","api_validated":False}); continue
        h=hits.iloc[0]; chrom,pos=parse_subject(h.subject_id); api_valid=False; api_error=""; ens_id=""; ens_gene=""; ens_chr=""; ens_pos=None
        data,status,err=ensembl_overlap(species,chrom,pos,1) if pos is not None else (None,None,"no_position")
        if isinstance(data,list):
            genes=[x for x in data if x.get("feature_type")=="gene" or x.get("object_type")=="Gene"]
            if genes:
                x=genes[0]; ens_id=x.get("id",""); ens_gene=x.get("external_name","") or x.get("display_name",""); ens_chr=x.get("seq_region_name",""); ens_pos=x.get("start"); api_valid=True
        else: api_error=err or f"HTTP {status}"
        same = api_valid and str(ens_chr)==str(chrom) and (pos is None or ens_pos is None or abs(int(ens_pos)-int(pos))<=1)
        rows.append({"SNP_ID":sid,"Best_Allele":h.get("Allele",""),"Local_Database":h.get("Database_Source",""),"Local_Subject":h.get("subject_id",""),"Local_Chromosome":chrom,"Local_Position":pos,"Local_Identity":h.get("Identity"),"Local_Coverage":h.get("Coverage"),"Ensembl_Chromosome":ens_chr,"Ensembl_Position":ens_pos,"Ensembl_ID":ens_id,"Ensembl_Gene":ens_gene,"api_validated":bool(same),"validation_status":"MATCH" if same else ("API_HIT_CONFLICT" if api_valid else "API_UNAVAILABLE"),"api_error":api_error})
    out=pd.DataFrame(rows); write_table(out,output_path); return out

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("blast"); ap.add_argument("output"); ap.add_argument("--species",default="triticum_aestivum"); ap.add_argument("--min-identity",type=float,default=90); ap.add_argument("--min-coverage",type=float,default=70); a=ap.parse_args(); run(a.blast,a.output,a.species,a.min_identity,a.min_coverage)
