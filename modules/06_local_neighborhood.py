"""Module 06: extract local gene neighborhood from GFF3/GTF."""
from __future__ import annotations
import argparse
import gzip
import pandas as pd
from pipeline_common import read_table, write_table

def opener(p): return gzip.open(p,"rt") if str(p).endswith(".gz") else open(p,"rt",encoding="utf8",errors="replace")
def attrs(s):
    d={}
    for x in s.strip().split(";"):
        if not x.strip(): continue
        if "=" in x: k,v=x.split("=",1)
        elif " " in x: k,v=x.split(None,1)
        else: continue
        d[k.strip()]=v.strip().strip('"')
    return d

def run(coords_path,gff_path,output_path,window=100000):
    c=read_table(coords_path); genes=[]
    with opener(gff_path) as h:
        for line in h:
            if line.startswith("#"): continue
            f=line.rstrip("\n").split("\t")
            if len(f)<9 or f[2].lower()!="gene": continue
            try: st,en=int(f[3]),int(f[4])
            except: continue
            a=attrs(f[8]); genes.append({"chrom":f[0],"start":st,"end":en,"strand":f[6],"gene_id":a.get("ID") or a.get("gene_id") or a.get("Name") or "","gene_name":a.get("Name") or a.get("gene_name") or a.get("gene") or ""})
    g=pd.DataFrame(genes); rows=[]
    for _,r in c.iterrows():
        chrom=str(r.get("Ensembl_Chromosome") or r.get("Local_Chromosome") or ""); pos=pd.to_numeric(r.get("Ensembl_Position") if pd.notna(r.get("Ensembl_Position")) else r.get("Local_Position"),errors="coerce")
        if not chrom or pd.isna(pos) or g.empty: rows.append({"SNP_ID":r.SNP_ID,"local_neighborhood_status":"NO_COORDINATE"}); continue
        x=g[g.chrom.astype(str)==chrom].copy(); x["distance"]=(x["start"]-int(pos)).abs().where(pos<x["start"],(int(pos)-x["end"]).abs().where(pos>x["end"],0)); x=x[x["distance"]<=window].sort_values("distance")
        for rank,(_,q) in enumerate(x.head(20).iterrows(),1): rows.append({"SNP_ID":r.SNP_ID,"Neighborhood_Rank":rank,"Local_Gene_ID":q.gene_id,"Local_Gene_Name":q.gene_name,"Local_Gene_Chromosome":q.chrom,"Local_Gene_Start":q.start,"Local_Gene_End":q.end,"Local_Gene_Strand":q.strand,"Distance_bp":int(q.distance),"local_neighborhood_status":"HIT"})
        if x.empty: rows.append({"SNP_ID":r.SNP_ID,"local_neighborhood_status":"NO_LOCAL_GENE"})
    out=pd.DataFrame(rows); write_table(out,output_path); return out
if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("coordinates"); ap.add_argument("gff"); ap.add_argument("output"); ap.add_argument("--window",type=int,default=100000); a=ap.parse_args(); run(a.coordinates,a.gff,a.output,a.window)
