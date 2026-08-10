"""Module 03: local BLASTn against configured NCBI and Ensembl databases.

Both REF and ALT records are queried. The output schema is stable even when
BLAST is unavailable, preventing downstream KeyErrors.
"""
from __future__ import annotations
import argparse, shutil, subprocess, tempfile
from pathlib import Path
import pandas as pd
from pipeline_common import read_table, write_table

BLAST_COLS=["SNP_ID","Allele","query_id","subject_id","Chromosome","Position","Identity","Alignment_Length","Mismatches","Gap_Openings","QStart","QEnd","SStart","SEnd","Evalue","Bitscore","Coverage","Database_Source","blast_status","blast_error"]


def empty(): return pd.DataFrame(columns=BLAST_COLS)

def run_one(fasta, db, source, max_targets=10):
    if not db or not shutil.which("blastn") or not Path(str(db)+".nhr").exists() and not Path(str(db)+".nin").exists() and not Path(str(db)+".ndb").exists():
        return empty().assign(Database_Source=source,blast_status="NOT_RUN",blast_error="BLAST database or blastn executable unavailable")
    with tempfile.TemporaryDirectory() as td:
        out=Path(td)/"blast.tsv"
        cmd=["blastn","-query",str(fasta),"-db",db,"-out",str(out),"-max_target_seqs",str(max_targets),"-evalue","1e-10","-outfmt","6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"]
        try: subprocess.run(cmd,check=True,capture_output=True,text=True)
        except Exception as e: return empty().assign(Database_Source=source,blast_status="ERROR",blast_error=str(e))
        if not out.exists() or not out.read_text().strip(): return empty().assign(Database_Source=source,blast_status="NO_HIT",blast_error="")
        cols=["query_id","subject_id","Identity","Alignment_Length","Mismatches","Gap_Openings","QStart","QEnd","SStart","SEnd","Evalue","Bitscore"]
        d=pd.read_csv(out,sep="\t",names=cols)
        d["SNP_ID"]=d.query_id.str.rsplit("_",n=1).str[0]; d["Allele"]=d.query_id.str.rsplit("_",n=1).str[-1]
        d["Coverage"]=(d["Alignment_Length"]/(d["QEnd"]-d["QStart"].abs()+1)*100).clip(0,100)
        d["Database_Source"]=source; d["blast_status"]="HIT"; d["blast_error"]=""
        d["Chromosome"]=d["subject_id"].astype(str).str.extract(r"(?:chr|chromosome|Chr)([^|:_\s]+)",expand=False).fillna("")
        d["Position"]=pd.NA
        return d.reindex(columns=BLAST_COLS)

def run(fasta_path, output_path, ncbi_db="", ensembl_db=""):
    frames=[run_one(fasta_path,ncbi_db,"LOCAL_NCBI"),run_one(fasta_path,ensembl_db,"LOCAL_ENSEMBL")]
    out=pd.concat(frames,ignore_index=True) if frames else empty(); out=out.reindex(columns=BLAST_COLS); write_table(out,output_path); return out

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("fasta"); ap.add_argument("output"); ap.add_argument("--ncbi-db",default=""); ap.add_argument("--ensembl-db",default=""); a=ap.parse_args(); run(a.fasta,a.output,a.ncbi_db,a.ensembl_db)
