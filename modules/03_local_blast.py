"""Module 03: local BLASTn plus optional online NCBI BLAST.

Both REF and ALT records are queried. Ensembl has no public BLAST REST endpoint;
its independent API sequence/lookup/overlap validation is performed in module 04.
"""
from __future__ import annotations
import argparse, shutil, subprocess, tempfile
from pathlib import Path
import pandas as pd
from pipeline_common import read_table, write_table
BLAST_COLS=["SNP_ID","Allele","query_id","subject_id","Chromosome","Position","Identity","Alignment_Length","Mismatches","Gap_Openings","QStart","QEnd","SStart","SEnd","Evalue","Bitscore","Coverage","Database_Source","blast_status","blast_error"]

def empty(source="",status="NOT_RUN",error=""):
    d=pd.DataFrame(columns=BLAST_COLS); d["Database_Source"]=source; d["blast_status"]=status; d["blast_error"]=error; return d

def db_exists(db): return bool(db) and any(Path(str(db)+ext).exists() for ext in (".nhr",".nin",".nsq",".ndb"))

def run_one(fasta,db,source,max_targets=10):
    if not db_exists(db): return empty(source,"NOT_RUN","BLAST database unavailable")
    if not shutil.which("blastn"): return empty(source,"NOT_RUN","blastn executable not found in PATH")
    with tempfile.TemporaryDirectory() as td:
        out=Path(td)/"blast.tsv"; cmd=["blastn","-query",str(fasta),"-db",db,"-out",str(out),"-max_target_seqs",str(max_targets),"-evalue","1e-10","-outfmt","6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"]
        try: subprocess.run(cmd,check=True,capture_output=True,text=True)
        except Exception as e:return empty(source,"ERROR",str(e))
        if not out.exists() or not out.read_text().strip():return empty(source,"NO_HIT","")
        cols=["query_id","subject_id","Identity","Alignment_Length","Mismatches","Gap_Openings","QStart","QEnd","SStart","SEnd","Evalue","Bitscore"]; d=pd.read_csv(out,sep="\t",names=cols)
        d["SNP_ID"]=d.query_id.str.rsplit("_",n=1).str[0]; d["Allele"]=d.query_id.str.rsplit("_",n=1).str[-1]; qlen=(d.QEnd-d.QStart).abs()+1; d["Coverage"]=(d.Alignment_Length/qlen*100).clip(0,100); d["Database_Source"]=source; d["blast_status"]="HIT"; d["blast_error"]=""
        d["Chromosome"]=d.subject_id.astype(str).str.extract(r"(?:chr|chromosome|Chr)[^0-9A-Za-z]*(\\w+)",expand=False).fillna(""); d["Position"]=pd.NA
        return d.reindex(columns=BLAST_COLS)

def online_ncbi(fasta,max_targets=5):
    try:
        from Bio.Blast import NCBIWWW, NCBIXML
        from Bio import SeqIO
    except Exception as e:return empty("API_NCBI","NOT_RUN",f"Biopython unavailable: {e}")
    rows=[]
    for rec in SeqIO.parse(fasta,"fasta"):
        try:
            handle=NCBIWWW.qblast("blastn","nt",str(rec.seq),hitlist_size=max_targets,expect=1e-10,format_type="XML"); result=NCBIXML.read(handle)
            for aln in result.alignments[:max_targets]:
                for hsp in aln.hsps[:1]:
                    qlen=len(rec.seq); cov=hsp.align_length/qlen*100
                    rows.append({"SNP_ID":rec.id.rsplit("_",1)[0],"Allele":rec.id.rsplit("_",1)[-1],"query_id":rec.id,"subject_id":aln.accession or aln.title,"Chromosome":"","Position":pd.NA,"Identity":hsp.identities/hsp.align_length*100,"Alignment_Length":hsp.align_length,"Mismatches":hsp.align_length-hsp.identities,"Gap_Openings":hsp.gaps,"QStart":hsp.query_start,"QEnd":hsp.query_end,"SStart":hsp.sbjct_start,"SEnd":hsp.sbjct_end,"Evalue":hsp.expect,"Bitscore":hsp.bits,"Coverage":cov,"Database_Source":"API_NCBI","blast_status":"HIT","blast_error":""})
        except Exception as e: rows.append({"SNP_ID":rec.id.rsplit("_",1)[0],"Allele":rec.id.rsplit("_",1)[-1],"query_id":rec.id,"Database_Source":"API_NCBI","blast_status":"ERROR","blast_error":str(e)})
    return pd.DataFrame(rows).reindex(columns=BLAST_COLS) if rows else empty("API_NCBI","NO_HIT","")

def run(fasta_path,output_path,ncbi_db="",ensembl_db="",online_ncbi=False):
    frames=[run_one(fasta_path,ncbi_db,"LOCAL_NCBI"),run_one(fasta_path,ensembl_db,"LOCAL_ENSEMBL")]
    if online_ncbi: frames.append(online_ncbi_fn(fasta_path))
    out=pd.concat(frames,ignore_index=True).reindex(columns=BLAST_COLS); write_table(out,output_path); return out

def online_ncbi_fn(fasta_path): return online_ncbi(fasta_path)

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("fasta"); ap.add_argument("output"); ap.add_argument("--ncbi-db",default=""); ap.add_argument("--ensembl-db",default=""); ap.add_argument("--online-ncbi",action="store_true"); a=ap.parse_args(); run(a.fasta,a.output,a.ncbi_db,a.ensembl_db,a.online_ncbi)
