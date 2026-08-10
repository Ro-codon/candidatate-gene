"""Module 02: generate separate REF and ALT sequence records for every SNP."""
from __future__ import annotations
import argparse, gzip
from pathlib import Path
import pandas as pd
from pipeline_common import clean_dna, read_table, write_table, reverse_complement, first_value
SEQ_COLS=["Clean_Sequence","Sequence","sequence","Flanking_Sequence","flank_sequence","Seq"]
OFFSET_COLS=["Variant_Offset","variant_offset","SNP_Offset","Position_In_Sequence","variant_position"]
CHR_COLS=["Chromosome","Chr","chromosome","CHROM","Local_Chromosome","Ensembl_Chromosome"]
POS_COLS=["Position","POS","position","BP","bp","Local_Position","Ensembl_Position"]

def read_fasta(path):
    rec={}; cur=None; chunks=[]; opener=gzip.open if str(path).endswith(".gz") else open
    with opener(path,"rt",encoding="utf8",errors="replace") as h:
        for line in h:
            line=line.strip()
            if line.startswith(">"):
                if cur: rec[cur]="".join(chunks)
                cur=line[1:].split()[0]; chunks=[]
            elif cur: chunks.append(line)
        if cur: rec[cur]="".join(chunks)
    return {k:clean_dna(v) for k,v in rec.items()}

def get_contig(refs,chrom):
    c=str(chrom); candidates=[c,c.replace("chr","",1) if c.lower().startswith("chr") else "chr"+c]
    for x in candidates:
        if x in refs:return x
    return next((k for k in refs if k.split("|")[0]==c),None)

def mutate(seq,ref,alt,offset):
    seq,ref,alt=clean_dna(seq),clean_dna(ref),clean_dna(alt); off=int(offset)
    if not seq or not ref or not alt: raise ValueError("Missing sequence or allele")
    if off<0 or off+len(ref)>len(seq): raise ValueError("Variant offset outside sequence")
    observed=seq[off:off+len(ref)]
    if observed!=ref: raise ValueError(f"REF mismatch: sequence has {observed}, expected {ref}")
    return seq,seq[:off]+alt+seq[off+len(ref):]

def run(input_path,fasta_path,output_table,allow_midpoint=False,reference_path="",flank=250):
    df=read_table(input_path).copy(); refs=read_fasta(reference_path) if reference_path else {}; rows=[]; fasta={}
    for _,row in df.iterrows():
        sid=str(row["SNP_ID"]); ref=clean_dna(row.get("REF","")); alt=clean_dna(row.get("ALT","")); seq=clean_dna(first_value(row,SEQ_COLS,"")); off=first_value(row,OFFSET_COLS,None); source="Clean_Sequence" if seq else "Reference_FASTA"; status="OK"; error=""; chrom=first_value(row,CHR_COLS,""); pos=first_value(row,POS_COLS,None)
        try:
            if not seq and refs:
                contig=get_contig(refs,chrom)
                if not contig or pos is None: raise ValueError("No usable chromosome/position for reference FASTA extraction")
                pos=int(float(pos)); start=max(1,pos-flank); end=min(len(refs[contig]),pos+flank); seq=refs[contig][start-1:end]; off=pos-start
            if off is None:
                if not allow_midpoint: raise ValueError("Variant_Offset missing; provide it or use --allow-midpoint")
                off=(len(seq)-len(ref))//2
            refseq,altseq=mutate(seq,ref,alt,int(off))
        except Exception as e: refseq=altseq=""; status="ERROR"; error=str(e)
        rows.append({"SNP_ID":sid,"REF":ref,"ALT":alt,"Variant_Offset":off,"Sequence_Source":source,"REF_Sequence":refseq,"ALT_Sequence":altseq,"REF_ReverseComplement":reverse_complement(refseq) if refseq else "","ALT_ReverseComplement":reverse_complement(altseq) if altseq else "","sequence_status":status,"sequence_error":error})
        if refseq: fasta[f"{sid}_REF"]=refseq; fasta[f"{sid}_ALT"]=altseq
    out=pd.DataFrame(rows); write_table(out,output_table); p=Path(fasta_path); p.parent.mkdir(parents=True,exist_ok=True)
    with p.open("w") as h:
        for k,s in fasta.items(): h.write(f">{k}\n{s}\n")
    return out
if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("fasta"); ap.add_argument("output_table"); ap.add_argument("--reference",default=""); ap.add_argument("--flank",type=int,default=250); ap.add_argument("--allow-midpoint",action="store_true"); a=ap.parse_args(); run(a.input,a.fasta,a.output_table,a.allow_midpoint,a.reference,a.flank)
