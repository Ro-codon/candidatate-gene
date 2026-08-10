"""Module 02: create separate REF and ALT FASTA records for every SNP.

If Clean_Sequence is supplied it is treated as the reference sequence and the
variant position is taken from Variant_Offset (0-based). If no offset is
provided, the sequence midpoint is used only when --allow-midpoint is set.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from pipeline_common import clean_dna, read_table, write_table, reverse_complement

SEQ_COLS=["Clean_Sequence","Sequence","sequence","Flanking_Sequence","flank_sequence","Seq"]
OFFSET_COLS=["Variant_Offset","variant_offset","SNP_Offset","Position_In_Sequence","variant_position"]

def find(row,cols,default=None):
    for c in cols:
        if c in row.index and pd.notna(row[c]) and str(row[c]).strip(): return row[c]
    return default

def mutate(seq, ref, alt, offset):
    seq=clean_dna(seq); ref=clean_dna(ref); alt=clean_dna(alt)
    if not seq or not ref or not alt: return "",""
    off=int(offset)
    if off<0 or off+len(ref)>len(seq): raise ValueError(f"Variant offset {off} outside sequence length {len(seq)}")
    observed=seq[off:off+len(ref)]
    if observed!=ref:
        raise ValueError(f"REF mismatch at offset {off}: sequence has {observed}, expected {ref}")
    return seq,seq[:off]+alt+seq[off+len(ref):]

def run(input_path, fasta_path, output_table, allow_midpoint=False):
    df=read_table(input_path).copy(); rows=[]; fasta={}
    for i,row in df.iterrows():
        sid=str(row["SNP_ID"]); ref=str(row.get("REF","")).upper(); alt=str(row.get("ALT","")).upper(); seq=find(row,SEQ_COLS,"")
        off=find(row,OFFSET_COLS,None)
        status="OK"; error=""
        try:
            if off is None:
                if not allow_midpoint: raise ValueError("Variant_Offset is missing; provide it or use --allow-midpoint")
                off=(len(clean_dna(seq))-len(ref))//2
            refseq,altseq=mutate(seq,ref,alt,int(off))
        except Exception as e:
            refseq=altseq=""; status="ERROR"; error=str(e)
        rec={"SNP_ID":sid,"REF":ref,"ALT":alt,"Variant_Offset":off,"REF_Sequence":refseq,"ALT_Sequence":altseq,"REF_ReverseComplement":reverse_complement(refseq) if refseq else "","ALT_ReverseComplement":reverse_complement(altseq) if altseq else "","sequence_status":status,"sequence_error":error}
        rows.append(rec)
        if refseq: fasta[f"{sid}_REF"]=refseq; fasta[f"{sid}_ALT"]=altseq
    out=pd.DataFrame(rows); write_table(out,output_table)
    p=Path(fasta_path); p.parent.mkdir(parents=True,exist_ok=True)
    with p.open("w") as h:
        for k,s in fasta.items(): h.write(f">{k}\n{s}\n")
    return out

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("fasta"); ap.add_argument("output_table"); ap.add_argument("--allow-midpoint",action="store_true"); a=ap.parse_args(); run(a.input,a.fasta,a.output_table,a.allow_midpoint)
