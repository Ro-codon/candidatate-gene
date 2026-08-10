"""Module 08: functional annotation from Ensembl, NCBI, UniProt, InterPro, KEGG and Reactome.

Unavailable services are explicitly flagged; no failed request is treated as a
biological negative result.
"""
from __future__ import annotations
import argparse, urllib.parse
import pandas as pd
from pipeline_common import cached_json, request_json, request_text, read_table, write_table


def ensembl(gene_id,species):
    if not gene_id: return {}
    d,_,_=request_json(f"https://rest.ensembl.org/lookup/id/{urllib.parse.quote(gene_id)}",params={"expand":1},headers={"Content-Type":"application/json","Accept":"application/json"})
    return d if isinstance(d,dict) else {}

def ncbi(gene_name, email=""):
    if not gene_name: return {}
    params={"db":"gene","term":f"{gene_name}[Gene Name] AND Triticum aestivum[Organism]","retmode":"json"}
    if email: params["email"]=email
    d,_,_=request_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",params=params)
    ids=(d or {}).get("esearchresult",{}).get("idlist",[]) if isinstance(d,dict) else []
    if not ids: return {}
    s,_,_=request_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",params={"db":"gene","id":ids[0],"retmode":"json"})
    return (s or {}).get("result",{}).get(ids[0],{}) if isinstance(s,dict) else {}

def uniprot(query):
    if not query: return {}
    d,_,_=request_json("https://rest.uniprot.org/uniprotkb/search",params={"query":query,"format":"json","size":1})
    r=(d or {}).get("results",[]) if isinstance(d,dict) else []
    return r[0] if r else {}

def interpro(gene_id):
    if not gene_id:return []
    d,_,_=request_json(f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{urllib.parse.quote(gene_id)}/",params={"page_size":100})
    return d.get("results",[]) if isinstance(d,dict) else []

def kegg(gene_name):
    if not gene_name:return ""
    text,_,_=request_text(f"https://rest.kegg.jp/find/genes/{urllib.parse.quote(gene_name)}")
    return text.strip()

def reactome(gene_id):
    if not gene_id:return []
    d,_,_=request_json(f"https://reactome.org/ContentService/data/mapping/UniProt/{urllib.parse.quote(gene_id)}/pathways")
    return d if isinstance(d,list) else []

def run(input_path,output_path,species="triticum_aestivum",ncbi_email=""):
    df=read_table(input_path); rows=[]
    for _,r in df.iterrows():
        eid=str(r.get("API_Gene_ID") or r.get("Ensembl_ID") or ""); en=ensembl(eid,species); gene=str(r.get("API_Gene_Name") or r.get("Local_Gene_Name") or en.get("display_name") or ""); n=ncbi(gene,ncbi_email)
        uni=uniprot(gene); uniacc=uni.get("primaryAccession","") if uni else ""; ip=interpro(uniacc) if uniacc else []; kp=kegg(gene); rp=reactome(uniacc) if uniacc else []
        gos=[]
        for x in en.get("GO",[]) if isinstance(en,dict) else []:
            if x.get("GO_id") or x.get("description"): gos.append(f"{x.get('GO_id','')}|{x.get('description','')}|{x.get('namespace','')}")
        rows.append({"SNP_ID":r.SNP_ID,"Gene_ID":eid,"Gene_Name":gene,"Ensembl_Biotype":en.get("biotype",""),"Ensembl_Description":en.get("description",""),"Transcript_ID":en.get("canonical_transcript",""),"NCBI_GeneID":n.get("uid","") if n else "","NCBI_Description":n.get("description","") if n else "","NCBI_Summary":n.get("summary","") if n else "","UniProt_Accession":uniacc,"UniProt_Protein_Name":(uni.get("proteinDescription",{}).get("recommendedName",{}).get("fullName",{}).get("value","") if uni else ""),"UniProt_Function":(uni.get("comments",[{}])[0].get("texts",[{}])[0].get("value","") if uni and uni.get("comments") else ""),"GO_Annotations":"; ".join(gos),"InterPro_Entries":"; ".join(str(x.get("metadata",{}).get("accession",x.get("metadata",{}).get("name",""))) for x in ip),"KEGG_Search":kp,"Reactome_Pathways":"; ".join(str(x.get("displayName",x.get("stId",""))) for x in rp),"annotation_status":"OK" if any([en,n,uni,ip,kp,rp]) else "NO_ANNOTATION"})
    out=pd.DataFrame(rows); write_table(out,output_path); return out
if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("output"); ap.add_argument("--species",default="triticum_aestivum"); ap.add_argument("--ncbi-email",default=""); a=ap.parse_args(); run(a.input,a.output,a.species,a.ncbi_email)
