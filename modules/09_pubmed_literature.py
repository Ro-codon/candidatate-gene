"""Module 09: PubMed literature mining for candidate genes."""
from __future__ import annotations
import argparse, xml.etree.ElementTree as ET
import pandas as pd
from pipeline_common import request_json, request_text, read_table, write_table

def run(input_path,output_path,email="",max_articles=20):
    df=read_table(input_path); rows=[]
    for _,r in df.drop_duplicates("SNP_ID").iterrows():
        gene=str(r.get("Gene_Name") or r.get("API_Gene_Name") or r.get("Local_Gene_Name") or "");
        if not gene: rows.append({"SNP_ID":r.SNP_ID,"PubMed_Count":0,"PubMed_Status":"NO_GENE"}); continue
        params={"db":"pubmed","term":f'"{gene}" AND (wheat OR Triticum) AND (stress OR yield OR heat OR drought OR disease)',"retmode":"json","retmax":max_articles}
        if email: params["email"]=email
        d,_,err=request_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",params=params); ids=(d or {}).get("esearchresult",{}).get("idlist",[]) if isinstance(d,dict) else []
        titles=[]; years=[]
        if ids:
            text,_,_=request_text("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",params={"db":"pubmed","id":",".join(ids),"retmode":"xml"})
            try:
                root=ET.fromstring(text)
                for art in root.findall(".//PubmedArticle"):
                    title="".join(art.findtext(".//ArticleTitle",default="")).strip(); year=art.findtext(".//PubDate/Year") or art.findtext(".//PubDate/MedlineDate",default="")
                    if title: titles.append(title); years.append(year)
            except ET.ParseError: pass
        rows.append({"SNP_ID":r.SNP_ID,"Gene_Name":gene,"PubMed_Count":len(ids),"PubMed_PMIDs":";".join(ids),"PubMed_Titles":" || ".join(titles),"PubMed_Years":";".join(years),"PubMed_Status":"OK" if ids else ("ERROR" if err else "NO_HITS")})
    out=pd.DataFrame(rows); write_table(out,output_path); return out
if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("output"); ap.add_argument("--email",default=""); ap.add_argument("--max-articles",type=int,default=20); a=ap.parse_args(); run(a.input,a.output,a.email,a.max_articles)
