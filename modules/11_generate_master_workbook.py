"""Module 11: generate the final single Excel workbook and preserve all evidence sheets."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

def load(p):
    q=Path(p)
    if not q.exists(): return pd.DataFrame()
    return pd.read_excel(q) if q.suffix.lower()==".xlsx" else pd.read_csv(q)

def merge_on_snp(base,other):
    if other.empty or "SNP_ID" not in other: return base
    x=other.drop_duplicates("SNP_ID")
    cols=[c for c in x.columns if c!="SNP_ID" and c not in base.columns]
    return base.merge(x[["SNP_ID"]+cols],on="SNP_ID",how="left")

def run(master_input,output_path,blast="",validation="",local_neighborhood="",api_neighborhood="",annotation="",literature="",scoring=""):
    base=load(master_input)
    for p in [blast,validation,local_neighborhood,api_neighborhood,annotation,literature,scoring]: base=merge_on_snp(base,load(p)) if p else base
    out=Path(output_path); out.parent.mkdir(parents=True,exist_ok=True)
    with pd.ExcelWriter(out,engine="openpyxl") as w:
        base.to_excel(w,index=False,sheet_name="Master_SNP_Results")
        for name,p in [("BLAST",blast),("Coordinate_Validation",validation),("Local_Neighborhood",local_neighborhood),("API_Neighborhood",api_neighborhood),("Functional_Annotation",annotation),("PubMed",literature),("Scoring",scoring)]:
            if p:
                d=load(p); d.to_excel(w,index=False,sheet_name=name[:31])
        pd.DataFrame({"Field":["Pipeline","Validation rule","Scoring"],"Value":["Corrected SNP candidate-gene pipeline","API validation is True only when an independent Ensembl API hit agrees with the resolved local coordinate.","Composite score is a prioritization heuristic, not a biological proof."]}).to_excel(w,index=False,sheet_name="README")
    return base
if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("master_input"); ap.add_argument("output");
    for x in ["blast","validation","local_neighborhood","api_neighborhood","annotation","literature","scoring"]: ap.add_argument("--"+x.replace("_","-"),default="")
    a=ap.parse_args(); run(a.master_input,a.output,a.blast,a.validation,a.local_neighborhood,a.api_neighborhood,a.annotation,a.literature,a.scoring)
