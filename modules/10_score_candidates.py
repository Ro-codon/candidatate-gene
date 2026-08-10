"""Module 10: reproducible composite SNP/candidate-gene scoring (0-100)."""
from __future__ import annotations
import argparse, math
import pandas as pd
from pipeline_common import read_table, write_table

def clamp(x): return max(0.0,min(1.0,float(x)))
def run(input_path,output_path):
    df=read_table(input_path).copy(); scores=[]
    for _,r in df.iterrows():
        p=pd.to_numeric(r.get("csv_GWAS_P",r.get("GWAS_P")),errors="coerce"); pscore=clamp((-math.log10(max(float(p),1e-300)))/10) if pd.notna(p) else 0
        dist=pd.to_numeric(r.get("API_Distance_bp",r.get("Distance_bp")),errors="coerce"); dscore=1/(1+float(dist)/10000) if pd.notna(dist) else 0
        api=1 if bool(r.get("api_validated",False)) else 0
        literature=clamp(math.log10(1+float(pd.to_numeric(r.get("PubMed_Count"),errors="coerce") or 0))/2)
        go=str(r.get("GO_Annotations","")).lower(); stress=1 if any(k in go for k in ["stress","heat","drought","oxidative","abiotic","biotic"]) else 0
        annotation=1 if str(r.get("annotation_status","")).upper()=="OK" else 0
        total=30*pscore+20*dscore+15*api+15*stress+10*annotation+10*literature
        scores.append((pscore,dscore,api,stress,annotation,literature,total))
    s=pd.DataFrame(scores,columns=["GWAS_Score","Proximity_Score","API_Validation_Score","Stress_GO_Score","Annotation_Score","Literature_Score","Composite_Score"])
    out=pd.concat([df.reset_index(drop=True),s],axis=1); out["Priority"]=pd.cut(out.Composite_Score,bins=[-1,40,70,101],labels=["C","B","A"]); write_table(out,output_path); return out
if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("output"); a=ap.parse_args(); run(a.input,a.output)
