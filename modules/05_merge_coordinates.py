"""Module 05: merge resolved coordinates back to every original SNP row."""
from __future__ import annotations
import argparse
import pandas as pd
from pipeline_common import read_table, write_table

def run(original_path, coordinates_path, output_path):
    original=read_table(original_path).copy(); coords=read_table(coordinates_path).copy()
    if "SNP_ID" not in original or "SNP_ID" not in coords: raise ValueError("Both tables require SNP_ID")
    original.columns=[c if str(c).startswith("csv_") or c=="SNP_ID" else f"csv_{c}" for c in original.columns]
    coords=coords.drop_duplicates("SNP_ID",keep="first")
    out=original.merge(coords,on="SNP_ID",how="left",suffixes=("","_coord"))
    write_table(out,output_path); return out
if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("original"); ap.add_argument("coordinates"); ap.add_argument("output"); a=ap.parse_args(); run(a.original,a.coordinates,a.output)
