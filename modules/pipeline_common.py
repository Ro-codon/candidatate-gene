from __future__ import annotations

import json, logging, os, re, time
from pathlib import Path
from typing import Any, Iterable
import pandas as pd
import requests

LOG = logging.getLogger("candidate_gene_pipeline")
DNA = set("ACGTRYSWKMBDHVN")


def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(level=getattr(logging, level.upper(), logging.INFO), format="%(asctime)s | %(levelname)s | %(message)s")


def clean_dna(seq: Any) -> str:
    s = re.sub(r"\s+", "", str(seq or "")).upper()
    return "".join(c for c in s if c in DNA)


def reverse_complement(seq: str) -> str:
    return clean_dna(seq).translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


def first_value(row: pd.Series, names: Iterable[str], default=""):
    for n in names:
        if n in row.index and pd.notna(row[n]) and str(row[n]).strip(): return row[n]
    return default


def normalize_snp_id(x: Any, fallback: int | None = None) -> str:
    s = str(x or "").strip()
    if not s or s.lower() in {"nan", "none", "na", "null"}:
        return f"SNP_{fallback}" if fallback is not None else "SNP_UNKNOWN"
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", s)


def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() in {".xlsx", ".xls"}: return pd.read_excel(p)
    return pd.read_csv(p, sep="\t" if p.suffix.lower() in {".tsv", ".txt"} else ",")


def write_table(df: pd.DataFrame, path: str) -> None:
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    if p.suffix.lower() == ".xlsx": df.to_excel(p, index=False)
    else: df.to_csv(p, index=False)


def request_json(url: str, params=None, headers=None, timeout=60, retries=3, sleep=1.5) -> tuple[Any, int | None, str]:
    last = ""
    for i in range(retries):
        try:
            r = requests.get(url, params=params, headers=headers or {"Accept":"application/json"}, timeout=timeout)
            last = r.text[:500]
            if r.ok: return r.json(), r.status_code, ""
            if r.status_code in {429, 500, 502, 503, 504}: time.sleep(sleep * (i + 1)); continue
            return None, r.status_code, last
        except Exception as e:
            last = str(e); time.sleep(sleep * (i + 1))
    return None, None, last


def request_text(url: str, params=None, headers=None, timeout=120, retries=3) -> tuple[str, int | None, str]:
    last = ""
    for i in range(retries):
        try:
            r = requests.get(url, params=params, headers=headers or {}, timeout=timeout)
            last = r.text[:500]
            if r.ok: return r.text, r.status_code, ""
            if r.status_code in {429,500,502,503,504}: time.sleep(1.5*(i+1)); continue
            return "", r.status_code, last
        except Exception as e:
            last = str(e); time.sleep(1.5*(i+1))
    return "", None, last


def post_json(url: str, payload: dict, headers=None, timeout=120, retries=3) -> tuple[Any, int | None, str]:
    last=""
    for i in range(retries):
        try:
            r=requests.post(url,json=payload,headers=headers or {"Accept":"application/json"},timeout=timeout)
            last=r.text[:500]
            if r.ok:
                try: return r.json(),r.status_code,""
                except Exception: return r.text,r.status_code,""
            if r.status_code in {429,500,502,503,504}: time.sleep(1.5*(i+1)); continue
            return None,r.status_code,last
        except Exception as e: last=str(e); time.sleep(1.5*(i+1))
    return None,None,last


def cached_json(path: str, key: str, fetcher):
    p=Path(path); p.parent.mkdir(parents=True,exist_ok=True)
    cache={}
    if p.exists():
        try: cache=json.loads(p.read_text())
        except Exception: cache={}
    if key in cache: return cache[key]
    value=fetcher(); cache[key]=value
    p.write_text(json.dumps(cache,indent=2,default=str))
    return value


def excel_safe(x: Any, max_len=32000):
    s="" if x is None else str(x)
    return s[:max_len]


def ensure_columns(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    for c in cols:
        if c not in df.columns: df[c] = pd.NA
    return df
