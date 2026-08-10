# Corrected SNP -> Candidate Gene Pipeline

Branch: `corrected-pipeline`

## Modules

1. `01_prepare_snps.py` — normalize filtered SNP input, preserve original columns, validate REF/ALT.
2. `02_generate_allele_sequences.py` — generate separate `SNP_ID_REF` and `SNP_ID_ALT` FASTA records; can extract the reference window from a genome FASTA.
3. `03_local_blast.py` — local BLASTn against NCBI and Ensembl databases with stable fallback schema; optional online NCBI BLAST.
4. `04_coordinate_validation.py` — resolve best local hit and independently validate with Ensembl REST overlap/lookup.
5. `05_merge_coordinates.py` — merge coordinates back to every original row using `SNP_ID` while preserving `csv_` columns.
6. `06_local_neighborhood.py` — local GFF3/GTF gene neighborhood extraction.
7. `07_api_neighborhood.py` — Ensembl Plants REST neighborhood cross-validation and `api_validated` flag.
8. `08_functional_annotation.py` — Ensembl, NCBI, UniProt, InterPro, GO, KEGG and Reactome annotation.
9. `09_pubmed_literature.py` — PubMed evidence mining.
10. `10_score_candidates.py` — reproducible 0–100 prioritization score.
11. `11_generate_master_workbook.py` — one Excel workbook plus evidence sheets.
12. `12_run_pipeline.py` — complete orchestrator.

## Important API correction

The current public Ensembl REST documentation exposes sequence, lookup, overlap and VEP endpoints, but not a public BLAST REST endpoint. Therefore this implementation does **not** fabricate an `Ensembl BLAST API` call. Ensembl API validation uses independent sequence/lookup/overlap evidence. NCBI online BLAST is available through Biopython's `NCBIWWW.qblast` and is optional because it can be slow and rate-limited.

## Run

```bash
pip install -r requirements.txt
python modules/12_run_pipeline.py \
  --input filtered_snps.csv \
  --reference wheat_reference.fa \
  --ncbi-db db/ncbi/wheat_dna_db \
  --ensembl-db db/ensembl/wheat_dna_db \
  --gff annotation.gff3 \
  --out results
```

For online NCBI BLAST add `--online-ncbi`. For explicit midpoint inference only, add `--allow-midpoint`; otherwise a missing variant offset is reported as an error rather than silently inventing the allele position.

Final file: `results/candidate_gene_master.xlsx`.
