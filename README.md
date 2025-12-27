# Overview

This pipeline performs downstream RNA-seq analysis starting from precomputed NCBI GEO count tables (GSE datasets).

The pipeline automates the following steps:

1. Downloading supplementary RNA-seq count data from NCBI GEO using a GSE accession ID
2. Extracting and parsing compressed count tables
3. Inferring sample metadata (condition and replicate) from file names
4. Structuring expression data and metadata for downstream analysis
5. Generating quality control visualization plots.
6. Running differential expression analysis and generating volcano plots for visualization.

# Project structure

nf-geo-rnaseq-qc-de/
├── main.nf              
├── nextflow.config     
├── scripts/            
│   ├── download_counts.py
│   ├── build_matrices.py
│   ├── qc_plots.py
│   └── differential_expression.py
├── README.md

# Parameters

| Parameter | Description             | Default     |
| --------- | ----------------------- | ----------- |
| `gse_id`  | GEO Series accession ID | `GSE116239` |
| `outdir`  | Output directory        | `results`   |
