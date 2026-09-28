# ARDS cell-state signature — reproducibility package

This repository contains the analysis code for the manuscript:

**Cell-state scoring outperforms conventional biomarkers for ARDS diagnosis: a multi-cohort and multi-omics validation study**

## Overview

A cell-state–based diagnostic signature was developed using an independent single-cell reference (GSE216009) to define ssGSEA gene sets for immature neutrophils and non-classical monocytes. The model was trained on GSE32707 (18 ARDS vs 34 untreated controls) and validated in three independent cohorts (GSE243066, GSE200847, GSE171524) and spatial transcriptomic data (GSE312053).

## Data sources

All raw data are publicly available from GEO. See `data_processed/ledger.csv` for the training cohort sample ledger.

| Dataset | Platform | Sample type | Use |
|:---|:---|:---|:---|
| GSE216009 | 10x Genomics | whole blood | Single-cell reference |
| GSE32707 | Illumina HT-12 v4 | whole blood | Training (18 ARDS vs 34 untreated) |
| GSE243066 | DNBSEQ-G400 | whole blood | Validation (fixed model) |
| GSE200847 | Illumina NovaSeq | tracheal aspirate | Validation (fixed model) |
| GSE171524 | 10x Genomics | lung tissue | Validation (fixed model) |
| GSE312053 | Visium HD | lung tissue | Spatial localization |

## Requirements

- R 4.3.0
- GSVA, pROC, dplyr
- See `sessionInfo.txt`

## How to run

```r
setwd("path/to/ARDS_repo")
source("scripts/run_all.R")
