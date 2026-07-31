# Single-nucleus RNA-seq analysis of Major Depressive Disorder

## Overview

This repository contains the scripts used for Mendelian randomization, colocalization , snRNA-seq preprocessing, quality control, cell annotation, differential expression analysis and downstream analyses in our study:

The workflow integrates publicly available snRNA-seq datasets and performs donor-aware differential expression analysis together with genetic causal inference.

---

## Repository structure
```
Project/

├── scripts/
│   ├── SMR.R
│   ├── TSMR.R
│   └── MDD_snRNA_analysis.R
│
├── data/
│   ├── raw/   ###original files
│   └── processed/  ###processed files
│
├── results/
│   ├── SMR/
│   ├── TSMR/
│   ├── colocalization/
│   ├── QC/
│   ├── CellAnnotation/
│   ├── DEanalysis/
│   │      ├── edgeR/
│   │      └── MASTRE/
│   ├── singlegeneGSEA/
│   ├── COQ8AEx/
│   └── DAanalysis/
│
└── README.md
```

---

## Software requirements

Analysis was performed using

- SMR v1.3.1
- R v4.4.3
- TwoSampleMR v0.5.10
- coloc v5.2.3
- locuscomparer v1.0.0
- Seurat v5.1.0 
- harmony v1.2.3
- decontX v1.4.1
- SCTransform v0.4.2
- MAST v1.32.0
- edgeR v4.4.2
- ggvenn v0.1.10
- ggplot2 v3.5.1
- clusterProfiler package v4.12.6
- enrichplot v1.24.4

Additional packages are listed in each script.

---

## Input datasets

### Single-cell RNA-seq

Public datasets

- GSE213982
- GSE144136

### eQTL

Bryois et al. human brain single-cell eQTL

### GWAS

PGC Major Depressive Disorder GWAS (without 23andMe)

---

## Analysis workflow

### Step 1. SMR

Run

```
01.SMR.R
```

Output

```
results/SMR/
```

---

### Step 2. TSMR and colocalization

Run

```
02.TSMR_coloc.R
```

Output

```
results/TSMR/
results/colocalization/
```

---

### Step 3. snRNA-seq

Run

```
3.MDD_snRNA_analysis.R
```

Output

```
results/QC/
results/CellAnnotation/
results/DEanalysis/
results/DEanalysis/edgeR/
results/DEanalysis/MASTRE/
results/singlegeneGSEA/
results/COQ8AEx/
results/DAanalysis/
```

---

## Reproducibility

All scripts use relative paths.

Users only need to define

```r
project_dir <- "."
```

at the beginning of each script.

Intermediate processed objects used in the manuscript are available at

10.5281/zenodo.21722584

---

## Citation

If you use these scripts, please cite

Cheng-Long Dong et al.
