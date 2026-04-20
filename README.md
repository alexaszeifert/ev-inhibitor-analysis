# Comparative analysis of extracellular vesicles released from resting neutrophils under distinct pharmacological inhibition

Code and figures for the paper: [Szeifert et al. (2026), *Scientific Reports*](https://doi.org/10.1038/s41598-026-48570-8)

## Repository structure

```
├── figs/        # figures from the paper
└── scripts/     # scripts to reproduce each figure
```

## Requirements

R version 4.3.1 was used. The following packages are required:

```r
install.packages(c("tidyverse", "openxlsx", "pheatmap", "UpSetR"))

# from Bioconductor:
BiocManager::install("org.Hs.eg.db")
```

## Citation

```bibtex
@article{szeifert2026,
  title   = {Comparative analysis of extracellular vesicles released from resting neutrophils under distinct pharmacological inhibition},
  author  = {Szeifert, Vikt{\'o}ria and Szeifert, Alexa and Arnold, William Randall and Terasaki, Azusa and Bhatnagar, Keshav and Okwan-Duodu, Derick},
  journal = {Scientific Reports},
  year    = {2026},
  doi     = {10.1038/s41598-026-48570-8}
}
```