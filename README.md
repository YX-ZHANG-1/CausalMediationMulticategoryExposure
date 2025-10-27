# Causal Mediation Analysis via Double Machine Learning for Multi-Category Exposure
<p style="font-size:18px;">
**License:** MIT License  
**Repository:** [https://github.com/YX-ZHANG-1/CausalMediationMulticategoryExposure](https://github.com/YX-ZHANG-1/CausalMediationMulticategoryExposure)

---

## 🧩 Overview

This repository provides reproducible R code for performing **causal mediation analysis** under **multi-category exposures** using the **Double Machine Learning (DML)** framework.  

The scripts generate realistic synthetic datasets, estimate direct and indirect causal pathways, and produce publication-ready results consistent with *Nature Communications* reporting standards.

The repository contains four main scripts:

| Script | Description |
|---------|--------------|
| **`generate_sample_data.R`** | Generates synthetic data mimicking real-world biomedical data (e.g., ADNI) for exposure (`Z`), mediators (`M`), confounders (`X`), and outcome (`Y`). |
| **`simulate_dataset.R`** | Simulates datasets under user-specified parameters, integrating multiple variable types and distributions. |
| **`medDML_multicategory_exposure.R`** | Implements the Double Machine Learning estimator for total, direct, and indirect effects under multi-category exposure variables. |
| **`conduct_causal_mediation_analysis.R`** | Main driver script that loads data, runs mediation analysis, summarizes results, and saves outputs for reproducibility. |

---

## 🧠 Scientific Context

This code supports the manuscript:

> Zhang, Y. *et al.* (2025). **Causal Mediation Analysis with Multi-category Exposures via Double Machine Learning.** *Nature Communications.*

The method estimates:
- **Total Effect (TE)**  
- **Natural Direct Effect (NDE)**  
- **Natural Indirect Effect (NIE)**  

while adjusting for high-dimensional confounders using regularized regression (Lasso / GLMNet).

---

## ⚙️ System Requirements

| Requirement | Recommended |
|--------------|-------------|
| **R version** | ≥ 4.2.0 |
| **Operating System** | macOS ≥ 12, Linux ≥ Ubuntu 20.04, or Windows ≥ 10 |
| **Memory** | ≥ 8 GB RAM |
| **R Packages** | `dplyr`, `VGAM`, `glmnet`, `hdm`, `Matrix` |
| **Optional** | `ggplot2`, `knitr` for visualization or reproducible notebooks |

Install dependencies:
```r
install.packages(c("dplyr", "VGAM", "glmnet", "hdm", "Matrix"))
