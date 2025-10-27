# Causal Mediation Analysis via Double Machine Learning for Multi-Category Exposure
<p style="font-size:18px;">
**License:** MIT License  
**Repository:** [https://github.com/YX-ZHANG-1/CausalMediationMulticategoryExposure](https://github.com/YX-ZHANG-1/CausalMediationMulticategoryExposure)

---

## 🧩 Overview

This repository provides reproducible R code for performing **causal mediation analysis** under **multi-category exposures** using the **Double Machine Learning (DML)** framework.  

The repository contains four main R scripts:

| R Script | Description |
|---------|--------------|
| **`generate_sample_data.R`** | Generates synthetic data mimicking real-world data (e.g., ADNI) for exposure (`Z`), potential mediators (`M`), potential confounders (`X`), and outcome (`Y`). |
| **`simulate_dataset.R`** | Simulates datasets under user-specified parameters, integrating multiple variable types and distributions. |
| **`medDML_multicategory_exposure.R`** | Implements the Double Machine Learning estimator for total, direct, and indirect effects under multi-category exposure variables. |
| **`conduct_causal_mediation_analysis.R`** | Main driver script that loads data, runs mediation analysis, summarizes results, and saves outputs for reproducibility. |

---

## 🧠 Scientific Context

This code supports the manuscript:

> Yuexia Zhang, Yubai Yuan, Fei Xue, Kecheng Wei, Jin Zhou, and Annie Qu (2025+). **Causal Mediation Analysis  of the Effect of Depression on Alzheimer's Disease Risk in Older Adults.** 

The method estimates:
- **Total Effect (TE)**  
- **Natural Direct Effect (NDE)**  
- **Natural Indirect Effect (NIE)**  

while adjusting for high-dimensional confounders using penalized regression (e.g., Lasso).

---

## ⚙️ System Requirements

| Requirement | Recommended |
|--------------|-------------|
| **R version** | ≥ 3.6.0 |
| **Operating System** | macOS ≥ 12, Linux ≥ Ubuntu 20.04, or Windows ≥ 10 |
| **Memory** | ≥ 8 GB RAM |
| **R Packages** | `truncnorm`, `glmnet`, `hdm`, `Matrix` |

Install dependencies:
```r
install.packages(c("truncnorm", "glmnet", "hdm", "Matrix"))
