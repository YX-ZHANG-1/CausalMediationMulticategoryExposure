# Causal Mediation Analysis via Double Machine Learning for Multi-Category Exposure Variables
---


## 📖 Overview

This repository provides R code implementing a **causal mediation analysis framework** for **multi-category exposure variables** using the **Double Machine Learning (DML)** approach.  
It includes all scripts necessary to **generate synthetic datasets**, **estimate total, direct, and indirect effects**, and **reproduce key analyses**.

---

## 📂 Repository Structure

| File | Description |
|------|--------------|
| **`generate_sample_data.R`** | Defines functions for simulating confounders, exposures, mediators, and outcomes based on realistic biomedical distributions (e.g., Beta-mixture methylation profiles). |
| **`simulate_dataset.R`** | Specifies simulation parameters and calls `generate_sample_data.R` to produce a complete dataset. |
| **`medDML_multicategory_exposure.R`** | Implements the DML estimation algorithm with cross-fitting and Lasso regularization for multi-category exposures. |
| **`conduct_causal_mediation_analysis.R`** | Conducts the causal mediation analysis using the above functions and summarizes total, direct, and indirect effects. |

All code is written in **base R**, using only **CRAN-available packages** for full reproducibility.

---

## 🧠 Scientific Context

This code accompanies the manuscript:

> **Zhang, Y.**, Yuan, Y., Xue, F., Wei, K., Zhou, J., & Qu, A. (2025).  
> *Causal Mediation Analysis of the Effect of Depression on Alzheimer’s Disease Risk in Older Adults.*  
> Submitted to *Nature Communications.*

The framework estimates:
- **Total Effect (TE)**  
- **Natural Direct Effect (NDE)**  
- **Natural Indirect Effect (NIE)**  

while adjusting for high-dimensional confounders using  **Double Machine Learning** with **Lasso-based regularization**.

---

## ⚙️ System Requirements

### 🖥️ Hardware Requirements

The R code runs on a standard computer capable of handling moderate matrix operations and cross-validation.

| Component | Recommended Specification |
|------------|----------------------------|
| **CPU** | ≥ 4 cores @ 2.5–3.5 GHz (Intel i7, AMD Ryzen 5, or Apple M2/M4) |
| **RAM** | Minimum 3 GB ; Recommended ≥ 6–8 GB |
| **Disk Space** | ≥ 500 MB free |
| **Internet Speed** | ≥ 25 Mbps for downloading R packages |

> 💡 *Example test system:* MacBook Pro (Apple M4 Max, 14 cores, 36 GB RAM)

---

### 💻 Software Requirements

| Component | Requirement |
|------------|-------------|
| **R version** | ≥ 3.6.0 (tested on R 4.4.2) |
| **Operating Systems Tested** | macOS 15.6.1 (Sequoia, ARM64), Windows 11 (23H2, x64) |
| **Required Packages** | `truncnorm`, `glmnet`, `Matrix`, `hdm` |

All packages are platform-independent and available via **CRAN**.

## 📦 Installation Guide

### 🧰 Package Installation

Before running the R scripts, install required packages from **CRAN**:

```r
install.packages(c("truncnorm", "glmnet", "Matrix", "hdm"))
```

> ⏱️ **Example installation runtime:**  
> On a **MacBook Pro (Apple M4 Max, 14 cores, 36 GB RAM)** with a stable **AT&T Fiber Internet**  
> (**Download:** 424.7 Mbps | **Upload:** 272.8 Mbps),  
> installation of all required packages completes in approximately **1.5 seconds**.


### Package Versions Tested (as of October 2025)

| Package | Version |
|------------------|----------|
| truncnorm| 1.0.9 | 
| glmnet| 4.1.10 |
| Matrix| 1.7.4 |
| hdm| 0.3.2 |

## 🧪 Demo

### 1️⃣ Run the Complete Workflow

From your R terminal or RStudio console, execute:

```r
# Step 1. Set working directory to this repository 
setwd("path/to/CausalMediationMulticategoryExposure")

# Step 2. Generate a simulated dataset
source("simulate_dataset.R")

# Step 3. Conduct the causal mediation analysis 
source("conduct_causal_mediation_analysis.R")
```

### 2️⃣ Expected Output

| Output File                        | Description                                                                       |
| ---------------------------------- | --------------------------------------------------------------------------------- |
| `simulated_data.csv`               | Generated dataset containing exposure, mediators, confounders, and outcome        |
| `mediation_analysis_results.RData` | Saved results of total, direct, and indirect effects                              |
| `df_total`, `df_direct`, `df_med`  | Summary data frames (in R environment) with estimates, standard errors, 95% CIs, and p-values |


### 3️⃣ Expected Runtime 
On a **MacBook Pro (Apple M4 Max, 14 cores, 36 GB RAM)**, the complete simulation and analysis run in approximately **50 seconds**.

## 🧭 Instructions for Use

### 🔹 Run on Your Own Data

To apply the DML-based causal mediation framework to your dataset:

#### 1. Prepare your data
Your dataset should include:
- **Outcome variable (`Y`)** — continuous or binary 
- **Exposure variable (`Z`)** — multi-category (control level coded as 0)  
- **Potential mediator(s) (`M`)** — one or more continuous/binary potentail mediators
- **Potential confounders (`X`)** — potential pre-treatment covariates

Ensure no missing values and save as **your_data.csv**.

#### 2. Load your dataset
```r
data <- read.csv("your_data.csv", header = TRUE)
```
#### 3. Format variables 
Example:
```r
y <- as.matrix(data$Y)
z <- as.matrix(data$Z)
m <- as.matrix(data[, c("Mediator1", "Mediator2")])
x <- as.matrix(data[, c("Confounder1", "Confouer2")])
```
#### 4. Run the analysis
Edit **conduct_causal_mediation_analysis.R** to load your dataset, then run:
```r
source("conduct_causal_mediation_analysis.R")
```

## 🪪 License

This repository is distributed under the **MIT License**.  
You are free to use, modify, and distribute this code with proper attribution.  
See the [LICENSE](LICENSE) file for full terms.

