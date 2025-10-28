# Causal Mediation Analysis via Double Machine Learning for a Multi-Category Exposure Variable
<p style="font-size:18px;">
**License:** MIT License  
**Repository:** [https://github.com/YX-ZHANG-1/CausalMediationMulticategoryExposure](https://github.com/YX-ZHANG-1/CausalMediationMulticategoryExposure)

---


## 📖 Overview

This repository provides R code for conducting **causal mediation analysis with a multi-category exposure variable** using a **Double Machine Learning (DML)** framework.  
It contains reproducible scripts for:

- Generating simulated datasets mimicking real biomedical data (e.g., Alzheimer’s Disease Neuroimaging Initiative [ADNI]);
- Estimating **total**, **direct**, and **indirect** effects via DML with cross-fitting and Lasso regularization;

---

## 📂 Repository Structure

| File | Description |
|------|--------------|
| **`generate_sample_data.R`** | Defines functions for simulating confounders, exposure, mediators, and outcomes based on realistic biomedical distributions (including Beta-mixture DNA methylation models). |
| **`simulate_dataset.R`** | Specifies model parameters and generates the full simulated dataset using `generate_sample_data.R`. |
| **`medDML_multicategory_exposure.R`** | Implements the Double Machine Learning estimation procedure for a multi-category exposure variable. |
| **`conduct_causal_mediation_analysis.R`**| Conducts the causal mediation analysis on simulated data using `medDML_multicategory_exposure.R` and produces summaries of total, direct, and indirect effects. |

All scripts are written in **base R** and depend only on **CRAN-available packages**.

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

### 🖥️ Hardware Requirements

The scripts require only a standard computer with sufficient RAM for matrix operations and cross-validation.

| Component | Recommended Specification |
|------------|----------------------------|
| **RAM** | ≥ 6 GB (minimum 3 GB) |
| **CPU** | ≥ 4 cores, 2.5–3.5 GHz per core (Intel i7, AMD Ryzen 5, Apple M2 or higher) |
| **Disk Space** | ≥ 500 MB free |
| **Internet Connection** | ≥ 25 Mbps (for CRAN package downloads) |

---

### 💻 Software Requirements

#### Operating Systems Tested

| Operating System | Version | Architecture |
|------------------|----------|---------------|
| **macOS Sequoia** | 15.6.1 | Apple Silicon (M4 Max) |
| **Windows 11** | 23H2 | x64 |

All dependencies are CRAN-available and cross-platform compatible. Before running the R scripts, users should have R version 3.6.0 or higher. The scripts have been tested on R version 4.4.2.

## 📦 Installation Guide

### 🧰 Package Dependencies

Before running the R scripts, please ensure that all required packages are installed from **CRAN**.  
You can install them directly from an R terminal:

```r
install.packages(c(
  "truncnorm", "glmnet", "Matrix", "hdm"
))
```

> ⏱️ **Example installation runtime:**  
> On a **MacBook Pro (Apple M4 Max, 14 cores, 36 GB RAM)** with a stable **AT&T Fiber connection**  
> (**Download:** 424.7 Mbps | **Upload:** 272.8 Mbps),  
> installation of all required packages completes in approximately **1.5 seconds**.


### Package Versions

All functions in this repository have been tested with the latest CRAN releases as of October 2025.
The specific versions used during testing are:

| Package | Version |
|------------------|----------|
| truncnorm| 1.0.9 | 
| glmnet| 4.1.10 |
| Matrix| 1.7.4 |
| hdm| 0.3.2 |

## 🧪 Demo

### 1️⃣ Instructions to Run the Demo

This demo illustrates the complete workflow for **causal mediation analysis with a multi-category exposure**, using the four main R scripts provided in this repository.

From your R terminal or RStudio console, run the following commands **in order**:

```r
# Step 1. Set working directory to the repository folder
setwd("path/to/CausalMediationMulticategoryExposure")

# Step 2. Generate a simulated dataset mimicking real biomedical data
source("simulate_dataset.R")

# Step 3. Conduct the causal mediation analysis and save results
source("conduct_causal_mediation_analysis.R")
```

### 2️⃣ Expected Output
After the demo completes, the following output files will appear in your working directory:

| File                                                 | Description                                                                                        |
| ---------------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| `simulated_data.csv`                                 | Contains the generated sample dataset used for causal mediation analysis                                            |
| `mediation_analysis_results.RData`                   | Stores estimated total, direct, and indirect effects and their related information                                                |
| `df_total`, `df_direct`, `df_med` (in R environment) | Data frames summarizing estimated effects, standard errors, 95% confidence intervals, and p-values |

### 3️⃣ Expected Run Time 
On a **MacBook Pro (Apple M4 Max, 14 cores, 36 GB RAM)**, the simulated data generation and causal mediation analysis completes in approximately **50 seconds**.

## 🧭 Instructions for Use

### 🔹 How to Run the Software on Your Own Data 

The R scripts in this repository can be easily adapted to analyze your own dataset  
for **causal mediation analysis with a multi-category exposure variable**.

#### 1. Prepare your data
Your dataset should include:
- **Outcome variable (`Y`)** — continuous or binary outcome.  
- **Exposure variable (`Z`)** — multi-category exposure, ensure the control level of the exposure variable is coded as 0.  
- **Potential mediator(s) (`M`)** — one or more potential mediators (e.g., biomarkers, imaging measures).  
- **Potential confounders (`X`)** — potential confounders (continuous or categorical) affecting exposure, outcome, and/or mediators.

The dataset must not contain missing values and should be stored as a .csv file. 

#### 2. Repalce simulated_data.csv in the conduct_causal_mediation_analysis. R with your own data file. 
Example:
```r
data <- read.csv("your_data.csv", header = TRUE)
```
#### 3. Format data for analysis 
Example:
```r
y <- as.matrix(data$Y)
z <- as.matrix(data$Z)
m <- as.matrix(data[, c("Mediator1", "Mediator2")])
x <- as.matrix(data[, c("Confounder1", "Confouer2")])
```
#### 4. Run the revised script conduct_causal_mediation_analysis. R for causal mediation analysis


## 🪪 License

This repository is distributed under the **MIT License**.  
You are free to use, modify, and distribute this code with proper attribution.  
See the [LICENSE](LICENSE) file for full terms.

