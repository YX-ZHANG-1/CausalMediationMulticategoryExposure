############################################################
##  Title: Preprocessing and Normalizing DNA Methylation (DNAm) Data
##  Description: This script reads Illumina IDAT files,
##               normalizes DNAm data using dasen, and 
##               saves M-values for downstream analysis.
############################################################

## ============================
## 1. Load or install packages
## ============================
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c(
  "minfi", "minfiData", "wateRmelon",
  "methylumi", "dplyr", "here"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    BiocManager::install(pkg, ask = FALSE, update = TRUE)
  }
  library(pkg, character.only = TRUE)
}

## ============================
## 2. Set working directory
## ============================
# Replace the path inside setwd() with the actual folder on your computer.
setwd("/Users/yourname/Documents/Project")

## ============================
## 3. Define paths
## ============================
# The base directory of your own DNAm IDAT files
# For example, the folder name may be called "IDAT_files"
baseDir <- "./IDAT_files"

# The name of your own DNAm sample sheet
# For example, the sample sheet name may be called "SampleSheet.csv"
sampleSheetFile <- "SampleSheet.csv"

## ============================
## 4. Read sample sheet
## ============================
targets <- read.metharray.sheet(
  base = baseDir,
  pattern = sampleSheetFile
)

## ============================
## 5. Read IDAT files
## ============================
RGset <- read.metharray.exp(targets = targets)

## ============================
## 6. Normalization using dasen
## ============================
RGset_dasen <- dasen(RGset)

## ============================
## 7. Calculate Beta and M values
## ============================
Beta_values <- getBeta(RGset_dasen, offset = 100, betaThreshold = 0.001)
M_values <- logit2(Beta_values)

## ============================
## 8. Save processed data
## ============================
save(M_values, file = "./Mvalues.RData")

