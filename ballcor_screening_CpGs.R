############################################################
##  Title: CpG Screening Using Ball Correlation
##  Description:
##    - Load DNA methylation M-values
##    - Apply Ball Correlation Screening (BcorSIS)
##    - Save selected CpGs for downstream mediation analysis
############################################################

## ============================
## 1. Load package
## ============================
# If MFSIS is not installed, uncomment the following line:
# install.packages("MFSIS")

library(MFSIS)

## ============================
## 2. Set working directory
## ============================
# Replace the path inside setwd() with the actual folder on your computer.
setwd("/Users/yourname/Documents/Project")

## ============================
## 3. Read data that include the outcome variable of interest
## ============================
# Specify the directory of your own data file
# For example, the data name may be called "data.csv"
data_path <- "./data.csv"
data <- read.csv(data_path, header = TRUE)

# Specify the outcome variable of interest
# For example, the outcome variable may be called "ADAS13"
outcome <- "ADAS13"
y <- data[, outcome]

## ============================
## 4. Load M-values
## ============================
mvalue_file <- "./Mvalues.RData"
load(mvalue_file)  

# Transpose so rows = samples, columns = CpGs
x0 <- as.matrix(t(M_values))

## ============================
## 5. CpG screening using Ball Correlation
## ============================
selected.cpgs <- BcorSIS(x0, y)
x.selected.cpgs <- x0[, selected.cpgs]

## ============================
## 6. Save results
## ============================
save(x.selected.cpgs, file = "./Selected_CpGs.RData")

