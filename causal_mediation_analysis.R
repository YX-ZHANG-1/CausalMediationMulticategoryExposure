####################################################################
##  Title: Conduct Causal Mediation Analysis
##  Description:
##    This script performs the full mediation analysis 
##    for examining the effect of depression on
##    cognitive outcomes (e.g., ADAS13) through potential mediators.
####################################################################

## ============================
## 1. Set working directory
## ============================
# Replace the path inside setwd() with the actual folder on your computer.
setwd("/Users/yourname/Documents/Project")

## ----------------------------
## 2. Set analysis parameters
## ----------------------------
outcome <- "ADAS13"
ncpgs <- 30
trim <- 0.00
fewsplits <- TRUE
normalized <- TRUE
seed <- 5
type.measure <- "deviance"
nfolds <- 5

## ----------------------------
## 3. Load packages and scripts
## ----------------------------
# If packages are not installed, uncomment the following lines:
# install.packages("dplyr")
# install.packages("LARF")
# install.packages("hdm")
# install.packages("VGAM")
# install.packages("glmnet")
# install.packages("Matrix")

library(dplyr)
library(LARF)
library(hdm)
library(VGAM)
library(glmnet)
library(Matrix)

source("./medDML_multicategory_exposure.R")

## ----------------------------
## 4. Load selected CpGs
## ----------------------------
load("./Selected_CpGs.RData")

## ----------------------------
## 5. Read data
## ----------------------------
# Specify the directory of your own data file
# For example, the data name may be called "data.csv"
data_path <- "./data.csv"
data <- read.csv(data_path, header = TRUE)
dim(data)

## ----------------------------
## 6. Depression variable processing
## ----------------------------
gd_vars <- c("GDSATIS", "GDDROP", "GDEMPTY", "GDBORED", "GDSPIRIT",
             "GDAFRAID", "GDHAPPY", "GDHELP", "GDHOME", "GDMEMORY", "GDALIVE",
             "GDWORTH", "GDENERGY", "GDHOPE", "GDBETTER")
data.GDSCALE <- data[, gd_vars]

# Reverses the coding for specific items: 1, 5, 7, 11, and 13.
data.GDSCALE[, c(1, 5, 7, 11, 13)] <- 1 - data.GDSCALE[, c(1, 5, 7, 11, 13)]

## ----------------------------
## 7. GDSCALE grouping
## ----------------------------
# Emotional domain
GDSCALE.group1 <- apply(data.GDSCALE[, c(1, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15)], 1, max)

# Social domain
GDSCALE.group2 <- apply(data.GDSCALE[, c(2, 9)], 1, max)

# Memory domain
GDSCALE.group3 <- data.GDSCALE[, 10]

## ------------------------------------------
## 8. Create a new variable about depression
## ------------------------------------------
data.GDSCALE.new <- data.frame(GDSCALE.group1, GDSCALE.group2, GDSCALE.group3)

data.GDSCALE.new %>%
  group_by(across(everything())) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

data.GDTOTAL <- rowSums(data.GDSCALE.new)

data$GDPROFILE <- ifelse(data.GDTOTAL == 0, 0,
                       ifelse((data.GDTOTAL==1)&(GDSCALE.group1==1), 1,
                              ifelse((data.GDTOTAL==1)&(GDSCALE.group2==1), 2,
                                     ifelse((data.GDTOTAL==1)&(GDSCALE.group3==1), 3,
                                            ifelse((data.GDTOTAL==2)&(GDSCALE.group1==1)&(GDSCALE.group2==1), 4,
                                                   ifelse((data.GDTOTAL==2)&(GDSCALE.group1==1)&(GDSCALE.group3==1), 5,
                                                          ifelse((data.GDTOTAL==2)&(GDSCALE.group2==1)&(GDSCALE.group3==1), 6, 
                                                                 7)))))))

## ----------------------------
## 9. Get a subset of data
## ----------------------------
vars_subset <- c(outcome, "GDPROFILE", "Hippocampus", "Ventricles", "ABETA", "TAU",
                 "FDG", "AGE", "PTGENDER", "PTEDUCAT", "PTMARRY", "PHS")
data.subset <- data[, vars_subset]
data.subset$PTGENDER <- ifelse(data.subset$PTGENDER == "M", 1, 0)
data.subset$PTMARRY <- ifelse(data.subset$PTMARRY == "Married", 1, 0)

## ----------------------------
## 10. Prepare data for mediation analysis
## ----------------------------
y <- as.matrix(data.subset[, 1])
z <- as.matrix(data.subset[, 2])
m <- as.matrix(data.subset[, 3:7])
x <- as.matrix(cbind(data.subset[, -c(1:7)], x.selected.cpgs[, 1:ncpgs]))

groups <- length(table(z)) - 1
med.effect <- direct.effect <- total.effect <- rep(0, groups)
med.se <- direct.se <- total.se <- rep(0, groups)
med.pvalue <- direct.pvalue <- total.pvalue <- rep(0, groups)
med.ci <- direct.ci <- total.ci <- matrix(0, nrow = groups, ncol = 2)
outputs <- vector("list", groups)

## ----------------------------
## 11. Mediation analysis
## ----------------------------
for (j in 1:groups) {
  output <- medDML(
    y = y,
    z = z,
    m = m,
    x = x,
    j = j,
    trim = trim,
    fewsplits = fewsplits,
    normalized = normalized,
    seed = seed,
    type.measure = type.measure,
    nfolds = nfolds
  )
  
  outputs[[j]] <- output
  results <- output$results
  
  med.effect[j] <- results[1, 5]
  med.se[j] <- results[2, 5]
  med.pvalue[j] <- results[3, 5]
  med.ci[j, ] <- c(med.effect[j] - 1.96 * med.se[j], med.effect[j] + 1.96 * med.se[j])
  
  direct.effect[j] <- results[1, 2]
  direct.se[j] <- results[2, 2]
  direct.pvalue[j] <- results[3, 2]
  direct.ci[j, ] <- c(direct.effect[j] - 1.96 * direct.se[j], direct.effect[j] + 1.96 * direct.se[j])
  
  total.effect[j] <- results[1, 1]
  total.se[j] <- results[2, 1]
  total.pvalue[j] <- results[3, 1]
  total.ci[j, ] <- c(total.effect[j] - 1.96 * total.se[j], total.effect[j] + 1.96 * total.se[j])
}

## ----------------------------
## 12. Summary tables
## ----------------------------
z.seq <- seq_len(groups)

df_total <- data.frame(
  depression = z.seq,
  estimate = round(total.effect, digits = 2),
  se = round(total.se, digits =2),
  CI = sprintf("[%.2f, %.2f]", total.ci[, 1], total.ci[, 2]),
  p.value = round(total.pvalue, digits = 3)
)
print(df_total)

df_direct <- data.frame(
  depression = z.seq,
  estimate = round(direct.effect, digits = 2),
  se = round(direct.se, digits =2),
  CI = sprintf("[%.2f, %.2f]", direct.ci[, 1], direct.ci[, 2]),
  p.value = round(direct.pvalue, digits = 3)
)
print(df_direct)

df_med <- data.frame(
  depression = z.seq,
  estimate = round(med.effect, digits = 2),
  se = round(med.se, digits =2),
  CI = sprintf("[%.2f, %.2f]", med.ci[, 1], med.ci[, 2]),
  p.value = round(med.pvalue, digits = 3)
)
print(df_med)

## ----------------------------
## 13. Save results
## ----------------------------
save(outputs, file = "./mediation_analysis_results.RData")

