###############################################################################
##  Title: Conduct Causal Mediation Analysis
## ----------------------------------------------------------------------------
##  Description:
##    This script performs a causal mediation analysis on simulated data.
##    It estimates total, direct, and indirect (mediated) effects of a
##    multi-category exposure variable (z) on a continuous outcome (y)
##    through potential mediators (m), adjusting for potential confounders (x).
###############################################################################


# ============================================================
# 1. Load required packages and scripts
# ============================================================
required_pkgs <- c("glmnet", "Matrix", "hdm")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

source("medDML_multicategory_exposure.R")

# ----------------------------
# 2. Set parameters
# ----------------------------
seed <- 5
fewsplits <- TRUE
type.measure <- "deviance"
nfolds <- 5
trim <- 0.00
normalized <- TRUE

# ----------------------------
# 3. Read simulated data
# ----------------------------
data <- read.csv("simulated_data.csv", header = TRUE)

# ----------------------------
# 4. Prepare variables
# ----------------------------
y <- as.matrix(data[, 1])          # Outcome
z <- as.matrix(data[, 2])          # Exposure (categorical)
m <- as.matrix(data[, 3:7])        # Potential mediators
x <- as.matrix(data[, -c(1:7)])    # Potential confounders

groups <- length(unique(z)) - 1
med.effect <- direct.effect <- total.effect <- rep(NA, groups)
med.se <- direct.se <- total.se <- rep(NA, groups)
med.pvalue <- direct.pvalue <- total.pvalue <- rep(NA, groups)
med.ci <- direct.ci <- total.ci <- matrix(NA, nrow = groups, ncol = 2)
outputs <- vector("list", groups)

# ----------------------------
# 5. Causal mediation analysis
# ----------------------------
for (j in 1:groups) {
  output <- medDML(
    y = y,
    z = z,
    m = m,
    x = x,
    j = j,
    seed = seed,
    fewsplits = fewsplits,
    type.measure = type.measure,
    nfolds = nfolds,
    trim = trim,
    normalized = normalized
  )
  
  outputs[[j]] <- output
  results <- output$results
  
  # Extract total, direct, and mediation effects, standard errors, p-values, and 95% confidence intervals
  total.effect[j] <- results[1, 1]
  total.se[j] <- results[2, 1]
  total.pvalue[j] <- results[3, 1]
  total.ci[j, ] <- total.effect[j]  + c(-1.96, 1.96) * total.se[j]
  
  direct.effect[j] <- results[1, 2]
  direct.se[j] <- results[2, 2]
  direct.pvalue[j] <- results[3, 2]
  direct.ci[j, ] <- direct.effect[j] + c(-1.96, 1.96) * direct.se[j]
  
  med.effect[j] <- results[1, 3]
  med.se[j] <- results[2, 3]
  med.pvalue[j] <- results[3, 3]
  med.ci[j, ] <- med.effect[j]    + c(-1.96, 1.96) * med.se[j]
}

# ----------------------------
# 6. Summary results
# ----------------------------
z.seq <- seq_len(groups)

df_total <- data.frame(
  z = z.seq,
  estimate = round(total.effect, 2),
  se = round(total.se, 2),
  CI = sprintf("[%.2f, %.2f]", total.ci[, 1], total.ci[, 2]),
  p.value = round(total.pvalue, 3)
)

df_direct <- data.frame(
  z = z.seq,
  estimate = round(direct.effect, 2),
  se = round(direct.se, 2),
  CI = sprintf("[%.2f, %.2f]", direct.ci[, 1], direct.ci[, 2]),
  p.value = round(direct.pvalue, 3)
)

df_med <- data.frame(
  z = z.seq,
  estimate = round(med.effect, 2),
  se = round(med.se, 2),
  CI = sprintf("[%.2f, %.2f]", med.ci[, 1], med.ci[, 2]),
  p.value = round(med.pvalue, 3)
)

# Nicely formatted output
cat("\n=== Total Effects ===\n");  print(df_total)
cat("\n=== Direct Effects ===\n"); print(df_direct)
cat("\n=== Indirect (Mediated) Effects ===\n"); print(df_med)

# ----------------------------
# 7. Save results
# ----------------------------
save(outputs, df_total, df_direct, df_med, file = "mediation_analysis_results.RData")

