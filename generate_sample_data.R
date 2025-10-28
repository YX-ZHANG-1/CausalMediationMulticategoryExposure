####################################################################
##  Title: Generate Sample Data for Causal Mediation Analysis
## ------------------------------------------------------------------
##  Description:
##    This script generates a synthetic dataset for causal mediation
##    analysis mimicking real-world data from Alzheimer’s 
##    Disease Neuroimaging Initiative (ADNI).
##
##    The dataset includes:
##      - Y: Outcome (mimic ADAS-Cog 13 score)
##      - Z: Multi-category exposure (mimic depressive symptom profile)
##      - M: Potential mediators (mimic hippocampal volume, brain ventricular volume, 
##           CSF Abeta42, CSF t-tau, FDG-PET)
##      - X: Potential confounders (mimic DNAm CpG M-values, age, gender, education, 
##           martial status, polygenectic hazard score)
##
##  Note:
##    This script defines the simulation function only.
##    Parameter initialization and execution should be done in a
##    separate driver script (i.e., simulate_dataset.R).
####################################################################


library(truncnorm)


# Logistic helper
expit <- function(x) 1 / (1 + exp(-x))

# ============================
# Main simulation function
# ============================
simulate_data <-  function(seed, n_subj, n_CpG, prop_bimodal, shape1_lower, shape1_upper, 
                           shape2_lower, shape2_upper, shape1a_lower, shape1a_upper, 
                           shape2a_lower, shape2a_upper, shape1b_lower, shape1b_upper, 
                           shape2b_lower, shape2b_upper, w1_lower, w1_upper, mean_X1, 
                           sd_X1, a_X1, b_X1, p_X2, p_X3, mean_X4, sd_X4, a_X4, b_X4, 
                           mean_X5, sd_X5, K, gamma, n_M, sd_e_M, theta_0, theta_X, 
                           theta_Z, beta_0, beta_X, beta_M, sd_e_Y){
  
  set.seed(seed)
  
  # ------------------------------------------------------------
  # 1. Generate potential confounders X
  # ------------------------------------------------------------
  CpG_mat <- matrix(NA, n_subj, n_CpG)
  CpG_type <- ifelse(runif(n_CpG) < prop_bimodal, "bimodal", "unimodal")
  colnames(CpG_mat) <- paste0("CpG", 1:n_CpG)
  
  for (j in 1:n_CpG) {
    if (CpG_type[j] == "unimodal") {
      shape1 <- runif(1, shape1_lower, shape1_upper)
      shape2 <- runif(1, shape2_lower, shape2_upper)
      beta_val <- rbeta(n_subj, shape1, shape2)
    } else {
      shape1a <- runif(1, shape1a_lower, shape1a_upper)
      shape2a <- runif(1, shape2a_lower, shape2a_upper)
      shape1b <- runif(1, shape1b_lower, shape1b_upper)
      shape2b <- runif(1, shape2b_lower, shape2b_upper)
      w1 <- runif(1, w1_lower, w1_upper)  
      comp <- rbinom(n_subj, 1, w1)
      beta_val <- ifelse(comp == 1,
                         rbeta(n_subj, shape1a, shape2a),
                         rbeta(n_subj, shape1b, shape2b))
    }
    beta_val <- pmin(pmax(beta_val, 1e-6), 1 - 1e-6) 
    CpG_mat[, j] <- log2(beta_val / (1 - beta_val))
  }
  
  X1 <- rtruncnorm(n_subj, a = a_X1, b = b_X1, mean = mean_X1, sd = sd_X1)
  X2 <- rbinom(n_subj, 1, p_X2) 
  X3 <- rbinom(n_subj, 1, p_X3) 
  X4 <- round(rtruncnorm(n_subj, a = a_X4, b = b_X4, mean = mean_X4, sd = sd_X4))
  X5 <- rnorm(n_subj, mean = mean_X5, sd = sd_X5)

  X <- as.matrix(cbind(CpG_mat, X1, X2, X3, X4, X5))
  X.scale <- scale(X)
  
  # ------------------------------------------------------------
  # 2. Generate exposure Z
  # ------------------------------------------------------------
  eta <- cbind(0, cbind(1, X.scale) %*% gamma)
  p_mat <- exp(eta) / rowSums(exp(eta)) 
  colnames(p_mat) <- paste0("cat", 0:(K-1))
  Z <- apply(p_mat, 1, function(prob) sample(0:(K-1), size = 1, prob = prob))
  Z_factor <- factor(Z, levels = 0:(K-1))
  Z_mat <-  model.matrix(~ Z_factor)[, -1, drop = FALSE]
  
  # ------------------------------------------------------------
  # 3. Generate potential mediators M
  # ------------------------------------------------------------
  M <- matrix(NA, n_subj, n_M)
  for (m in 1:n_M){
    M[, m] <- theta_0[m] + X.scale %*% theta_X[, m] + Z_mat %*% theta_Z[, m] + rnorm(n_subj, 0, sd_e_M[m])
  }
  colnames(M) <- paste0("M", 1:n_M)
  M.scale <- scale(M)
  
  # ------------------------------------------------------------
  # 4. Generate outcome Y
  # ------------------------------------------------------------
  Y <- numeric(n_subj)
  for (i in 1:n_subj) {
    j <- Z[i] + 1 
    Y[i] <- beta_0[j] + X.scale[i, ] %*% beta_X[, j] + M.scale[i, ] %*% beta_M[, j] + rnorm(1, 0, sd_e_Y)
  }
  
  # ------------------------------------------------------------
  # 5. Return final dataset
  # ------------------------------------------------------------
  fulldf <- as.data.frame(cbind(Y = Y, Z=Z, M, X))
  
  return(fulldf)
}

