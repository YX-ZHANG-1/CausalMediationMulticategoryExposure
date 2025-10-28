####################################################################
##  Title: Generate a Simulated Dataset
## -----------------------------------------------------------------
##  Description:
##    Driver script that sources `generate_sample_data.R`,
##    sets parameters, calls `simulate_data()`, and saves the
##    output dataset for reproducible experiments.
####################################################################


# ----------------------------
# 1. Load required packages and scripts
# ----------------------------
required_pkgs <- c("truncnorm")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

source("generate_sample_data.R")

# ----------------------------
# 2. Global settings
# ----------------------------
seed <- 123
set.seed(seed)
n_subj <- 430
n_CpG <- 30
p <- n_CpG + 5
K <- 8
n_M <- 5

# ----------------------------
# 3. Parameter setup
# ----------------------------
## Parameters for potentail confounders X
prop_bimodal <- 0.5
shape1_lower <- 1
shape1_upper <- 280
shape2_lower <- 5
shape2_upper <- 60
shape1a_lower <- 3
shape1a_upper <- 100
shape2a_lower <- 13
shape2a_upper <- 175
shape1b_lower <- 28
shape1b_upper <- 180
shape2b_lower <- 20
shape2b_upper <- 90
w1_lower <- 0.4
w1_upper <- 0.6
  
mean_X1 <- 72.63
sd_X1 <- 7.10
a_X1 <- 55
b_X1 <- 91.4
p_X2 <- 0.53 
p_X3 <- 0.75 
mean_X4 <- 16.29 
sd_X4 <- 2.63
a_X4 <- 8
b_X4 <- 20
mean_X5 <- 0.27 
sd_X5 <- 0.76

## Parameters for exposure Z
prop_zero_Z <- 0.3
prop_partial_Z <- 0.4
n_zero <- floor(prop_zero_Z * p)
n_partial <- floor(prop_partial_Z * p)
n_full <- p - n_zero - n_partial
var_types <- rep(NA, p)
var_types[sample(1:p, n_zero)] <- "zero"
remaining <- setdiff(1:p, which(var_types == "zero"))
var_types[sample(remaining, n_partial)] <- "partial"
var_types[is.na(var_types)] <- "full"
names(var_types) <- paste0("X", 1:p)

gamma <- matrix(0, nrow = p + 1, ncol = K - 1)
target_prob <- c(0.31, 0.13, 0.06, 0.13, 0.10, 0.13, 0.05, 0.09)
gamma[1, ] <- log(target_prob[-1] / target_prob[1])
for (j in 1:p) {
  if (var_types[j] == "full") {
    gamma[j + 1, ] <- rnorm(K - 1, 0, 0.3)
  } else if (var_types[j] == "partial") {
    active_cats <- sample(1:(K - 1), size = sample(1:(K - 2), 1))
    gamma[j + 1, active_cats] <- rnorm(length(active_cats), 0, 0.3)
  }
}

## Parameters for potential mediators M
prop_zero_x_M <- 0.3
prop_zero_z_M <- 0.3
sd_x_M <- 0.1
min_z_M <- 2
max_z_M <- 10
sd_e_M <- c(1100, 20300, 635, 119, 0.1)
theta_0 <- c(6900, 380000, 1900, 400, 1.5)
theta_X <- matrix(0, nrow = p, ncol = n_M)
theta_Z <- matrix(0, nrow = K - 1, ncol = n_M)
for (m in 1:n_M) {
  active_x <- sample(1:p, size = floor((1 - prop_zero_x_M) * p))
  theta_X[active_x, m] <- rnorm(length(active_x), 0, sd_x_M)
  active_z <- sample(1:(K-1), size = floor((1 - prop_zero_z_M) * (K - 1)))
  theta_Z[active_z, m] <- runif(length(active_z), min_z_M, max_z_M)
}

## Parameters for outcome Y
prop_zero_x_Y <- 0.3
prop_zero_m_Y <- 0.3
sd_x_Y <- 0.1
min_m_Y <- 2
max_m_Y <- 10
sd_e_Y <- 6
beta_0 <- runif(K, min = 40, max = 50)
beta_X <- matrix(0, nrow = p, ncol = K)
beta_M <- matrix(0, nrow = n_M, ncol = K)
for (j in 1:K) {
  active_x <- sample(1:p, size = floor((1 - prop_zero_x_Y) * p))
  beta_X[active_x, j] <- rnorm(length(active_x), 0, sd_x_Y)
  active_m <- sample(1:n_M, size = floor((1 - prop_zero_m_Y) * n_M))
  beta_M[active_m, j] <- runif(length(active_m), min_m_Y, max_m_Y)
}

# ----------------------------
# 4. Generate a dataset
# ----------------------------
simu_df <- simulate_data(seed, n_subj, n_CpG, prop_bimodal, shape1_lower, shape1_upper, 
                         shape2_lower, shape2_upper, shape1a_lower, shape1a_upper, 
                         shape2a_lower, shape2a_upper, shape1b_lower, shape1b_upper, 
                         shape2b_lower, shape2b_upper, w1_lower, w1_upper, mean_X1, 
                         sd_X1, a_X1, b_X1, p_X2, p_X3, mean_X4, sd_X4, a_X4, b_X4, 
                         mean_X5, sd_X5, K, gamma, n_M, sd_e_M, theta_0, theta_X, 
                         theta_Z, beta_0, beta_X, beta_M, sd_e_Y)

# ----------------------------
# 5. Save the generated dataset
# ----------------------------
write.csv(simu_df, file = "simulated_data.csv", row.names = FALSE)

