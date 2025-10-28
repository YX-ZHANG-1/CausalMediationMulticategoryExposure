########################################################################
##  Title: Double Machine Learning for Causal Mediation Analysis
##
##  Description:
##    This script implements the Double Machine Learning (DML)
##    algorithm for estimating total, direct, and indirect
##    (mediated) effects under a multi-category exposure.
##
##  This file defines two main functions:
##    - medDML():  Wrapper function for estimation and inference
##    - hdmedalt(): Core estimation routine with cross-fitting
########################################################################


# ============================================================
# Function: medDML
# ============================================================
#' @title medDML
#' @description
#' Estimates total, direct, and indirect effects using
#' double machine learning (DML) and computes standard
#' errors and p-values for inference.
#'
#' @param y Outcome variable.
#' @param z Multi-category exposure variable.
#' @param m Potential mediators.
#' @param x Potential pre-treatment confounders of the exposure, mediator, and/or outcome.
#' @param j Active exposure level to be analyzed (integer, e.g., 1 for Z=1 vs 0 baseline).
#' @param seed Random seed for reproducibility. 
#' @param fewsplits Logical; if TRUE, reuses the same training data for nested nuisance models.
#' @param type.measure Loss function used in cv.glmnet (e.g., "deviance").
#' @param nfolds Number of cross-validation folds for cv.glmnet.
#' @param trim Trimming threshold for small denominators.
#' @param normalized Logical; if TRUE, normalizes inverse-probability weights within exposure levels.
#'
#' @return A list containing:
#'   \describe{
#'     \item{results}{Matrix with effect estimates, SEs, and p-values.}
#'     \item{ntrimmed}{Number of observations trimmed due to small denominators.}
#'   }
#' @export
########################################################################

medDML <- function(y, z, m, x, j, seed,
                   fewsplits, type.measure,
                   nfolds, trim, normalized) {
  
  # Run main estimation routine
  temp <- hdmedalt(
    y = y, z = z, m = m, x = x, j = j,
    seed = seed, fewsplits = fewsplits,
    type.measure = type.measure, nfolds = nfolds,
    trim = trim, normalized = normalized
  )
  
  # Extract point estimates and variances
  eff <- temp[1:3]
  se <- sqrt(temp[4:6] / temp[7])
  
  # Construct summary table
  results <- rbind(
    Estimate = eff,
    SE = se,
    P_value = 2 * pnorm(-abs(eff / se))
  )
  colnames(results) <- c("total", "direct", "indirect")
  
  # Number of trimmed observations
  ntrimmed <- length(z) - temp[13]
  
  list(results = results, ntrimmed = ntrimmed)
}


# ============================================================
# Function: hdmedalt
# ============================================================
#' @title hdmedalt
#' @description
#' Core cross-fitting procedure implementing double machine
#' learning for mediation analysis with multi-category exposures.
#'
#' This function is internal and is typically called by `medDML()`.
########################################################################

hdmedalt <- function(y, z, m, x, j, seed,
                     fewsplits, type.measure,
                     nfolds, trim, normalized) {
  
  # Identify binary outcome
  ybin <- 1 * (length(unique(y)) == 2 & min(y) == 0 & max(y) == 1)
  
  # Combine x and m for joint modeling
  xm <- cbind(x, m)
  
  # Determine sample splitting sizes
  stepsize <- ceiling(length(z) / 3)
  nobs <- min(3 * stepsize, length(z))
  
  # Randomly split samples
  set.seed(seed)
  idx <- sample(nobs)
  sample1 <- idx[1:stepsize]
  sample2 <- idx[(stepsize + 1):(2 * stepsize)]
  sample3 <- idx[(2 * stepsize + 1):nobs]
  
  score <- c()
  selall <- c()
  
  # ------------------------------------------------------------
  # Cross-fitting procedure that splits sample into
  # training and testing data
  # ------------------------------------------------------------
  for (i in 1:3) {
    if (i == 1) { tesample <- sample1; musample <- sample2; deltasample <- sample3 }
    if (i == 2) { tesample <- sample3; musample <- sample1; deltasample <- sample2 }
    if (i == 3) { tesample <- sample2; musample <- sample3; deltasample <- sample1 }
    
    trsample <- c(musample, deltasample)
    zte <- z[tesample]
    yte <- y[tesample]
    
    if (fewsplits == 1) {
      musample <- c(musample, deltasample)
      deltasample <- musample
    }
    
    x <- as.matrix(x)
    xm <- as.matrix(xm)
    
    ## ------------------------------------------------------------
    ## Exposure model: Pr(Z = z | M, X)
    ## ------------------------------------------------------------
    # Fit Pr(Z = z | M, X) in total training data
    pmx <- suppressWarnings(cv.glmnet(xm[trsample, ], z[trsample], family = "multinomial",
                                      type.measure = type.measure, nfolds = nfolds))
    # Predict Pr(Z = j | M, X) and Pr(Z = 0 | M, X) in test data
    pmx_pred <- as.data.frame(predict(pmx, newx = xm[tesample, ], s = "lambda.min", type = "response"))
    pmxtej <- pmx_pred[, 1 + j]
    pmxte0 <- pmx_pred[, 1]
    
    ## ------------------------------------------------------------
    ## Exposure model: Pr(Z = z | X)
    ## ------------------------------------------------------------
    # Fit Pr(Z = z | X) in total training data
    px <- suppressWarnings(cv.glmnet(x[trsample, ], z[trsample], family = "multinomial",
                                     type.measure = type.measure, nfolds = nfolds))
    # Predict Pr(Z = j | X) and Pr(Z = 0 | X) in test data
    px_pred <- as.data.frame(predict(px, newx = x[tesample, ], s = "lambda.min", type = "response"))
    pxtej <- px_pred[, 1 + j]
    pxte0 <- px_pred[, 1]
    
    ## ------------------------------------------------------------
    ## Outcome model: E(Y | M, X, Z = 0)
    ## ------------------------------------------------------------
    if (ybin != 1) {
      # Fit E(Y | M, X, Z = 0) in first training data
      eymx0 <- rlasso(y[musample[z[musample] == 0]] ~ xm[musample[z[musample] == 0], ])
      # Predict E(Y | M, X, Z = 0) in test data
      eymx0te <- predict(eymx0, xm[tesample, ])
      # Predict E(Y | M, X, Z = 0) in delta sample
      eymx0trte <- predict(eymx0, xm[deltasample, ])
    } else {
      # Fit E(Y | M, X, Z = 0) in first training data
      eymx0 <- rlassologit(y[musample[z[musample] == 0]] ~ xm[musample[z[musample] == 0], ])
      # Predict E(Y | M, X, Z = 0) in test data
      eymx0te <- predict(eymx0, xm[tesample, ], type = "response")
      # Predict E(Y | M, X, Z = 0) in delta sample
      eymx0trte <- predict(eymx0, xm[deltasample, ], type = "response")
    }
    
    # Fit E[E(Y | M, X, Z = 0) | Z = j, X] in delta sample
    ztrte <- z[deltasample]
    xtrte <- x[deltasample, ]
    regweymx0j <- rlasso(eymx0trte[ztrte == j] ~ xtrte[ztrte == j, ])
    # Predict E[E(Y | M, X, Z = 0) | Z = j, X] in test data
    regweymx0jte <- predict(regweymx0j, x[tesample, ])
    
    ## ------------------------------------------------------------
    ## Outcome models: E(Y | X, Z = 0) and E(Y | M, X, Z = j)
    ## ------------------------------------------------------------
    if (ybin != 1) {
      # Fit E(Y | X, Z = 0) in total training data 
      eyx0 <- rlasso(y[trsample[z[trsample] == 0]] ~ x[trsample[z[trsample] == 0], ])
      # Predict E(Y | X, Z = 0) in test data
      eyx0te <- predict(eyx0, x[tesample, ])
      
      # Fit E(Y | M, X, Z = j) in first training data
      eymxj <- rlasso(y[musample[z[musample] == j]] ~ xm[musample[z[musample] == j], ])
      # Predict E(Y | M, X, Z = j) in delta sample
      eymxjtrte <- predict(eymxj, xm[deltasample, ])
    } else {
      # Fit E(Y | X, Z = 0) in total training data 
      eyx0 <- rlassologit(y[trsample[z[trsample] == 0]] ~ x[trsample[z[trsample] == 0], ])
      # Predict E(Y | X, Z = 0) in test data
      eyx0te <- predict(eyx0, x[tesample, ], type = "response")
      
      # Fit E(Y | M, X, Z = j) in first training data
      eymxj <- rlassologit(y[musample[z[musample] == j]] ~ xm[musample[z[musample] == j], ])
      # Predict E(Y | M, X, Z = j) in delta sample
      eymxjtrte <- predict(eymxj, xm[deltasample, ], type = "response")
    }
    
    # Fit E[E(Y | M, X, Z = j) | Z = 0, X] in delta sample
    regweymxj0 <- rlasso(eymxjtrte[ztrte == 0] ~ xtrte[ztrte == 0, ])
    
    ## ------------------------------------------------------------
    ## Outcome model: E(Y | X, Z = j)
    ## ------------------------------------------------------------
    if (ybin != 1) {
      # Fit E(Y | X, Z = j) in total training data 
      eyxj <- rlasso(y[trsample[z[trsample] == j]] ~ x[trsample[z[trsample] == j], ])
      # Predict E(Y | X, Z = j) in test data
      eyxjte <- predict(eyxj, x[tesample, ])
    } else {
      # Fit E(Y | X, Z = j) in total training data
      eyxj <- rlassologit(y[trsample[z[trsample] == j]] ~ x[trsample[z[trsample] == j], ])
      # Predict E(Y | X, Z = j) in test data
      eyxjte <- predict(eyxj, x[tesample, ], type = "response")
    }
    
    ## ------------------------------------------------------------
    ## Select observations satisfying trimming restriction
    ## ------------------------------------------------------------
    sel <- 1 * (((pmxte0 * pxtej) >= trim) &
                  (pxtej >= trim) &
                  (pxte0 >= trim) &
                  ((pmxtej * pxte0) >= trim))
    
    # Select elements of the score functions
    score <- rbind(score,
                   cbind(zte, pmxte0, pmxtej, pxte0, pxtej,
                         yte, eymx0te, regweymx0jte, eyx0te,
                          eyxjte)[sel == 1, ])
    
    # Collect selection dummies
    selall <- c(selall, sel)
  }
  
  # ------------------------------------------------------------
  # Compute scores for potential outcomes
  # ------------------------------------------------------------
  if (!normalized) {
    yjmj <- (score[, 1] == j) * (score[, 6] - score[, 10]) / score[, 5] + score[, 10]
    y0mj <- (score[, 1] == 0) * score[, 3] / (score[, 2] * score[, 5]) * (score[, 6] - score[, 7]) +
      (score[, 1] == j) / score[, 5] * (score[, 7] - score[, 8]) + score[, 8]
    y0m0 <- (score[, 1] == 0) * (score[, 6] - score[, 9]) / score[, 4] + score[, 9]
  } else {
    nobs <- nrow(score)
    sumscore1 <- sum((score[, 1] == 0) / score[, 4])
    sumscore2 <- sum((score[, 1] == j) / score[, 5])
    sumscore3 <- sum((score[, 1] == 0) * score[, 3] / (score[, 2] * score[, 5]))
    
    
    yjmj <- (nobs * (score[, 1] == j) * (score[, 6] - score[, 10]) /
               score[, 5]) / sumscore2 + score[, 10]
    
    y0mj <- (nobs * (score[, 1] == 0) * score[, 3] / (score[, 2] * score[, 5]) *
               (score[, 6] - score[, 7])) / sumscore3 +
      (nobs * (score[, 1] == j) / score[, 5] *
         (score[, 7] - score[, 8])) / sumscore2 + score[, 8]
    
    y0m0 <- (nobs * (score[, 1] == 0) * (score[, 6] - score[, 9]) /
               score[, 4]) / sumscore1 + score[, 9]
  }
  
  # ------------------------------------------------------------
  # Compute mean potential outcomes and effects
  # ------------------------------------------------------------
  myjmj <- mean(yjmj)
  my0mj <- mean(y0mj)
  my0m0 <- mean(y0m0)
  
  tot <- myjmj - my0m0
  dir <- myjmj - my0mj
  indir <- my0mj - my0m0
  
  # ------------------------------------------------------------
  # Compute variances
  # ------------------------------------------------------------
  vtot <- mean((yjmj - y0m0 - tot)^2)
  vdir <- mean((yjmj - y0mj - dir)^2)
  vindir <- mean((y0mj - y0m0 - indir)^2)

  # ------------------------------------------------------------
  # Return results
  # ------------------------------------------------------------ 
  return(c(tot, dir, indir, vtot, vdir, vindir, sum(selall)))
}
