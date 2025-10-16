## ========================================================
## Function: medDML
## Description: Estimates total, direct, and indirect effects
##              along with standard errors and p-values using
##              double machine learning
## ========================================================

medDML <- function(y, z, m, x, j, trim, fewsplits, normalized, seed, type.measure, nfolds) {
  temp <- hdmedalt(
    y = y, z = z, m = m, x = x, j = j,
    trim = trim, fewsplits = fewsplits,
    normalized = normalized, seed = seed,
    type.measure = type.measure, nfolds = nfolds
  )
  
  eff <- temp[1:6]
  se <- sqrt(temp[7:12] / temp[13])
  
  results <- rbind(
    effect = eff,
    se = se,
    pval = 2 * pnorm(-abs(eff / se))
  )
  
  colnames(results) <- c("total", "dir.treat", "dir.control", 
                         "indir.treat", "indir.control", "Y(0,M(0))")
  ntrimmed <- length(z) - temp[13]
  
  list(results = results, ntrimmed = ntrimmed)
}



## ============================================================
## Function: hdmedalt
## Description: Implements the underlying estimation procedure
##              with cross-fitting.
## =============================================================

hdmedalt <- function(y, z, m, x, j, trim, fewsplits,
                     normalized, seed, type.measure, nfolds) {
  
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
  
  ## ------------------------------------------------------------
  ## Cross-fitting procedure that splits sample into
  ## training and testing data
  ## ------------------------------------------------------------
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
    pmx <- cv.glmnet(xm[trsample, ], z[trsample], family = "multinomial",
                     type.measure = type.measure, nfolds = nfolds)
    # Predict Pr(Z = j | M, X) and Pr(Z = 0 | M, X) in test data
    pmx_pred <- as.data.frame(predict(pmx, newx = xm[tesample, ], s = "lambda.min", type = "response"))
    pmxtej <- pmx_pred[, 1 + j]
    pmxte0 <- pmx_pred[, 1]
    
    ## ------------------------------------------------------------
    ## Exposure model: Pr(Z = z | X)
    ## ------------------------------------------------------------
    # Fit Pr(Z = z | X) in total training data
    px <- cv.glmnet(x[trsample, ], z[trsample], family = "multinomial",
                    type.measure = type.measure, nfolds = nfolds)
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
      # Predict E(Y | M, X, Z = j) in test data
      eymxjte <- predict(eymxj, xm[tesample, ])
      # Predict E(Y | M, X, Z = j) in delta sample
      eymxjtrte <- predict(eymxj, xm[deltasample, ])
    } else {
      # Fit E(Y | X, Z = 0) in total training data 
      eyx0 <- rlassologit(y[trsample[z[trsample] == 0]] ~ x[trsample[z[trsample] == 0], ])
      # Predict E(Y | X, Z = 0) in test data
      eyx0te <- predict(eyx0, x[tesample, ], type = "response")
      
      # Fit E(Y | M, X, Z = j) in first training data
      eymxj <- rlassologit(y[musample[z[musample] == j]] ~ xm[musample[z[musample] == j], ])
      # Predict E(Y | M, X, Z = j) in test data
      eymxjte <- predict(eymxj, xm[tesample, ], type = "response")
      # Predict E(Y | M, X, Z = j) in delta sample
      eymxjtrte <- predict(eymxj, xm[deltasample, ], type = "response")
    }
    
    # Fit E[E(Y | M, X, Z = j) | Z = 0, X] in delta sample
    regweymxj0 <- rlasso(eymxjtrte[ztrte == 0] ~ xtrte[ztrte == 0, ])
    # Predict E[E(Y | M, X, Z = j) | Z = 0, X] in test data
    regweymxj0te <- predict(regweymxj0, x[tesample, ])
    
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
                         eymxjte, regweymxj0te, eyxjte)[sel == 1, ])
    
    # Collect selection dummies
    selall <- c(selall, sel)
  }
  
  ## ------------------------------------------------------------
  ## Compute scores for potential outcomes
  ## ------------------------------------------------------------
  if (!normalized) {
    yjm0 <- (score[, 1] == j) * score[, 2] / (score[, 3] * score[, 4]) * (score[, 6] - score[, 10]) +
      (score[, 1] == 0) / score[, 4] * (score[, 10] - score[, 11]) + score[, 11]
    yjmj <- (score[, 1] == j) * (score[, 6] - score[, 12]) / score[, 5] + score[, 12]
    y0mj <- (score[, 1] == 0) * score[, 3] / (score[, 2] * score[, 5]) * (score[, 6] - score[, 7]) +
      (score[, 1] == j) / score[, 5] * (score[, 7] - score[, 8]) + score[, 8]
    y0m0 <- (score[, 1] == 0) * (score[, 6] - score[, 9]) / score[, 4] + score[, 9]
  } else {
    nobs <- nrow(score)
    sumscore1 <- sum((score[, 1] == j) * score[, 2] / (score[, 3] * score[, 4]))
    sumscore2 <- sum((score[, 1] == 0) / score[, 4])
    sumscore3 <- sum((score[, 1] == j) / score[, 5])
    sumscore4 <- sum((score[, 1] == 0) * score[, 3] / (score[, 2] * score[, 5]))
    
    yjm0 <- (nobs * (score[, 1] == j) * score[, 2] / (score[, 3] * score[, 4]) *
               (score[, 6] - score[, 10])) / sumscore1 +
      (nobs * (score[, 1] == 0) / score[, 4] *
         (score[, 10] - score[, 11])) / sumscore2 + score[, 11]
    
    yjmj <- (nobs * (score[, 1] == j) * (score[, 6] - score[, 12]) /
               score[, 5]) / sumscore3 + score[, 12]
    
    y0mj <- (nobs * (score[, 1] == 0) * score[, 3] / (score[, 2] * score[, 5]) *
               (score[, 6] - score[, 7])) / sumscore4 +
      (nobs * (score[, 1] == j) / score[, 5] *
         (score[, 7] - score[, 8])) / sumscore3 + score[, 8]
    
    y0m0 <- (nobs * (score[, 1] == 0) * (score[, 6] - score[, 9]) /
               score[, 4]) / sumscore2 + score[, 9]
  }
  
  ## ------------------------------------------------------------
  ## Compute mean potential outcomes and effects
  ## ------------------------------------------------------------
  myjm0 <- mean(yjm0)
  myjmj <- mean(yjmj)
  my0mj <- mean(y0mj)
  my0m0 <- mean(y0m0)
  
  tot <- myjmj - my0m0
  dir1 <- myjmj - my0mj
  dir0 <- myjm0 - my0m0
  indir1 <- myjmj - myjm0
  indir0 <- my0mj - my0m0
  
  ## ------------------------------------------------------------
  ## Compute variances
  ## ------------------------------------------------------------
  vtot <- mean((yjmj - y0m0 - tot)^2)
  vdir1 <- mean((yjmj - y0mj - dir1)^2)
  vdir0 <- mean((yjm0 - y0m0 - dir0)^2)
  vindir1 <- mean((yjmj - yjm0 - indir1)^2)
  vindir0 <- mean((y0mj - y0m0 - indir0)^2)
  vcontrol <- mean((y0m0 - my0m0)^2)
  
  ## ------------------------------------------------------------
  ## Report effects, mean of Y(0,M(0)), variances,
  ## and number of non-trimmed observations
  ## ------------------------------------------------------------ 
  c(tot, dir1, dir0, indir1, indir0, my0m0,
    vtot, vdir1, vdir0, vindir1, vindir0, vcontrol, sum(selall))
}
