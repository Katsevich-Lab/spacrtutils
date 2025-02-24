spaCRT_old <- function(data, X_on_Z_fam, Y_on_Z_fam,
                       normalize = FALSE, return_cdf = FALSE, R = 5) {

  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  # fit X on Z regression
  X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z, family = X_on_Z_fam))

  # fit Y on Z regression
  if(Y_on_Z_fam == "negative.binomial"){
    aux_info_Y_on_Z <- nb_precomp(list(Y = Y, Z = Z))

    Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z,
                                              family = MASS::negative.binomial(aux_info_Y_on_Z$theta_hat),
                                              mustart = aux_info_Y_on_Z$fitted_values))
    NB.disp.param <- aux_info_Y_on_Z$theta_hat
  }else{
    Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z, family = Y_on_Z_fam))
    NB.disp.param <- "Invalid request"
  }

  W <- Y - Y_on_Z_fit$fitted.values
  P <- X_on_Z_fit$fitted.values

  # compute the products of residuals for each observation
  prod_resids <- (X - X_on_Z_fit$fitted.values) * W

  # compute the test statistic
  test_stat <- 1/sqrt(n) * sum(prod_resids)

  # perform saddlepoint approximation
  p_value_opp <- suppressWarnings(spa_cdf_old(t = test_stat + 1/sqrt(n) * sum(P*W),
                                              P = P, W = W,
                                              fam = X_on_Z_fam,
                                              R = abs(R),
                                              max_expansions = 6))

  if(is.nan(p_value_opp) == TRUE | p_value_opp < 0 | p_value_opp > 1){
    temp.gcm <- GCM(data, X_on_Z_fam, Y_on_Z_fam)

    # return test statistic, GCM p-values, and related quantities
    return(list(test_stat = temp.gcm$test_stat,
                p.left = temp.gcm$p.left,
                p.right = temp.gcm$p.right,
                p.both = temp.gcm$p.both,
                NB.disp.param = NB.disp.param,
                gcm.default = TRUE,
                nan.spacrt = is.nan(p_value_opp)))
  }else{
    # return test statistic, spaCRT p-values, and related quantities
    return(list(test_stat = test_stat,
                p.left = p_value_opp,
                p.right = 1 - p_value_opp,
                p.both = 2*min(c(p_value_opp, 1 - p_value_opp)),
                NB.disp.param = NB.disp.param,
                gcm.default = FALSE,
                nan.spacrt = is.nan(p_value_opp)))
  }
}
