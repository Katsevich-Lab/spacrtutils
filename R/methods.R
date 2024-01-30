#####################################################################################
#' The generalized covariance measure test of Shah and Peters.
#'
#' \code{GCM} is a function carrying out the GCM test based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param fit_glm_X A logical variable indicating whether to use GLM for estimating `E[X|Z]`.
#' Default is TRUE. If set to FALSE, mean of X is used to estimate `E[X|Z]`.
#' @param fit_glm_Y A logical variable indicating whether to use GLM for estimating `E[Y|Z]`.
#' Default is TRUE. If set to FALSE, mean of Y is used to estimate `E[Y|Z]`.
#' @param aux_info_X_on_Z The auxiliary information that may be used for complex GLM regression
#' (For instance, when X_on_Z_fam = negative.binomial, the dispersion parameter should be provided).
#' @param aux_info_Y_on_Z The auxiliary information that may be used for complex GLM regression
#' (For instance, when Y_on_Z_fam = negative.binomial, the dispersion parameter should be provided).
#'
#' @return A named list with fields \code{test_stat} and \code{left_side_p_value},
#' \code{right_side_p_value}, \code{both_side_p_value}.
#'
#' @examples
#' n <- 20; p <- 2
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' results <- GCM(data, X_on_Z_fam, Y_on_Z_fam)
#' results$test_stat
#' results$p_value
#' @export
GCM <- function(data, X_on_Z_fam = NULL, Y_on_Z_fam = NULL,
                fit_glm_X = TRUE, fit_glm_Y = TRUE,
                aux_info_X_on_Z = NULL, aux_info_Y_on_Z = NULL) {

  if(is.null(X_on_Z_fam) | is.null(Y_on_Z_fam)){
    stop("X_on_Z_fam and Y_on_Z_fam can't be empty!")
  }

  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z

  # fit X on Z regression
  if(!fit_glm_X){
    X_on_Z_fit <- list(fitted.values = rep(mean(X), length(X)))
  }else{
    if(X_on_Z_fam == "negative.binomial"){
      X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z,
                          family = MASS::negative.binomial(aux_info_X_on_Z$theta_hat),
                          mustart = aux_info_X_on_Z$fitted_values))
    }else{
      X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z, family = X_on_Z_fam))
    }
  }

  # fit Y on Z regression
  if(!fit_glm_Y){
    Y_on_Z_fit <- list(fitted.values = rep(mean(Y), length(Y)))
  }else{
    if(Y_on_Z_fam == "negative.binomial"){
      Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z,
                          family = MASS::negative.binomial(aux_info_Y_on_Z$theta_hat),
                          mustart = aux_info_Y_on_Z$fitted_values))
      NB.disp.param <- aux_info_Y_on_Z$theta_hat
    }else{
      Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z, family = Y_on_Z_fam))
      NB.disp.param <- "Invalid request"
    }
  }

  # compute the products of residuals for each observation
  prod_resids <- (X - X_on_Z_fit$fitted.values)*(Y - Y_on_Z_fit$fitted.values)

  # compute the test statistic
  n <- length(X)
  test_stat <- 1/sqrt(n)*sum(prod_resids)/stats::sd(prod_resids) * sqrt(n/(n-1))

  # return test statistic, GCM p-values, and related quantities
  return(list(test_stat = test_stat,
              left_side_p_value = stats::pnorm(test_stat, lower.tail = TRUE),
              right_side_p_value = stats::pnorm(test_stat, lower.tail = FALSE),
              both_side_p_value = 2*stats::pnorm(abs(test_stat), lower.tail = FALSE),
              NB.disp.param = NB.disp.param,
              unnormalized_test_stat = 1/sqrt(n)*sum(prod_resids)))
}


#####################################################################################
#' The distilled conditional randomization test.
#'
#' \code{dCRT} is a function carrying out the dCRT based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param B The number of resamples to draw.
#' @param normalize A logical variable indicating whether the dCRT should be
#' normalized.
#' @param return_resamples A logical variable indicating whether to return the
#' resampled test statistics.
#' @param fit_glm_X A logical variable indicating whether to use GLM for estimating `E[X|Z]`.
#' Default is TRUE. If set to FALSE, mean of X is used to estimate `E[X|Z]`.
#' @param fit_glm_Y A logical variable indicating whether to use GLM for estimating `E[Y|Z]`.
#' Default is TRUE. If set to FALSE, mean of Y is used to estimate `E[Y|Z]`.
#' @param aux_info_X_on_Z The auxiliary information that may be used for complex GLM regression
#' (For instance, when X_on_Z_fam = negative.binomial, the dispersion parameter should be provided).
#' @param aux_info_Y_on_Z The auxiliary information that may be used for complex GLM regression
#' (For instance, when Y_on_Z_fam = negative.binomial, the dispersion parameter should be provided).
#'
#' @return A named list with fields \code{test_stat}, \code{left_side_p_value},
#' \code{right_side_p_value}, \code{both_side_p_value} and
#' \code{resamples}. Here, \code{resamples} is a vector of length \code{B}. It
#' is returned only if \code{return_resamples == TRUE}.
#'
#' @examples
#' n <- 20; p <- 2; B <- 100; normalize <- FALSE; return_resamples <- FALSE
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' results <- dCRT(data, X_on_Z_fam, Y_on_Z_fam, B, normalize, return_resamples)
#' results$test_stat
#' results$p_value
#' @export
dCRT <- function(data, X_on_Z_fam = NULL, Y_on_Z_fam = NULL, B = 2000,
                 normalize = FALSE, return_resamples = FALSE,
                 fit_glm_X = TRUE, fit_glm_Y = TRUE,
                 aux_info_X_on_Z = NULL, aux_info_Y_on_Z = NULL) {

  if(is.null(X_on_Z_fam) | is.null(Y_on_Z_fam)){
    stop("X_on_Z_fam and Y_on_Z_fam can't be empty!")
  }

  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  # fit X on Z regression
  if(!fit_glm_X){
    X_on_Z_fit <- list(fitted.values = rep(mean(X), length(X)))
  }else{
    if(X_on_Z_fam == "negative.binomial"){
      X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z,
                            family = MASS::negative.binomial(aux_info_X_on_Z$theta_hat),
                            mustart = aux_info_X_on_Z$fitted_values))
    }else{
      X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z, family = X_on_Z_fam))
    }
  }

  # fit Y on Z regression
  if(!fit_glm_Y){
    Y_on_Z_fit <- list(fitted.values = rep(mean(Y), length(Y)))
  }else{
    if(Y_on_Z_fam == "negative.binomial"){
      Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z,
                            family = MASS::negative.binomial(aux_info_Y_on_Z$theta_hat),
                            mustart = aux_info_Y_on_Z$fitted_values))
      NB.disp.param <- aux_info_Y_on_Z$theta_hat
    }else{
      Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z, family = Y_on_Z_fam))
      NB.disp.param <- "Invalid request"
    }
  }

  # compute the products of residuals for each observation
  prod_resids <- (X - X_on_Z_fit$fitted.values)*(Y - Y_on_Z_fit$fitted.values)
  # compute the test statistic
  test_stat <- 1/sqrt(n) * sum(prod_resids)

  prod_resid_resamp <- c()

  for(b in 1:B){
    # resampling X from X|Z
    resamp_X <- spacrt::dCRT_dist(n = n, fitted.val = X_on_Z_fit$fitted.values,
                                  fam = X_on_Z_fam)

    # compute the products of residuals for each resampled observation
    prod_resid_resamp[b] <- 1/sqrt(n) * sum((resamp_X - X_on_Z_fit$fitted.values)*
                                              (Y - Y_on_Z_fit$fitted.values))
  }

  # compute p-values
  left_side_p_value <- 1/(B+1) * (1 + sum(prod_resid_resamp <= test_stat))
  right_side_p_value <- 1/(B+1) * (1 + sum(prod_resid_resamp >= test_stat))

  # return test statistic, dCRT p-values, and related quantities
  return(list(test_stat = test_stat,
              left_side_p_value = left_side_p_value,
              right_side_p_value = right_side_p_value,
              both_side_p_value = 2*min(left_side_p_value, right_side_p_value),
              NB.disp.param = NB.disp.param))
}


#####################################################################################
#' The saddlepoint approximation to the dCRT.
#'
#' \code{spaCRT} is a function carrying out the saddlepoint approximation to the
#' dCRT based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param normalize A logical variable indicating whether the spaCRT is based on
#' the normalized test statistic.
#' @param return_cdf A logical variable indicating whether to return the CDF
#' @param fit_glm_X A logical variable indicating whether to use GLM for estimating `E[X|Z]`.
#' Default is TRUE. If set to FALSE, mean of X is used to estimate `E[X|Z]`.
#' @param fit_glm_Y A logical variable indicating whether to use GLM for estimating `E[Y|Z]`.
#' Default is TRUE. If set to FALSE, mean of Y is used to estimate `E[Y|Z]`.
#' @param aux_info_X_on_Z The auxiliary information that may be used for complex GLM regression
#' (For instance, when X_on_Z_fam = negative.binomial, the dispersion parameter should be provided).
#' @param aux_info_Y_on_Z The auxiliary information that may be used for complex GLM regression
#' (For instance, when Y_on_Z_fam = negative.binomial, the dispersion parameter should be provided).
#' @param R Upper bound of search space for the saddlepoint
#'
#' @examples
#' n <- 20; p <- 2; normalize <- FALSE; return_cdf <- FALSE
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' results <- spaCRT(data, X_on_Z_fam, Y_on_Z_fam, normalize)
#' results$test_stat
#' results$p_value
#'
#' @return A named list with fields \code{test_stat}, \code{left_side_p_value},
#' \code{right_side_p_value}, \code{both_side_p_value}, \code{cdf}, and
#' \code{gcm.default}.
#' Here, cdf is a function that takes in a value t and returns the
#' saddlepoint approximation to the CDF of the resampling distribution of the
#' test statistic evaluated at t. This function is returned only if
#' return_cdf == TRUE.
#' gcm.default returns TRUE if spacrt::GCM was employed due to the failure of spaCRT.
#'
#' @export
spaCRT <- function(data, X_on_Z_fam = NULL, Y_on_Z_fam = NULL,
                   normalize = FALSE, return_cdf = FALSE,
                   fit_glm_X = TRUE, fit_glm_Y = TRUE,
                   aux_info_X_on_Z = NULL, aux_info_Y_on_Z = NULL, R = 5) {

  if(is.null(X_on_Z_fam) | is.null(Y_on_Z_fam)){
    stop("X_on_Z_fam and Y_on_Z_fam can't be empty!")
  }

  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  # fit X on Z regression
  if(!fit_glm_X){
    X_on_Z_fit <- list(fitted.values = rep(mean(X), length(X)))
  }else{
    if(X_on_Z_fam == "negative.binomial"){
      X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z,
                              family = MASS::negative.binomial(aux_info_X_on_Z$theta_hat),
                              mustart = aux_info_X_on_Z$fitted_values))
    }else{
      X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z, family = X_on_Z_fam))
    }
  }

  # fit Y on Z regression
  if(!fit_glm_Y){
    Y_on_Z_fit <- list(fitted.values = rep(mean(Y), length(Y)))
  }else{
    if(Y_on_Z_fam == "negative.binomial"){
      Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z,
                              family = MASS::negative.binomial(aux_info_Y_on_Z$theta_hat),
                              mustart = aux_info_Y_on_Z$fitted_values))
      NB.disp.param <- aux_info_Y_on_Z$theta_hat
    }else{
      Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z, family = Y_on_Z_fam))
      NB.disp.param <- "Invalid request"
    }
  }

  W <- Y - Y_on_Z_fit$fitted.values
  P <- X_on_Z_fit$fitted.values

  # compute the products of residuals for each observation
  prod_resids <- (X - X_on_Z_fit$fitted.values) * W

  # compute the test statistic
  test_stat <- 1/sqrt(n) * sum(prod_resids)

  ##### SPA to CDF of T_n = S_n / sqrt(n)
  spa.cdf <- function(t, P = P, W = W, fam = X_on_Z_fam, R){
    n <- length(P)

    temp.gcm <- "NO"

    if(tryCatch(s.hat <- stats::uniroot(function(s){
      spacrt::d1.wcgf(s, P = P, W = W, fam) - sqrt(n)*t},
      lower = -R, upper = R, tol = .Machine$double.eps)$root,
      error = function(e) FALSE) == FALSE){

      if(tryCatch(s.hat <- stats::uniroot(function(s){
        spacrt::d1.wcgf(s, P = P, W = W, fam) - sqrt(n)*t},
        lower = -10*R, upper = 10*R, tol = .Machine$double.eps)$root,
        error = function(e) FALSE) == FALSE){

        if(tryCatch(s.hat <- stats::uniroot(function(s){
          spacrt::d1.wcgf(s, P = P, W = W, fam) - sqrt(n)*t},
          lower = -100*R, upper = 100*R, tol = .Machine$double.eps)$root,
          error = function(e) FALSE) == FALSE){

          if(tryCatch(s.hat <- stats::uniroot(function(s){
            spacrt::d1.wcgf(s, P = P, W = W, fam) - sqrt(n)*t},
            lower = -1000*R, upper = 1000*R, tol = .Machine$double.eps)$root,
            error = function(e) FALSE) == FALSE){

              temp.gcm <- "YES"
          }
        }
      }
    }

    if(temp.gcm == "NO"){
      r.hat <- sign(s.hat) * sqrt(2 * (n*s.hat*t/sqrt(n) -
                                         spacrt::wcgf(s = s.hat, P = P, W = W, fam)))

      F.hat <- stats::pnorm(r.hat) + stats::dnorm(r.hat) *
        (1/r.hat - 1/(s.hat*sqrt(spacrt::d2.wcgf(s = s.hat, P = P, W = W, fam))))

      return(F.hat)
    }else{
      return(NaN)
    }
  }

  # perform saddlepoint approximation
  p_value_opp <- suppressWarnings(spa.cdf(test_stat + 1/sqrt(n) * sum(P*W),
                                          P = P, W = W,
                                          fam = X_on_Z_fam,
                                          R = abs(R)))

  if(is.nan(p_value_opp) == TRUE){
    temp.gcm <- spacrt::GCM(data, X_on_Z_fam, Y_on_Z_fam,
                            fit_glm_X = fit_glm_X, fit_glm_Y = fit_glm_Y,
                            aux_info_X_on_Z = aux_info_X_on_Z,
                            aux_info_Y_on_Z = aux_info_Y_on_Z)

    # return test statistic, GCM p-values, and related quantities
    return(list(test_stat = temp.gcm$test_stat,
                left_side_p_value = temp.gcm$left_side_p_value,
                right_side_p_value = temp.gcm$right_side_p_value,
                both_side_p_value = temp.gcm$both_side_p_value,
                NB.disp.param = NB.disp.param,
                cdf = NULL,
                gcm.default = TRUE,
                nan.spacrt = is.nan(p_value_opp)))
  }else{
    if(p_value_opp < 0 | p_value_opp > 1){
      temp.gcm <- spacrt::GCM(data, X_on_Z_fam, Y_on_Z_fam,
                              fit_glm_X = fit_glm_X, fit_glm_Y = fit_glm_Y,
                              aux_info_X_on_Z = aux_info_X_on_Z,
                              aux_info_Y_on_Z = aux_info_Y_on_Z)

      # return test statistic, GCM p-values, and related quantities
      return(list(test_stat = temp.gcm$test_stat,
                  left_side_p_value = temp.gcm$left_side_p_value,
                  right_side_p_value = temp.gcm$right_side_p_value,
                  both_side_p_value = temp.gcm$both_side_p_value,
                  NB.disp.param = NB.disp.param,
                  cdf = NULL,
                  gcm.default = TRUE,
                  nan.spacrt = is.nan(p_value_opp)))
    }else{
      # return test statistic, spaCRT p-values, and related quantities
      return(list(test_stat = test_stat,
                  left_side_p_value = p_value_opp,
                  right_side_p_value = 1 - p_value_opp,
                  both_side_p_value = 2*min(c(p_value_opp, 1 - p_value_opp)),
                  NB.disp.param = NB.disp.param,
                  cdf = spa.cdf,
                  gcm.default = FALSE,
                  nan.spacrt = is.nan(p_value_opp)))
    }
  }
}


#####################################################################################
#' Score Test
#'
#' \code{score.test} is a function carrying out the saddlepoint approximation to the
#' dCRT based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param fit_glm_X A logical variable indicating whether to use GLM for estimating `E[X|Z]`.
#' Default is TRUE. If set to FALSE, mean of X is used to estimate `E[X|Z]`.
#' @param fit_glm_Y A logical variable indicating whether to use GLM for estimating `E[Y|Z]`.
#' Default is TRUE. If set to FALSE, mean of Y is used to estimate `E[Y|Z]`.
#' @param aux_info_X_on_Z The auxiliary information that may be used for complex GLM regression
#' (For instance, when X_on_Z_fam = negative.binomial, the dispersion parameter should be provided).
#' @param aux_info_Y_on_Z The auxiliary information that may be used for complex GLM regression
#' (For instance, when Y_on_Z_fam = negative.binomial, the dispersion parameter should be provided).
#'
#' @examples
#' n <- 50; p <- 2; normalize <- FALSE; return_cdf <- FALSE
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "negative.binomial"
#' results <- score.test(data, X_on_Z_fam, Y_on_Z_fam)
#' results$test_stat
#'
#' @return A named list with fields \code{test_stat} and \code{left_side_p_value},
#' \code{right_side_p_value} and \code{both_side_p_value}.
#'
#' @export
score.test <- function(data, X_on_Z_fam = NULL, Y_on_Z_fam = NULL,
                       fit_glm_X = TRUE, fit_glm_Y = TRUE,
                       aux_info_X_on_Z = NULL, aux_info_Y_on_Z = NULL){

  if(is.null(X_on_Z_fam) | is.null(Y_on_Z_fam)){
    stop("X_on_Z_fam and Y_on_Z_fam can't be empty!")
  }

  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  # fit X on Z regression
  if(!fit_glm_X){
    X_on_Z_fit <- list(fitted.values = rep(mean(X), length(X)))
  }else{
    if(X_on_Z_fam == "negative.binomial"){
      X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z,
                                  family = MASS::negative.binomial(aux_info_X_on_Z$theta_hat),
                                  mustart = aux_info_X_on_Z$fitted_values))
    }else{
      X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z, family = X_on_Z_fam))
    }
  }

  # fit Y on Z regression
  if(!fit_glm_Y){
    Y_on_Z_fit <- list(fitted.values = rep(mean(Y), length(Y)))
  }else{
    if(Y_on_Z_fam == "negative.binomial"){
      # First try to fit the model using glm.nb
      temp.result <- tryCatch({
        Y_on_Z_fit <- suppressWarnings(MASS::glm.nb(Y ~ Z))
        NB.disp.param <- Y_on_Z_fit$theta
        list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
      },
      error = function(e) {
        if(is.null(aux_info_Y_on_Z) == TRUE){
          aux_info_Y_on_Z <- spacrt::nb_precomp(list(Y = Y, Z = Z))
        }

        Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z,
                                family = MASS::negative.binomial(aux_info_Y_on_Z$theta_hat),
                                mustart = aux_info_Y_on_Z$fitted_values))
        NB.disp.param <- aux_info_Y_on_Z$theta_hat

        list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
      })
    }else if(Y_on_Z_fam == 'poisson'){
      if(is.null(aux_info_Y_on_Z) == TRUE){
        aux_info_Y_on_Z <- spacrt::nb_precomp(list(Y = Y, Z = Z))
      }

      Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z,
                                                family = stats::poisson(),
                                                mustart = aux_info_Y_on_Z$fitted_values))
      NB.disp.param <- "Invalid request"

      temp.result <- list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
    }else{
      Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z, family = Y_on_Z_fam))
      NB.disp.param <- "Invalid request"

      temp.result <- list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
    }
  }

  # perform score test
  test_stat <- statmod::glm.scoretest(fit = temp.result$Y_on_Z_fit, x2 = X)

  # return test statistic, score test p-values, and related quantities
  return(list(test_stat = test_stat,
              left_side_p_value = stats::pnorm(test_stat, lower.tail = TRUE),
              right_side_p_value = stats::pnorm(test_stat, lower.tail = FALSE),
              both_side_p_value = 2*stats::pnorm(abs(test_stat), lower.tail = FALSE),
              NB.disp.param = temp.result$NB.disp.param))
}




# spacrt - methods.R
