#####################################################################################
#' The generalized covariance measure test of Shah and Peters.
#'
#' \code{GCM} is a function carrying out the GCM test based on GLMs for X|Z and Y|Z.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#'
#' @return A named list with fields \code{test_stat} and \code{p_value}.
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
GCM <- function(data, X_on_Z_fam, Y_on_Z_fam) {
  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z
  # fit X on Z and Y on Z regressions
  X_on_Z_fit <- stats::glm(X ~ Z, family = X_on_Z_fam)
  Y_on_Z_fit <- stats::glm(Y ~ Z, family = Y_on_Z_fam)
  # compute the products of residuals for each observation
  prod_resids <- (X-X_on_Z_fit$fitted.values)*(Y-Y_on_Z_fit$fitted.values)
  # compute the test statistic
  n <- length(X)
  test_stat <- 1/sqrt(n)*sum(prod_resids)/stats::sd(prod_resids)
  # compute the p-value by comparing test statistic to normal distribution
  p_value <- 2*stats::pnorm(abs(test_stat), lower.tail = FALSE)
  # return test statistic and p-value
  list(test_stat = test_stat, p_value = p_value)
}


#####################################################################################
#' The distilled conditional randomization test.
#'
#' \code{dCRT} is a function carrying out the dCRT based on GLMs for X|Z and Y|Z.
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
#'
#' @return A named list with fields \code{test_stat}, \code{p_value}, and
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
dCRT <- function(data, X_on_Z_fam, Y_on_Z_fam, B, normalize = FALSE, return_resamples = FALSE) {

  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  X_on_Z_fit <- stats::glm(X ~ Z, family = X_on_Z_fam)
  Y_on_Z_fit <- stats::glm(Y ~ Z, family = Y_on_Z_fam)

  prod_resids <- (X - X_on_Z_fit$fitted.values)*(Y - Y_on_Z_fit$fitted.values)

  test_stat <- 1/sqrt(n) * sum(prod_resids) #/stats::sd(prod_resids)

  resamp_X <- list()
  prod_resid_resamp <- c()

  for(b in 1:B){
    resamp_X[[b]] <- cbind(1,Z) %*% as.matrix(X_on_Z_fit$coefficients) %>%
                        sapply(function(val) rbinom(1, 1, expit(val)))

    prod_resid_resamp[b] <- 1/sqrt(n) * sum((resamp_X[[b]] - X_on_Z_fit$fitted.values)*
                                                  (Y - Y_on_Z_fit$fitted.values))
  }

  p_value <- 1/(B+1) * (1 + sum(prod_resid_resamp >= test_stat))

  return(list(test_stat = test_stat, p_value = p_value))
}


#####################################################################################
#' The saddlepoint approximation to the dCRT.
#'
#' \code{spaCRT} is a function carrying out the saddlepoint approximation to the
#' dCRT based on GLMs for X|Z and Y|Z.
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
#' @return A named list with fields \code{test_stat}, \code{p_value}, and
#' \code{cdf}. Here, cdf is a function that takes in a value t and returns the
#' saddlepoint approximation to the CDF of the resampling distribution of the
#' test statistic evaluated at t. This function is returned only if
#' return_cdf == TRUE.
#'
#' @export
spaCRT <- function(data, X_on_Z_fam, Y_on_Z_fam, normalize, return_cdf) {

  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  X_on_Z_fit <- stats::glm(X ~ Z, family = X_on_Z_fam)
  Y_on_Z_fit <- stats::glm(Y ~ Z, family = Y_on_Z_fam)


}
