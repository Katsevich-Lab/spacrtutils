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
}

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
dCRT <- function(data, X_on_Z_fam, Y_on_Z_fam, B, normalize, return_resamples) {
}

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
#'
#' @examples
#' n <- 20; p <- 2; normalize <- FALSE
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' results <- spaCRT(data, X_on_Z_fam, Y_on_Z_fam, normalize)
#' results$test_stat
#' results$p_value
#'
#' @return A named list with fields \code{test_stat} and \code{p_value}.
#'
#' @export
spaCRT <- function(data, X_on_Z_fam, Y_on_Z_fam, normalize) {
}
