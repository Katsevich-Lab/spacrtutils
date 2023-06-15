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
#' @return A named list with fields \code{test_stat} and \code{p_value}.
#'
#' @export
spaCRT <- function(data, X_on_Z_fam, Y_on_Z_fam, normalize) {
}
