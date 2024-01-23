#' \code{GCM_f} is a simulation helper function carrying out the GCM test from the spacrt package.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @return A named list with fields \code{left_side_p_value}, \code{right_side_p_value}, and \code{both_side_p_value}.
#'
#' @examples
#' n <- 100; p <- 5; gamma_1 = 1; beta_0 = -3; beta_1 = 1; theta = 1; gamma_0 = -2
#' expit <- function(theta)(exp(theta)/(1+exp(theta)))
#' Z <- as.matrix(rnorm(n = n, mean = 0, sd = 1))
#' X <- rbinom(n = n, size = 1, prob = expit(gamma_0 + gamma_1*Z))
#' Y <- MASS::rnegbin(n = n, mu = exp(beta_0 + beta_1*Z), theta = theta)
#' data <- list(X = X, Y = Y, Z = Z)
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' results <- GCM_f(data, X_on_Z_fam, Y_on_Z_fam)
#' @export
GCM_f <- function(data,
                  X_on_Z_fam = 'binomial',
                  Y_on_Z_fam = 'negative.binomial'){

  aux_info_Y_on_Z <- spacrt::nb_precomp(list(Y = data$Y, Z = data$Z))

  results.GCM <- spacrt::GCM(data, X_on_Z_fam, Y_on_Z_fam,
                             aux_info_Y_on_Z = aux_info_Y_on_Z)

  return(list(p.left = results.GCM$left_side_p_value,
              p.right = results.GCM$right_side_p_value,
              p.both = results.GCM$both_side_p_value))
}


#' \code{dCRT_f} is a simulation helper function carrying out the dCRT test from the spacrt package.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param B Number of resamples
#' @return A named list with fields \code{left_side_p_value}, \code{right_side_p_value}, and \code{both_side_p_value}.
#'
#' @examples
#' n <- 100; p <- 5; gamma_1 = 1; beta_0 = -3; beta_1 = 1; theta = 1; gamma_0 = -2
#' B <- 1000
#' expit <- function(theta)(exp(theta)/(1+exp(theta)))
#' Z <- as.matrix(rnorm(n = n, mean = 0, sd = 1))
#' X <- rbinom(n = n, size = 1, prob = expit(gamma_0 + gamma_1*Z))
#' Y <- MASS::rnegbin(n = n, mu = exp(beta_0 + beta_1*Z), theta = theta)
#' data <- list(X = X, Y = Y, Z = Z)
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' results <- dCRT_f(data, X_on_Z_fam, Y_on_Z_fam, B)
#' @export
dCRT_f <- function(data,
                   X_on_Z_fam = 'binomial',
                   Y_on_Z_fam = 'negative.binomial', B = 10000){

  aux_info_Y_on_Z <- spacrt::nb_precomp(list(Y = data$Y, Z = data$Z))

  results.dCRT <- spacrt::dCRT(data, X_on_Z_fam, Y_on_Z_fam,
                               B = 10000,
                               normalize = FALSE, return_resamples = FALSE,
                               # test_side = 'right',
                               aux_info_Y_on_Z = aux_info_Y_on_Z)

  return(list(p.left = results.dCRT$left_side_p_value,
              p.right = results.dCRT$right_side_p_value,
              p.both = results.dCRT$both_side_p_value))
}


#' \code{spaCRT_f} is a simulation helper function carrying out the spaCRT test from the spacrt package.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @return A named list with fields \code{left_side_p_value}, \code{right_side_p_value}, and \code{both_side_p_value}.
#'
#' @examples
#' n <- 100; p <- 5; gamma_1 = 1; beta_0 = -3; beta_1 = 1; theta = 1; gamma_0 = -2
#' expit <- function(theta)(exp(theta)/(1+exp(theta)))
#' Z <- as.matrix(rnorm(n = n, mean = 0, sd = 1))
#' X <- rbinom(n = n, size = 1, prob = expit(gamma_0 + gamma_1*Z))
#' Y <- MASS::rnegbin(n = n, mu = exp(beta_0 + beta_1*Z), theta = theta)
#' data <- list(X = X, Y = Y, Z = Z)
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' results <- spaCRT_f(data, X_on_Z_fam, Y_on_Z_fam)
#' @export
spaCRT_f <- function(data,
                     X_on_Z_fam = 'binomial',
                     Y_on_Z_fam = 'negative.binomial'){

  aux_info_Y_on_Z <- spacrt::nb_precomp(list(Y = data$Y, Z = data$Z))

  results.spaCRT <- spacrt::spaCRT(data, X_on_Z_fam, Y_on_Z_fam,
                                   normalize = F, return_cdf = F,
                                   fit_glm_X = TRUE, fit_glm_Y = TRUE,
                                   aux_info_Y_on_Z = aux_info_Y_on_Z)

  return(list(p.left = results.spaCRT$left_side_p_value,
              p.right = results.spaCRT$right_side_p_value,
              p.both = results.spaCRT$both_side_p_value))
}


#' \code{scoretest_glm_nb_f} is a simulation helper function carrying out the score test from the spacrt package.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @return A named list with fields \code{left_side_p_value}, \code{right_side_p_value}, and \code{both_side_p_value}.
#'
#' @examples
#' n <- 100; p <- 5; gamma_1 = 1; beta_0 = -3; beta_1 = 1; theta = 1; gamma_0 = -2
#' expit <- function(theta)(exp(theta)/(1+exp(theta)))
#' Z <- as.matrix(rnorm(n = n, mean = 0, sd = 1))
#' X <- rbinom(n = n, size = 1, prob = expit(gamma_0 + gamma_1*Z))
#' Y <- MASS::rnegbin(n = n, mu = exp(beta_0 + beta_1*Z), theta = theta)
#' data <- list(X = X, Y = Y, Z = Z)
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' results <- scoretest_glm_nb_f(data, X_on_Z_fam, Y_on_Z_fam)
#' @export
scoretest_glm_nb_f <- function(data,
                               X_on_Z_fam = 'binomial',
                               Y_on_Z_fam = 'negative.binomial'){

  aux_info_Y_on_Z <- spacrt::nb_precomp(list(Y = data$Y, Z = data$Z))

  results.scoretest <- suppressWarnings(
    spacrt::score.test(data, X_on_Z_fam, Y_on_Z_fam,
                       fit_glm_X = TRUE, fit_glm_Y = TRUE,
                       aux_info_Y_on_Z = aux_info_Y_on_Z)
  )

  return(list(p.left = results.scoretest$left_side_p_value,
              p.right = results.scoretest$right_side_p_value,
              p.both = results.scoretest$both_side_p_value))
}
