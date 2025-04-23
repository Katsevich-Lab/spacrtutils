#####################################################################################
#' The generalized covariance measure test of Shah and Peters.
#'
#' \code{GCM_internal} is a function carrying out the GCM test based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param data A (non-empty) named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z} (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param fitting_X_on_Z The fitting method for the regression X on Z.
#' (values can be \code{glm} (default), \code{rf}, \code{prob_forest}, or \code{own})
#' @param fitting_Y_on_Z The fitting method for the regression Y on Z.
#' (values can be \code{glm} (default), \code{rf}, \code{prob_forest}, or \code{own})
#' @param fit_vals_X_on_Z_own Vector of fitted values for X on Z in case the user's custom method.
#' Works only if fitting_X_on_Z = 'own'.
#' @param fit_vals_Y_on_Z_own Vector of fitted values for Y on Z in case the user's custom method.
#' Works only if fitting_Y_on_Z = 'own'.
#'
#' @return A named list with fields \code{test_stat}, \code{p.left} (Left-sided p-value),
#' \code{p.right} (Right-sided p-value), \code{p.both} (Two-sided p-value).
#'
#' @examples
#' n <- 20; p <- 2
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' GCM_internal(data, X_on_Z_fam, Y_on_Z_fam,
#'              fitting_X_on_Z = 'rf',
#'              fitting_Y_on_Z = 'glm')
#'
#' @export
GCM_internal <- function(data, X_on_Z_fam, Y_on_Z_fam,
                         fitting_X_on_Z = 'glm',
                         fitting_Y_on_Z = 'glm',
                         fit_vals_X_on_Z_own = NULL,
                         fit_vals_Y_on_Z_own = NULL){

   results <- spacrt::GCM(data, X_on_Z_fam, Y_on_Z_fam,
                          fitting_X_on_Z = 'glm',
                          fitting_Y_on_Z = 'glm',
                          fit_vals_X_on_Z_own = NULL,
                          fit_vals_Y_on_Z_own = NULL,
                          alternative = 'less')

   test_stat <- results$test_stat

   # return test statistic, GCM p-values, and related quantities
   return(list(test_stat = test_stat,
               p.left = stats::pnorm(test_stat, lower.tail = TRUE),
               p.right = stats::pnorm(test_stat, lower.tail = FALSE),
               p.both = 2*stats::pnorm(abs(test_stat), lower.tail = FALSE)))
}



#####################################################################################
#' The distilled conditional randomization test.
#'
#' \code{dCRT_internal} is a function carrying out the dCRT based on GLMs for `X|Z` and `Y|Z`.
#'
#' @inheritParams GCM_internal
#' @param B The number of resamples to draw (Default value is 2000).
#'
#' @return A named list with fields \code{test_stat}, \code{p.left} (Left-sided p-value),
#' \code{p.right} (Right-sided p-value), \code{p.both} (Two-sided p-value).
#'
#' @examples
#' n <- 80; p <- 2
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' dCRT_internal(data, X_on_Z_fam, Y_on_Z_fam,
#'               fitting_X_on_Z = 'rf',
#'               fitting_Y_on_Z = 'glm',
#'               B = 2000)
#'
#' @export
dCRT_internal <- function(data, X_on_Z_fam, Y_on_Z_fam,
                          fitting_X_on_Z = 'glm',
                          fitting_Y_on_Z = 'glm',
                          fit_vals_X_on_Z_own = NULL,
                          fit_vals_Y_on_Z_own = NULL,
                          B = 2000) {

   results <- spacrt::dCRT(data, X_on_Z_fam, Y_on_Z_fam,
                           fitting_X_on_Z = 'glm',
                           fitting_Y_on_Z = 'glm',
                           fit_vals_X_on_Z_own = NULL,
                           fit_vals_Y_on_Z_own = NULL,
                           alternative = 'less',
                           B = B)

   test_stat <- results$test_stat

   # compute p-values
   p.left <- results$p_value |> unname()
   p.right <- 1 - p.left + 1/(B+1)
   p.both <- 2 * min(c(p.left, p.right))

   # return test statistic, dCRT p-values, and related quantities
   return(list(test_stat = test_stat,
               p.left = p.left,
               p.right = p.right,
               p.both = p.both))
}


#####################################################################################
#' The saddlepoint approximation to the dCRT.
#'
#' \code{spaCRT_internal} is a function carrying out the saddlepoint approximation to the
#' dCRT based on GLMs for `X|Z` and `Y|Z`.
#'
#' @inheritParams GCM_internal
#'
#' @return A named list with fields \code{test_stat}, \code{p.left} (Left-sided p-value),
#' \code{p.right} (Right-sided p-value), \code{p.both} (Two-sided p-value), and \code{spa.success}.
#' \code{spa.success} returns TRUE if the saddlepoint equation could be solved; otherwise,
#' the backup method (GCM) was employed.
#' .
#'
#' @examples
#' n <- 50; p <- 4
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rbinom(n = n, size = 1, prob = 0.7),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "binomial"
#' spaCRT_internal(data, X_on_Z_fam, Y_on_Z_fam,
#'                 fitting_X_on_Z = 'rf',
#'                 fitting_Y_on_Z = 'glm')
#'
#' @export
spaCRT_internal <- function(data, X_on_Z_fam, Y_on_Z_fam,
                            fitting_X_on_Z = 'glm',
                            fitting_Y_on_Z = 'glm',
                            fit_vals_X_on_Z_own = NULL,
                            fit_vals_Y_on_Z_own = NULL) {

   results <- spacrt::spaCRT(data, X_on_Z_fam, Y_on_Z_fam,
                             fitting_X_on_Z = 'glm',
                             fitting_Y_on_Z = 'glm',
                             fit_vals_X_on_Z_own = NULL,
                             fit_vals_Y_on_Z_own = NULL,
                             alternative = 'less')

   test_stat <- results$test_stat

   # compute p-values
   p.left <- results$p_value |> unname()
   p.right <- 1 - p.left
   p.both <- 2 * min(c(p.left, p.right))

   return(list(test_stat = test_stat,
               p.left = p.left,
               p.right = p.right,
               p.both = p.both,
               spa.success = results$spa.success))
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
#'
#' @examples
#' n <- 50; p <- 2; normalize <- FALSE; return_cdf <- FALSE
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "negative.binomial"
#' score.test(data, X_on_Z_fam, Y_on_Z_fam)
#'
#' @return A named list with fields \code{test_stat}, \code{p.left} (Left-sided p-value),
#' \code{p.right} (Right-sided p-value), \code{p.both} (Two-sided p-value), and
#' \code{NB.disp.param}.
#'
#' @export
score.test <- function(data, X_on_Z_fam, Y_on_Z_fam){

   # extract (X,Y,Z) from inputted data
   X <- data$X; Y <- data$Y; Z <- data$Z
   n <- length(X)

   # fit X on Z regression
   X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z, family = X_on_Z_fam))

   # fit Y on Z regression
   if(Y_on_Z_fam == "negative.binomial"){
      # First try to fit the model using glm.nb
      temp.result <- tryCatch({
         Y_on_Z_fit <- suppressWarnings(MASS::glm.nb(Y ~ Z))
         NB.disp.param <- Y_on_Z_fit$theta
         list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
      },
      error = function(e) {
         aux_info_Y_on_Z <- nb_precomp(list(Y = Y, Z = Z))

         Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z,
                                                   family = MASS::negative.binomial(aux_info_Y_on_Z$theta_hat),
                                                   mustart = aux_info_Y_on_Z$fitted_values))
         NB.disp.param <- aux_info_Y_on_Z$theta_hat

         list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
      })
   }else if(Y_on_Z_fam == 'poisson'){
      aux_info_Y_on_Z <- nb_precomp(list(Y = Y, Z = Z))

      Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z,
                                                family = stats::poisson(),
                                                mustart = aux_info_Y_on_Z$fitted_values))
      NB.disp.param <- NA

      temp.result <- list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
   }else{
      Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z, family = Y_on_Z_fam))
      NB.disp.param <- NA

      temp.result <- list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
   }

   # perform score test
   test_stat <- statmod::glm.scoretest(fit = temp.result$Y_on_Z_fit, x2 = X)

   # return test statistic, score test p-values, and related quantities
   return(list(test_stat = test_stat,
               p.left = stats::pnorm(test_stat, lower.tail = TRUE),
               p.right = stats::pnorm(test_stat, lower.tail = FALSE),
               p.both = 2*stats::pnorm(abs(test_stat), lower.tail = FALSE),
               NB.disp.param = temp.result$NB.disp.param))
}

