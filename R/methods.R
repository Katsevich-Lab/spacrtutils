#####################################################################################
#' The generalized covariance measure test of Shah and Peters.
#'
#' \code{GCM_internal} is a function carrying out the GCM test based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param X
#'   Numeric vector of length \eqn{n}, representing the predictor variable.
#' @param Y
#'   Numeric vector of length \eqn{n}, representing the response variable.
#' @param Z
#'   Numeric matrix with \eqn{n} rows and \eqn{p} columns, representing covariates.
#' @param family
#'   Named list with elements `XZ` and `YZ` specifying the model family for \eqn{X \mid Z} and
#'   \eqn{Y \mid Z} for each model. Each list element must be a string (e.g. `"binomial"`,
#'   `"poisson"`). Ignored for any model where you supply your own fitted values via `fitted`.
#' @param method
#'   Named list with elements `XZ` and `YZ` that selects the model-fitting method to use
#'   for each model. Each element must be a string (e.g. `"glm"`, `"random_forest"`).
#'   Ignored for any model where you supply your own fitted values via `fitted`.
#' @param fitted
#'   Named list with elements `XZ` and `YZ` of user-supplied fitted values
#'   (numeric vectors of length n). For each non-NULL element, that model is
#'   treated as custom and neither `method` nor `family` is used.
#'
#' @return A named list with fields \code{test_stat}, \code{p.left} (Left-sided p-value),
#' \code{p.right} (Right-sided p-value), \code{p.both} (Two-sided p-value).
#'
#' @examples
#' n <- 200; p <- 4
#' set.seed(1234)
#' X <- rbinom(n = n, size = 1, prob = 0.3)
#' Y <- rpois(n = n, lambda = 1)
#' Z <- matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p)
#'
#' # fit both models via GLM
#' res.GCM.1 <- GCM_internal(X = X, Y = Y, Z = Z,
#'                           family = list(XZ = "binomial", YZ = "poisson"),
#'                           method = list(XZ = "glm", YZ = "glm"))
#' res.GCM.1
#'
#' # custom fit for Y|Z, explicit method/family only for X|Z
#' user_fit_Y <- glm(Y ~ Z,
#'                   family = poisson(),
#'                   data = data.frame(Y = Y, Z = Z))$fitted.values |> unname()
#'
#' res.GCM.2 <- GCM_internal(X = X, Y = Y, Z = Z,
#'                           family = list(XZ = "binomial"),
#'                           method = list(XZ = "random_forest"),
#'                           fitted = list(XZ = NULL, YZ = user_fit_Y))
#' res.GCM.2
#'
#' @export
GCM_internal <- function(X, Y, Z,
                         family,
                         method,
                         fitted = list(XZ = NULL, YZ = NULL)){

   results <- spacrt::GCM(X = X, Y = Y, Z = Z,
                          family, method, fitted,
                          alternative = "less")

   test_stat <- results$test_stat

   # compute p-values
   p.left <- stats::pnorm(test_stat, lower.tail = TRUE)
   p.right <- stats::pnorm(test_stat, lower.tail = FALSE)
   p.both <- 2*stats::pnorm(abs(test_stat), lower.tail = FALSE)

   # return test statistic and GCM p-values
   return(list(test_stat = test_stat,
               p.left = p.left,
               p.right = p.right,
               p.both = p.both))
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
#' n <- 200; p <- 4
#' set.seed(1234)
#' X <- rbinom(n = n, size = 1, prob = 0.3)
#' Y <- rpois(n = n, lambda = 1)
#' Z <- matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p)
#'
#' # fit both models via GLM
#' res.dCRT.1 <- dCRT_internal(X = X, Y = Y, Z = Z,
#'                             family = list(XZ = "binomial", YZ = "poisson"),
#'                             method = list(XZ = "glm", YZ = "glm"))
#' res.dCRT.1
#'
#' # custom fit for Y|Z, explicit method/family only for X|Z
#' user_fit_Y <- glm(Y ~ Z,
#'                   family = poisson(),
#'                   data = data.frame(Y = Y, Z = Z))$fitted.values |> unname()
#'
#' res.dCRT.2 <- dCRT_internal(X = X, Y = Y, Z = Z,
#'                             family = list(XZ = "binomial"),
#'                             method = list(XZ = "random_forest"),
#'                             fitted = list(XZ = NULL, YZ = user_fit_Y),
#'                             B = 10000)
#' res.dCRT.2
#'
#' @export
dCRT_internal <- function(X, Y, Z,
                          family,
                          method,
                          fitted = list(XZ = NULL, YZ = NULL),
                          B = 5000) {

   results <- spacrt::dCRT(X = X, Y = Y, Z = Z,
                           family, method, fitted,
                           alternative = 'less',
                           B = B)

   test_stat <- results$test_stat

   # compute p-values
   p.left <- results$p_value |> unname()
   p.right <- 1 - p.left + 1/(B+1)
   p.both <- 2 * min(c(p.left, p.right))

   # return test statistic, dCRT p-values
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
#' n <- 200; p <- 4
#' set.seed(1234)
#' X <- rbinom(n = n, size = 1, prob = 0.3)
#' Y <- rpois(n = n, lambda = 1)
#' Z <- matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p)
#'
#' # fit both models via GLM
#' res.spaCRT.1 <- spaCRT_internal(X = X, Y = Y, Z = Z,
#'                                 family = list(XZ = "binomial", YZ = "poisson"),
#'                                 method = list(XZ = "glm", YZ = "glm"))
#' res.spaCRT.1
#'
#' # custom fit for Y|Z, explicit method/family only for X|Z
#' dtrain <- xgboost::xgb.DMatrix(data = Z, label = Y)
#' model.Y <- xgboost::xgboost(data = dtrain,
#'                             objective = "count:poisson",
#'                             nrounds = 100, verbose = 0)
#' user_fit_Y <- stats::predict(model.Y, newdata = Z)
#'
#' res.spaCRT.2 <- spaCRT_internal(X = X, Y = Y, Z = Z,
#'                                 family = list(XZ = "binomial"),
#'                                 method = list(XZ = "random_forest"),
#'                                 fitted = list(XZ = NULL, YZ = user_fit_Y))
#' res.spaCRT.2
#'
#' @export
spaCRT_internal <- function(X, Y, Z,
                            family,
                            method,
                            fitted = list(XZ = NULL, YZ = NULL)) {

   results <- spacrt::spaCRT(X = X, Y = Y, Z = Z,
                             family, method, fitted,
                             alternative = "less")

   test_stat <- results$test_stat

   # compute p-values
   p.left <- results$p_value |> unname()
   p.right <- 1 - p.left
   p.both <- 2 * min(c(p.left, p.right))

   # return test statistic, spaCRT p-values, and spa.success
   return(list(test_stat = test_stat,
               p.left = p.left,
               p.right = p.right,
               p.both = p.both,
               spa.success = results$spa.success))
}


#####################################################################################
#' Score Test
#'
#' \code{score.test} is a function carrying out the score test based on GLMs for
#' `X|Z` and `Y|Z`.
#'
#' @param X
#'   Numeric vector of length \eqn{n}, representing the predictor variable.
#' @param Y
#'   Numeric vector of length \eqn{n}, representing the response variable.
#' @param Z
#'   Numeric matrix with \eqn{n} rows and \eqn{p} columns, representing covariates.
#' @param family
#'   Named list with elements `XZ` and `YZ` specifying the model family for \eqn{X \mid Z} and
#'   \eqn{Y \mid Z} for each model. Each list element must be a string (e.g. `"binomial"`,
#'   `"poisson"`).
#'
#' @examples
#' n <- 200; p <- 4
#' set.seed(1234)
#' X <- rbinom(n = n, size = 1, prob = 0.3)
#' Y <- rpois(n = n, lambda = 1)
#' Z <- matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p)
#'
#' res.score.1 <- score.test(X = X, Y = Y, Z = Z,
#'                           family = list(XZ = "binomial", YZ = "negative.binomial"))
#' res.score.1
#'
#' @return A named list with fields \code{test_stat}, \code{p.left} (Left-sided p-value),
#' \code{p.right} (Right-sided p-value), \code{p.both} (Two-sided p-value), and
#' \code{NB.disp.param}.
#'
#' @export
score.test <- function(X, Y, Z,
                       family){

   # extract (X,Y,Z) from inputted data
   n <- length(X)

   X_on_Z_fam <- family$XZ; Y_on_Z_fam <- family$YZ

   # fit X on Z regression
   X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z, family = X_on_Z_fam))

   # fit Y on Z regression
   if(Y_on_Z_fam == "negative.binomial"){
      # First try to fit the model using glm.nb
      temp.res <- tryCatch({
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

      temp.res <- list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
   }else{
      Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z, family = Y_on_Z_fam))
      NB.disp.param <- NA

      temp.res <- list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
   }

   # perform score test
   test_stat <- statmod::glm.scoretest(fit = temp.res$Y_on_Z_fit, x2 = X)

   # return test statistic, score test p-values, and related quantities
   return(list(test_stat = test_stat,
               p.left = stats::pnorm(test_stat, lower.tail = TRUE),
               p.right = stats::pnorm(test_stat, lower.tail = FALSE),
               p.both = 2*stats::pnorm(abs(test_stat), lower.tail = FALSE),
               NB.disp.param = temp.res$NB.disp.param))
}

