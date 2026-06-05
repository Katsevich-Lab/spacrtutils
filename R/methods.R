#####################################################################################
#' The generalized covariance measure test of Shah and Peters.
#'
#' \code{GCM_internal} is a function carrying out the GCM test based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param data
#'    A (non-empty) named list with fields \code{X} (an nx1 vector for the predictor
#'    variable of interest), \code{Y} (an nx1 response vector), and \code{Z} (an nxp matrix
#'    of covariates).
#' @param X_on_Z_fam
#'    The GLM family for the regression of X on Z (values can be \code{binomial},
#'    \code{poisson}, etc).
#' @param Y_on_Z_fam
#'    The GLM family for the regression of Y on Z (values can be \code{binomial},
#'    \code{poisson}, \code{negative.binomial}, etc).
#' @param fitting_X_on_Z
#'    The fitting method for the regression X on Z (values can be \code{glm} (default),
#'    \code{random_forest}, or \code{NA}).
#' @param fitting_Y_on_Z
#'    The fitting method for the regression Y on Z (values can be \code{glm} (default),
#'    \code{random_forest}, or \code{NA}).
#' @param fit_vals_X_on_Z_own
#'    Vector of fitted values for X on Z in case the user's custom method.
#' @param fit_vals_Y_on_Z_own
#'    Vector of fitted values for Y on Z in case the user's custom method.
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
#' data <- list(X = X, Y = Y, Z = Z)
#'
#' # fit both models via GLM
#' res.GCM.1 <- GCM_internal(data,
#'                           X_on_Z_fam = "binomial", Y_on_Z_fam = "poisson",
#'                           fitting_X_on_Z = 'glm', fitting_Y_on_Z = 'glm',
#'                           fit_vals_X_on_Z_own = NULL, fit_vals_Y_on_Z_own = NULL)
#' res.GCM.1
#'
#' # custom fit for Y|Z, explicit method/family only for X|Z
#' user_fit_Y <- glm(Y ~ Z,
#'                   family = poisson(),
#'                   data = data.frame(Y = Y, Z = Z))$fitted.values |> unname()
#'
#' res.GCM.2 <- GCM_internal(data,
#'                           X_on_Z_fam = "binomial", Y_on_Z_fam = NULL,
#'                           fitting_X_on_Z = 'random_forest', fitting_Y_on_Z = NULL,
#'                           fit_vals_X_on_Z_own = NULL, fit_vals_Y_on_Z_own = user_fit_Y)
#' res.GCM.2
#'
#' @export
GCM_internal <- function(data, X_on_Z_fam, Y_on_Z_fam,
                         fitting_X_on_Z = 'glm',
                         fitting_Y_on_Z = 'glm',
                         fit_vals_X_on_Z_own = NULL,
                         fit_vals_Y_on_Z_own = NULL){

   results <- spacrt::GCM(X = data$X, Y = data$Y, Z = data$Z,
                          family = list(XZ = X_on_Z_fam,
                                        YZ = Y_on_Z_fam),
                          method = list(XZ = fitting_X_on_Z,
                                        YZ = fitting_Y_on_Z),
                          fitted = list(XZ = fit_vals_X_on_Z_own,
                                        YZ = fit_vals_Y_on_Z_own),
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
#' data <- list(X = X, Y = Y, Z = Z)
#'
#' # fit both models via GLM
#' res.dCRT.1 <- dCRT_internal(data,
#'                             X_on_Z_fam = "binomial", Y_on_Z_fam = "poisson",
#'                             fitting_X_on_Z = 'glm', fitting_Y_on_Z = 'glm',
#'                             fit_vals_X_on_Z_own = NULL, fit_vals_Y_on_Z_own = NULL)
#' res.dCRT.1
#'
#' # custom fit for Y|Z, explicit method/family only for X|Z
#' user_fit_Y <- glm(Y ~ Z,
#'                   family = poisson(),
#'                   data = data.frame(Y = Y, Z = Z))$fitted.values |> unname()
#'
#' res.dCRT.2 <- dCRT_internal(data,
#'                             X_on_Z_fam = "binomial", Y_on_Z_fam = NULL,
#'                             fitting_X_on_Z = 'random_forest', fitting_Y_on_Z = NULL,
#'                             fit_vals_X_on_Z_own = NULL, fit_vals_Y_on_Z_own = user_fit_Y,
#'                             B = 10000)
#' res.dCRT.2
#'
#' @export
dCRT_internal <- function(data, X_on_Z_fam, Y_on_Z_fam,
                          fitting_X_on_Z = 'glm',
                          fitting_Y_on_Z = 'glm',
                          fit_vals_X_on_Z_own = NULL,
                          fit_vals_Y_on_Z_own = NULL,
                          B = 5000) {

   results <- spacrt::dCRT(X = data$X, Y = data$Y, Z = data$Z,
                           family = list(XZ = X_on_Z_fam,
                                         YZ = Y_on_Z_fam),
                           method = list(XZ = fitting_X_on_Z,
                                         YZ = fitting_Y_on_Z),
                           fitted = list(XZ = fit_vals_X_on_Z_own,
                                         YZ = fit_vals_Y_on_Z_own),
                           alternative = "less",
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
#' data <- list(X = X, Y = Y, Z = Z)
#'
#' # fit both models via GLM
#' res.spaCRT.1 <- spaCRT_internal(data,
#'                                 X_on_Z_fam = "binomial", Y_on_Z_fam = "poisson",
#'                                 fitting_X_on_Z = 'glm', fitting_Y_on_Z = 'glm',
#'                                 fit_vals_X_on_Z_own = NULL, fit_vals_Y_on_Z_own = NULL)
#' res.spaCRT.1
#'
#' # custom fit for Y|Z, explicit method/family only for X|Z
#' dtrain <- xgboost::xgb.DMatrix(data = Z, label = Y)
#' model.Y <- xgboost::xgboost(data = dtrain,
#'                             objective = "count:poisson",
#'                             nrounds = 100, verbose = 0)
#' user_fit_Y <- stats::predict(model.Y, newdata = Z)
#'
#' res.spaCRT.2 <- spaCRT_internal(data,
#'                                 X_on_Z_fam = "binomial", Y_on_Z_fam = NULL,
#'                                 fitting_X_on_Z = 'random_forest', fitting_Y_on_Z = NULL,
#'                                 fit_vals_X_on_Z_own = NULL, fit_vals_Y_on_Z_own = user_fit_Y)
#' res.spaCRT.2
#'
#' @export
spaCRT_internal <- function(data, X_on_Z_fam, Y_on_Z_fam,
                            fitting_X_on_Z = 'glm',
                            fitting_Y_on_Z = 'glm',
                            fit_vals_X_on_Z_own = NULL,
                            fit_vals_Y_on_Z_own = NULL) {

   results <- spacrt::spaCRT(X = data$X, Y = data$Y, Z = data$Z,
                             family = list(XZ = X_on_Z_fam,
                                           YZ = Y_on_Z_fam),
                             method = list(XZ = fitting_X_on_Z,
                                           YZ = fitting_Y_on_Z),
                             fitted = list(XZ = fit_vals_X_on_Z_own,
                                           YZ = fit_vals_Y_on_Z_own),
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
#' @param data
#'    A (non-empty) named list with fields \code{X} (an nx1 vector for the predictor
#'    variable of interest), \code{Y} (an nx1 response vector), and \code{Z} (an nxp matrix
#'    of covariates).
#' @param X_on_Z_fam
#'    The GLM family for the regression of X on Z (values can be \code{binomial},
#'    \code{poisson}, etc).
#' @param Y_on_Z_fam
#'    The GLM family for the regression of Y on Z (values can be \code{binomial},
#'    \code{poisson}, \code{negative.binomial}, etc).
#'
#' @examples
#' n <- 200; p <- 4
#' set.seed(1234)
#' X <- rbinom(n = n, size = 1, prob = 0.3)
#' Y <- rpois(n = n, lambda = 1)
#' Z <- matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p)
#' data <- list(X = X, Y = Y, Z = Z)
#'
#' res.score.1 <- score.test(data,
#'                           X_on_Z_fam = "binomial", Y_on_Z_fam = "negative.binomial")
#' res.score.1
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
      temp.res <- tryCatch({
         Y_on_Z_fit <- suppressWarnings(MASS::glm.nb(Y ~ Z))
         NB.disp.param <- Y_on_Z_fit$theta
         list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
      },
      error = function(e) {
         aux_info_Y_on_Z <- nb_precomp(Y,Z)

         # switch to Poisson if negative binomial regression fails
         Y_on_Z_fit <- tryCatch({
            suppressWarnings(stats::glm(Y ~ Z,
                                        family = MASS::negative.binomial(aux_info_Y_on_Z$theta_hat),
                                        mustart = aux_info_Y_on_Z$fitted_values))
         },
         error = function(e) {
            fit <- suppressWarnings(stats::glm(Y ~ Z,
                                               family = stats::poisson(),
                                               mustart = aux_info_Y_on_Z$fitted_values))
            fit$family <- MASS::negative.binomial(aux_info_Y_on_Z$theta_hat)
            fit
         })

         NB.disp.param <- aux_info_Y_on_Z$theta_hat

         list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
      })
   }else if(Y_on_Z_fam == 'poisson'){
      aux_info_Y_on_Z <- nb_precomp(Y,Z)

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


#####################################################################################
#' SAIGE-Bernoulli randomize-X saddlepoint test
#'
#' \code{SAIGE_Bern_internal} carries out the unified-vision SAIGE variant that
#' randomizes \code{X} from a logistic \code{X|Z} fit and residualizes \code{Y}
#' linearly on \code{Z} via GLS with weights \eqn{W = \hat\mu_x (1 - \hat\mu_x)}
#' (the Bernoulli conditional variance, i.e. the X-model IRLS metric of
#' Section 3 of the gls-validity note). Lugannani-Rice is applied to the
#' Bernoulli CGF on weights \eqn{a_i = Y_i - Z_d \hat\gamma_y}.
#'
#' This is the canonical "SAIGE in the unified randomize-X framework" when
#' \code{X} is binary: the observed statistic is invariant to the linear
#' Y-fit (Section 3, Proposition 1), and the SPA reference variance matches
#' the sampling variance under the logistic-GLM score equation. The Y outcome
#' need not be binary -- it enters only as a fixed weight.
#'
#' @param data
#'    A (non-empty) named list with fields \code{X} (an nx1 binary vector),
#'    \code{Y} (an nx1 response vector, possibly count), and \code{Z} (an nxp
#'    matrix of covariates).
#' @param fitting_X_on_Z
#'    Fitting method for the logistic X|Z model. Currently \code{"glm"} (default).
#' @param fit_vals_X_on_Z_own
#'    Optional vector of user-supplied fitted values \eqn{\hat\mu_x} (e.g. from
#'    a GAM or random forest). If supplied, \code{fitting_X_on_Z} is ignored.
#'
#' @return A named list with fields \code{test_stat}, \code{p.left}
#' (left-sided p-value), \code{p.right} (right-sided p-value), and \code{p.both}
#' (two-sided p-value).
#'
#' @examples
#' set.seed(1)
#' n <- 1000
#' Z <- matrix(stats::rnorm(2 * n), n, 2)
#' X <- stats::rbinom(n, 1, plogis(-3 + Z[, 1] + Z[, 2]))
#' Y <- MASS::rnegbin(n, mu = exp(-2 + 0.5 * Z[, 1] + 0.5 * Z[, 2]), theta = 0.5)
#' res <- SAIGE_Bern_internal(list(X = X, Y = Y, Z = Z))
#' res
#'
#' @export
SAIGE_Bern_internal <- function(data,
                                fitting_X_on_Z = "glm",
                                fit_vals_X_on_Z_own = NULL) {

   if (is.null(fit_vals_X_on_Z_own) && !identical(fitting_X_on_Z, "glm")) {
      stop("Unsupported fitting_X_on_Z: ", fitting_X_on_Z,
           ". Supply fit_vals_X_on_Z_own instead.")
   }

   X <- data$X; Y <- data$Y; Z <- data$Z

   res <- tryCatch(withCallingHandlers({
      if (!is.null(fit_vals_X_on_Z_own)) {
         mu_x <- as.numeric(fit_vals_X_on_Z_own)
      } else {
         fit_x <- suppressWarnings(stats::glm(X ~ Z, family = stats::binomial()))
         mu_x  <- as.numeric(stats::fitted(fit_x))
      }

      Omega   <- mu_x * (1 - mu_x)
      Zd      <- cbind(1, Z)
      WZ      <- Zd * Omega
      gamma_y <- solve(crossprod(Zd, WZ), crossprod(WZ, Y))
      mu_y    <- as.numeric(Zd %*% gamma_y)
      a       <- Y - mu_y
      T_      <- sum(a * (X - mu_x))
      p       <- .bin_spa_pvalues(w = a, mu = mu_x, S_obs = T_)
      list(test_stat = T_, p.left = p$p.left, p.right = p$p.right,
           p.both = 2 * min(p$p.left, p$p.right))
   }, warning = function(w) invokeRestart("muffleWarning")),
   error = function(e) list(test_stat = NA_real_, p.left = NA_real_,
                            p.right = NA_real_, p.both = NA_real_))
   res
}


#####################################################################################
#' SAIGE-NB randomize-X saddlepoint test
#'
#' \code{SAIGE_NB_internal} models \code{X|Z} as negative binomial and
#' randomizes \eqn{\tilde X_i \sim \mathrm{NegBin}(\hat\mu_{x,i}, \hat r)} for
#' the SPA reference. \code{Y} is residualized linearly on \code{Z} via GLS
#' with the NB-GLM IRLS Fisher-information weight
#' \eqn{W = \hat\mu_x/(1+\hat\mu_x/\hat r) = \hat\mu_x^2/V}, where
#' \eqn{V = \hat\mu_x + \hat\mu_x^2/\hat r} is the NB conditional variance.
#' The same Fisher-information factor \eqn{1/(1+\hat\mu_x/\hat r) = \hat\mu_x/V}
#' is applied per observation to the Y-residual, giving the test statistic
#' \deqn{T = \sum_i \frac{(Y_i - Z_d \hat\gamma_y)\,(X_i - \hat\mu_{x,i})}{1 + \hat\mu_{x,i}/\hat r}.}
#' Lugannani-Rice is applied to the NB CGF with the same reweighted weights.
#'
#' This is the NB analog of \code{\link{SAIGE_Bern_internal}}. By
#' construction (X residualization in the NB Fisher-information metric, plus
#' per-observation Fisher-information reweighting of the Y residual) the
#' statistic and its SPA studentization are bit-for-bit the NB Rao score
#' test of \code{glm.nb} (\code{statmod::glm.scoretest}) with the asymptotic
#' \eqn{N(0,1)} reference replaced by the NB SPA. For Bernoulli+logit
#' (canonical link) the Fisher-information factor collapses to \eqn{1}, and
#' one recovers the unweighted residual product of
#' \code{\link{SAIGE_Bern_internal}}.
#' Like the other randomize-X methods (\code{\link{spaCRT_internal}},
#' \code{\link{SAIGE_Bern_internal}}), only \code{X|Z}-side parameters are
#' exposed; \code{Y} enters as a fixed weight.
#'
#' @param data
#'    A (non-empty) named list with fields \code{X} (an nx1 count vector),
#'    \code{Y} (an nx1 response vector), and \code{Z} (an nxp matrix of
#'    covariates).
#' @param fitting_X_on_Z
#'    Fitting method for the NB X|Z model. Currently \code{"glm"} (default),
#'    via \code{MASS::glm.nb}.
#' @param fit_vals_X_on_Z_own
#'    Optional vector of user-supplied fitted values \eqn{\hat\mu_x}. If
#'    supplied, \code{size_hat_own} must also be supplied.
#' @param size_hat_own
#'    Optional scalar with the NB dispersion / size estimate \eqn{\hat r} for
#'    the X|Z model. Required when \code{fit_vals_X_on_Z_own} is supplied.
#'
#' @return A named list with fields \code{test_stat}, \code{p.left}, \code{p.right},
#' \code{p.both}, and \code{NB.disp.param} (the NB \eqn{\hat r} used).
#'
#' @examples
#' set.seed(1)
#' n <- 1000
#' Z <- matrix(stats::rnorm(2 * n), n, 2)
#' X <- MASS::rnegbin(n, mu = exp(-1 + Z[, 1] + Z[, 2]), theta = 5)
#' Y <- stats::rbinom(n, 1, plogis(-3 + Z[, 1] + Z[, 2]))
#' res <- SAIGE_NB_internal(list(X = X, Y = Y, Z = Z))
#' res
#'
#' @export
SAIGE_NB_internal <- function(data,
                              fitting_X_on_Z = "glm",
                              fit_vals_X_on_Z_own = NULL,
                              size_hat_own = NULL) {

   # argument validation (loud, outside the numerical tryCatch)
   if (!is.null(fit_vals_X_on_Z_own) && is.null(size_hat_own)) {
      stop("size_hat_own must be supplied when fit_vals_X_on_Z_own is given.")
   }
   if (is.null(fit_vals_X_on_Z_own) && !identical(fitting_X_on_Z, "glm")) {
      stop("Unsupported fitting_X_on_Z: ", fitting_X_on_Z,
           ". Supply fit_vals_X_on_Z_own + size_hat_own instead.")
   }

   X <- data$X; Y <- data$Y; Z <- data$Z

   res <- tryCatch(withCallingHandlers({
      if (!is.null(fit_vals_X_on_Z_own)) {
         mu_x     <- as.numeric(fit_vals_X_on_Z_own)
         size_hat <- as.numeric(size_hat_own)
      } else {
         fit_x    <- suppressWarnings(MASS::glm.nb(X ~ Z))
         mu_x     <- as.numeric(stats::fitted(fit_x))
         size_hat <- fit_x$theta
      }

      # NB-GLM IRLS info weight: w = mu_x / (1 + mu_x/size_hat) = mu_x^2 / V_x,
      # the Fisher information per obs for the linear predictor under NB+log.
      # Residualizing X on Z in this metric makes the X residual coincide with
      # the NB Rao-score projection (Schur complement). The same factor then
      # reweights (Y - mu_y) to produce the Rao score numerator; combined,
      # the SAIGE-NB statistic and its SPA studentization are bit-for-bit the
      # NB Rao score with SPA-corrected tails.
      info_w  <- 1 / (1 + mu_x / size_hat)
      Omega   <- mu_x * info_w
      Zd      <- cbind(1, Z)
      WZ      <- Zd * Omega
      gamma_y <- solve(crossprod(Zd, WZ), crossprod(WZ, Y))
      mu_y    <- as.numeric(Zd %*% gamma_y)
      a       <- (Y - mu_y) * info_w
      T_      <- sum(a * (X - mu_x))
      p       <- .nb_spa_pvalues(w = a, mu = mu_x, size_ = size_hat, S_obs = T_)
      list(test_stat = T_, p.left = p$p.left, p.right = p$p.right,
           p.both = 2 * min(p$p.left, p$p.right),
           NB.disp.param = size_hat)
   }, warning = function(w) invokeRestart("muffleWarning")),
   error = function(e) list(test_stat = NA_real_, p.left = NA_real_,
                            p.right = NA_real_, p.both = NA_real_,
                            NB.disp.param = NA_real_))
   res
}


#####################################################################################
#' SAIGE-NB randomize-X saddlepoint test (simple Poisson + Pearson-theta variant)
#'
#' Identical in structure to \code{\link{SAIGE_NB_internal}}, but the
#' nuisance fit of \code{X|Z} uses the lightweight Poisson + Pearson-theta
#' estimator \code{nb_precomp} from \pkg{spacrt} rather than the iterative
#' alternating IRLS / profile-likelihood routine in \code{MASS::glm.nb}.
#' Concretely \code{nb_precomp} fits a Poisson GLM of \code{X} on \code{Z}
#' and then derives \eqn{\hat r} from the Pearson chi-square statistic of
#' the Poisson fit. The Poisson fitted values play the role of
#' \eqn{\hat\mu_x} (they are the MLE of the conditional mean under both
#' Poisson and NB log-linear models), and \eqn{\hat r} is the
#' method-of-moments NB dispersion. All downstream steps (NB-IRLS GLS
#' projection of \code{Y} on \code{Z}, Fisher-information reweighting,
#' Lugannani--Rice on the NB CGF) match \code{\link{SAIGE_NB_internal}}.
#'
#' This variant trades a small efficiency loss (Pearson-moments \eqn{\hat r}
#' is slightly noisier than the NB MLE) for substantially faster fitting
#' and no risk of \code{glm.nb} convergence pathologies.
#'
#' @inheritParams SAIGE_NB_internal
#'
#' @return A named list with the same fields as
#' \code{\link{SAIGE_NB_internal}}.
#'
#' @examples
#' set.seed(1)
#' n <- 1000
#' Z <- matrix(stats::rnorm(2 * n), n, 2)
#' X <- MASS::rnegbin(n, mu = exp(-1 + Z[, 1] + Z[, 2]), theta = 5)
#' Y <- stats::rbinom(n, 1, plogis(-3 + Z[, 1] + Z[, 2]))
#' SAIGE_NB_simple_internal(list(X = X, Y = Y, Z = Z))
#'
#' @export
SAIGE_NB_simple_internal <- function(data,
                                     fitting_X_on_Z = "glm",
                                     fit_vals_X_on_Z_own = NULL,
                                     size_hat_own = NULL) {

   if (!is.null(fit_vals_X_on_Z_own) && is.null(size_hat_own)) {
      stop("size_hat_own must be supplied when fit_vals_X_on_Z_own is given.")
   }
   if (is.null(fit_vals_X_on_Z_own) && !identical(fitting_X_on_Z, "glm")) {
      stop("Unsupported fitting_X_on_Z: ", fitting_X_on_Z,
           ". Supply fit_vals_X_on_Z_own + size_hat_own instead.")
   }

   X <- data$X; Y <- data$Y; Z <- data$Z

   res <- tryCatch(withCallingHandlers({
      if (!is.null(fit_vals_X_on_Z_own)) {
         mu_x     <- as.numeric(fit_vals_X_on_Z_own)
         size_hat <- as.numeric(size_hat_own)
      } else {
         pre      <- nb_precomp(X, Z)
         mu_x     <- as.numeric(pre$fitted_values)
         size_hat <- as.numeric(pre$theta_hat)
      }

      info_w  <- 1 / (1 + mu_x / size_hat)
      Omega   <- mu_x * info_w
      Zd      <- cbind(1, Z)
      WZ      <- Zd * Omega
      gamma_y <- solve(crossprod(Zd, WZ), crossprod(WZ, Y))
      mu_y    <- as.numeric(Zd %*% gamma_y)
      a       <- (Y - mu_y) * info_w
      T_      <- sum(a * (X - mu_x))
      p       <- .nb_spa_pvalues(w = a, mu = mu_x, size_ = size_hat, S_obs = T_)
      list(test_stat = T_, p.left = p$p.left, p.right = p$p.right,
           p.both = 2 * min(p$p.left, p$p.right),
           NB.disp.param = size_hat)
   }, warning = function(w) invokeRestart("muffleWarning")),
   error = function(e) list(test_stat = NA_real_, p.left = NA_real_,
                            p.right = NA_real_, p.both = NA_real_,
                            NB.disp.param = NA_real_))
   res
}


#####################################################################################
#' Generic GLM sign-flipping score test (Hemerik-Goeman-Finos, 2020)
#'
#' \code{signflip_score_internal} carries out the Rao-score test for
#' \eqn{H_0:\beta=0} in the GLM
#' \deqn{g(\mathbb{E}[Y_i \mid X_i, Z_i]) = \beta X_i + Z_i^\top \gamma,}
#' using the sign-flipping reference distribution of
#' Hemerik--Goeman--Finos (2020, JRSSB) and the Lugannani--Rice saddlepoint
#' approximation of the conditional CGF
#' \eqn{K(t) = \sum_i \log\cosh(t \nu^*_i)} (see
#' \code{.signflip_spa_pvalues}). The implementation is generic across GLM
#' families exposed through the \pkg{stats} family interface (Poisson,
#' Binomial, Gaussian, Gamma, inverse Gaussian, quasi-families), and handles
#' negative binomial through \code{MASS::glm.nb} / \code{MASS::negative.binomial}.
#'
#' For a GLM with link \eqn{g}, variance function \eqn{V} and IRLS Fisher
#' weight \eqn{W_i = (d\mu_i/d\eta_i)^2 / V(\mu_i)}, the effective score for
#' \eqn{\beta} after profiling out \eqn{\gamma} is
#' \deqn{\nu^*_i = \tilde X_i \, H_i \, (Y_i - \hat\mu_i), \quad
#'       H_i = (d\mu_i/d\eta_i)/V(\hat\mu_i), \quad
#'       \tilde X_i = X_i - Z_i^\top \hat\delta,}
#' with \eqn{\hat\delta = (Z^\top W Z)^{-1} Z^\top W X}. For canonical links
#' (Poisson log, Binomial logit, Gaussian identity) \eqn{H_i = 1} and the
#' score collapses to \eqn{\tilde X_i (Y_i - \hat\mu_i)}; for the negative
#' binomial with log link, \eqn{H_i = 1/(1 + \hat\mu_i/\hat r)}.
#'
#' The observed statistic is \eqn{T = \sum_i \nu^*_i}. Tail probabilities
#' come from Lugannani--Rice applied to the exact sign-flipping CGF
#' \eqn{K(t) = \sum_i \log\cosh(t\,\nu^*_i)}; see
#' \code{.signflip_spa_pvalues}.
#'
#' Implementation reference:
#' \code{writeups/sign-flipping-implementaiton/sign_flipping_spa.tex}.
#'
#' @param data
#'    A (non-empty) named list with fields \code{X} (an nx1 vector for the
#'    predictor variable of interest), \code{Y} (an nx1 response vector),
#'    and \code{Z} (an nxp matrix of covariates).
#' @param Y_on_Z_fam
#'    String giving the GLM family for the regression of Y on Z under
#'    \eqn{H_0:\beta=0}. One of \code{"poisson"}, \code{"binomial"},
#'    \code{"gaussian"}, \code{"Gamma"}, \code{"inverse.gaussian"}, or
#'    \code{"negative.binomial"}. Quasi-families (e.g.\ \code{"quasipoisson"})
#'    are supported via \code{do.call}.
#' @param fitting_Y_on_Z
#'    Fitting method for the Y on Z regression. Currently only \code{"glm"}
#'    (default) is supported when fitting is auto-performed; custom fits
#'    should be passed through \code{fit_vals_Y_on_Z_own}.
#' @param fit_vals_Y_on_Z_own
#'    Optional vector of user-supplied fitted means \eqn{\hat\mu_i}. When
#'    supplied, the function does not refit \eqn{Y} on \eqn{Z}. For
#'    \code{Y_on_Z_fam = "negative.binomial"}, \code{size_hat_own} must also
#'    be supplied.
#' @param size_hat_own
#'    Optional scalar NB dispersion / size \eqn{\hat r}. Required when
#'    \code{fit_vals_Y_on_Z_own} is supplied together with
#'    \code{Y_on_Z_fam = "negative.binomial"}.
#'
#' @return A named list with fields \code{test_stat}, \code{p.left}
#' (left-sided), \code{p.right} (right-sided), \code{p.both} (two-sided),
#' \code{spa.success} (TRUE if the saddlepoint equation was solved), and
#' \code{NB.disp.param} (the NB \eqn{\hat r} used, or \code{NA} for non-NB
#' families).
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' Z <- matrix(stats::rnorm(2 * n), n, 2)
#' X <- stats::rnorm(n)
#' Y <- stats::rpois(n, lambda = exp(0.5 + 0.3 * Z[, 1] + 0.2 * Z[, 2]))
#' data <- list(X = X, Y = Y, Z = Z)
#'
#' # Poisson coefficient test, auto-fit Y ~ Z.
#' signflip_score_internal(data, Y_on_Z_fam = "poisson")
#'
#' # User-supplied fitted means.
#' mu_y <- stats::fitted(stats::glm(Y ~ Z, family = stats::poisson()))
#' signflip_score_internal(data, Y_on_Z_fam = "poisson",
#'                         fit_vals_Y_on_Z_own = mu_y)
#'
#' @export
signflip_score_internal <- function(data,
                                    Y_on_Z_fam,
                                    fitting_Y_on_Z = "glm",
                                    fit_vals_Y_on_Z_own = NULL,
                                    size_hat_own = NULL) {

   is_nb <- identical(Y_on_Z_fam, "negative.binomial")
   if (!is.null(fit_vals_Y_on_Z_own) && is_nb && is.null(size_hat_own)) {
      stop("size_hat_own must be supplied when fit_vals_Y_on_Z_own is given ",
           "for Y_on_Z_fam = 'negative.binomial'.")
   }
   if (is.null(fit_vals_Y_on_Z_own) && !identical(fitting_Y_on_Z, "glm")) {
      stop("Unsupported fitting_Y_on_Z: ", fitting_Y_on_Z,
           ". Supply fit_vals_Y_on_Z_own instead.")
   }

   X <- data$X; Y <- data$Y; Z <- data$Z
   Zd <- cbind(1, Z)

   res <- tryCatch(withCallingHandlers({

      # --- Fit Y ~ Z under H0 to obtain (mu_y, eta_y, family object). ---
      if (!is.null(fit_vals_Y_on_Z_own)) {
         mu_y     <- as.numeric(fit_vals_Y_on_Z_own)
         size_hat <- if (is_nb) as.numeric(size_hat_own) else NA_real_
         fam      <- if (is_nb) MASS::negative.binomial(size_hat)
                     else       do.call(Y_on_Z_fam, list())
         eta_y    <- fam$linkfun(mu_y)
      } else if (is_nb) {
         fit_y    <- suppressWarnings(MASS::glm.nb(Y ~ Z))
         mu_y     <- as.numeric(stats::fitted(fit_y))
         eta_y    <- as.numeric(stats::predict(fit_y, type = "link"))
         size_hat <- fit_y$theta
         fam      <- MASS::negative.binomial(size_hat)
      } else {
         fam      <- do.call(Y_on_Z_fam, list())
         fit_y    <- suppressWarnings(stats::glm(Y ~ Z, family = fam))
         mu_y     <- as.numeric(stats::fitted(fit_y))
         eta_y    <- as.numeric(stats::predict(fit_y, type = "link"))
         size_hat <- NA_real_
      }

      # --- GLM score machinery: H = (dmu/deta)/V(mu), W = (dmu/deta)^2/V(mu). ---
      mu_eta <- fam$mu.eta(eta_y)
      V_mu   <- fam$variance(mu_y)
      H      <- mu_eta / V_mu                 # canonical link => H = 1
      W      <- mu_eta * H                    # = (dmu/deta)^2 / V(mu)

      # --- Residualize X on Z in the W-metric (profile out gamma). ---
      WZ      <- Zd * W
      delta   <- solve(crossprod(Zd, WZ), crossprod(WZ, X))
      X_tilde <- as.numeric(X - Zd %*% delta)

      # --- Effective score and observed sign-flip statistic. ---
      nu_star <- X_tilde * H * (Y - mu_y)
      T_      <- sum(nu_star)

      # --- Lugannani-Rice SPA on K(t) = sum log cosh(t nu_i). ---
      p       <- .signflip_spa_pvalues(nu = nu_star, S_obs = T_)
      ok      <- is.finite(p$p.left) && is.finite(p$p.right)
      list(test_stat     = T_,
           p.left        = p$p.left,
           p.right       = p$p.right,
           p.both        = 2 * min(p$p.left, p$p.right),
           spa.success   = ok,
           NB.disp.param = size_hat)
   }, warning = function(w) invokeRestart("muffleWarning")),
   error = function(e) list(test_stat     = NA_real_,
                            p.left        = NA_real_,
                            p.right       = NA_real_,
                            p.both        = NA_real_,
                            spa.success   = FALSE,
                            NB.disp.param = NA_real_))
   res
}
