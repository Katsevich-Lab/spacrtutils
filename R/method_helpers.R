#####################################################################################
#' \code{spa.cdf} SPA to CDF of T_n = S_n / sqrt(n)
#'
#' @param n The point where the CGF will be computed.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param W <- Y - Y_on_Z_fit$fitted.values
#' @param P <- X_on_Z_fit$fitted.values

#' @return Simulated data from an appropriate distribution.
#'
#' @examples
#' n <- 20; p <- 2; normalize <- FALSE; return_cdf <- FALSE
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X <- data$X; Y <- data$Y; Z <- data$Z
#' X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z, family = "binomial"))
#' Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z, family = "poisson"))
#' W <- Y - Y_on_Z_fit$fitted.values
#' P <- X_on_Z_fit$fitted.values
#' prod_resids <- (X - X_on_Z_fit$fitted.values) * W
#' test_stat <- 1/sqrt(n) * sum(prod_resids)
#' spa.cdf(t = test_stat + sum(P*W)/sqrt(n), P = P, W = W, fam = X_on_Z_fam, R = 10)
#'
#' @export
spa.cdf <- function(t, P = P, W = W, fam = X_on_Z_fam, R, max_expansions = 10){
  n <- length(P)

  # temp.gcm <- "NO"
  current_lower <- -abs(R)
  current_upper <- abs(R)
  success_uniroot <- FALSE

  for (i in seq_len(max_expansions)) {
    tryCatch({
      s.hat <- uniroot(
        f = function(s) spacrt::d1.wcgf(s, P = P, W = W, fam) - sqrt(n)*t,
        lower = current_lower, upper = current_upper, tol = .Machine$double.eps)$root

      success_uniroot <<- TRUE
      break
    }, error = function(e) {
                  message(sprintf("Attempt %d failed, expanding interval...", i))
                  expansion_factor <- ifelse(i <= max_expansions/2, 2, 10)
                  current_lower <<- current_lower * expansion_factor
                  current_upper <<- current_upper * expansion_factor
                }
    )

    if(success_uniroot == TRUE) break
  }

  if(success_uniroot == TRUE){
    r.hat <- sign(s.hat) * sqrt(2 * (n * s.hat * t / sqrt(n) -
                                       spacrt::wcgf(s = s.hat, P = P, W = W, fam)))

    F.hat <- stats::pnorm(r.hat) + stats::dnorm(r.hat) *
      (1/r.hat - 1/(s.hat*sqrt(spacrt::d2.wcgf(s = s.hat, P = P, W = W, fam))))

    return(F.hat)
  }else{
    return(NaN)
  }
}



#####################################################################################
#' \code{dCRT_dist} is a function that returns simulated data from an appropriate
#' distribution depending on a specified  GLM family
#'
#' @param n The point where the CGF will be computed.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param fitted.val A vector containing the fitted parameter values by
#' fitting a GLM to X on Z.
#' @return Simulated data from an appropriate distribution.
#'
#' @examples
#' dCRT_dist(n = 5, fitted.val = c(1,2,1,1,1), fam = 'poisson')
#'
#' @export
dCRT_dist <- function(n, fitted.val, fam){

  if(fam == 'binomial') return(stats::rbinom(n = n, size = 1, prob = fitted.val))
  if(fam == 'gaussian') return()
  if(fam == 'Gamma') return()
  if(fam == 'inverse.gaussian') return()
  if(fam == 'poisson') return(stats::rpois(n = n, lambda = fitted.val))
  if(fam == 'quasi') return()
  if(fam == 'quasibinomial') return()
  if(fam == 'quasipoisson') return()
}


#####################################################################################
#' \code{wcgf} is a function computing the cumulant generating function (CGF) of
#' distributions, multiplied by a weight function, from the GLM family
#'
#' @param s The point where the CGF will be computed.
#' @param P A vector containing the parameter values of the family of distributions.
#' @param W A vector containing the weights.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @return CGF of the weighted distribution evaluated at \code{s}.
#'
#' @examples
#' wcgf(s = 5, P = c(1,1,2), W = c(1,1,1), fam = 'poisson')
#'
#' @export
wcgf <- function(s, P, W, fam){

  if(fam == 'binomial') return(sum(log(exp(s*W)*P + 1 - P)))
  if(fam == 'gaussian') return()
  if(fam == 'Gamma') return()
  if(fam == 'inverse.gaussian') return()
  if(fam == 'poisson') return(sum(P*(exp(s*W) - 1)))
  if(fam == 'quasi') return()
  if(fam == 'quasibinomial') return()
  if(fam == 'quasipoisson') return()
}


#####################################################################################
#' \code{d1.wcgf} is a function computing the derivative of the cumulant generating
#' function (CGF) of distributions, multiplied by a weight function, from GLM family
#'
#' @param s The point where the CGF will be computed.
#' @param P A vector containing the parameter values of the family of distributions.
#' @param W A vector containing the weights.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @return CGF of the weighted distribution evaluated at \code{s}.
#'
#' @examples
#' d1.wcgf(s = 5, P = c(1,1,2), W = c(1,1,1), fam = 'poisson')
#'
#' @export
d1.wcgf <- function(s, P, W, fam){

  # if(fam == 'binomial') return(sum((W*P*exp(s*W)) / (exp(s*W)*P + 1 - P)))
  if(fam == 'binomial') return(sum((W*P) / (P + (1 - P) * exp(-s*W))))
  if(fam == 'gaussian') return()
  if(fam == 'Gamma') return()
  if(fam == 'inverse.gaussian') return()
  if(fam == 'poisson') return(sum(P*exp(s*W)*W))
  if(fam == 'quasi') return()
  if(fam == 'quasibinomial') return()
  if(fam == 'quasipoisson') return()
}


#####################################################################################
#' \code{d2.wcgf} is a function computing the hessian of the cumulant generating
#' function (CGF) of distributions, multiplied by a weight function, from GLM family
#'
#' @param s The point where the CGF will be computed.
#' @param P A vector containing the parameter values of the family of distributions.
#' @param W A vector containing the weights.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @return CGF of the weighted distribution evaluated at \code{s}.
#'
#' @examples
#' d2.wcgf(s = 5, P = c(1,1,2), W = c(1,1,1), fam = 'poisson')
#'
#' @export
d2.wcgf <- function(s, P, W, fam){

  if(fam == 'binomial'){
    Q <- 1 - P
    # return(sum((W^2*P*Q*exp(s*W)) / (exp(s*W)*P + Q)^2))
    return(sum((W^2*P*Q) / (exp(s*W)*P^2 + 2 * P * Q + Q^2 * exp(-s*W))))
  }

  if(fam == 'gaussian'){
    Q <- 1 - P
    return()
  }

  if(fam == 'Gamma'){
    Q <- 1 - P
    return()
  }

  if(fam == 'inverse.gaussian'){
    Q <- 1 - P
    return()
  }

  if(fam == 'poisson'){
    return(sum(P*exp(s*W)*W^2))
  }

  if(fam == 'quasi'){
    Q <- 1 - P
    return()
  }

  if(fam == 'quasibinomial'){
    Q <- 1 - P
    return()
  }

  if(fam == 'quasipoisson'){
    Q <- 1 - P
    return()
  }
}


#####################################################################################
#' A function computing the dispersion parameter in negative binomial regression
#'
#' @param data A list containing the response Y and covariate Z
#'
#' @return a list containing the Poisson model fitted values and estimate for dispersion
#' @export
nb_precomp <- function(data){

  Y <- data$Y; Z <- data$Z

  pois_fit <- stats::glm.fit(y = Y, x = Z, family = stats::poisson())

  theta_hat <- sceptre:::estimate_theta(
    y = Y,
    mu = pois_fit$fitted.values,
    dfr = pois_fit$df.residual,
    limit = 50,
    eps = (.Machine$double.eps)^(1/4)
  )[[1]]

  return(list(fitted_values = pois_fit$fitted.values, theta_hat = theta_hat))
}




# spacrt - method_helpers.R
