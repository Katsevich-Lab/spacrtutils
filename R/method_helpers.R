
# method helpers will go here (e.g. solving saddlepoint equations)

#' \code{dCRT_dist} is a function computing the cumulant generating function (CGF) of
#' distributions, multiplied by a weight function, from the GLM family
#'
#' @param n The point where the CGF will be computed.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param fitted.val A vector containing the fitted parameter values by
#' fitting a GLM to X on Z.
#' @return Simulated data from an appropriate distribution.
#'
#' @examples
#' dCRT_dist(n = n, fitted.val = X_on_Z_fit$fitted.values, fam = 'poisson')
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



# method helpers will go here (e.g. solving saddlepoint equations)

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

#' \code{d1.wcgf} is a function computing the derivative of the cumulant generating
#' function (CGF) of distributions, multiplied by a weight function, from the GLM family
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

  if(fam == 'binomial') return(sum((W*P*exp(s*W)) / (exp(s*W)*P + 1 - P)))
  if(fam == 'gaussian') return()
  if(fam == 'Gamma') return()
  if(fam == 'inverse.gaussian') return()
  if(fam == 'poisson') return(sum(P*exp(s*W)*W))
  if(fam == 'quasi') return()
  if(fam == 'quasibinomial') return()
  if(fam == 'quasipoisson') return()


}

#' \code{d2.wcgf} is a function computing the hessian of the cumulant generating
#' function (CGF) of distributions, multiplied by a weight function, from the GLM family
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
    return(sum((W^2*P*Q*exp(s*W)) / (exp(s*W)*P + Q)^2))
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














