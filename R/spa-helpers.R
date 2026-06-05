#' Bernoulli Lugannani-Rice saddlepoint p-values
#'
#' Computes one-sided saddlepoint tail probabilities for
#'   \eqn{S = \sum_i w_i (X_i - \mu_i)},  with \eqn{X_i \sim \mathrm{Bernoulli}(\mu_i)}
#' independently, evaluated at observed \eqn{S = S_{obs}}.
#'
#' @param w  numeric vector of (fixed) weights.
#' @param mu numeric vector of Bernoulli means in \eqn{(0,1)}.
#' @param S_obs scalar, the observed centered sum.
#' @return list with \code{p.left} and \code{p.right}.
#' @keywords internal
.bin_spa_pvalues <- function(w, mu, S_obs) {
   K  <- function(t) sum(log1p(mu * (exp(w * t) - 1)) - t * w * mu)
   K1 <- function(t) { e <- exp(w * t); p <- mu * e / (1 - mu + mu * e); sum(w * (p - mu)) }
   K2 <- function(t) { e <- exp(w * t); p <- mu * e / (1 - mu + mu * e); sum(w^2 * p * (1 - p)) }
   var_S <- K2(0)
   if (!is.finite(var_S) || var_S <= 0) {
      return(list(p.left = NA_real_, p.right = NA_real_))
   }
   fallback <- function() {
      z <- S_obs / sqrt(var_S)
      list(p.left = stats::pnorm(z), p.right = stats::pnorm(z, lower.tail = FALSE))
   }
   if (abs(S_obs) < 1e-8 * sqrt(var_S)) return(fallback())
   tb <- 0.9 * 700 / max(abs(w))
   t_hat <- tryCatch(
      stats::uniroot(function(t) K1(t) - S_obs, lower = -tb, upper = tb, tol = 1e-10)$root,
      error = function(e) NA_real_)
   if (!is.finite(t_hat)) return(fallback())
   r_lr <- sign(t_hat) * sqrt(pmax(2 * (t_hat * S_obs - K(t_hat)), 0))
   v_lr <- t_hat * sqrt(K2(t_hat))
   p.right <- stats::pnorm(r_lr, lower.tail = FALSE) + stats::dnorm(r_lr) * (1 / r_lr - 1 / v_lr)
   p.left  <- stats::pnorm(r_lr, lower.tail = TRUE)  - stats::dnorm(r_lr) * (1 / r_lr - 1 / v_lr)
   list(p.left  = min(max(p.left,  .Machine$double.xmin), 1),
        p.right = min(max(p.right, .Machine$double.xmin), 1))
}


#' Sign-flipping Lugannani-Rice saddlepoint p-values
#'
#' Computes one-sided saddlepoint tail probabilities for
#'   \eqn{T = \sum_i g_i \nu_i}, with \eqn{g_i \stackrel{\text{iid}}{\sim}
#'   \mathrm{Unif}\{-1,+1\}} (Rademacher signs) and \eqn{\nu_i \in \mathbb{R}}
#' fixed. This is the reference distribution of the Hemerik-Goeman-Finos
#' (2020) sign-flipping score test conditional on observed scores.
#'
#' The conditional cumulant generating function is
#' \deqn{K(t) = \sum_{i=1}^n \log\cosh(t \nu_i),}
#' an entire, strictly convex, even function of \eqn{t}. The Lugannani--Rice
#' formula is applied with saddle \eqn{\hat t} solving \eqn{K'(\hat t) =
#' S_{\mathrm{obs}}}.
#'
#' Algebraically this is the \eqn{\pi_i \equiv 1/2} slice of the Bernoulli
#' CGF used by \code{.bin_spa_pvalues}: with \eqn{w_i = 2\nu_i} and
#' \eqn{\mu_i = 1/2}, \code{.bin_spa_pvalues(w = 2*nu, mu = rep(0.5, n),
#' S_obs = T)} returns identical tail probabilities (used as a cross-check in
#' the testthat suite).
#'
#' @param nu numeric vector of (fixed) sign-flip scores. Typically the
#'   effective GLM score contributions \eqn{\nu^*_i = \tilde x_i H_i (Y_i -
#'   \hat\mu_i)} from \code{\link{signflip_score_internal}}.
#' @param S_obs scalar, the observed value of \eqn{T = \sum_i \nu_i}.
#' @return list with \code{p.left} and \code{p.right}.
#' @keywords internal
.signflip_spa_pvalues <- function(nu, S_obs) {
   # Numerically stable elementary functions.
   logcosh <- function(x) { ax <- abs(x); ax + log1p(exp(-2 * ax)) - log(2) }
   sech2   <- function(x) { e <- exp(-2 * abs(x)); 4 * e / (1 + e)^2 }

   K  <- function(t) sum(logcosh(t * nu))
   K1 <- function(t) sum(nu * tanh(t * nu))
   K2 <- function(t) sum(nu^2 * sech2(t * nu))

   var_S <- K2(0)
   if (!is.finite(var_S) || var_S <= 0) {
      return(list(p.left = NA_real_, p.right = NA_real_))
   }
   fallback <- function() {
      z <- S_obs / sqrt(var_S)
      list(p.left = stats::pnorm(z), p.right = stats::pnorm(z, lower.tail = FALSE))
   }
   if (abs(S_obs) < 1e-8 * sqrt(var_S)) return(fallback())
   max_nu <- max(abs(nu))
   if (!is.finite(max_nu) || max_nu <= 0) return(fallback())
   tb <- 0.9 * 700 / max_nu
   t_hat <- tryCatch(
      stats::uniroot(function(t) K1(t) - S_obs, lower = -tb, upper = tb, tol = 1e-10)$root,
      error = function(e) NA_real_)
   if (!is.finite(t_hat)) return(fallback())
   r_lr <- sign(t_hat) * sqrt(pmax(2 * (t_hat * S_obs - K(t_hat)), 0))
   v_lr <- t_hat * sqrt(K2(t_hat))
   p.right <- stats::pnorm(r_lr, lower.tail = FALSE) + stats::dnorm(r_lr) * (1 / r_lr - 1 / v_lr)
   p.left  <- stats::pnorm(r_lr, lower.tail = TRUE)  - stats::dnorm(r_lr) * (1 / r_lr - 1 / v_lr)
   list(p.left  = min(max(p.left,  .Machine$double.xmin), 1),
        p.right = min(max(p.right, .Machine$double.xmin), 1))
}


#' Negative-binomial Lugannani-Rice saddlepoint p-values
#'
#' Computes one-sided saddlepoint tail probabilities for
#'   \eqn{S = \sum_i w_i (Y_i - \mu_i)},  with \eqn{Y_i \sim \mathrm{NegBin}(\mu_i, r)}
#' independently (mean \eqn{\mu_i}, dispersion \eqn{r}; \eqn{\mathrm{Var}(Y_i) = \mu_i + \mu_i^2/r}).
#' The NB CGF in centered form is
#' \deqn{K(t) = -r \sum_i \log\!\bigl(1 - \mu_i (e^{w_i t} - 1)/r\bigr) - t \sum_i w_i \mu_i,}
#' and \eqn{K''(0) = \sum_i w_i^2 (\mu_i + \mu_i^2/r)}, matching the per-i variance.
#'
#' @param w  numeric vector of (fixed) weights.
#' @param mu numeric vector of NB means (positive).
#' @param size_ numeric, NB dispersion / size parameter (positive).
#' @param S_obs scalar, the observed centered sum.
#' @return list with \code{p.left} and \code{p.right}.
#' @keywords internal
.nb_spa_pvalues <- function(w, mu, size_, S_obs) {
   K  <- function(t) {
      fi <- 1 - mu * (exp(w * t) - 1) / size_
      if (any(fi <= 0) || any(!is.finite(fi))) return(NA_real_)
      -size_ * sum(log(fi)) - t * sum(w * mu)
   }
   K1 <- function(t) {
      e <- exp(w * t); fi <- 1 - mu * (e - 1) / size_
      if (any(fi <= 0) || any(!is.finite(fi))) return(NA_real_)
      sum(w * mu * e / fi) - sum(w * mu)
   }
   K2 <- function(t) {
      e <- exp(w * t); fi <- 1 - mu * (e - 1) / size_
      if (any(fi <= 0) || any(!is.finite(fi))) return(NA_real_)
      sum(w^2 * mu * e * (fi + mu * e / size_) / fi^2)
   }
   var_S <- K2(0)
   if (!is.finite(var_S) || var_S <= 0) {
      return(list(p.left = NA_real_, p.right = NA_real_))
   }
   fallback <- function() {
      z <- S_obs / sqrt(var_S)
      list(p.left = stats::pnorm(z), p.right = stats::pnorm(z, lower.tail = FALSE))
   }
   if (abs(S_obs) < 1e-8 * sqrt(var_S)) return(fallback())

   # Convergence radius: |w_i t| < log(1 + size/mu_i) for all i.
   bnd <- log1p(size_ / pmax(mu, 1e-300)) / pmax(abs(w), 1e-300)
   bnd <- bnd[is.finite(bnd) & bnd > 0]
   if (!length(bnd)) return(fallback())
   tb <- 0.9 * min(bnd)

   t_hat <- tryCatch(
      stats::uniroot(function(t) {
         v <- K1(t); if (is.na(v)) return(NA_real_); v - S_obs
      }, lower = -tb, upper = tb, tol = 1e-10)$root,
      error = function(e) NA_real_)
   if (!is.finite(t_hat)) return(fallback())

   kt <- K(t_hat); k2t <- K2(t_hat)
   if (!is.finite(kt) || !is.finite(k2t) || k2t <= 0) return(fallback())
   r_lr <- sign(t_hat) * sqrt(pmax(2 * (t_hat * S_obs - kt), 0))
   v_lr <- t_hat * sqrt(k2t)
   p.right <- stats::pnorm(r_lr, lower.tail = FALSE) + stats::dnorm(r_lr) * (1 / r_lr - 1 / v_lr)
   p.left  <- stats::pnorm(r_lr, lower.tail = TRUE)  - stats::dnorm(r_lr) * (1 / r_lr - 1 / v_lr)
   list(p.left  = min(max(p.left,  .Machine$double.xmin), 1),
        p.right = min(max(p.right, .Machine$double.xmin), 1))
}
