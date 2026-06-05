library(spacrtutils)

# Reference null DGP. X _||_ Y | Z under H0: beta = 0 in the GLM
#   g(E[Y | X, Z]) = beta * X + Z gamma.
make_data <- function(n = 800, kind = c("pois", "bern", "nb", "gauss"),
                      seed = 42, beta = 0) {
   kind <- match.arg(kind)
   set.seed(seed)
   Z <- matrix(rnorm(2 * n), n, 2)
   X <- rnorm(n)
   lp <- 0.3 * Z[, 1] + 0.2 * Z[, 2] + beta * X
   if (kind == "pois") {
      Y <- rpois(n, exp(0.5 + lp))
   } else if (kind == "bern") {
      Y <- rbinom(n, 1, plogis(-0.5 + lp))
   } else if (kind == "nb") {
      Y <- MASS::rnegbin(n, mu = exp(0.5 + lp), theta = 2)
   } else {  # gauss
      Y <- 0.5 + lp + rnorm(n, sd = 0.7)
   }
   list(X = X, Y = Y, Z = Z)
}


# ---- structural tests ------------------------------------------------------
test_that("signflip_score_internal: Poisson returns valid output", {
   d <- make_data(n = 500, kind = "pois", seed = 1)
   r <- signflip_score_internal(d, Y_on_Z_fam = "poisson")
   expect_named(r, c("test_stat", "p.left", "p.right", "p.both",
                     "spa.success", "NB.disp.param"))
   expect_true(is.finite(r$test_stat))
   expect_true(r$spa.success)
   expect_gte(r$p.left, 0);  expect_lte(r$p.left, 1)
   expect_gte(r$p.right, 0); expect_lte(r$p.right, 1)
   expect_lte(abs(r$p.left + r$p.right - 1), 1e-8)        # LR identity
   expect_equal(r$p.both, 2 * min(r$p.left, r$p.right))
   expect_true(is.na(r$NB.disp.param))
})

test_that("signflip_score_internal works on Binomial / Gaussian / NB", {
   for (kind in c("bern", "gauss", "nb")) {
      fam <- switch(kind, bern = "binomial", gauss = "gaussian",
                          nb = "negative.binomial")
      d <- make_data(n = 600, kind = kind, seed = 7)
      r <- signflip_score_internal(d, Y_on_Z_fam = fam)
      expect_true(r$spa.success, info = paste("family:", fam))
      expect_true(is.finite(r$test_stat))
      expect_gte(r$p.both, 0); expect_lte(r$p.both, 1)
      expect_lte(abs(r$p.left + r$p.right - 1), 1e-8)
      if (kind == "nb") expect_true(r$NB.disp.param > 0)
      else              expect_true(is.na(r$NB.disp.param))
   }
})


# ---- user-supplied fitted values reproduce auto-fit ------------------------
test_that("Poisson user-supplied mu_y_hat equals auto-fit", {
   d <- make_data(n = 800, kind = "pois", seed = 5)
   auto <- signflip_score_internal(d, Y_on_Z_fam = "poisson")
   mu_y <- as.numeric(stats::fitted(
      stats::glm(d$Y ~ d$Z, family = stats::poisson())))
   man  <- signflip_score_internal(d, Y_on_Z_fam = "poisson",
                                   fit_vals_Y_on_Z_own = mu_y)
   expect_equal(man$test_stat, auto$test_stat, tolerance = 1e-10)
   expect_equal(man$p.left,    auto$p.left,    tolerance = 1e-10)
   expect_equal(man$p.right,   auto$p.right,   tolerance = 1e-10)
})

test_that("NB user-supplied (mu_y_hat, size_hat) equals auto-fit", {
   d <- make_data(n = 800, kind = "nb", seed = 5)
   auto <- signflip_score_internal(d, Y_on_Z_fam = "negative.binomial")
   fit  <- suppressWarnings(MASS::glm.nb(d$Y ~ d$Z))
   man  <- signflip_score_internal(d, Y_on_Z_fam = "negative.binomial",
                                   fit_vals_Y_on_Z_own = as.numeric(stats::fitted(fit)),
                                   size_hat_own = fit$theta)
   expect_equal(man$test_stat,    auto$test_stat,    tolerance = 1e-10)
   expect_equal(man$p.left,       auto$p.left,       tolerance = 1e-10)
   expect_equal(man$p.right,      auto$p.right,      tolerance = 1e-10)
   expect_equal(man$NB.disp.param, auto$NB.disp.param, tolerance = 1e-10)
})

test_that("NB user-supplied mu without size_hat errors", {
   d <- make_data(n = 200, kind = "nb", seed = 1)
   expect_error(
      signflip_score_internal(d, Y_on_Z_fam = "negative.binomial",
                              fit_vals_Y_on_Z_own = rep(1, length(d$Y))),
      "size_hat_own"
   )
})


# ---- SPA primitive: mathematical identities --------------------------------
test_that("sign-flip K''(0) equals sum(nu^2)", {
   set.seed(2)
   nu <- rnorm(200)
   sech2 <- function(x) { e <- exp(-2*abs(x)); 4*e/(1+e)^2 }
   K2_0 <- sum(nu^2 * sech2(0))
   expect_equal(K2_0, sum(nu^2), tolerance = 1e-12)
})

test_that("sign-flip CGF is even: right tail at +T == left tail at -T", {
   set.seed(3)
   nu  <- rnorm(300)
   T_  <- 1.5 * sqrt(sum(nu^2))             # ~ 1.5 sd tail
   p_p <- spacrtutils:::.signflip_spa_pvalues(nu,  T_)
   p_n <- spacrtutils:::.signflip_spa_pvalues(nu, -T_)
   expect_equal(p_p$p.right, p_n$p.left,  tolerance = 1e-10)
   expect_equal(p_p$p.left,  p_n$p.right, tolerance = 1e-10)
})

test_that("sign-flip SPA == Bernoulli SPA at mu = 1/2 with w = 2 nu", {
   set.seed(4)
   nu   <- rnorm(250)
   T_   <- 0.8 * sqrt(sum(nu^2))
   p_sf <- spacrtutils:::.signflip_spa_pvalues(nu, T_)
   p_bn <- spacrtutils:::.bin_spa_pvalues(w     = 2 * nu,
                                          mu    = rep(0.5, length(nu)),
                                          S_obs = T_)
   expect_equal(p_sf$p.left,  p_bn$p.left,  tolerance = 1e-9)
   expect_equal(p_sf$p.right, p_bn$p.right, tolerance = 1e-9)
})


# ---- SPA primitive: Monte-Carlo cross-check (smoke) ------------------------
test_that("sign-flip SPA agrees with Monte-Carlo sign-flip reference", {
   skip_on_cran()
   set.seed(11)
   n  <- 300
   nu <- rnorm(n, sd = 1.2) * runif(n, 0.5, 1.5)
   T_obs <- sum(nu) + 0.7 * sqrt(sum(nu^2))     # ~ 0.7 sd to the right of mean

   p_spa <- spacrtutils:::.signflip_spa_pvalues(nu, T_obs)$p.right

   B <- 5e4
   G <- matrix(sample(c(-1, 1), n * B, replace = TRUE), nrow = n)
   T_b <- as.numeric(crossprod(nu, G))
   p_mc <- (sum(T_b >= T_obs) + 1) / (B + 1)

   # MC SE ~ sqrt(p(1-p)/B) ~ 2e-3 at p~0.25, SPA bias ~ O(1/n) = 3e-3.
   expect_lt(abs(p_spa - p_mc), 1e-2)
})


# ---- canonical-link identity: H == 1 ---------------------------------------
test_that("canonical Poisson and Binomial give H = 1 (score = X_tilde * (Y - mu))", {
   d  <- make_data(n = 400, kind = "pois", seed = 21)
   fam <- stats::poisson()
   fit <- stats::glm(d$Y ~ d$Z, family = fam)
   mu  <- as.numeric(stats::fitted(fit))
   eta <- as.numeric(stats::predict(fit, type = "link"))
   H   <- fam$mu.eta(eta) / fam$variance(mu)
   expect_equal(max(abs(H - 1)), 0, tolerance = 1e-10)

   fam <- stats::binomial()
   db  <- make_data(n = 400, kind = "bern", seed = 21)
   fit <- stats::glm(db$Y ~ db$Z, family = fam)
   mu  <- as.numeric(stats::fitted(fit))
   eta <- as.numeric(stats::predict(fit, type = "link"))
   H   <- fam$mu.eta(eta) / fam$variance(mu)
   expect_equal(max(abs(H - 1)), 0, tolerance = 1e-10)
})

test_that("NB log-link has H = 1/(1 + mu/r)", {
   d  <- make_data(n = 400, kind = "nb", seed = 31)
   fit <- suppressWarnings(MASS::glm.nb(d$Y ~ d$Z))
   r   <- fit$theta
   mu  <- as.numeric(stats::fitted(fit))
   eta <- as.numeric(stats::predict(fit, type = "link"))
   fam <- MASS::negative.binomial(r)
   H   <- fam$mu.eta(eta) / fam$variance(mu)
   expect_equal(H, 1 / (1 + mu / r), tolerance = 1e-10)
})


###############################################################################
# NB GLM score test: empirical correctness (calibration / MC / power)
###############################################################################

test_that("NB null calibration: Type-I rate near nominal at n=500, R=300", {
   skip_on_cran()
   set.seed(101)
   R <- 300L
   pvals <- vapply(seq_len(R), function(i) {
      d <- make_data(n = 500, kind = "nb", seed = 100 + i)
      r <- signflip_score_internal(d, Y_on_Z_fam = "negative.binomial")
      r$p.both
   }, numeric(1))
   pvals <- pvals[!is.na(pvals)]
   expect_gt(length(pvals), 0.95 * R)        # >=95% of reps succeed
   N <- length(pvals)
   # MC SE at alpha=0.10 is sqrt(0.09/N). Allow 4 SE for two-sided slack.
   t1_10 <- mean(pvals < 0.10)
   t1_05 <- mean(pvals < 0.05)
   expect_lt(abs(t1_10 - 0.10), 4 * sqrt(0.10 * 0.90 / N))
   expect_lt(abs(t1_05 - 0.05), 4 * sqrt(0.05 * 0.95 / N))
   # KS test on full p-value distribution.
   ks <- suppressWarnings(ks.test(pvals, "punif"))
   expect_gt(ks$p.value, 0.005)              # only reject at very small level
})

test_that("NB SPA matches Monte-Carlo sign-flip (n=300, B=2e4)", {
   skip_on_cran()
   set.seed(202)
   R <- 25L
   diffs <- numeric(R)
   for (i in seq_len(R)) {
      d <- make_data(n = 300, kind = "nb", seed = 200 + i)
      # Package SPA p-value
      r <- signflip_score_internal(d, Y_on_Z_fam = "negative.binomial")
      if (is.na(r$p.both)) { diffs[i] <- NA_real_; next }
      # Reconstruct effective score (identical code path to the package).
      fit <- tryCatch(suppressWarnings(MASS::glm.nb(d$Y ~ d$Z)),
                      error = function(e) NULL)
      if (is.null(fit)) { diffs[i] <- NA_real_; next }
      mu <- as.numeric(stats::fitted(fit))
      eta <- as.numeric(stats::predict(fit, type = "link"))
      fam <- MASS::negative.binomial(fit$theta)
      mu_eta <- fam$mu.eta(eta); V_mu <- fam$variance(mu)
      H <- mu_eta / V_mu; W <- mu_eta * H
      Zd <- cbind(1, d$Z)
      delta <- solve(crossprod(Zd, Zd * W), crossprod(Zd * W, d$X))
      X_tilde <- as.numeric(d$X - Zd %*% delta)
      nu <- X_tilde * H * (d$Y - mu)
      t0 <- sum(nu); n <- length(nu); B <- 2e4L
      G <- matrix(sample(c(-1, 1), n * B, replace = TRUE), nrow = n)
      T_b <- as.numeric(crossprod(nu, G))
      p_mc_R <- (sum(T_b >= abs(t0)) + 1) / (B + 1)
      p_mc_L <- (sum(T_b <= -abs(t0)) + 1) / (B + 1)
      p_mc   <- min(1, p_mc_R + p_mc_L)
      diffs[i] <- abs(r$p.both - p_mc)
   }
   diffs <- diffs[!is.na(diffs)]
   expect_gt(length(diffs), 0.9 * R)
   # Median diff should be small (SPA O(1/n)~3e-3 + MC SE ~7e-3 at p~0.5).
   expect_lt(median(diffs), 0.02)
   expect_lt(quantile(diffs, 0.9), 0.05)
})

test_that("NB sign-flip has power under alternative (beta=0.5, n=500)", {
   skip_on_cran()
   set.seed(303)
   R <- 200L
   pvals <- vapply(seq_len(R), function(i) {
      d <- make_data(n = 500, kind = "nb", seed = 300 + i, beta = 0.5)
      r <- signflip_score_internal(d, Y_on_Z_fam = "negative.binomial")
      r$p.both
   }, numeric(1))
   pvals <- pvals[!is.na(pvals)]
   power <- mean(pvals < 0.05)
   # Under beta=0.5 with theta=2, n=500, power at alpha=0.05 should be large.
   expect_gt(power, 0.5)
})

test_that("NB sign-flip is robust to mild Poisson misspecification", {
   skip_on_cran()
   # True DGP is Poisson; we fit NB. The NB MLE should converge to theta
   # large (-> Poisson), and the test should still calibrate.
   set.seed(404)
   R <- 200L
   pvals <- vapply(seq_len(R), function(i) {
      d <- make_data(n = 500, kind = "pois", seed = 400 + i)
      r <- signflip_score_internal(d, Y_on_Z_fam = "negative.binomial")
      r$p.both
   }, numeric(1))
   pvals <- pvals[!is.na(pvals)]
   N <- length(pvals)
   expect_gt(N, 0.95 * R)
   expect_lt(abs(mean(pvals < 0.05) - 0.05),
             4 * sqrt(0.05 * 0.95 / N))
})


###############################################################################
# Multinomial regression via one-vs-rest binary reduction
###############################################################################

# Generate a K-category multinomial Y with softmax(eta + beta_k X + Z gamma_k).
# Under H0: beta_k = 0 for all k=1..K-1, X is independent of Y given Z, so each
# binary indicator I(Y=k) | Z is independent of X. The signflip score test on
# any such indicator should calibrate.
make_multinom_data <- function(n, K = 3L, seed = 1L, beta_k = NULL) {
   if (is.null(beta_k)) beta_k <- rep(0, K - 1L)
   stopifnot(length(beta_k) == K - 1L)
   set.seed(seed)
   Z <- matrix(rnorm(2 * n), n, 2)
   X <- rnorm(n)
   eta_nonref <- matrix(0, n, K - 1L)
   for (k in seq_len(K - 1L)) {
      eta_nonref[, k] <- 0.2 * (k - 1) + 0.4 * Z[, 1] - 0.3 * Z[, 2] +
                         beta_k[k] * X
   }
   eta_full <- cbind(eta_nonref, 0)                          # reference = K
   probs    <- exp(eta_full); probs <- probs / rowSums(probs)
   Y <- vapply(seq_len(n), function(i)
      sample.int(K, 1L, prob = probs[i, ]), integer(1))
   list(X = X, Y = Y, Z = Z)
}

test_that("multinomial DGP yields balanced categories under H0", {
   set.seed(1)
   dm <- make_multinom_data(n = 5000, K = 3L)
   tb <- as.numeric(table(dm$Y) / 5000)
   # Marginal proportions depend on Z marginal; check none are degenerate.
   expect_true(all(tb > 0.10))
   expect_true(all(tb < 0.80))
})

test_that("multinomial one-vs-rest: each category's binary null calibrates", {
   skip_on_cran()
   set.seed(505)
   K <- 3L; R <- 250L
   pvals_by_k <- vector("list", K)
   for (i in seq_len(R)) {
      dm <- make_multinom_data(n = 500, K = K, seed = 500 + i)
      for (k in seq_len(K)) {
         dk <- list(X = dm$X, Y = as.integer(dm$Y == k), Z = dm$Z)
         r  <- signflip_score_internal(dk, Y_on_Z_fam = "binomial")
         pvals_by_k[[k]] <- c(pvals_by_k[[k]], r$p.both)
      }
   }
   for (k in seq_len(K)) {
      pv <- pvals_by_k[[k]]; pv <- pv[!is.na(pv)]
      N <- length(pv)
      expect_gt(N, 0.95 * R)
      t1 <- mean(pv < 0.05)
      expect_lt(abs(t1 - 0.05), 4 * sqrt(0.05 * 0.95 / N),
                label = sprintf("|T1(k=%d) - 0.05| (T1=%.4f, N=%d)",
                                k, t1, N))
   }
})

test_that("multinomial one-vs-rest: power detects per-category X effect", {
   skip_on_cran()
   set.seed(606)
   K <- 3L; R <- 150L
   beta_k <- c(0.6, 0)            # X drives category 1 only (vs ref=3)
   power_by_k <- numeric(K)
   for (k in seq_len(K)) {
      pv <- vapply(seq_len(R), function(i) {
         dm <- make_multinom_data(n = 600, K = K, seed = 600 + i,
                                  beta_k = beta_k)
         dk <- list(X = dm$X, Y = as.integer(dm$Y == k), Z = dm$Z)
         r  <- signflip_score_internal(dk, Y_on_Z_fam = "binomial")
         r$p.both
      }, numeric(1))
      pv <- pv[!is.na(pv)]
      power_by_k[k] <- mean(pv < 0.05)
   }
   # In the multinomial-via-binary reduction, ALL categories' indicators
   # I(Y=k) acquire X-dependence when any one beta_k != 0, because
   #   d/dX log P(Y=j) = beta_j - sum_k pi_k beta_k
   # is generically nonzero for all j (with magnitude proportional to
   # beta_1 * (1 - pi_1) for j=1 and -beta_1 * pi_j for j != 1). All three
   # binary tests should therefore have non-trivial power; we just require
   # the direct-effect category to clear a reasonable bar.
   expect_gt(power_by_k[1], 0.30,
             label = sprintf("power_by_k = (%.3f, %.3f, %.3f)",
                             power_by_k[1], power_by_k[2], power_by_k[3]))
})

test_that("multinomial one-vs-rest matches MC sign-flip on each category", {
   skip_on_cran()
   set.seed(707)
   K <- 3L
   diffs <- c()
   for (i in seq_len(15L)) {
      dm <- make_multinom_data(n = 300, K = K, seed = 700 + i)
      for (k in seq_len(K)) {
         dk <- list(X = dm$X, Y = as.integer(dm$Y == k), Z = dm$Z)
         r  <- signflip_score_internal(dk, Y_on_Z_fam = "binomial")
         if (is.na(r$p.both)) next
         # Manual MC effective-score reconstruction (Bernoulli canonical: H=1).
         fit <- suppressWarnings(stats::glm(dk$Y ~ dk$Z,
                                            family = stats::binomial()))
         mu <- as.numeric(stats::fitted(fit))
         W  <- mu * (1 - mu)
         Zd <- cbind(1, dk$Z)
         delta <- solve(crossprod(Zd, Zd * W), crossprod(Zd * W, dk$X))
         X_tilde <- as.numeric(dk$X - Zd %*% delta)
         nu <- X_tilde * (dk$Y - mu)
         t0 <- sum(nu); B <- 2e4L; n <- length(nu)
         G <- matrix(sample(c(-1, 1), n * B, replace = TRUE), nrow = n)
         T_b <- as.numeric(crossprod(nu, G))
         p_mc_R <- (sum(T_b >= abs(t0)) + 1) / (B + 1)
         p_mc_L <- (sum(T_b <= -abs(t0)) + 1) / (B + 1)
         p_mc   <- min(1, p_mc_R + p_mc_L)
         diffs <- c(diffs, abs(r$p.both - p_mc))
      }
   }
   expect_lt(median(diffs), 0.02)
   expect_lt(quantile(diffs, 0.9), 0.05)
})
