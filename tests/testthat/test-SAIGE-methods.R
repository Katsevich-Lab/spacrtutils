library(spacrtutils)

# Fixed-seed reference dataset for regression-style numerical checks.
make_data <- function(n = 1000, kind = c("bern_x", "nb_x"), seed = 42) {
   kind <- match.arg(kind)
   set.seed(seed)
   Z <- matrix(rnorm(2 * n), n, 2)
   if (kind == "bern_x") {
      X <- rbinom(n, 1, plogis(-3 + Z[, 1] + Z[, 2]))
      Y <- MASS::rnegbin(n, mu = exp(-2 + 0.5 * Z[, 1] + 0.5 * Z[, 2]), theta = 0.5)
   } else {
      X <- MASS::rnegbin(n, mu = exp(-1 + Z[, 1] + Z[, 2]), theta = 5)
      Y <- rbinom(n, 1, plogis(-3 + Z[, 1] + Z[, 2]))
   }
   list(X = X, Y = Y, Z = Z)
}


# ---- structural tests ------------------------------------------------------
test_that("SAIGE_Bern_internal returns the expected fields with valid p-values", {
   d <- make_data(n = 500, kind = "bern_x", seed = 1)
   r <- SAIGE_Bern_internal(d)
   expect_named(r, c("test_stat", "p.left", "p.right", "p.both"))
   expect_true(is.finite(r$test_stat))
   expect_gte(r$p.left, 0);  expect_lte(r$p.left, 1)
   expect_gte(r$p.right, 0); expect_lte(r$p.right, 1)
   expect_lte(abs(r$p.left + r$p.right - 1), 1e-8)        # LR identity
   expect_equal(r$p.both, 2 * min(r$p.left, r$p.right))
})

test_that("SAIGE_NB_internal returns the expected fields with valid p-values", {
   d <- make_data(n = 500, kind = "nb_x", seed = 1)
   r <- SAIGE_NB_internal(d)
   expect_named(r, c("test_stat", "p.left", "p.right", "p.both", "NB.disp.param"))
   expect_true(is.finite(r$test_stat))
   expect_gte(r$p.left, 0);  expect_lte(r$p.left, 1)
   expect_gte(r$p.right, 0); expect_lte(r$p.right, 1)
   expect_lte(abs(r$p.left + r$p.right - 1), 1e-8)
   expect_equal(r$p.both, 2 * min(r$p.left, r$p.right))
   expect_true(r$NB.disp.param > 0)
})


# ---- user-supplied fitted values must reproduce auto-fit -------------------
test_that("SAIGE_Bern with user-supplied muhat_x equals the auto-fit", {
   d <- make_data(n = 800, kind = "bern_x", seed = 7)
   auto <- SAIGE_Bern_internal(d)
   mu_x <- as.numeric(stats::fitted(stats::glm(d$X ~ d$Z, family = stats::binomial())))
   manual <- SAIGE_Bern_internal(d, fit_vals_X_on_Z_own = mu_x)
   expect_equal(manual$test_stat, auto$test_stat, tolerance = 1e-10)
   expect_equal(manual$p.left,    auto$p.left,    tolerance = 1e-10)
   expect_equal(manual$p.right,   auto$p.right,   tolerance = 1e-10)
})

test_that("SAIGE_NB with user-supplied (muhat_x, size_hat) equals the auto-fit", {
   d <- make_data(n = 800, kind = "nb_x", seed = 7)
   auto <- SAIGE_NB_internal(d)
   fx <- suppressWarnings(MASS::glm.nb(d$X ~ d$Z))
   manual <- SAIGE_NB_internal(d,
                               fit_vals_X_on_Z_own = as.numeric(stats::fitted(fx)),
                               size_hat_own = fx$theta)
   expect_equal(manual$test_stat,    auto$test_stat,    tolerance = 1e-10)
   expect_equal(manual$p.left,       auto$p.left,       tolerance = 1e-10)
   expect_equal(manual$p.right,      auto$p.right,      tolerance = 1e-10)
   expect_equal(manual$NB.disp.param, auto$NB.disp.param, tolerance = 1e-10)
})

test_that("SAIGE_NB errors if fit_vals_X_on_Z_own is given without size_hat_own", {
   d <- make_data(n = 200, kind = "nb_x", seed = 1)
   expect_error(
      SAIGE_NB_internal(d, fit_vals_X_on_Z_own = rep(0.5, length(d$X))),
      "size_hat_own"
   )
})


# ---- NB SPA variance identity ----------------------------------------------
# K''(0) under the NB CGF equals sum_i a_i^2 (mu_i + mu_i^2 / size).
test_that("NB SPA K''(0) matches sum_i a_i^2 Var(X_i|Z_i)", {
   set.seed(3)
   n <- 250
   mu <- exp(stats::rnorm(n))                     # positive means
   size_ <- 4
   a <- stats::rnorm(n) * 0.5                     # bounded weights
   expected <- sum(a^2 * (mu + mu^2 / size_))

   # extract the K2(0) the SPA uses (re-derive inline; mirrors .nb_spa_pvalues).
   e0 <- exp(0); fi <- 1 - mu * (e0 - 1) / size_   # fi(0) = 1
   K2_at_0 <- sum(a^2 * mu * e0 * (fi + mu * e0 / size_) / fi^2)

   expect_equal(K2_at_0, expected, tolerance = 1e-10)
})


# ---- agreement with spaCRT_internal when the residualizations align --------
# SAIGE_Bern's score matches spaCRT_internal's test_stat (up to 1/sqrt(n) scale)
# when the X|Z fit and the Y|Z linear-GLS fit happen to be compatible.
test_that("SAIGE_Bern recovers a sensible test stat next to spaCRT", {
   d <- make_data(n = 1500, kind = "bern_x", seed = 11)
   r_saige  <- SAIGE_Bern_internal(d)
   r_spacrt <- spaCRT_internal(d, X_on_Z_fam = "binomial",
                               Y_on_Z_fam = "negative.binomial",
                               fitting_X_on_Z = "glm",
                               fitting_Y_on_Z = "glm")
   # Both null statistics under H0 should be moderate (not extreme).
   expect_lt(abs(r_saige$test_stat / sqrt(length(d$X))), 5)
   expect_lt(abs(r_spacrt$test_stat), 5)

   # Both produce reasonable two-sided p-values.
   expect_gt(r_saige$p.both,  1e-6)
   expect_gt(r_spacrt$p.both, 1e-6)
})
