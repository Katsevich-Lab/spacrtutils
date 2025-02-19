# set.seed(237)
#
# n <- 50; p <- 2; normalize <- FALSE; return_cdf <- FALSE
#
# data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#              Y = rpois(n = n, lambda = 1),
#              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#
# X_on_Z_fam <- "binomial"
# Y_on_Z_fam <- "poisson"
#
# (results_old <- spaCRT_old(data, X_on_Z_fam, Y_on_Z_fam, normalize))
# (results_new <- spaCRT(data, X_on_Z_fam, Y_on_Z_fam, normalize))
# (results_gcm <- GCM(data, X_on_Z_fam, Y_on_Z_fam))
