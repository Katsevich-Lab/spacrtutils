library(testthat)
library(glmnet)

# generate pEmit, Q and pInit using spacrt package
p <- 500
K <- 5
M <- 2
n <- 1000
aux_info <- withr::with_seed(1, spacrt::help_dgp(p = p, K = K, M = M,
                                                 stay_prob = .9, beta_prior = TRUE,
                                                 alpha = 0.5, beta = 1.5))

# sample from HMM
X <- SNPknock::sampleHMM(aux_info$pInit, aux_info$Q, aux_info$pEmit, n = n)

# create effect size
signal_prop <- 0.1
eta <- 1
gamma_0 <- -2
num_signal <- round(signal_prop * p)
positive_signal <- withr::with_seed(1, runif(n = num_signal / 2, min = eta, max = eta))
negative_signal <- withr::with_seed(1, runif(n = num_signal / 2, min = -eta, max = -eta))
effect_size <- c(positive_signal, negative_signal)

# sample high dimensional logistic model for Y
Y <- rbinom(n = n, size = 1, prob = 1 / (1 + exp(-gamma_0 - colSums(t(X[, 1:num_signal]) * effect_size))))


# compute the conditional probability with oracle Q, pEmit and pInit
conditional_mean <- t(
  apply(X, 1, function(x){
    spacrt::compute_conditional_prob_Rcpp(x = x, pInit = aux_info$pInit,
                                          pEmit = aux_info$pEmit, Q = aux_info$Q)[, 2]

  })
)

# model Y on Z via lasso
fitted_lasso <- glmnet::cv.glmnet(x = X, y = Y, family = "binomial")

# test new implementation
test_that("compute_all_means_efficient works!", {

  # use lambda.min
  system.time({
    leave_one_fit <- t(
      sapply(1:n, function(i){
        spacrt::compute_all_means(fitted_model = fitted_lasso, lambda = "min", x = X[i, ],
                                  support_x = c(0, 1),
                                  conditional_prob = matrix(c(1 - conditional_mean[i, ],
                                                              conditional_mean[i, ]), ncol = 2))
      })
    )
  })["elapsed"]

  # more efficient implementation
  system.time({
    leave_one_fit_efficient <- spacrt::compute_all_means_efficient(fitted_model = fitted_lasso,
                                                                   X = X, conditional_prob_mat = cbind(1 - conditional_mean, conditional_mean),
                                                                   support_x = c(0, 1), lambda = "min")
  })["elapsed"]

  # test if the computation difference is only within the machine error
  expect_lt(max(abs(leave_one_fit - leave_one_fit_efficient)), 1e-10)
})


