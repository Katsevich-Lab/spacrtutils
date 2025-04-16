library(testthat)
library(spacrtutils)
library(SNPknock)

# necessary parameters
p <- 500
M <- 2
K <- 5
stay_prob <- .9
aux_info <- help_dgp(K = K, M = M, p = p,
                     stay_prob = stay_prob, beta_prior = TRUE,
                     alpha = .5, beta = 1.5)

# extract information
Q <- aux_info$Q
pEmit <- aux_info$pEmit
pInit <- aux_info$pInit

# test the validity of new function
test_that("compute_conditional_prob_efficient works!", {

  # sample a x
  n <- 1000
  X <- SNPknock::sampleHMM(pInit = pInit, Q = Q, pEmit = pEmit, n = n)

  # use the old function
  system.time(old_result <- apply(X, 1, function(x){
    compute_conditional_prob(x = x, pEmit = pEmit, Q = Q, pInit = pInit)[, 1]
  }))["elapsed"]

  # use the new function
  system.time(new_result <- apply(X, 1, function(x){
    compute_conditional_prob_Rcpp(x = x, pEmit = pEmit, Q = Q, pInit = pInit)[, 1]
  }))["elapsed"]

  # test if the computation difference is only within the machine error
  expect_lt(max(abs(old_result - new_result)), 1e-10)

})
