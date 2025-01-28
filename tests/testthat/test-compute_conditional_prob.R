library(testthat)

# necessary parameters
p <- 100
M <- 3
K <- 10
gamma <- .8
stay_prob <- .9
aux_info <- help_dgp(K = K, M = M, p = p, gamma = gamma, stay_prob = stay_prob)

# extract information
Q <- aux_info$Q
pEmit <- aux_info$pEmit
pInit <- aux_info$pInit

# test the validity of new function
test_that("compute_conditional_prob_efficient works!", {

  # sample a x
  n <- 500
  X <- matrix(stats::rbinom(n = p * n, size = 1, prob = rep(rep(.5, p), n)),
              nrow = n, byrow = TRUE)

  # use the old function
  system.time(old_result <- apply(X, 1, function(x){
    compute_conditional_prob(x = x,
                             pEmit = pEmit,
                             Q = Q, pInit = pInit)[, 1]
  }))["elapsed"]

  # use the new function
  system.time(new_result <- apply(X, 1, function(x){
    compute_conditional_prob_efficient(x = x,
                                       pEmit = pEmit,
                                       Q = Q, pInit = pInit)[, 1]
  }))["elapsed"]

  # test the equality
  expect_equal(max(abs(old_result - new_result)), 0)

})
