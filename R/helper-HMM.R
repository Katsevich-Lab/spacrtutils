#' Forward algorithm computing matrix A
#'
#' @param x Observation from HMM
#' @param pInit Initial probability vector of length K
#' @param pEmit Emission array, of dimension p-by-M-K. p is number of SNPs, M is support size of observation and K is support size of latent variable from a Markov Chain
#' @param Q Transition array, same dimension with pEmit
#'
#' @return A matrix of dimension p times K
#' @export
compute_forward_prob <- function(x, pInit, pEmit, Q){

  # extract the number of SNPs and number of latent states
  num_snp <- length(x)
  num_state <- length(pInit)

  ###################### compute the forward prob ##############################
  # define empty matrix
  A <- matrix(0, nrow = num_snp, ncol = num_state,
              dimnames = list(
                SNP = sprintf("SNP_%d", 1:num_snp),
                state = sprintf("state_%d", 1:num_state)
              ))

  # impute the first row with initial probability
  A[1, ] <- pInit[colnames(A)]

  # start from 1
  cur_snp <- 1
  while (cur_snp < num_snp) {

    # current emission distribution
    cur_emission <- pEmit[cur_snp, , ]

    # current transition distribution
    cur_transition <- Q[cur_snp, , ]

    # compute the current A value with recursive formula
    A[cur_snp + 1, ] <- sapply(1:num_state, function(state){
      sum(A[cur_snp, ] * cur_emission[sprintf("obs_%d", x[cur_snp] + 1), ] * cur_transition[, sprintf("transit_%d", state)])
    })

    # update the cur_snp
    cur_snp <- cur_snp + 1
  }

  # return matrix of interest A
  return(A)
}

#' Backward algorithm computing matrix B
#'
#' @inheritParams compute_forward_prob
#'
#' @return A matrix of dimension p times K
#' @export
compute_backward_prob <- function(x, pInit, pEmit, Q){

  # extract the number of SNPs and number of latent states
  num_snp <- length(x)
  num_state <- length(pInit)

  ###################### compute the backward prob ##############################
  # define empty matrix
  B <- matrix(0, nrow = num_snp, ncol = num_state,
              dimnames = list(
                SNP = sprintf("SNP_%d", 1:num_snp),
                state = sprintf("state_%d", 1:num_state)
              ))

  # impute the first row with initial probability
  B[num_snp, ] <- 1

  # start from num_snp
  cur_snp <- num_snp
  while (cur_snp > 1) {

    # current emission distribution
    cur_emission <- pEmit[cur_snp, , ]

    # current transition distribution
    cur_transition <- Q[cur_snp - 1, , ]

    # compute the current A value with recursive formula
    B[cur_snp - 1, ] <- sapply(1:num_state, function(state){
      sum(B[cur_snp, ] * cur_emission[sprintf("obs_%d", x[cur_snp] + 1), ] * cur_transition[sprintf("current_%d", state), ])
    })

    # update the cur_snp
    cur_snp <- cur_snp - 1
  }

  # return matrix of interest B
  return(B)
}

#' Compute conditional probability for observed SNPs
#'
#' @inheritParams compute_forward_prob
#'
#' @return A vector of length of total number of SNP
#' @export
compute_conditional_prob <- function(x, pInit, pEmit, Q){

  # compute the total number of x
  num_snp <- length(x)

  # compute A and B matrix
  forward_mat <- compute_forward_prob(x = x, pInit = pInit, pEmit = pEmit, Q = Q)
  backward_mat <- compute_backward_prob(x = x, pInit = pInit, pEmit = pEmit, Q = Q)

  # compute the conditional probability
  conditional_prob <- sapply(1:num_snp, function(SNP_id){

    # compute the emission prob for jth SNP
    emission_prob <- pEmit[SNP_id, sprintf("obs_%d", x[SNP_id] + 1), ]

    # extract the row id for forward and backward probability
    row_id <- sprintf("SNP_%d", SNP_id)

    # compute the numerator of the conditional prob
    numerator <- sum(forward_mat[row_id, ] * emission_prob * backward_mat[row_id, ])

    # compute the denominator of the conditional prob
    denominator <- sum(forward_mat[row_id, ] * backward_mat[row_id, ])

    # return the final prob
    return(numerator / denominator)
  })

  # return the conditional probability vector
  return(conditional_prob)
}
