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
    emission_prob <- pEmit[SNP_id, , ]

    # extract the row id for forward and backward probability
    row_id <- sprintf("SNP_%d", SNP_id)

    # compute the numerator of the conditional prob
    numerator <- t(emission_prob) *  (forward_mat[row_id, ] * backward_mat[row_id, ])

    # compute the denominator of the conditional prob
    denominator <- sum(forward_mat[row_id, ] * backward_mat[row_id, ])

    # return the final prob
    return(colSums(numerator) / denominator)
  })

  # return the conditional probability vector
  return(t(conditional_prob))
}

#' Compile fastPhase from R
#'
#' @param fp_path Path to fastPhase software
#' @param X_file Path to inp file
#' @param out_path Path to output file (can be either relative or absolute)
#' @param K Number of latent states in Markov chain
#' @param numit Number of iterations in EM algorithm
#' @param phased Haplotyp (0, 1) or SNP (0,1,2)
#' @param seed Seed for random start in fastPhase
#'
#' @return Output NULL if there is no error from fastPhase
#' @export
fastPhase_new <- function (fp_path, X_file, out_path = NULL, K = 12, numit = 25,
                           phased = FALSE, seed = 1)
{

  # regulate/check the inputs
  K = as.integer(K)
  numit = as.integer(numit)
  seed = as.integer(seed)
  stopifnot(is.character(fp_path))
  stopifnot(is.character(X_file))
  stopifnot(is.null(out_path) | is.character(out_path))
  stopifnot(is.integer(K))
  stopifnot(is.integer(numit))
  stopifnot(is.logical(phased))
  stopifnot(is.integer(seed))
  if (!file.exists(fp_path)) {
    message(paste("SNPknock could find the fastPHASE executable: '",
                  fp_path, "' does not exist.\nIf you have not downloaded it yet, you can obtain fastPHASE from: http://scheet.org/software.html",
                  sep = ""))
    return(NULL)
  }
  if (is.null(out_path)) {
    out_path = tempfile(pattern = "file", tmpdir = tempdir(),
                        fileext = "")
  }

  # compose the inputs for fastPhase
  out_path_dirname = tools::file_path_as_absolute(dirname(out_path))
  out_path_basename = basename(out_path)
  out_path_abs = paste(out_path_dirname, sprintf("%s/", out_path_basename),
                       sep = "/")
  command = fp_path
  command = paste(command, " -Pp -T1 -K", K, sep = "")
  command = paste(command, " -M1 -g -H-4 -C", numit, sep = "")
  if (phased) {
    command = paste(command, " -B", sep = "")
  }
  command = paste(command, " -S", seed, sep = "")
  command = paste(command, " -o'", out_path_abs, "' ", X_file,
                  sep = "")
  cat(command)

  # run fastPhase
  tryCatch(system(command), error = function(e) 1)

  return(NULL)
}

#' Generate emission, transition and initial probabilites
#'
#' @param K Number of latent states
#' @param M Support size of haplotype
#' @param p Number of markers
#' @param gamma Emission probability for X being 0
#' @param stay_prob Stay probability for the Markov jump process
#' @param initial_prob Initial probability vector
#'
#' @return List including final pEmit, Q and pInit
#' @export
help_dgp <- function(K, M, p,
                     gamma,
                     stay_prob,
                     initial_prob = rep(1 / K, K)){

  # initial state distribution
  pInit <- stats::setNames(initial_prob, sprintf("state_%d", 1:K))

  # define transition matrix
  transition_mat <- matrix(0, nrow = K, ncol = K)
  diag(transition_mat) <- stay_prob
  transition_mat[cbind(1:K, (1:K) %% K + 1)] <- 1 - stay_prob

  # create emission distribution
  emission_mat <- matrix(1 / M, nrow = M, ncol = K)
  emission_mat[, 2 : K] <- c(gamma, rep((1 - gamma) / (M - 1), M - 1))

  # Create an array for storing transition and emission distributions
  Q <- array(0, dim = c(p - 1, K, K),
             dimnames = list(
               transition = 1:(p-1),
               cur_state = sprintf("current_%d", 1:K),
               tran_state = sprintf("transit_%d", 1:K)
             ))
  pEmit <- array(0, dim = c(p, M, K),
                 dimnames = list(
                   SNP = 1:p,
                   obs_state = sprintf("obs_%d", 1:M),
                   latent_state = sprintf("latent_%d", 1:K)
                 ))

  # construct array for emission and transition probabilities
  for (SNP in 1:p) {
    pEmit[SNP,,] <- emission_mat
    if(SNP != p){
      Q[SNP,,] <- transition_mat
    }
  }

  # return the outputs
  return(list(pEmit = pEmit, pInit = pInit, Q = Q))
}

#' rename the output from fastPhase
#'
#' @param Q Transition probability array
#' @param pEmit Emission probability array
#' @param pInit Initial probability vector
#'
#' @return A list including the renamed inputs
#' @export
name_output <- function(Q, pEmit, pInit){

  # name transition probability array
  dimnames(Q) <- list(
    transition = 1:(dim(Q)[1]),
    cur_state = sprintf("current_%d", 1:(dim(Q)[2])),
    tran_state = sprintf("transit_%d", 1:(dim(Q)[3]))
  )

  # name emission probability array
  dimnames(pEmit) <- list(
    SNP = 1:(dim(pEmit)[1]),
    obs_state = sprintf("obs_%d", 1:(dim(pEmit)[2])),
    latent_state = sprintf("latent_%d", 1:(dim(pEmit)[3]))
  )

  # return the final results
  return(list(
    Q = Q,
    pEmit = pEmit,
    pInit = stats::setNames(pInit, sprintf("state_%d", 1:length(pInit)))
  ))
}


compute_conditional_prob_efficient <- function(x, pInit, pEmit, Q){

  # compute the total number of x
  num_snp <- length(x)

  # compute A and B matrix
  forward_mat <- compute_forward_prob(x = x, pInit = pInit, pEmit = pEmit, Q = Q)
  backward_mat <- compute_backward_prob(x = x, pInit = pInit, pEmit = pEmit, Q = Q)

  # transform forward_mat and backward_mat to arrays
  forward_array <- aperm(array(forward_mat,
                               c(nrow(forward_mat),
                                 ncol(forward_mat), dim(pEmit)[2])), c(1, 3, 2))
  backward_array <- aperm(array(backward_mat,
                                c(nrow(backward_mat),
                                  ncol(backward_mat), dim(pEmit)[2])), c(1, 3, 2))

  # compute the final probability array
  numerator_array <- forward_array * backward_array * pEmit
  denominator_mat <- forward_mat * backward_mat
  conditional_prob <- apply(numerator_array, c(1, 2), sum) / apply(denominator_mat, 1, sum)

  # return the conditional probability vector
  return(conditional_prob)
}
