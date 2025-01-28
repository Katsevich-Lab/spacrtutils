compute_all_means <- function(fitted_model, x, conditional_prob, support_x){

  # extract the dimnesion of the porblem
  p <- length(x)
  M <- length(support_x)

  # compute the all the conditioinal expectations Y|X_{-j}
  # obtain the fitted model on (x_1,..,x_{j-1}, \{0,\ldots,M\},x_{j+1},...,x_p)
  x_impute <- matrix(rep(x, M), nrow = p, ncol = M)
  diag(x_impute) <- 1

}
