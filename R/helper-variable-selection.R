#' Compute all conditional means for Y|X_{-j}
#'
#' @param fitted_model Fitted regression model
#' @param x Query point
#' @param conditional_prob Conditional probability matrix (p-by-M)
#' @param support_x Support set of X
#' @param lambda Lambda regularization choice in cv.glmnet
#'
#' @return Conditional probability vector E(Y|X_{-j} = x_{-j})
#' @export
compute_all_means <- function(fitted_model, x, conditional_prob, support_x, lambda = "min"){

  # extract the dimension of the problem
  p <- length(x)

  # obtain a square matrix
  x_impute <- matrix(rep(x, p), nrow = p, ncol = p, byrow = TRUE)

  # compute the all the conditional expectations Y|X_{-j}
  fix_one_fitted <- sapply(sort(support_x), function(s){

    # obtain the fitted model on (x_1,..,x_{j-1}, \{0,\ldots,M\},x_{j+1},...,x_p)
    diag(x_impute) <- s

    # use predict function to obtain the leave-one-out fitted values
    stats::predict(fitted_model, newx = x_impute, s = sprintf("lambda.%s", lambda),
                   type = "response")

  })

  # compute the
  integrate_one_fitted <- rowSums(fix_one_fitted * conditional_prob)

  # return the output
  return(integrate_one_fitted)
}
