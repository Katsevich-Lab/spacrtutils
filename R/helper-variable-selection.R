#' Compute all conditional means for Y|X_{-j}
#'
#' @param fitted_model Fitted regression model
#' @param x Query point
#' @param conditional_prob Conditional probability matrix (p-by-M)
#' @param support_x Support set of X
#' @param lambda Lambda regularization choice in cv.glmnet
#' @param post_lasso A logical value. If TRUE, post_lasso should be used and the
#' `fitted_model` should be a low-dimensional GLM
#'
#' @return Conditional probability vector E(Y|X_{-j} = x_{-j})
#' @export
compute_all_means <- function(fitted_model, x, conditional_prob, support_x, lambda = "min",
                              post_lasso = FALSE){

  # extract the dimension of the problem
  p <- length(x)

  # obtain a square matrix
  x_impute <- matrix(rep(x, p), nrow = p, ncol = p, byrow = TRUE,
                     dimnames = list(
                       sample = 1:p,
                       predictor = 1:p
                     ))

  # separate the analysis for different post_lasso parameters
  if(post_lasso){

    # compute the all the conditional expectations Y|X_{-j}
    fix_one_fitted <- sapply(sort(support_x), function(s){

      # obtain the fitted model on (x_1,..,x_{j-1}, \{0,\ldots,M\},x_{j+1},...,x_p)
      diag(x_impute) <- s

      # use predict function to obtain the leave-one-out fitted values
      stats::predict(fitted_model$model,
                     newx = x_impute[, sort(fitted_model$act_coef)],
                     type = "response")

    })

  }else{

    # compute the all the conditional expectations Y|X_{-j}
    fix_one_fitted <- sapply(sort(support_x), function(s){

      # obtain the fitted model on (x_1,..,x_{j-1}, \{0,\ldots,M\},x_{j+1},...,x_p)
      diag(x_impute) <- s

      # use predict function to obtain the leave-one-out fitted values
      stats::predict(fitted_model, newx = x_impute, s = sprintf("lambda.%s", lambda),
                     type = "response")

    })
  }

  # compute the
  integrate_one_fitted <- rowSums(fix_one_fitted * conditional_prob)

  # return the output
  return(integrate_one_fitted)
}


#' This is a post-lasso fitting function
#'
#' @param X A matrix including n rows and p columns
#' @param Y A vector of length n
#' @param family A GLM family
#' @param lambda_param_vec A vector including either "lambda.min" or "lambda.1se"
#'
#' @return A list of fitted models and active set selected from lasso algorithm
#' @export
post_lasso <- function(X, Y, family = "binomial",
                       lambda_param_vec = c("lambda.min", "lambda.1se")){

  # fit Y on X using lasso
  lasso_model <- glmnet::glmnet(x = X, y = Y, family = family)

  # transform the data to data frame
  dimnames(X) <- list(
    sample = 1:nrow(X),
    predictor = 1:ncol(X)
  )

  # extract the lambda.min and lambda.1se parameters
  fitted_model <- list(
    lambda.1se = list(NULL),
    lambda.min = list(NULL)
  )
  for (lambda_param in lambda_param_vec) {

    # extract the non-zero coefficient
    act_set <- which(as.vector(stats::coef(lasso_model, s = lambda_param))[-1] != 0)
    fitted_mode[[lambda_param]]$act_set <- act_set

    # perform the glm.fit
    glm_fitted <- stats::glm.fit(x = X[, sort(act_set)], y = Y, family = get(family))

    # extract the fitted model
    fitted_model[[lambda_param]]$model <- glm_fitted
  }

  # return the output
  return(fitted_model)
}
