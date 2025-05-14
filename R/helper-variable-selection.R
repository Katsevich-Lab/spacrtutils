#' Compute all conditional means for Y given X_-j
#'
#' @param fitted_model Fitted regression model
#' @param x Query point
#' @param conditional_prob Conditional probability matrix (p-by-M)
#' @param support_x Support set of X
#' @param lambda Lambda regularization choice in cv.glmnet
#' @param post_lasso A logical value. If TRUE, post_lasso should be used and the
#' `fitted_model` should be a low-dimensional GLM
#'
#' @return Conditional probability vector
#' @export
compute_all_means <- function(fitted_model, x, conditional_prob, support_x, lambda = "min",
                              post_lasso = FALSE){

  # extract the dimension of the problem
  p <- length(x)

  # obtain a square matrix
  x_impute <- matrix(rep(x, p), nrow = p, ncol = p, byrow = TRUE)

  # separate the analysis for different post_lasso parameters
  if(post_lasso){

    # compute the all the conditional expectations Y given X_-j
    fix_one_fitted <- sapply(sort(support_x), function(s){

      # obtain the fitted model on (x_1,..,x_{j-1}, \{0,\ldots,M\},x_{j+1},...,x_p)
      diag(x_impute) <- s

      # separate the case when act_set is NULL
      if(length(fitted_model$act_set) == 0){
        # use just intercept term
        stats::predict(fitted_model$model,
                       newdata = data.frame(rep(1, nrow(x_impute))),
                       type = "response")
      }else{
        # use predict function to obtain the leave-one-out fitted values
        stats::predict(fitted_model$model,
                       newdata = data.frame(x_impute)[, sort(fitted_model$act_set), drop = FALSE],
                       type = "response")
      }

    })

  }else{

    # compute the all the conditional expectations Y given X_-j
    fix_one_fitted <- sapply(sort(support_x), function(s){

      # obtain the fitted model on (x_1,..,x_{j-1}, \{0,\ldots,M\},x_{j+1},...,x_p)
      diag(x_impute) <- s

      # use predict function to obtain the leave-one-out fitted values
      stats::predict(fitted_model, newx = x_impute, s = sprintf("lambda.%s", lambda),
                     type = "response")

    })
  }

  # compute the
  integrate_one_fitted <- rowSums(as.matrix(fix_one_fitted) * conditional_prob)

  # return the output
  return(integrate_one_fitted)
}


#' Compute the conditional mean using power trick (efficient version)
#'
#' @inheritParams compute_all_means
#' @param X data matrix of dim n-by-p
#' @param conditional_prob_mat Conditional probability matrix of dim n-by-(p*M) and M is the size of support_x
#'
#' @return Matrix of dim n-by-p
#' @export
compute_all_means_efficient <- function(fitted_model, X, conditional_prob_mat,
                                        support_x, lambda = "min"){

  # extract the dimension of the problem
  p <- ncol(X)
  n <- nrow(X)

  # predicted matrix to be imputed
  predicted_mat <- matrix(0, nrow = n, ncol = p)

  # extract the nonzero component in the model
  act_set <- sort(which(as.vector(stats::coef(fitted_model, s = sprintf("lambda.%s", lambda)))[-1] != 0))
  predicted_mat[, dplyr::setdiff(1:p, act_set)] <- stats::predict(fitted_model,
                                                                  s = sprintf("lambda.%s", lambda),
                                                                  newx = X, type = "response")

  # separate the case when act_set is of length zero or not
  if (length(act_set) > 0){

    # loop over the act_set
    integrate_one_fitted <- sapply(act_set, function(act_beta){

      # compute the all the conditional expectations Y|X_{-j} for j in act_set
      fix_one_fitted <- sapply(sort(support_x), function(s){

        # impute the X to be the value in support_x
        x_impute <- X
        x_impute[, act_beta] <- s

        # use predict function to obtain the leave-one-out fitted values
        stats::predict(fitted_model, newx = x_impute, s = sprintf("lambda.%s", lambda),
                       type = "response")

      })

      # compute the integrated prediction
      rowSums(fix_one_fitted * conditional_prob_mat[, act_beta + p * (0 : (length(support_x) - 1))])
    })

    # finish the imputation
    predicted_mat[, act_set] <- integrate_one_fitted
  }

  # return the output
  return(predicted_mat)
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
  lasso_model <- glmnet::cv.glmnet(x = X, y = Y, family = family)

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
    fitted_model[[lambda_param]]$act_set <- act_set

    # perform the glm.fit
    if(length(act_set) == 0){
      glm_fitted <- stats::glm(Y ~ 1, data = data.frame(Y),
                               family = family)
    }else{
      X_act <- X[, sort(act_set), drop = FALSE]
      glm_fitted <- stats::glm(Y ~ ., data = data.frame(Y, X_act),
                               family = family)
    }

    # extract the fitted model
    fitted_model[[lambda_param]]$model <- glm_fitted
  }

  # return the output
  return(fitted_model)
}


#' Compute GCM p-values
#'
#' @param data A list consisting of data, the necessary HMM parameters, file paths and hashing_id
#'
#' @return A list with two resulting p-values obtained using different regularization parameters
#' @export
GCM_HMM <- function(data){

   # extract necessary from data
   X <- data$X
   Y <- data$Y
   n <- nrow(X)
   p <- ncol(X)
   K <- length(data$pInit)
   M <- dim(data$pEmit)[2]
   X_file <- data$X_file
   hashing_id <- data$hashing_id
   fp_path <- data$fp_path
   out_path <- data$out_path

   ########################### estimate hmm parameter ##########################
   # fit the HMM model
   HMM_parameter <- fit_HMM(fp_path = fp_path, out_path = out_path,
                            hashing_id = hashing_id, X_file = X_file, K = K)

   ########################### perform GCM test #################################
   # compute the conditional probability with oracle Q, pEmit and pInit
   conditional_mean <- t(
      apply(X, 1, function(x){
         compute_conditional_prob_Rcpp(x = x, pInit = HMM_parameter$pInit,
                                               pEmit = HMM_parameter$pEmit, Q = HMM_parameter$Q)[, 2]

      })
   )

   # model Y on Z via lasso
   fitted_lasso <- glmnet::cv.glmnet(x = X, y = Y, family = "binomial")

   # loop over lasso_model vector
   pvalue_list <- list(
      min_pvalue = list(NULL),
      ose_pvalue = list(NULL)
   )
   for (lambda_param in data$lasso_model) {
      # use lambada.model for inference
      leave_one_fit <- compute_all_means_efficient(fitted_model = fitted_lasso,
                                                           X = X, conditional_prob_mat = cbind(1 - conditional_mean, conditional_mean),
                                                           support_x = c(0, 1), lambda = lambda_param)

      # compute the p-value
      test_stat_numerator <- (X - conditional_mean)  * (matrix(rep(Y, p), ncol = p) - leave_one_fit)
      test_stat_denominator <- (X - conditional_mean)^2  * (matrix(rep(Y, p), ncol = p) - leave_one_fit)^2

      # compute test_stat
      test_stat <- sqrt(n) * colMeans(test_stat_numerator) / sqrt(colMeans(test_stat_denominator))

      # compute p-value
      left_pvalue <- stats::pnorm(q = test_stat, lower.tail = TRUE)
      right_pvalue <- stats::pnorm(q = test_stat, lower.tail = FALSE)
      both_pvalue <- 2 * apply(cbind(left_pvalue, right_pvalue), 1, function(x) min(x))

      # list to be saved
      list_to_save <- dplyr::if_else(lambda_param == "min", "min_pvalue", "ose_pvalue")

      # save the resutls
      pvalue_list[[list_to_save]] <- list(
         left_pvalue = left_pvalue,
         right_pvalue = right_pvalue,
         both_pvalue = both_pvalue
      )
   }

   # return the output
   return(pvalue_list)
}


#' Compute dCRT p-values
#'
#' @inheritParams GCM_HMM
#'
#' @return A list with two resulting p-values obtained using different regularization parameters
#' @export
#' @importFrom stats rbinom
dCRT_HMM <- function(data){
   # extract necessary from data
   X <- data$X
   Y <- data$Y
   n <- nrow(X)
   p <- ncol(X)
   K <- length(data$pInit)
   M <- dim(data$pEmit)[2]
   X_file <- data$X_file
   hashing_id <- data$hashing_id
   fp_path <- data$fp_path
   out_path <- data$out_path
   B <- data$B

   ########################### estimate hmm parameter ##########################
   # fit the HMM model
   HMM_parameter <- fit_HMM(fp_path = fp_path, out_path = out_path,
                            hashing_id = hashing_id, X_file = X_file, K = K)

   ########################### perform dCRT test ################################
   # compute the conditional probability with oracle Q, pEmit and pInit
   conditional_mean <- t(
      apply(X, 1, function(x){
         compute_conditional_prob_Rcpp(x = x, pInit = HMM_parameter$pInit,
                                               pEmit = HMM_parameter$pEmit,
                                               Q = HMM_parameter$Q)[, 2]

      })
   )

   # model Y on Z via lasso
   fitted_lasso <- glmnet::cv.glmnet(x = X, y = Y, family = "binomial")

   # loop over lasso_model vector
   pvalue_list <- list(
      min_pvalue = list(NULL),
      ose_pvalue = list(NULL)
   )
   for (lambda_param in data$lasso_model) {
      # use lambda_param for inference
      leave_one_fit <- compute_all_means_efficient(fitted_model = fitted_lasso,
                                                           X = X, conditional_prob_mat = cbind(1 - conditional_mean, conditional_mean),
                                                           support_x = c(0, 1), lambda = lambda_param)

      # compute the p-value
      test_stat_numerator <- (X - conditional_mean)  * (matrix(rep(Y, p), ncol = p) - leave_one_fit)

      # compute test_stat
      test_stat <- sqrt(n) * colMeans(test_stat_numerator)

      # do the following resampling
      resampled_test <- matrix(0, nrow = p, ncol = B)
      for (b in 1:B) {

         # resample from conditional_mean
         X_resampled <- matrix(rbinom(n = n * p,
                                      size = 1,
                                      prob = as.vector(conditional_mean)),
                               nrow = n, ncol = p)

         # recompute the test statistic
         resampled_numerator <- (X_resampled - conditional_mean)  * (matrix(rep(Y, p), ncol = p) - leave_one_fit)

         # compute test_stat
         resampled_test[, b] <- sqrt(n) * colMeans(resampled_numerator)

      }

      # compute the p-value
      left_pvalue <- apply(resampled_test - test_stat, 1, function(x) (sum(x <= 0) + 1) / (B + 1) )
      right_pvalue <- apply(resampled_test - test_stat, 1, function(x) (sum(x >= 0) + 1) / (B + 1) )
      both_pvalue <- 2 * apply(cbind(left_pvalue, right_pvalue), 1, function(x) min(x))

      # list to be saved
      list_to_save <- dplyr::if_else(lambda_param == "min", "min_pvalue", "ose_pvalue")

      # save the resutls
      pvalue_list[[list_to_save]] <- list(
         left_pvalue = left_pvalue,
         right_pvalue = right_pvalue,
         both_pvalue = both_pvalue
      )
   }

   # return the output
   return(pvalue_list)
}

#' Compute spaCRT p-values
#'
#' @inheritParams GCM_HMM
#'
#' @return A list with two resulting p-values obtained using different regularization parameters
#' @export
spaCRT_HMM <- function(data){
   # extract necessary from data
   X <- data$X
   Y <- data$Y
   n <- nrow(X)
   p <- ncol(X)
   K <- length(data$pInit)
   M <- dim(data$pEmit)[2]
   X_file <- data$X_file
   hashing_id <- data$hashing_id
   fp_path <- data$fp_path
   out_path <- data$out_path

   ########################### estimate hmm parameter ##########################
   # fit the HMM model
   HMM_parameter <- fit_HMM(fp_path = fp_path, out_path = out_path,
                            hashing_id = hashing_id, X_file = X_file, K = K)

   ########################### perform GCM test #################################
   # compute the conditional probability with oracle Q, pEmit and pInit
   conditional_mean <- t(
      apply(X, 1, function(x){
         compute_conditional_prob_Rcpp(x = x, pInit = HMM_parameter$pInit,
                                               pEmit = HMM_parameter$pEmit, Q = HMM_parameter$Q)[, 2]

      })
   )

   # model Y on Z via lasso
   fitted_lasso <- glmnet::cv.glmnet(x = X, y = Y, family = "binomial")

   # loop over lasso_model vector
   pvalue_list <- list(
      min_pvalue = list(NULL),
      ose_pvalue = list(NULL)
   )
   for (lambda_param in data$lasso_model) {
      # use lambda_param for inference
      leave_one_fit <- compute_all_means_efficient(fitted_model = fitted_lasso,
                                                           X = X, conditional_prob_mat = cbind(1 - conditional_mean, conditional_mean),
                                                           support_x = c(0, 1), lambda = lambda_param)

      # pass the argument to spa_cdf in spacrt package
      p_values_output <- sapply(1:p, function(j){
         suppressWarnings(spacrt:::spa_cdf(X = X[, j],
                                           Y = Y,
                                           X_on_Z_fit_vals = conditional_mean[, j],
                                           Y_on_Z_fit_vals = leave_one_fit[, j],
                                           fam = "binomial",
                                           R = 5))
      })

      # list to be saved
      list_to_save <- dplyr::if_else(lambda_param == "min", "min_pvalue", "ose_pvalue")

      # save the resutls
      pvalue_list[[list_to_save]] <- list(
         left_pvalue = sapply(1:p, function(j) p_values_output[, j]$p.left),
         right_pvalue = sapply(1:p, function(j) p_values_output[, j]$p.right),
         both_pvalue = sapply(1:p, function(j) p_values_output[, j]$p.both),
         GCM_default = sapply(1:p, function(j) !p_values_output[, j]$spa.success)
      )
   }

   # return the output
   return(pvalue_list)
}


#' Compute knockoff rejection sets
#'
#' @inheritParams GCM_HMM
#'
#' @return A list of different rejection sets subject to different choice of regularization parameters
#' @export
knockoff_HMM <- function(data){
   # extract necessary from data
   X <- data$X
   Y <- data$Y
   n <- nrow(X)
   p <- ncol(X)
   K <- length(data$pInit)
   X_file <- data$X_file
   hashing_id <- data$hashing_id
   fp_path <- data$fp_path
   out_path <- data$out_path

   ########################### estimate hmm parameter ###########################
   # create the output directory
   dir.create(out_path, recursive = TRUE)

   # run fastPhase
   fastPhase_new(fp_path = fp_path,
                         X_file = X_file,
                         out_path = out_path,
                         K = K, phased = TRUE)

   # transform the results to R
   estimated_parameter <- SNPknock::loadHMM(r_file = paste0(out_path, "/_rhat.txt"),
                                            alpha_file = paste0(out_path, "/_alphahat.txt"),
                                            theta_file = paste0(out_path, "/_thetahat.txt"),
                                            char_file = paste0(out_path, "/_origchars")
   )

   ########################## perform knockoff ##################################
   # restructure the parameter
   estimated_parameter$r <- estimated_parameter$r[1:p]
   estimated_parameter$alpha <- estimated_parameter$alpha[1:p,]
   estimated_parameter$theta <- estimated_parameter$theta[1:p,]

   # Generate knockoffs
   Xk = SNPknock::knockoffHaplotypes(X, estimated_parameter$r,
                                     estimated_parameter$alpha, estimated_parameter$theta)

   # stack X and Xk
   X_full <- cbind(X, Xk)  # Combine original features and knockoffs

   # Fit Lasso Logistic Regression
   fit <- glmnet::cv.glmnet(x = X_full, y = Y, family = "binomial")

   # loop over lasso_model
   rejection_list <- list(
      min_rejection = list(NULL),
      ose_rejection = list(NULL)
   )
   for (lambda_param in data$lasso_model) {
      # Extract coefficients
      coef_selected <- as.vector(stats::coef(fit,
                                             s = sprintf("lambda.%s", lambda_param)))[-1]

      # Identify which features (original vs knockoff) are selected
      selected_features <- which(coef_selected != 0)

      # Compute importance statistics
      W_stat <- abs(coef_selected[1:p]) - abs(coef_selected[(p + 1):(2 * p)])

      # Apply Knockoff Filter (assuming FDR control at alpha)
      T_threshold <- knockoff::knockoff.threshold(W_stat, fdr = data$alpha)

      # Select features passing threshold
      final_selected <- which(W_stat >= T_threshold)

      # list to be saved
      list_to_save <- dplyr::if_else(lambda_param == "min", "min_rejection", "ose_rejection")

      # save the rejection set
      rejection_list[[list_to_save]] <- list(
         rejection_set = final_selected
      )
   }

   # return the rejection set
   return(rejection_list)
}
