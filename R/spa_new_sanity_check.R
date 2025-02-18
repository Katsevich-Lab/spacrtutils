set.seed(125)

n <- 40; p <- 1; normalize <- FALSE; return_cdf <- FALSE

gamma_0 = -6              # intercept in X on Z model
gamma_1 = 1               # slope in X on Z model
beta_0 = -5               # intercept in Y on Z model
beta_1 = 1                # slope in Y on Z model
theta = 5

expit <- function(theta)(exp(theta)/(1+exp(theta)))

generate_data_nb <- function(n,
                             gamma_0 = gamma_0, gamma_1 = gamma_1,
                             beta_0 = beta_0, beta_1 = beta_1,
                             theta = theta){

  Z <- as.matrix(rnorm(n = n, mean = 0, sd = 1))
  X <- rbinom(n = n, size = 1, prob = expit(gamma_0 + gamma_1*Z))
  Y <- MASS::rnegbin(n = n, mu = exp(beta_0 + beta_1*Z), theta = theta)

  return(list(X = X, Y = Y, Z = Z))
}


# data <- list(X = rbinom(n = n, size = 1, prob = 0.1),
#              Y = MASS::negative.binomial(n = n, ),
#              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))

data <- generate_data_nb(n = n, gamma_0 = gamma_0, gamma_1 = gamma_1,
                         beta_0 = beta_0, beta_1 = beta_1,
                         theta = theta)

X_on_Z_fam <- "binomial"
Y_on_Z_fam <- "negative.binomial"
(results_old <- spaCRT_old(data, X_on_Z_fam, Y_on_Z_fam, normalize))
(results_new <- spaCRT(data, X_on_Z_fam, Y_on_Z_fam, normalize))
