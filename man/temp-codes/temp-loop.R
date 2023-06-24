
library(magrittr)

n <- 20; p <- 2; B <- 100; normalize <- FALSE; return_resamples <- FALSE
X_on_Z_fam <- "binomial"
Y_on_Z_fam <- "poisson"

expit <- function(x) return(exp(x)/(1+exp(x)))
results <- list()

for(k in 1:100){
  data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
               Y = rpois(n = n, lambda = 1),
               Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))

  results[[k]] <- dCRT(data, X_on_Z_fam, Y_on_Z_fam, B, normalize, return_resamples)

  print(k)
}

results %>% lapply(function(val) val$p_value) %>% unlist() %>% hist()
