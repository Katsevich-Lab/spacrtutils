library(spacrtutils)

n <- 80; p <- 2; iterations <- 100

X_on_Z_fam <- "binomial"; Y_on_Z_fam <- "poisson"
fitting_X_on_Z <- 'rf'; fitting_Y_on_Z <- 'glm'

results_GCM <- results_dCRT <- results_spaCRT <- list()

for(i in 1:iterations){
   set.seed(i)

   data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
                Y = rpois(n = n, lambda = 1),
                Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))

   results_GCM[[i]] <-
      GCM_internal(data, X_on_Z_fam, Y_on_Z_fam,
                   fitting_X_on_Z = fitting_X_on_Z,
                   fitting_Y_on_Z = fitting_Y_on_Z)

   results_dCRT[[i]] <-
      dCRT_internal(data, X_on_Z_fam, Y_on_Z_fam,
                    fitting_X_on_Z = fitting_X_on_Z,
                    fitting_Y_on_Z = fitting_Y_on_Z,
                    B = 2000)

   results_spaCRT[[i]] <-
      spaCRT_internal(data, X_on_Z_fam, Y_on_Z_fam,
                      fitting_X_on_Z = fitting_X_on_Z,
                      fitting_Y_on_Z = fitting_Y_on_Z)
}

df.pvals <- data.frame(
   id = rep(1:iterations, times = 3),
   method = rep(c("GCM", "dCRT", "spaCRT"), each = iterations),
   p_value = c(results_GCM |> lapply(function(lst) lst$p.left |> unname()) |> unlist(),
               results_dCRT |> lapply(function(lst) lst$p.left |> unname()) |> unlist(),
               results_spaCRT |> lapply(function(lst) lst$p.left |> unname()) |> unlist())
)

df.pvals |>
   dplyr::filter(method == 'GCM' & p_value <= 0.05) |>
   dplyr::summarize(n = dplyr::n()) |>
   dplyr::pull(n) / iterations

df.pvals |>
   dplyr::filter(method == 'dCRT' & p_value <= 0.05) |>
   dplyr::summarize(n = dplyr::n()) |>
   dplyr::pull(n) / iterations

df.pvals |>
   dplyr::filter(method == 'spaCRT' & p_value <= 0.05) |>
   dplyr::summarize(n = dplyr::n()) |>
   dplyr::pull(n) / iterations









