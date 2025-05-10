library(spacrtutils)

n <- 80; p <- 2; iterations <- 100

X_on_Z_fam <- "binomial"; Y_on_Z_fam <- "poisson"
fitting_X_on_Z <- 'random_forest'; fitting_Y_on_Z <- 'glm'

results_GCM <- results_dCRT <- results_spaCRT <- results_score <- list()
pb <- txtProgressBar(min = 0, max = n, style = 3)

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

   results_score[[i]] <- score.test(data, X_on_Z_fam, Y_on_Z_fam)

   Sys.sleep(0.005)
   setTxtProgressBar(pb, i)
}

df.pvals <- data.frame(
   id = rep(1:iterations, times = 3),
   method = rep(c("GCM", "dCRT", "spaCRT"), each = iterations),
   p_value = c(results_GCM |> lapply(function(lst) lst$p.left |> unname()) |> unlist(),
               results_dCRT |> lapply(function(lst) lst$p.left |> unname()) |> unlist(),
               results_spaCRT |> lapply(function(lst) lst$p.left |> unname()) |> unlist())
)

res.GCM <- df.pvals |>
               dplyr::filter(method == 'GCM' & p_value <= 0.05) |>
               dplyr::summarize(n = dplyr::n()) |>
               dplyr::pull(n) / iterations

res.dCRT <- df.pvals |>
               dplyr::filter(method == 'dCRT' & p_value <= 0.05) |>
               dplyr::summarize(n = dplyr::n()) |>
               dplyr::pull(n) / iterations

res.spaCRT <- df.pvals |>
                  dplyr::filter(method == 'spaCRT' & p_value <= 0.05) |>
                  dplyr::summarize(n = dplyr::n()) |>
                  dplyr::pull(n) / iterations

res.score <- results_score |>
                  lapply(function(elem) elem$p.left) |>
                  unlist() |>
                  {\(x) x[x <= 0.05]}() %>%
                  mean() |>
                  round(3)

writeLines("")

test_that("Hello",{
   expect_equal(res.GCM, 0.08)
   expect_equal(res.dCRT, 0.02)
   expect_equal(res.spaCRT, 0.02)
   expect_equal(res.score, 0.023)
})








