# for suppressing dependency startup messages
.onLoad <- function(libname, pkgname) {
   ns <- asNamespace(pkgname)

   invisible(suppressMessages(
      sapply(c("MASS", "statmod"),
             requireNamespace, quietly = TRUE)
   ))

   # Dynamically retrieve the unexported functions from the spacrt namespace
   if (requireNamespace("spacrt", quietly = TRUE)) {
      estimate_theta <- get("estimate_theta", envir = asNamespace("spacrt"))
      assign("estimate_theta", estimate_theta, envir = ns)

      spa_cdf <- get("spa_cdf", envir = asNamespace("spacrt"))
      assign("spa_cdf", spa_cdf, envir = ns)

      nb_precomp <- get("nb_precomp", envir = asNamespace("spacrt"))
      assign("nb_precomp", nb_precomp, envir = ns)
   }
}
