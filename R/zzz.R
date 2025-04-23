# for suppressing dependency startup messages
.onLoad <- function(libname, pkgname) {
   ns <- asNamespace(pkgname)

   invisible(suppressMessages(
      sapply(c("MASS", "statmod"),
             requireNamespace, quietly = TRUE)
   ))

   # Dynamically retrieve the unexported function from the sceptre namespace
   if (requireNamespace("sceptre", quietly = TRUE)) {
      estimate_theta <- get("estimate_theta", envir = asNamespace("sceptre"))
      assign("estimate_theta", estimate_theta, envir = ns)
   }

   # Dynamically retrieve the unexported functions from the spacrt namespace
   if (requireNamespace("spacrt", quietly = TRUE)) {
      spa_cdf <- get("spa_cdf", envir = asNamespace("spacrt"))
      assign("spa_cdf", spa_cdf, envir = ns)

      nb_precomp <- get("nb_precomp", envir = asNamespace("spacrt"))
      assign("nb_precomp", nb_precomp, envir = ns)
   }
}
