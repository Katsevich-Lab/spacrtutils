# for suppressing dependency startup messages
.onLoad <- function(libname, pkgname) {
   ns <- asNamespace(pkgname)

   invisible(suppressMessages(
      sapply(c("MASS", "statmod"),
             requireNamespace, quietly = TRUE)
   ))

   # Dynamically retrieve the unexported functions from the spacrt namespace
   if (requireNamespace("spacrt", quietly = TRUE)) {
      spa_cdf <- get("spa_cdf", envir = asNamespace("spacrt"))
      assign("spa_cdf", spa_cdf, envir = ns)

      nb_precomp <- get("nb_precomp", envir = asNamespace("spacrt"))
      assign("nb_precomp", nb_precomp, envir = ns)
   }
}
