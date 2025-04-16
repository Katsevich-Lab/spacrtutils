# for suppressing dependency startup messages
.onLoad <- function(libname, pkgname) {
   ns <- asNamespace(pkgname)

   invisible(suppressMessages(
      sapply(c("MASS", "statmod"),
             requireNamespace, quietly = TRUE)
   ))

   if (requireNamespace("sceptre", quietly = TRUE)) {
      estimate_theta <- get("estimate_theta", envir = asNamespace("sceptre"))
      assign("estimate_theta", estimate_theta, envir = ns)
   }
}
