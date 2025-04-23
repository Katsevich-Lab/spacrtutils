# for suppressing dependency startup messages
.onLoad <- function(libname, pkgname) {
   ns <- asNamespace(pkgname)

   invisible(suppressMessages(
      sapply(c("MASS", "statmod"),
             requireNamespace, quietly = TRUE)
   ))

   # call estimate_theta() from sceptre
   if (requireNamespace("sceptre", quietly = TRUE)) {
      estimate_theta <- get("estimate_theta", envir = asNamespace("sceptre"))
      assign("estimate_theta", estimate_theta, envir = ns)
   }
}
