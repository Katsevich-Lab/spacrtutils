# for suppressing dependency startup messages
.onLoad <- function(libname, pkgname) {
  invisible(suppressPackageStartupMessages(
    sapply(c("MASS", "stats", "sceptre", "statmod"),
           requireNamespace, quietly = TRUE)
  ))
}
