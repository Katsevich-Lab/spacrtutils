# for suppressing dependency startup messages
.onLoad <- function(libname, pkgname) {
  invisible(suppressMessages(
    sapply(c("MASS", "stats", "sceptre", "statmod"),
           requireNamespace, quietly = TRUE)
  ))
}
