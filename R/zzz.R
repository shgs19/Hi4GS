.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("mlr3learners", quietly = TRUE)) {
    warning("Package 'mlr3learners' is not installed.")
  }
  if (!requireNamespace("mlr3extralearners", quietly = TRUE)) {
    warning("Package 'mlr3extralearners' is not installed.")
  }
}

# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("√ mlr3learners")
#   packageStartupMessage("√ mlr3extralearners")
#   packageStartupMessage("√ Hi4GS")
# }
