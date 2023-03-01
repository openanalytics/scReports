# To configure the active Python environment and 
# install the required Python packages
.onLoad <- function(libname, pkgname) {
  reticulate::configure_environment(pkgname)
}