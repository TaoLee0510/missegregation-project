ensure_packages <- function(pkgs) {
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(to_install) > 0) {
    message("Installing missing packages: ", paste(to_install, collapse = ", "))
    install.packages(to_install)
  }
  
  failed <- pkgs[!sapply(pkgs, require, character.only = TRUE, quietly = TRUE)]
  if (length(failed) > 0) {
    stop("The following packages could not be loaded: ", paste(failed, collapse = ", "))
  }
}
