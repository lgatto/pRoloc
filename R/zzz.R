.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("\nThis is pRoloc version",packageVersion("pRoloc"),"\n"))
}

