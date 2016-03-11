.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
      paste("\nThis is pRoloc version", packageVersion("pRoloc"), "\n",
            " Read '?pRoloc' and references therein for information\n",
            " about the package and how to get started.\n"))

    if (interactive() && .Platform$OS.type == "windows" &&
        .Platform$GUI == "Rgui") {
        Biobase::addVigs2WinMenu("pRoloc")
    }
}
