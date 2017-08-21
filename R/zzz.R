.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
      paste("\nThis is pRoloc version", packageVersion("pRoloc"), "\n",
            " Visit https://lgatto.github.io/pRoloc/ to get started.\n"))

    if (interactive() && .Platform$OS.type == "windows" &&
        .Platform$GUI == "Rgui") {
        Biobase::addVigs2WinMenu("pRoloc")
    }
}
