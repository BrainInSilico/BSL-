# install BSL package
installBSL <- function() {
  load.fun("devtools")
  library(devtools)
  devtools::install_github("BrainInSilico/BreedingSchemeLanguage/")  
}

# install a package
# adapted from http://r.789695.n4.nabble.com/Install-package-automatically-if-not-there-td2267532.html
load.fun <- function(x) {
  print(paste("*** installing package", x, "****"))
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text=paste("require(", x, ")", sep="")))
  } else {
    install.packages(x, repos="https://cloud.r-project.org/")
    eval(parse(text=paste("require(", x, ")", sep="")))
  }
}