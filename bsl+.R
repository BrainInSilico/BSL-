######################## bsl+ ######################
#
# creation date : 09/06/2017
# last update : 09/06/2017
# author : Nicolas Beaume (nicolasbeaume.consultancy@gmail.com)
#
# description : run a simulation of population using BSL (Shiori Yabe, Hiroyoshi Iwata, and Jean-Luc Jannink)
# with extra features
#
#####################################################

#******************** functions *********************

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

installBSL <- function() {
  load.fun("devtools")
  install_github("syabe/BreedingSchemeLanguage")  
}
#******************** main **************************
# test if the package is installed, in case not, install it
if(!isTRUE("BreedingSchemeLanguage" %in% .packages(all.available=TRUE))) {
  installBSL()
}
library(BreedingSchemeLanguage)
#** set some parameters **
initPop <- 100
percentQTL <- 0.2
nMarkers <- 1002
nQTL <- ceiling(percentQTL*nMarker)
nselectedParents <- 50
#** run simulation **

# simulate phenotipic selection
defineSpecies(nSim = 1, nMarkers = nMarkers, nQTL = nQTL)
initializePopulation(nPop = initPop)
phenotype()
select(nSelect = nselectedParents)
cross()
phenotype()
select()
cross()
phenotype()
select()
cross()
plotData()
