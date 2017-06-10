######################### testParents ##################################
# author : Nicolas Beaume
# created : 30/05/2017
# last update : 09/06/2017
#
# this file contains method to access the object create by defineSpecies()
# from BSL package, as well as some function to test simulation
#
#
#########################################################################
library(BreedingSchemeLanguage)
#********************** functions *******************************
# access to the mapData field of the ith simulation
getMapData <- function(i=1){
  return(lists[[i]]$mapData)
}

# access to the breedingData field of the ith simulation
getBreedingData <- function(i=1){
  lists[[i]]$breedingData
}
#*********************** simulations ****************************

defineSpecies(nSim = 1, effPopSize = 20)
initializePopulation(20)
phenotype()
beforeSelect <- lists
select(nSelect = 10)
afterSelect <- lists
cross()
phenotype()
select(nSelect = 5)
cross()
plot(afterSelect[[1]]$breedingData$GID, afterSelect[[1]]$breedingData$pValue,xlab = "lines index", ylab = "pValue")
sel <- which(afterSelect[[1]]$breedingData$popID == 1)
points(afterSelect[[1]]$breedingData$GID[sel], afterSelect[[1]]$breedingData$pValue[sel], col="red")
points(afterSelect[[1]]$breedingData$GID[sel], afterSelect[[1]]$breedingData$pValue[sel], col="red", pch=16)