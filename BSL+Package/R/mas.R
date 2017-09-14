######################## MAS ######################
#
# creation date : 01/08/2017
# last update : 23/08/2017
# author : Nicolas Beaume (nicolasbeaume.consultancy@gmail.com)
#
# description : fonction to do a Marker Assisted Selection
#
#####################################################
mas <- function(sEnv=simEnv, popID=NULL) {
  parent.env(sEnv) <- environment()
  mas.func <- function(bsl, popID=NULL, masObject=NULL) {
    # get all elements
    if(is.null(popID)){
      popID <- max(bsl$genoRec$popID)
    }
    pop <- which(bsl$genoRec$popID %in% popID)
    pop <- sort(c(pop*2-1,pop*2))
    # check if ref marker are fixed
    refMarkerID <- getReferenceMarkersID(masObject)
    MAS <<- update(masObject, bsl$geno[pop,refMarkerID])
    if(!is.null(masObject)) {
      # select individuals
      refMarkers <- getReferenceMarkers(masObject)
      selectedGID <- getIndividualWithGoodMarkers(bsl$geno[pop,refMarkerID], refMarkers)
      if(!is.null( selectedGID)) {
        popID.new <- max(bsl$genoRec$popID) + 1
        bsl$genoRec$popID[bsl$genoRec$GID %in% selectedGID] <- popID.new
      } else{
        warning(paste("for population", popID," no individuals where selected by MAS, no selection was performed"))
      }
    } else {
      warning("no more reference marker. MAS was not executed")
    }
    return(bsl)
  }
  with(sEnv, {
    if(!is.null(MAS)) {
      if(nCore > 1){
        sfInit(parallel=T, cpus=nCore)
        sims <- sfLapply(sims, mas.func, popID=popID, masObject=MAS)
        sfStop()
      }else{
        sims <- lapply(sims, mas.func, popID=popID, masObject=MAS)
      }
    }
  })
}

compareToRef <- function(geno, refMarkers) {
  # direct but restrictive comparison of geno vs marker
  isEqual <- all(geno == refMarkers)
  # if it is FALSE then test the opposite
  if(!isEqual) {
    temp <- geno[1,]
    geno[1,] <- geno[2,]
    geno[2,] <- temp
    isEqual <- all(geno == refMarkers)
  }
  return(isEqual)
}

getIndividualWithGoodMarkers <- function(genotypes, refMarkers) {
  # transform the matrix into a list of genotypes
  genoList <- NULL
  for(i in seq(2,nrow(genotypes),2)) {
    genoList <- c(genoList, list(genotypes[(i-1):i,]))
  }
  goodMarkers <- unlist(lapply(genoList, compareToRef, refMarkers))
  if(any(goodMarkers)){
    goodMarkers <- which(goodMarkers)
  } else {
    goodMarkers <- NULL
  }
  return(goodMarkers)
}