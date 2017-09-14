# return index of the fixed markers, -1 if no marker is fixed
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param popID population ID to evaluate (default: the last population generated)
getFixedMarkers <- function(sEnv=simEnv, popID=NULL){
  parent.env(sEnv) <- environment()
  if (is.null(popID)){
    popID <- max(sEnv$sims[[1]]$genoRec$popID)  
  }    
  getFixedMarkers.func <- function(data,popID){
    # select population to check for fixation
    index <- which(data$genoRec$popID %in% popID)
    # get index of gametes
    genoIndex <- sort(c(index*2-1,index*2))
    # check for fixation
    fixed <- apply(data$geno[genoIndex,], 2, isFixed)
    if(!any(fixed)){
      fixed <- -1
    } else {
      fixed <- which(fixed)
    }
    return(fixed)
  }
  with(sEnv, {
    # This is too fast to want to parallelize
    return(lapply(sims, getFixedMarkers.func, popID=popID))
  })
}

# determine whether a marker is fixed or not
#'@param marker a vector containg allele
#'@return TRUE if all the values of the vector are the same, false if not
isFixed <- function(marker) {
  return(all(marker == marker[1]))
}