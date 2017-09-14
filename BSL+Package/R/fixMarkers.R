# fix some markers with value provided by the user
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param makerIDs the ID of the markers to fix. The IDs refer to mapData$markerPos
#'@param values maqueur genotype to set to geno
fixMarkers <- function(sEnv=simEnv, markerIDs, values=NULL) {
  parent.env(sEnv) <- environment()
  fixMarkers.func <- function(data, markerIDs, values=NULL) {
    # draw random values
    if(is.null(values)){
      sample(c(-1,1),length(markerIDs), replace = T)
    }
    # deal with issue of length(values) != length(markerIDs) 
    if(length(values) < length(markerIDs)){
      interval <- (length(values)+1):length(markerIDs)
      values[interval] <- sample(c(-1,1),length(interval), replace = T)
      warnings("less values than markers IDs, random values have been assigned to markers without defined ones")
    } else if(length(values) > length(markerIDs)) {
      values <- values[1:length(markerIDs)]
      warnings(paste("more values than markers IDs, only the",length(markerIDs, "first values have beeb considered")))
    }
    values <- lapply(values, rep, nrow(data$geno))
    for(i in 1:length(markerIDs)) {
      data$geno[,markerIDs[i]] <- values[[i]]
    }
    return(data)
  }
  with(sEnv, {
    # This is too fast to want to parallelize
    sims <- lapply(sims, fixMarkers.func, markerIDs=markerIDs, values=values)
  })
}