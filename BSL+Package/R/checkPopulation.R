checkPopulation <- function(popIDs, allPopulations=NULL){
  if (is.null(popIDs)) {
    # popIDs <- c(max(allPopulations), (max(allPopulations)-1))
    popIDs <- unique(allPopulations)
  }
  # check for population that does not exist in the data
  toRemove <- NULL
  for(i in 1:length(popIDs)) {
    if (!any(unique(allPopulations) == popIDs[i])) {
      toRemove <- c(toRemove, i)
      warning(paste("cannot use", popIDs[i], "as population of comparison, 
                    it has been removed from the comparison set"))
    }
    }
  if(!is.null(toRemove)) {
    popIDs <- popIDs[-toRemove]
  }
  # case of trying to test only one population
  if(length(popIDs) <= 1) {
    warning("one or less population to test, please check your IDs of population")
    return(NULL)
  }
  # sort to ensure to always compare the most recent population to the oldest one
  return(sort(popIDs))
  }