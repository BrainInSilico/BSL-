# compute the size of a population
computePopulationSize <- function(data, popID=NULL){
  if(is.null(popID)) {
    popID <- max(data$genoRec$popID)
  }
  tf <- data$genoRec$popID %in% popID
  GIDpar <- data$genoRec$GID[tf]
  return(length(GIDpar))
}