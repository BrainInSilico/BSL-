# compute selection rate for all simulations
selectionRate <- function(sEnv = simEnv, pops=NULL) {
  parent.env(sEnv) <- environment()
  # compute selection rate for each simulations
  with(sEnv, {
    if (nCore > 1) {
      sfInit(parallel = T, cpus = simEnv$nCore)
      selRate <<- sfLapply(simEnv$sims, selectionRate.fun, pops)
      sfStop()
    }
    else {
      selRate <<- lapply(simEnv$sims, selectionRate.fun, pops)
    }
  })
}

selectionRate.fun <- function(data, popIDs=NULL) {
  popIDs <- checkPopulation(popIDs, data$genoRec$popID)
  if(!is.null(popIDs)) {
    # case of only 2 populations to compare
    if(length(popIDs) ==2) {
      sr <- length(which(data$genoRec$popID == popIDs[2]))/length(which(data$genoRec$popID == popIDs[1]))
      names(sr) <- paste(popIDs[2], "_", popIDs[1], sep="")
    } else {
      # find combination to test
      toTest <- combn(popIDs, 2)
      # create a matrix to store results
      sr <- matrix(rep(NA, length(popIDs)^2), ncol=length(popIDs))
      for(i in 1:ncol(toTest)) {
        # locate populations to consider
        # in toTest, combination are in columns and the second line
        # is the most recent population
        res <- length(which(data$genoRec$popID == toTest[2,i]))/length(which(data$genoRec$popID == toTest[1,i]))
        # find the indexes of the populations implied in the comparison
        index1 <- which(popIDs == toTest[1,i])
        index2 <- which(popIDs == toTest[2,i])
        sr[index2, index1] <- res
        sr[index1,index2] <- 1/res # reverse comparison to allow the reader to find back the order of comparison
      }
      colnames(sr) <- paste("p", popIDs, sep = "")
      rownames(sr) <- paste("p", popIDs, sep = "")
      diag(sr) <- 0
    }
    return(sr)
  } else {
    warning("bad popIDs, selection rate cannot be computed")
    return(NULL)
  }
}