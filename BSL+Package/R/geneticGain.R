# compute mean genetic value of a population for all current simulations
# by default, the population used is the last population created
geneticGain <- function(sEnv=simEnv, pops=NULL){
  parent.env(sEnv) <- environment()
  # compute genetic gain for each simulations
  with(sEnv, {
    if (nCore > 1) {
      sfInit(parallel = T, cpus = simEnv$nCore)
      genGain <<- sfLapply(simEnv$sims, geneticGain.fun, pops)
      sfStop()
    }
    else {
      genGain <<- lapply(simEnv$sims, geneticGain.fun, pops)
    }
  })
  
}

# compute the genetic gain between 2 populations
# by default, the last population created and the previous populations
geneticGain.fun <- function(data, popIDs=NULL) {
  popIDs <- checkPopulation(popIDs, data$genoRec$popID)
  if(!is.null(popIDs)) {
    # case of only 2 populations to compare
    if(length(popIDs) ==2) {
      # locate populations to consider
      pop1 <- rep(F, length(data$genoRec$GID))
      pop2 <- rep(F, length(data$genoRec$GID))
      pop1[data$genoRec$popID == popIDs[2]] <- T  # in toTest, combination are in columns and the second line
      # is the most recent population
      pop2[data$genoRec$popID == popIDs[1]] <- T
      gg <- mean(data$gValue[pop1]) - mean(data$gValue[pop2])
      names(gg) <- paste(popIDs[2], "_", popIDs[1], sep="")
    } else {
      # find combination to test
      toTest <- combn(popIDs, 2)
      # create a matrix to store results
      gg <- matrix(rep(NA, length(popIDs)^2), ncol=length(popIDs))
      for(i in 1:ncol(toTest)) {
        # locate populations to consider
        pop1 <- rep(F, length(data$genoRec$GID))
        pop2 <- rep(F, length(data$genoRec$GID))
        pop1[data$genoRec$popID == toTest[2,i]] <- T  # in toTest, combination are in columns and the second line
        # is the most recent population
        pop2[data$genoRec$popID == toTest[1,i]] <- T
        res <- mean(data$gValue[pop1]) - mean(data$gValue[pop2])
        # find the indexes of the populations implied in the comparison
        index1 <- which(popIDs == toTest[1,i])
        index2 <- which(popIDs == toTest[2,i])
        gg[index2, index1] <- res
        gg[index1,index2] <- -res # reverse comparison to allow the reader to find back the order of comparison
      }
      colnames(gg) <- paste("p", popIDs, sep = "")
      rownames(gg) <- paste("p", popIDs, sep = "")
      diag(gg) <- 0
    }
    return(gg)
  } else {
    warning("bad popIDs, genetic gain cannot be computed")
    return(NULL)
  }
}
