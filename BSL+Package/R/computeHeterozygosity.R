# compute heterozygositie inside a single genome
#'@param IDs the geno field of a simulation
#'@param genome a vector with 2 elements : the indexes of the gamete1 and gamete2 of an individual
computeHeterozygosity.core <- function(IDs, genome) {
  # extarct gamete
  gamete1 <- genome[IDs[1],]
  gamete2 <- genome[IDs[2],]
  # compute heterozygosity
  return(length(which(gamete1 != gamete2))/length(gamete1))
}

# compute heterozygositie between all genome of a population
#'@param genomes the "geno" field of a simulation
computeHeterozygosity.fun <- function(bsl, popID) {
  # extract genomes
  genomes <- bsl$geno
  if(is.null(popID)){
    popID <- max(bsl$genoRec$popID)
  }
  id <- bsl$genoRec$popID %in% popID
  GID <- bsl$genoRec$GID[id]
  nb <- length(GID)
  geno <- bsl$geno[rep(GID*2, each=2) + rep(-1:0, nb), ]
  # paste gamete1 and gamete2 of a single genome :
  # gamete1 are odd row indexes while gamete2 are even row indexes
  indexes <- rbind(1:nb*2-1, 1:nb*2)
  # add a "heterozygosityRate" field in genoRec
  bsl$genoRec <- data.frame(bsl$genoRec, heterozygosityRate=rep(NA, nrow(bsl$genoRec)))
  # apply computeHeterozygosity.core to all genomes
  bsl$genoRec$heterozygosityRate[id] <- unlist(apply(indexes, 2, computeHeterozygosity.core, genomes))
  bsl$selCriterion <- list(popID=popID, criterion="hetz")
  return(bsl)
}

# run computeHeterozygosity.fun over all simulation and a selected population
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param popID population ID to which heterozygosity must be computed. Default: the last population

computeHeterozygosity <- function(sEnv=simEnv, popID=NULL) {
  parent.env(sEnv) <- environment()
  with(sEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, computeHeterozygosity.fun, popID)
      sfStop()
    } else{
      sims <- lapply(sims, computeHeterozygosity.fun, popID)
    }
  })
}
