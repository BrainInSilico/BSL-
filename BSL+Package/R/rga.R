rga <- function(sEnv = simEnv, nGeneration=4) {
  parent.env(sEnv) <- environment()
  rga.fun <- function(data, nGeneration=4) {
    for(i in 1:nGeneration) {
      # copy paste of selfFertilize.func() but with a change in cost
      locPos <- data$mapData$map$Pos
      popID <- max(data$genoRec$popID)
      tf <- data$genoRec$popID %in% popID
      GIDpar <- data$genoRec$GID[tf]
      nPar <- length(GIDpar)
      geno <- data$geno[rep(GIDpar * 2, each = 2) + rep(-1:0, 
                                                        nPar), ]
      geno <- makeSelfs(popSize = nPar, geno = geno, pos = locPos)
      pedigree <- cbind(matrix(GIDpar[geno$pedigree], nPar), 
                        0)
      geno <- geno$progenies
      data <- addProgenyData(data, geno, pedigree)
    }
    return(data)
  }
  with(sEnv, {
    if (nCore > 1) {
      sfInit(parallel = T, cpus = nCore)
      sims <- sfLapply(sims, rga.fun, nGeneration)
      sfStop()
    }
    else {
      sims <- lapply(sims, rga.fun, nGeneration)
    }
  })
}