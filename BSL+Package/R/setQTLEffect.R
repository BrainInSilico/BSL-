setQTLeffect <- function(sEnv=simEnv, effects){
  parent.env(sEnv) <- environment()
  setQTLeffect.fun <- function(bsl, effects) {
    nQTL <- length(bsl$mapData$effects)
    nEffects <- length(effects) 
    # if more effect than number of QTL provided, just set the first n and send a warning
    if(nEffects > nQTL) {
      warning(paste("more effect define than QTLs, only the", nEffects, "first effects will be used"))
      effects <- effects[1:nQTL]
    }
    bsl$mapData$effects[1:nEffects] <- effects
    return(bsl)
  }
  # manage different simulation with potential multiple cores
  with(sEnv, {  
    if (nCore > 1) {
      sfInit(parallel = T, cpus = nCore)
      sims <- sfLapply(sims, setQTLeffect.fun, effects=effects)
      sfStop()
    } else {
      sims <- lapply(sims, setQTLeffect.fun, effects=effects)
    }
  })
}