inBredObs <- function(env=simEnv, nSelectedByMAS=NULL, nSelectedByPheno=NULL, finalProgenySize=NULL) {
  costList <- NULL
  # keep track of the initial population size for setiing default finalProgenySize
  parentSize <- computePopulationSize(env$sims[[1]])
  # keep current total cost for further computations
  for(i in 1:env$nSim) {
    costList <- c(costList, list(env$sims[[i]]$totalCost))
  }
  # set default values of nSelectedByMAS
  # default is 80% of the population
  if(is.null(nSelectedByMAS)) {
    curPopSize <- computePopulationSize(env$sims[[1]])
    nSelectedByMAS <- ceiling(curPopSize*0.8)
  }
  # genotype(sEnv = env)
  # select(sEnv = env, nSelect=nSelectedByMAS)
  # set default values of nSelectedByPheno
  # default is 80% of the population
  if(is.null(nSelectedByPheno)) {
    curPopSize <- computePopulationSize(env$sims[[1]])
    nSelectedByPheno <- ceiling(curPopSize*0.8)
  }
  phenotype(sEnv = env)
  select(sEnv = env, nSelect=nSelectedByPheno)
  # set default values of finalProgenySize
  # default is 100 times th initial propulation
  if(is.null(finalProgenySize)) {
    finalProgenySize <- parentSize*100
  }
  selfFertilize(sEnv = env, nProgeny = finalProgenySize)
  for(i in 1:env$nSim) {
    env$sims[[i]]$totalCost <- costList[[i]]+(env$sims[[i]]$costs$lstCost*parentSize)
  }
}