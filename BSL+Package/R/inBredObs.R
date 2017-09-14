inBredObs <- function(env=simEnv, ratioSelectedByMAS=0.8, ratioSelectedByPheno=0.8, progenySize) {
  ## MAS
  mas(sEnv=env)
  ## phenotypic selection
  # set a number of individual to select
  # this number is the size of the smallest parent population x ratioProgenySize
  parentSize <- NULL
  for(i in env$nSim) {
      parentSize <- c(parentSize, computePopulationSize(env$sims[[i]]))
  }
  nSelectedByPheno <- min(parentSize)*ratioSelectedByPheno
  phenotype(sEnv = env)
  select(sEnv = env, nSelect=nSelectedByPheno)
  selfFertilize(sEnv = env, nProgeny = progenySize)
}