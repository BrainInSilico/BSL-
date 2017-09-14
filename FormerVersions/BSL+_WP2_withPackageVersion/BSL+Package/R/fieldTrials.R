fieldTrials <- function(env=simEnv, type="oyt") {
  costList <- NULL
  # keep track of the initial population size for setiing default finalProgenySize
  parentSize <- computePopulationSize(env$sims[[1]])
  # keep current total cost for further computations
  for(i in 1:env$nSim) {
    costList <- c(costList, list(env$sims[[i]]$totalCost))
  }
  # do something to reflect effect of OYT on the population
  #### coming soon ####

  # update total cost
  cost <- switch(type,
                 oyt = env$sims[[1]]$costs$oytCost,
                 pyt = env$sims[[1]]$costs$pytCost,
                 ayt = env$sims[[1]]$costs$aytCost,
                 0)
  for(i in 1:env$nSim) {
    env$sims[[i]]$totalCost <- costList[[i]]+(cost*parentSize)
  }
}