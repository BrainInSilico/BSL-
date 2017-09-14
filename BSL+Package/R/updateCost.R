#'update costs for each simulations
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param cost a vector of costs for the current breeding action. There must be as many costs as simulations.
updateCost <- function(sEnv= simEnv, cost=0) {
  parent.env(sEnv) <- environment()
  with(sEnv, {
    for(i in 1:nSim) {
      sims[[i]]$totalCost <- simEnv$sims[[i]]$totalCost + cost[i];
      sims[[i]]$costs[currentStep] <- cost;
      currentStep <- currentStep + 1 
    }
  })
}