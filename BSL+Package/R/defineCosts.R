#' Define the costs that go into breeding
#' Default for some costs is zero because they probably belong to fixed costs
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param phenoCost named vector: names are the plot types and contents are the cost for one plot of that type (default: Standard=1; if plotTypes have been defined in defineVariances and phenoCost is not specified, all plotTypes are given a cost of 1)
#'@param genoCost scalar: cost to genotype one individual default (0.25)
#'@param crossCost scalar: cost of creating a new individual from a cross (1)
#'@param selfCost scalar: the cost of creating a selfed seed (1)
#'@param doubHapCost scalar: the cost of creating a doubled haploid seed (5)
#'@param predCost scalar: the cost of running the analysis to make predictions (0)
#'@param selectCost scalar: the cost of running the analysis to do selection (0)
#'@param locCost scalar: the cost of maintaining a location for a year (0)
#'@param yearCost scalar: the cost of program upkeep for a year (0)
#'@param rgaCost scalar: cost of creating a new individual from one generation of RGA (1)
#'@param lstCost scalar: cost of inbred observation (3)
#'@param masCost scalar: cost for MAS (2)
#'
#'@return Breeding scheme simulation object supplemented with a list of costs
#'
#'@export
defineCosts <- function(sEnv=simEnv, phenoCost=NULL, genoCost=0.25, crossCost=1, selfCost=1,
                        doubHapCost=5, predCost=0, selectCost=0, locCost=0, yearCost=0,
                        rgaCost=0.5, lstCost=3, oytCost=0, pytCost=0, aytCost=0, masCost=2, nSteps=1){
  parent.env(sEnv) <- environment()
  cost.func <- function(bsl, phenoCost, genoCost, crossCost, selfCost, doubHapCost, predCost, selectCost,
                        locCost, yearCost, rgaCost, lstCost){
    if (is.null(phenoCost)){
      phenoCost <- rep(1, length(bsl$varParms$plotTypeErrVars))
      names(phenoCost) <- names(bsl$varParms$plotTypeErrVars)
    }
    stepCosts <- rep(0, nSteps)
    bsl$costs <- list(phenoCost=phenoCost, genoCost=genoCost, crossCost=crossCost, selfCost=selfCost,
                      doubHapCost=doubHapCost, predCost=predCost, selectCost=selectCost, locCost=locCost,
                      yearCost=yearCost, rgaCost=rgaCost, lstCost=lstCost, oytCost=oytCost, pytCost=pytCost,
                      aytCost=aytCost, masCost=masCost, stepCosts=stepCosts, currentStep=1)
    if (bsl$varParms$randLoc) bsl$totalCost <- 0 else bsl$totalCost <- locCost * ncol(bsl$varParms$locCov)
    return(bsl)
  }
  
  with(sEnv, {
    # This is too fast to want to parallelize
    sims <- lapply(sims, cost.func, phenoCost=phenoCost, genoCost=genoCost, crossCost=crossCost,
                   selfCost=selfCost, doubHapCost=doubHapCost, predCost=predCost, selectCost=selectCost,
                   locCost=locCost, yearCost=yearCost, rgaCost=rgaCost, lstCost)
  })
}
  