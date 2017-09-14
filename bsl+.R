######################## bsl+ ######################
#
# creation date : 03/07/2017
# last update : 11/09/2017
# author : Nicolas Beaume (nicolasbeaume.consultancy@gmail.com)
#
# description : run a simulation of population using BSL (Shiori Yabe, Hiroyoshi Iwata, and Jean-Luc Jannink)
# with extra features
#
#####################################################


#******************** main **************************
# test if the package is installed, in case not, install it
if(!isTRUE("BreedingSchemeLanguage" %in% .packages(all.available=TRUE))) {
  installBSL()
}
library(BreedingSchemeLanguage)
print(paste("**** this is BSL+, version ***", packageVersion("BreedingSchemeLanguage")))
library(miscTools)
#** set some parameters **
args <- "bslParameters.txt"
args <- readParameters(args)
#** run simulation **

# initiate simulation
simEnv <- defineSpecies(nSim = args$nSim, nMarkers = args$nMarkers, nQTL = args$nQTL, MAS=args$mas)
# set optional QTL effects
if(any("qtlEffects"%in%names(args))) {
  setQTLeffect(effects = args$qtlEffects)
}
defineVariances()
#defineCosts(phenoCost = args$phenoCost, crossCost = args$crossCost, genoCost = args$genoCost,
            # rgaCost = args$rgaCost, lstCost=args$lstCost, oytCost=args$oytCost, pytCost=args$pytCost,
            # aytCost=args$aytCost, masCost=args$masCost)
initializePopulation(nInd = args$initPop)
# marker fixation
if(any("fixedMarkerIDs"%in%names(args))) {
  fixMarkers(markerIDs=args$fixedMarkerIDs, values=args$fixedMarkerValues)
}
# read breeding schem
schemFile <- "scheme.txt"
schem <- read.table(schemFile, stringsAsFactors = F)
for(i in 1:nrow(schem)) {
  switch(tolower(schem[i,1]),
         ayt={
           print("AYT");
           # compute cost
           for(i in 1:simEnv$nSim) {
             cost <- args$aytCost * computePopulationSize(simEnv$sims[[i]]);
           }
           updateCost(cost=cost);
           # do AYT
           fieldTrials(type="ayt")
           },
         cross={
           print("crossing");
           for(i in 1:simEnv$nSim) {
             cost <- args$crossCost * computePopulationSize(simEnv$sims[[i]]);
           }
           updateCost(cost=cost);
           cross()},
         ebv={
           print("calculating breeding values");
           predictValue()},
         gselection={
           print("genotypic selection");
           for(i in 1:simEnv$nSim) {
             cost <- (args$genoCost+args$parentSelected*args$selectCost) * computePopulationSize(simEnv$sims[[i]]);
           }
           updateCost(cost=cost);
           genotype();
           select(nSelect=args$parentSelected)},
         lst={
           print("LST");
           for(i in 1:simEnv$nSim) {
             cost <- args$lstCost * computePopulationSize(simEnv$sims[[i]]);
           }
           updateCost(cost=cost);
           inBredObs(args$lstPop)},
         mas={
           print("MAS");
           for(i in 1:simEnv$nSim) {
             cost <- args$masCost * computePopulationSize(simEnv$sims[[i]]);
           }
           updateCost(cost=cost);
           mas()},
         oyt={
           print("OYT");
           for(i in 1:simEnv$nSim) {
             cost <- args$oytCost * computePopulationSize(simEnv$sims[[i]]);
           }
           updateCost(cost=cost);
           fieldTrials(type="oyt")},
         pyt={
           print("PYT");
           for(i in 1:simEnv$nSim) {
             cost <- args$pytCost * computePopulationSize(simEnv$sims[[i]]);
           }
           updateCost(cost=cost);
           fieldTrials(type="pyt")},
         pselection={
           print("phenotypic selection");
           for(i in 1:simEnv$nSim) {
             cost <- (args$phenoCost+args$parentSelected*args$selectCost) * computePopulationSize(simEnv$sims[[i]]);
           }
           updateCost(cost=cost);
           phenotype();
           select(nSelect = args$parentSelected)},
         rga={
           print("RGA");
           for(i in 1:simEnv$nSim) {
             cost <- args$rgaCost * computePopulationSize(simEnv$sims[[i]]);
           }
           updateCost(cost=cost);
           rga()},
         warning(paste("unknown option line", i, "IGNORED"))
  )
}
# write outputs
pdf("testBSL.pdf")
plotData()
dev.off()
selectionRate()
geneticGain()
# write total cost, selection rate, genetic gate, heterozygosity and EBV
for(i in 1:simEnv$nSim) {
  computeHeterozygosity(popID = unique(simEnv$sims[[i]]$genoRec$popID))
  heterozygosity <- data.frame(lineID=simEnv$sims[[i]]$genoRec$GID)
  if(exists("predRec", simEnv$sims[[i]])) {
    ebv <- data.frame(lineID=simEnv$sims[[i]]$predRec$predGID,
                      generationID=simEnv$sims[[i]]$predRec$predNo)
  }
  if(exists("predRec", simEnv$sims[[i]])) {
    ebv <- data.frame(ebv, simEnv$sims[[i]]$predRec$predict)
  }
  heterozygosity <- data.frame(heterozygosity, simEnv$sims[[i]]$genoRec$heterozygosityRate)
  write.csv(selRate[[i]], file = paste("selectionRate_sim",i,".csv", sep=""), row.names = T)
  write.csv(genGain[[i]], file = paste("geneticGain_sim",i,".csv", sep=""), row.names = T)
  if(exists("ebv")) {
    colnames(ebv) <- c("lineID", "generationID", "EBV")
    write.csv(ebv, file = paste("breedingValues_sim",i,".csv"), row.names = F)
  }
  colnames(heterozygosity) <- c("lineID", "heterozygosity")
  write.csv(heterozygosity, file = paste("heterozygosity_sim",i,".csv"), row.names = F)
}