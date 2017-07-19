######################## bsl+ ######################
#
# creation date : 03/07/2017
# last update : 10/07/2017
# author : Nicolas Beaume (nicolasbeaume.consultancy@gmail.com)
#
# description : run a simulation of population using BSL (Shiori Yabe, Hiroyoshi Iwata, and Jean-Luc Jannink)
# with extra features
#
#####################################################

installBSL <- function() {
  load.fun("devtools")
  library(devtools)
  devtools::install_github("BrainInSilico/BreedingSchemeLanguage/tree/nicolas")  
}

# install a package
# adapted from http://r.789695.n4.nabble.com/Install-package-automatically-if-not-there-td2267532.html
load.fun <- function(x) {
  print(paste("*** installing package", x, "****"))
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text=paste("require(", x, ")", sep="")))
  } else {
    install.packages(x, repos="https://cloud.r-project.org/")
    eval(parse(text=paste("require(", x, ")", sep="")))
  }
}

#******************** main **************************
# test if the package is installed, in case not, install it
if(!isTRUE("BreedingSchemeLanguage" %in% .packages(all.available=TRUE))) {
  installBSL()
}
library(BreedingSchemeLanguage)
library(miscTools)
#** set some parameters **
args <- "bslParameters.txt"
args <- readParameters(args)
#** run simulation **

# initiate simulation
simEnv <- defineSpecies(nSim = args$nSim, nMarkers = args$nMarkers, nQTL = args$nQTL)
# set optional QTL effects
if(any("qtlEffects"%in%names(args))) {
  setQTLeffect(effects = args$qtlEffects)
}
defineVariances()
defineCosts(phenoCost = args$phenoCost, crossCost = args$crossCost, genoCost = args$genoCost,
            rgaCost = args$rgaCost, lstCost=args$lstCost, oytCost=args$oytCost, pytCost=args$pytCost,
            aytCost=args$aytCost)
initializePopulation(nInd = args$initPop)
# read breeding schem
schemFile <- "scheme.txt"
schem <- read.table(schemFile, stringsAsFactors = F)
for(i in 1:nrow(schem)) {
  switch(tolower(schem[i,1]),
         cross=cross(),
         rga=rga(),
         lst=inBredObs(),
         oyt=fieldTrials(type="oyt"),
         pyt=fieldTrials(type="pyt"),
         ayt=fieldTrials(type="ayt"),
         pselection={phenotype();
           select(nSelect = args$parentSelected)},
         gselection={genotype();
           select(nSelect=args$parentSelected)},
         ebv=predictValue(),
         warning(paste("unknown option line", i, "IGNORED"))
         )
}
# write outputs
pdf("testBSL.pdf")
plotData()
dev.off()
selectionRate()
geneticGain()
computeHeterozygosity(popID = unique(simEnv$sims[[1]]$genoRec$popID))
# write total cost, selection rate, genetic gate, heterozygosity and EBV
totalCost <- data.frame(simulation=paste("simulation", 1:simEnv$nSim), totalCost=rep(-1, simEnv$nSim))
if(exists("predRec", simEnv$sims[[1]])) {
  ebv <- data.frame(lineID=simEnv$sims[[1]]$predRec$predGID,
                    generationID=simEnv$sims[[1]]$predRec$predNo)
}
heterozygosity <- data.frame(lineID=simEnv$sims[[1]]$genoRec$GID)
for(i in 1:simEnv$nSim) {
  totalCost[i,2] <- simEnv$sims[[i]]$totalCost
  if(exists("predRec", simEnv$sims[[i]])) {
    ebv <- data.frame(ebv, simEnv$sims[[i]]$predRec$predict)
  }
  heterozygosity <- data.frame(heterozygosity, simEnv$sims[[i]]$genoRec$heterozygosityRate)
  write.csv(selRate[[i]], file = paste("selectionRate_sim",i,".csv", sep=""), row.names = T)
  write.csv(genGain[[i]], file = paste("geneticGain_sim",i,".csv", sep=""), row.names = T)
}
write.csv(totalCost, file = "totalCost.csv", row.names = F)
if(exists("ebv")) {
  colnames(ebv) <- c("lineID", "generationID", paste("sim", 1:simEnv$nSim, "_EBV", sep=""))
  write.csv(ebv, file = "breedingValues.csv", row.names = F)
}
colnames(heterozygosity) <- c("lineID", paste("sim", 1:simEnv$nSim, "_heterozygosity", sep=""))
write.csv(heterozygosity, file = "heterozygosity.csv", row.names = F)
