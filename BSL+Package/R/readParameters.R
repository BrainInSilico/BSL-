# parameter file structure
# Column 1 : description, column 2 : parameter value 
# row 1 : number of simulations
# row 2 : initial population size
# row 3 : number of markers
# row 4 : ratio of QTL
# row 5 set QTL effects (y/n)
# row 6 QTLs effects (optional)
# row 7 : number of parent to select based on phenotype
# row 8 progeny size (each time a progeny is generated)
# row 9 cross cost
# row 10 phenotyping cost
# row 11 genotyping cost
# row 12 selection cost
# row 13 number of RGA generation
# row 14 RGA cost per generation
# row 15 LST cost (total)
# row 16 LST final population
# row 17 OYT cost
# row 18 PYT cost
# row 19 AYT cost
# row 20 heritability
# row 21 fix some markers (y/n)
# row 22 ids of markers to fix
# row 23 values of fixed markers
# row 24 number of markers to use in MAS simultaneously
# row 25 indexes of reference markers for MAS
# row 26 values of reference markers for MAS
# row 27 mas cost

readParameters <- function(paramFile) {
  data <- read.table(paramFile, sep = ":", stringsAsFactors = F)
  nSim <- as.numeric(data[1,2])
  initPop <- as.numeric(data[2,2])
  nMarkers <- as.numeric(data[3,2])
  ratioQTL <- as.numeric(data[4,2])
  nQTL <- ceiling(ratioQTL*nMarkers)
  parentSelected <- as.numeric(data[7,2])
  progenySize <- as.numeric(data[8,2])
  crossCost <- as.numeric(data[9,2])
  phenoCost <- c(Standard=as.numeric(data[10,2]))
  genoCost <- as.numeric(data[11,2])
  selectCost <- as.numeric(data[12,2])
  rgaGeneration <- as.numeric(data[13,2])
  rgaCost <- as.numeric(data[14,2])
  lstCost <- as.numeric(data[15,2])
  lstPop <- as.numeric(data[16,2])
  oytCost <- as.numeric(data[17,2])
  pytCost <- as.numeric(data[18,2])
  aytCost <- as.numeric(data[19,2])
  # heritability
  h2 <- as.numeric(data[20,2])
  # MAS
  refMASid <- as.numeric(unlist(strsplit(data[25,2], split = ";")))
  refMASgeno <- unlist(strsplit(data[26,2], split = ";")) # genotype are 2 values separated by a space
  refMASgeno <- gsub("^\\s?", "", refMASgeno, perl = T)
  refMASgenoMat <- matrix(rep(0,length(refMASgeno)*2), nrow = 2)
  for(i in 1:length(refMASgeno)) {
    refMASgenoMat[,i] <- as.numeric(unlist(strsplit(refMASgeno[i], split=" ")))
  }
  mas <- refMAS(nbRef=as.numeric(data[24,2]), id=refMASid, genoRef=as.data.frame(refMASgenoMat))
  masCost <- as.numeric(data[27,2])
  args <- list(nSim, initPop, nMarkers, nQTL, parentSelected, progenySize, crossCost, phenoCost, genoCost,
               selectCost, rgaGeneration, rgaCost, lstCost, lstPop, oytCost, pytCost, aytCost, h2, mas, masCost)
  names(args) <- c("nSim","initPop", "nMarkers", "nQTL", "parentSelected", "progenySize", "crossCost", "phenoCost",
                   "genoCost", "selectCost", "rgaGeneration", "rgaCost", "lstCost", "lstPop", "oytCost", "pytCost", 
                   "aytCost", "h2", "mas", "masCost")
  # optional QTL effects
  data[5,2] <- tolower(gsub(" ", "", data[5,2], fixed = T))
  if(tolower(data[5,2]) == "y") {
    qtlEffects <- as.numeric(unlist(strsplit(data[6,2], split = ";")))
    n <- names(args)
    args <- c(args, list(qtlEffects))
    names(args) <- c(n, "qtlEffects")
  }
  # optional marker fixation
  data[21,2] <- tolower(gsub(" ", "", data[21,2], fixed = T))
  if(tolower(data[20,2]) == "y") {
    fixedMarkersIDs <- as.numeric(unlist(strsplit(data[22,2], split = ";")))
    if(gsub(" ", "", data[23,2]) == "") {
      fixedMarkersValues <- NULL
    } else {
      fixedMarkersValues <- as.numeric(unlist(strsplit(data[24,2], split = ";")))
      if(length(fixedMarkersValues) != length(fixedMarkersIDs)) {
        warning("number of values for fixed markers is different than number IDs for fixed markers,
                please check these options")
        stop()
      }
    }
    n <- names(args)
    args <- c(args, list(fixedMarkersIDs), list(fixedMarkersValues))
    names(args) <- c(n, "fixedMarkerIDs", "fixedMarkerValues")
  }
  return(args)
}