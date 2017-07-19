# parameter file structure
# Column 1 : description, column 2 : parameter value 
# row 1 : number of simulations
# row 2 : initial population size
# row 3 : number of markers
# row 4 : ratio of QTL
# row 5 : number of parent to select based on phenotype
# row 6 progeny size (each time a progeny is generated)
# row 7 cross cost
# row 8 phenotyping cost
# row 9 genotyping cost
# row 10 set QTL effects (y/n)
# row 11 QTLs effects (optional)
# row 12 number of RGA generation
# row 13 RGA cost per generation
# row 14 LST cost (total)
# row 15 OYT cost
# row 16 PYT cost
# row 17 AYT cost
# row 19 heritability

readParameters <- function(paramFile) {
  data <- read.table(paramFile, sep = ":", stringsAsFactors = F)
  nSim <- as.numeric(data[1,2])
  initPop <- as.numeric(data[2,2])
  nMarkers <- as.numeric(data[3,2])
  ratioQTL <- as.numeric(data[4,2])
  parentSelected <- as.numeric(data[5,2])
  progenySize <- as.numeric(data[6,2])
  crossCost <- as.numeric(data[7,2])
  phenoCost <- c(Standard=as.numeric(data[8,2]))
  genoCost <- as.numeric(data[9,2])
  nQTL <- ceiling(ratioQTL*nMarkers)
  rgaGeneration <- as.numeric(data[12,2])
  rgaCost <- as.numeric(data[13,2])
  lstCost <- as.numeric(data[14,2])
  oytCost <- as.numeric(data[15,2])
  pytCost <- as.numeric(data[16,2])
  aytCost <- as.numeric(data[17,2])
  # heritability
  h2 <- 1-as.numeric(unlist(data[18,2]))
  args <- list(nSim, initPop, nMarkers, nQTL, parentSelected, progenySize, crossCost, phenoCost, genoCost,
               rgaGeneration, rgaCost, lstCost, oytCost, pytCost, aytCost, h2)
  names(args) <- c("nSim","initPop", "nMarkers", "nQTL", "parentSelected", "progenySize", "crossCost", "phenoCost",
                   "genoCost", "rgaGeneration", "rgaCost", "lstCost", "oytCost", "pytCost", "aytCost", "h2")
  # optional QTL effects
  data[10,2] <- tolower(gsub(" ", "", data[10,2], fixed = T))
  if(tolower(data[10,2]) == "y") {
    qtlEffects <- as.numeric(unlist(strsplit(data[11,2], split = ";")))
    n <- names(args)
    args <- c(args, list(qtlEffects))
    names(args) <- c(n, "qtlEffects")
  }
  return(args)
}