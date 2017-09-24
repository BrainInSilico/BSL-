# parameter file structure
# Column 1 : description, column 2 : parameter value 
# row 1 : ratio of QTL
# row 2 set QTL effects (y/n)
# row 3 QTLs effects (optional)
# row 4 cross cost
# row 5 phenotyping cost
# row 6 genotyping cost
# row 7 selection cost
# row 8 RGA cost per generation
# row 9 LST cost (total)
# row 10 OYT cost
# row 11 PYT cost
# row 12 AYT cost
# row 13 heritability
# row 14 fix some markers (y/n)
# row 15 ids of markers to fix
# row 16 values of fixed markers
# row 17 number of markers to use in MAS simultaneously
# row 18 indexes of reference markers for MAS
# row 19 values of reference markers for MAS
# row 20 mas cost

readParameters <- function(paramFile) {
  data <- read.table(paramFile, sep = ":", stringsAsFactors = F)
  nQTL <- as.numeric(data[1,2])
  crossCost <- as.numeric(data[4,2])
  phenoCost <- c(Standard=as.numeric(data[5,2]))
  genoCost <- as.numeric(data[6,2])
  selectCost <- as.numeric(data[7,2])
  rgaCost <- as.numeric(data[8,2])
  lstCost <- as.numeric(data[9,2])
  oytCost <- as.numeric(data[10,2])
  pytCost <- as.numeric(data[11,2])
  aytCost <- as.numeric(data[12,2])
  # heritability
  h2 <- as.numeric(data[13,2])
  # MAS
  refMASid <- as.numeric(unlist(strsplit(data[18,2], split = ";")))
  refMASgeno <- unlist(strsplit(data[19,2], split = ";")) # genotype are 2 values separated by a space
  refMASgeno <- gsub("^\\s?", "", refMASgeno, perl = T)
  refMASgenoMat <- matrix(rep(0,length(refMASgeno)*2), nrow = 2)
  for(i in 1:length(refMASgeno)) {
    refMASgenoMat[,i] <- as.numeric(unlist(strsplit(refMASgeno[i], split=" ")))
  }
  mas <- refMAS(nbRef=as.numeric(data[17,2]), id=refMASid, genoRef=as.data.frame(refMASgenoMat))
  masCost <- as.numeric(data[20,2])
  args <- list(nQTL, crossCost, phenoCost, genoCost, selectCost, rgaCost, lstCost, oytCost, pytCost,
               aytCost, h2, mas, masCost)
  names(args) <- c("nQTL", "crossCost", "phenoCost", "genoCost", "selectCost","rgaCost", "lstCost",
                   "oytCost", "pytCost", "aytCost", "h2", "mas", "masCost")
  # optional QTL effects
  data[2,2] <- tolower(gsub(" ", "", data[2,2], fixed = T))
  if(tolower(data[2,2]) == "y") {
    qtlEffects <- as.numeric(unlist(strsplit(data[3,2], split = ";")))
    n <- names(args)
    args <- c(args, list(qtlEffects))
    names(args) <- c(n, "qtlEffects")
  }
  # optional marker fixation
  data[14,2] <- tolower(gsub(" ", "", data[14,2], fixed = T))
  if(tolower(data[14,2]) == "y") {
    fixedMarkersIDs <- as.numeric(unlist(strsplit(data[15,2], split = ";")))
    if(gsub(" ", "", data[16,2]) == "") {
      fixedMarkersValues <- NULL
    } else {
      fixedMarkersValues <- as.numeric(unlist(strsplit(data[15,2], split = ";")))
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