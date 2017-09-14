######################## refMAS ######################
#
# creation date : 12/08/2017
# last update : 12/08/2017
# author : Nicolas Beaume (nicolasbeaume.consultancy@gmail.com)
#
# description : R object for reference marker used in MAS
#
#####################################################
library(methods)

#**************** class *********************
refMAS <- setClass("refMAS",
  slots = c(nbRef = "numeric", id = "numeric", genoRef = "data.frame")
)
#**************** validity ******************
validRefMAS <- function(object) {
  valid <- TRUE
  if(length(object@id) != ncol(object@genoRef)){
    valid <- FALSE
    warning("number of id is not equal to the number of possible markers")
  }
  if(length(object@nbRef) > length(object@id)) {
    object@nbRef <- object@nbRef[1:length(object@id)]
  }
  if(nrow(object@genoRef) != 2) {
    # case of more than two alleles
    if(nrow(object@genoRef) > 2) {
      warning("too much allele for markers, only two alleles are allowed, extra allele will be ignored")
      object@genoRef <- object@genoRef[,1:2]
    # case of less than 2 alleles
    } else {
      warning("not enough allele for markers, two alleles are needed, please modify reference markers accordingly")
      valid <- FALSE
    }
  }
  return(valid)
}

setValidity("refMAS", validRefMAS)

#**************** methods ******************

setGeneric(
  "getReferenceMarkers",
  function(obj){return(obj@genoRef)}
)

setMethod(
  "getReferenceMarkers",
  "refMAS",
  function(obj){
    return(obj@genoRef[,1:obj@nbRef])
  }
)

setGeneric(
  "getReferenceMarkersID",
  function(obj){return(obj@id)}
)

setMethod(
  "getReferenceMarkersID",
  "refMAS",
  function(obj){
    return(obj@id[1:obj@nbRef])
  }
)


setGeneric("removeMarkers",
  function(obj, indexes){return(obj)}
)

setMethod("removeMarkers",
          "refMAS",
          function(obj, indexes) {
            # deal with the number of ref marker to remove
            # if more than one marker remains, everything is fine
            if(length(indexes) < (ncol(obj@genoRef)-1)) {
              obj@id <- obj@id[-which(obj@id%in%indexes)]
              obj@genoRef <- obj@genoRef[,-indexes]
            # only one marker remains : create a new dataframe with one column
            } else if(length(indexes) == (ncol(obj@genoRef)-1)) {
              obj@id <- obj@id[-which(obj@id%in%indexes)]
              obj@genoRef <- data.frame("m"=obj@genoRef[,-indexes])
            # no marker remains
            } else {
              obj <- NULL
            }
            # test if the number of remaining marker is > number of reference
            if(!is.null(obj) & length(obj@id) < obj@nbRef) {
              obj@nbRef <- length(obj@id)
            }
            return(obj)
          }
)

setGeneric("update",
  function(obj, genotypes){return(obj)}
)

setMethod("update",
  "refMAS",
  function(obj, genotypes){
    fixed <- unlist(apply(genotypes,2,isFixed))
    if(any(fixed)) {
      fixed <- which(fixed)
      # update refMarker list
      obj <- removeMarkers(obj, fixed)
    }
    return(obj)
  }
)

