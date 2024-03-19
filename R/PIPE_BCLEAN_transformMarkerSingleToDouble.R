transMarkerSingle <- function(
    markerDTfile= NULL,
    badCall=NULL,
    genoColumn=NULL,
    firstColum= NULL,
    lastColumn=NULL,
    verbose=FALSE
){
  ## THIS FUNCTION TRANSFORMS MARKER CALLS CODED IN A SINGLE LETTER (E.G., A, T, G, C) INTO A DOUBLE LETTER CODE (E.G., AA, TT, ...) FOR POSTERIOR ANALYSES
  ## IS USED IN THE BCLEAN APP UNDER THE TRANSFORMATION MODULES
  if(is.null(markerDTfile)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(genoColumn)){stop("Please provide the name of the column indicating the genotype id", call. = FALSE)}
  if(is.null(firstColum)){stop("Please provide the name of the column indicating the first marker", call. = FALSE)}
  if(is.null(lastColumn)){stop("Please provide the name of the column indicating the last marker", call. = FALSE)}
  ###################################
  # loading the dataset
  mydata <- markerDTfile

  v1 <- which(colnames(mydata) == firstColum)
  v2 <- which(colnames(mydata) == lastColumn)
  mydata[,v1:v2] <- apply(mydata[,v1:v2],2,function(x){gsub("/","",x)})
  mydata[,v1:v2] <- apply(mydata[,v1:v2],2,function(x){gsub(" ","",x)})

  markerDouble <- apply(mydata[,v1:v2],2,function(x){
    nDigits <- nchar(x)
    toTransform <- which(nDigits == 1)
    if(length(toTransform) > 0){
      x[toTransform] <- paste0(x[toTransform],x[toTransform])
    }
    toDelete <- which(nDigits > 2)
    if(length(toDelete) > 0){
      x[toDelete] <- NA
    }
    return(x)
  })
  if(!is.null(badCall)){
    for(k in badCall){ # k <- badCall[1]
      badOnes <- which(markerDouble == paste0(k,k), arr.ind = TRUE)
      if (nrow(badOnes) > 0){
        markerDouble[badOnes] = NA
      }
    }
  }
  v3 <- which(colnames(mydata) == genoColumn)
  newData <- cbind(mydata[[genoColumn]], markerDouble)

  ##########################################

  return(newData)
}
