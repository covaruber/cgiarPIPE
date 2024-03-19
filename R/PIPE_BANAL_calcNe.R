calcNe <- function(
    markerDTfile= NULL,
    maxNe=100, 
    maxMarker=1000, 
    nSamples=5,
    verbose=FALSE
){
  ## FUNCTION TO CALCULATE EFFECTIVE POPULATION SIZE USING MARKER INFORMATION
  ## IT IS USED IN THE NE MODULE IN BANAL
  id <- paste( paste("nec",cgiarPIPE::idGenerator(5,5),sep=""), markerDTfile$idOriginal, sep = "_")
  type <- "nec"
  ############################
  # loading the dataset
  if (is.null(markerDTfile)) stop("No input marker data file specified.")
  M0 <- markerDTfile$cleaned 
  
  ############################
  # calculate the relationship matrix
  M <- M0$M
  missingData <- apply(M,2,function(x){length(which(is.na(x)))/length(x)}) # check percentage of missing data
  M <- M[,which(missingData < .50)] # keep markers with at least 50 % of information
  ne <- neMarker(M, maxNe=maxNe, maxMarker=maxMarker, nSamples=nSamples)
  
  metrics <- data.frame(value=ne$allelesCovered,  stdError=ne$allelesCoveredSe,
                   fieldinst="across",  trait="Ne",
                   analysisId=id, method="ratio",
                   traitUnits=ne$Ne, parameter="alleleCoverage",
                   pipeline=NA,
                   stage =NA
  )
  #########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	markerDTfile$id,
    markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = NA
  )
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the files_cleaned folder under such id \n"))
  }
  result <- list(metrics=metrics, predictions=NA, modeling=NA, metadata=metadata,
                 cleaned=NA, outliers=NA, desire=NA, id=id, idOriginal=markerDTfile$idOriginal,
                 metadataFieldinst=NA)
  return(result)
}
