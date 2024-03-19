grm <- function(
    markerDTfile= NULL,#wd=NULL,
    verbose=FALSE
){
  ## THIS FUNCTION CALCULATES A GENOMIC RELATIONSHIP MATRIX USING GENETIC MARKERS
  ## IS USED IN THE BANAL APP UNDER THE STRUCTURE MODULES
  id <- paste( paste("grm",idGenerator(5,5),sep=""), markerDTfile$idOriginal, sep = "_")
  type <- "grm"

  ############################
  # loading the dataset
  if (is.null(markerDTfile)) stop("No input marker data file specified.")
  M0 <- markerDTfile$cleaned #readRDS(file.path(wd,"files_cleaned",paste0(markerDTfile)))

  ############################
  # calculate the relationship matrix
  M <- M0$M
  missingData <- apply(M,2,function(x){length(which(is.na(x)))/length(x)}) # check percentage of missing data
  M <- M[,which(missingData < .50)] # keep markers with at least 50 % of information
  M <- apply(M,2,function(x){x-(mean(c(min(x),max(x))))}) # center the matrix at 0
  A <- sommer::A.mat(M)

  A[lower.tri(A)] <- NA
  #########################################
  ## update databases: write the parameters to the parameter database
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
  result <- list(metrics=NA, predictions=NA, modeling=NA, metadata=metadata,
                 cleaned=A, outliers=NA, desire=NA, id=id, idOriginal=markerDTfile$idOriginal,
                 metadataFieldinst=NA)
  return(result)
}
