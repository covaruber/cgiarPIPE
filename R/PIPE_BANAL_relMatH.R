hrm <- function(
    nrmDTfile = NULL,grmDTfile = NULL, tau=1,
    omega=1, tolparinv=1e-6
){
  ## THIS FUNCTION CALCULATES THE H MATRIX OR SINGLE STEP
  ## IS USED IN THE BANAL APP UNDER THE STRUCTURE MODULES
  id <- paste( paste("hrm",idGenerator(5,5),sep=""), nrmDTfile$idOriginal, sep = "_")
  type <- "hrm"
  
  ############################
  # loading the dataset
  A <- nrmDTfile$cleaned
  G <- grmDTfile$cleaned
  H <- H.mat(A,G, tau=tau,
             omega=omega, tolparinv=tolparinv)
  
  #########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	pedigreeDTfile$id,
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
                 cleaned=H, outliers=NA, desire=NA, id=id, idOriginal=nrmDTfile$idOriginal,
                 metadataFieldinst=NA)
  return(result)
}
