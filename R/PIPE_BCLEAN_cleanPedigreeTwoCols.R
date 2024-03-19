cleanr <- function(
    pedigreeDTfile= NULL, indiv=NULL, dam=NULL, sire=NULL,
    verbose=FALSE, idOriginal=NULL, cleanCharacters=TRUE
){
  ## THIS FUNCTION MATCHES COLUMNS FOR A RAW PEDIGREE FILE WHEN BOTH PARENTAL COLUMNS ARE PROVIDED TO PUT IT IN A STANDARD FORMAT FOR POSTERIOR USE
  ## IT IS USED IN THE MATCH COLUMNS PEDIGREE MODULE IN BCLEAN
  id <- paste( paste("clr",idGenerator(5,5),sep=""), gsub(".csv","",gsub(" ","",idOriginal)), sep = "_")
  type <- "cleaningr"
  ############################
  # loading the dataset
  if (is.null(pedigreeDTfile)) stop("No input pedigree data file specified.")
  if (is.null(indiv)) stop("No column specifying the name of the individual specified.")
  if (is.null(sire)) stop("No column specifying the name of the sire/male specified.")
  if (is.null(dam)) stop("No column specifying the name of the dam/female specified.")
  
  mydata <- as.data.frame(pedigreeDTfile) 
  if(!indiv %in% colnames(mydata) ){
    stop("Please make sure that a column 'geno' with genotype names is included in the file and the rest of the columns is only markers",call. = FALSE)
  }
  if(!sire %in% colnames(mydata) ){
    stop("Please make sure that a column 'geno' with genotype names is included in the file and the rest of the columns is only markers",call. = FALSE)
  }
  if(!dam %in% colnames(mydata) ){
    stop("Please make sure that a column 'geno' with genotype names is included in the file and the rest of the columns is only markers",call. = FALSE)
  }

  ############################
  ## save with standard names and unique records
  orPed <- unique(mydata[,c(indiv,dam,sire)])
  colnames(orPed) <- c("indiv","dam","sire")

  if(cleanCharacters){
    for(iCol in c("indiv","dam","sire")){
      orPed[,iCol] <- cgiarBase::cleanChar(as.character(orPed[,iCol]))
    }
  }
  #########################################
  ## update databases: write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	NA,
    markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = NA
  )

  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  result <- list(metrics=NA, predictions=NA, modeling=NA, metadata=metadata,
                 cleaned=orPed, outliers=NA, desire=NA, id=id, idOriginal=idOriginal,
                 metadataFieldinst=NA)
  return(result)
}
