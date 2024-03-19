cleanm <- function(
    markerDTfile= NULL, geno="geno", useColumns=NULL,
    verbose=FALSE,
    missingData=c("NN","FAIL","FAILED","Uncallable","Unused",""), 
    specialToRemove=c("[:]","[|]"),
    idOriginal=NULL,
    refAlleles=NULL,
    cleanCharacters=TRUE,
    ploidy=2
){
  ## FUNCTION TO MATCH MARKER FILES TO DESIRED COLUMNS (GENO & USECOLUMNS) AND CONVERT TO NUMERIC CODE (-1,0,1) AND DO SOME BASIC CLEANING
  ## IT IS USED IN THE MATCH COLUMNS MARKERS MODULE IN BCLEAN
  id <- paste( paste("clm",cgiarPIPE::idGenerator(5,5),sep=""), gsub(".csv","",gsub(" ","",idOriginal)), sep = "_")
  type <- "cleaningm"
  ############################
  # loading the dataset
  if (is.null(markerDTfile)) stop("No input marker data file specified.")
  mydata <- markerDTfile #
  
  if(!geno %in% colnames(mydata) ){
    stop("Please make sure that a column 'geno' with genotype names is included in the file and the rest of the columns is only markers",call. = FALSE)
  }else{
    colnames(mydata)[which(colnames(mydata) == geno)] <- "geno"
  }
  if(!is.null(useColumns)){
    mydata<- mydata[,useColumns]
  }
  ############################
  # transform to numbers
  genos <- mydata[,"geno"]
  mydata <- mydata[,-which(colnames(mydata) %in% "geno")]
  # replace missing data foraneous call
  if(!is.null(missingData)){
    for(u in 1:length(missingData)){
      miss <- which(mydata == missingData[u], arr.ind = TRUE)
      if(nrow(miss) > 0){
        mydata[miss] = NA
      }
    }
  }
  # replace special characters with empty space
  if(!is.null(specialToRemove)){
    cn <- colnames(mydata)
    rn <- rownames(mydata)
    for(u in 1:length(specialToRemove)){
      mydata <- apply(mydata,2,function(x){gsub(specialToRemove[u],"",x)})
    }
    rownames(mydata) <- rn
  }
  shouldBeTransformed=FALSE # is input letters or numbers
  isLetter <- which(unlist(lapply(mydata, class)) == "character")
  if(length(isLetter) > 1){shouldBeTransformed=TRUE}  # if we found these are nor numbers but letters we transform
  if(shouldBeTransformed){
    
    ## remove markers with too much missing indivisuals
    missingDataPerc <- apply(mydata,2,function(x){length(which(is.na(x)))/length(x)}) # check percentage of missing data
    mydata <- mydata[,which(missingDataPerc < .50)] # keep markers with at least 50 % of information
    ## remove individuals with too many markers missing
    missingDataPerc <- apply(mydata,1,function(x){length(which(is.na(x)))/length(x)}) # check percentage of missing data
    mydata <- mydata[which(missingDataPerc < .50),] # keep individuals with at least 50 % of information
    
    if(median(nchar(mydata[,1]), na.rm=TRUE) <= 1){
      stop("This function only allows double-letter coded markers. Please look at other function to transform it first", call. = FALSE)
    }
    # identify if markers have the expected number of alleles
    genotypesPloidyN <- ploidy + 1 # accepted number of alleles
    genotypesN <- apply(mydata,2,function(x){length(table(x))}) # number of alleles per marker
    goodMarkersWithGenotypesN <- which(genotypesN > 0 & genotypesN <= genotypesPloidyN) # no monomorphic and beyond accepted #  # removed to keep monomorphic
    mydata <- mydata[,goodMarkersWithGenotypesN] # subset
    M0 <- cgiarBase::atcg1234(mydata, ref.alleles = refAlleles, silent = TRUE) # transform
    genos <- genos[which(missingDataPerc < .50)] 
    rownames(M0$M) <- genos
  }else{
    ## remove markers with too much missing indivisuals
    missingDataPerc <- apply(mydata,2,function(x){length(which(is.na(x)))/length(x)}) # check percentage of missing data
    mydata <- mydata[,which(missingDataPerc < .50)] # keep markers with at least 50 % of information
    ## remove individuals with too many markers missing
    missingDataPerc <- apply(mydata,1,function(x){length(which(is.na(x)))/length(x)}) # check percentage of missing data
    mydata <- mydata[which(missingDataPerc < .50),] # keep individuals with at least 50 % of information
    genos <- genos[which(missingDataPerc < .50)]
    #  impute the rest
    mydata <- apply(mydata,2,sommer::imputev) # impute
    ## center
    mydata <- apply(mydata,2,function(x){y <- x - (max(x)+min(x))/2; return(y)}) # make sure is centered (if monomorphic it will be zero)
    M0 <- list(M=mydata, ref.alleles=NA, Moriginal=mydata)
    rownames(M0$M) <- genos
  }
  # ensure no weird names in individuals
  if(cleanCharacters){
    rownames(M0$M) <- cgiarBase::cleanChar(as.character(rownames(M0$M))) # remove special characters
  }
  #########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,  analysisType =	type,
    fieldbooks	= NA,  phenoDataFile =	NA,
    markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,  stage = NA
  )
  #
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the files_cleaned folder under such id \n"))
  }
  result <- list(metrics=NA, predictions=NA, modeling=NA, metadata=metadata,
                 cleaned=M0, outliers=NA, desire=NA, id=id, idOriginal=idOriginal,
                 metadataFieldinst=NA)
  return(result)
}
