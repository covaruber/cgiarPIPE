gdd <- function(phenoDTfile= NULL, verbose=FALSE){
  ## THIS FUNCTION CALCULATES GROWING DEGREE DAYS FOR A WEATHER DATA FILE EXTRACTED IN BCLEAN
  ## THIS IS USED UNDER THE WEATHER MODULES ANALYSIS IN BANAL 
  id <- paste( paste("gdd",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "gdd"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(length(grep("wea",phenoDTfile$id)) == 0){stop("This function can only be used in weather files",call. = FALSE)}
  ############################
  # loading the dataset
  if(is.null(phenoDTfile$cleaned)){stop("There's no data for this file. Likely belongs to an older version of the cgiarPIPE package.", call. = FALSE)}
  mydata <- phenoDTfile$cleaned # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  ############################
  ## gdd calculation
  fields <- na.omit(unique(mydata$fieldinst))
  toAdd <- list()
  for(k in 1:length(fields)){ # k=1
    ## subset to dates
    mydataSubset <- mydata[which(mydata$fieldinst == fields[k]),]
    # growing degree days
    xmax <- ifelse(mydataSubset$T2M_MAX > 30, 30, mydataSubset$T2M_MAX)
    xmin <- ifelse(mydataSubset$T2M_MIN < 10, 10, mydataSubset$T2M_MIN)
    gdd <- ((xmax + xmin)/2) - 10
    mydataSubset$gdd <- gdd * 30 # average degree days times 30 days
    mydataSubset$cumgdd <- cumsum(mydataSubset$gdd)
    toAdd[[k]] <- mydataSubset
  }
  wdataDf <- do.call(rbind, toAdd)
  
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  phenoDTfile$cleaned <- wdataDf
  return(phenoDTfile)
}
