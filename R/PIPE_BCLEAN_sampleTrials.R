sampTrials <- function(
    phenoDTfile= NULL,
    trait=NULL, # per trait
    samplingVariable="genoYearTesting",
    maxNumSamples=200,
    verbose=TRUE
){
  ## THIS FUNCTION TAKES AN OUTPUT FILE FROM A SINGLE TRIAL ANALYSIS AND SAMPLES A MAXIMUM NUMBER OF GENOTYPES TO REDUCE THE SIZE OF GENOTYPES IN AN ANALYSIS (E.G., REALIZED GAIN)
  ## IS USED IN THE BANAL APP UNDER THE TRANSFORMATIONS MODULES
  id <- paste( paste("sta",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, "SAMP", sep = "_")
  type <- "sta"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  
  ############################
  # loading the dataset
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  parameters <- phenoDTfile$metadata # readRDS(file.path(wd,"metadata",paste0(phenoDTfile)))
  mydata <- phenoDTfile$predictions #readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  if(length(which(colnames(mydata) == samplingVariable)) == 0){stop("sampling variable not found in the dataset. Please double check. ")}
  ############################
  ## 
  traitPreds <- list()
  for(iTrait in trait){ # iTrait=trait[1]
    if(verbose){
      cat(paste("Analyzing trait", iTrait,"\n"))
    }
    # subset data
    mydataSub <- droplevels(mydata[which(mydata$trait == iTrait),])
    sampLevels <- sort(unique(mydataSub[,samplingVariable]))
    # table
    tt <- table(mydataSub[,"geno"],mydataSub[,samplingVariable])
    tt <- tt/tt
    tt[which(is.nan(tt), arr.ind = TRUE)] <- 0
    ttsum <- apply(tt,1,sum)
    tt <- cbind(tt,ttsum)
    
    levelPreds <- list()
    for(iLevel in sampLevels){ # iLevel = sampLevels[1]
      ttI <- tt[which(tt[,as.character(iLevel)] > 0),]
      ttI <- ttI[ order(-ttI[,"ttsum"]), ]
      indsToKeep <- rownames(ttI[1:min(maxNumSamples,nrow(ttI)),])
      levelPreds[[iLevel]] <- mydataSub[which(mydataSub[,"geno"] %in% indsToKeep),]
    }
    traitPreds[[iTrait]] <- do.call(rbind,levelPreds)
    
  }
  phenoDTfile$predictions <- do.call(rbind,traitPreds)
  phenoDTfile$id <- id
 
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  return(phenoDTfile)
}
