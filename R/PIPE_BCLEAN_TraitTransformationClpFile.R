trans <- function(
    phenoDTfile= NULL,
    trait=NULL, # per trait
    transformation = NULL,
    verbose=FALSE
){
  ## THIS FUNCTION TRANFORMS TRAITS USING MATHEMATICAL FUNCTIONS TO CREATE A NEW TRAIT. IT STARTS FROM A CLP FILE 
  ## IS USED IN THE BCLEAN APP UNDER THE TRANSFORMATION MODULES
  # transformation <- c(transformation,rep("I",100)) # in case the user forgets to provide all the transformations
  if(length(trait) != length(transformation)){stop("The vector of transformations required need to be as many as the number of traits and viceversa.", call. = FALSE)}
  cbrt <- function(x){return(x^(1/3))}
  mpr <- function(x){ return(cgiarBase::replaceValues(x, c("[-]","[+]","?","[Het]","[Het]-[+/Bold]","[-]-Bold","[-]-Medium"),c(-1,1,NA,0,0,-1,-1)))}

  traitOrig <- trait
  id <- paste("clpt",idGenerator(5,5),sep="")
  type <- "clpt"
  if(is.null(phenoDTfile)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}

  ###################################
  # loading the dataset
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  mydata <- phenoDTfile$cleaned #readRDS(file.path(wd,"files_cleaned",paste0(phenoDTfile)))
  outliers <- phenoDTfile$outliers

  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% colnames(mydata)){
      if(verbose){
        cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))
      }
      traitToRemove <- c(traitToRemove,trait[k])
      # stop(paste0("'", trait[k], "' is not a column in the given dataset"))
    }
  }
  trait <- setdiff(trait,traitToRemove)
  transformation <- transformation[which(traitOrig %in% trait)]
  if(!is.null(outliers)){
    outliers <- outliers[which(outliers$traitName %in% trait),]
    if(nrow(outliers) > 0){outliers$analysisId <- id}
  }
  #####################################
  # transformation
  counter <- 1
  for(iTrait in trait){ # iTrait=trait[1]
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    mydata[,paste(iTrait,transformation[counter], sep = "_")] <- do.call(transformation[counter],list(mydata[,iTrait]))
    identOut <- which(outliers$traitName %in% iTrait)
    if(length(identOut) > 0) {outliers[identOut,"traitName"] <- paste(iTrait,transformation[counter], sep = "_")}
    counter <- counter + 1
  }
  phenoDTfile$cleaned <- mydata
  if(!is.null(outliers)){
    if(nrow(outliers) > 0){
      phenoDTfile$outliers <- rbind(phenoDTfile$outliers, outliers)
    }
  }
  ##########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,analysisType =	type,fieldbooks	= NA,
    phenoDataFile =	phenoDTfile$id, markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,rowcoord =	NA,  colcoord = NA,
    stage =NA
  )
  phenoDTfile$metadata <- metadata
  ## write the values used for cleaning to the modeling database
  modeling <- data.frame(
    trait = trait, traitLb = NA,traitUb = NA,outlierCoef = NA,
    analysisId = id,analysisType = type,fixedModel = transformation,
    randomModel = NA,residualModel = NA,h2Threshold = NA
  )
  phenoDTfile$modeling <- modeling
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  phenoDTfile$id <- id
  return(phenoDTfile)
}
