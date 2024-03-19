transRaw <- function(
    phenoDTfile= NULL,
    trait=NULL, # per trait
    transformation = NULL,
    verbose=FALSE
){
  ## THIS FUNCTION TRANFORMS TRAITS USING MATHEMATICAL FUNCTIONS TO CREATE A NEW TRAIT. IT STARTS FROM A CSV FILE
  ## IS USED IN THE BCLEAN APP UNDER THE TRANSFORMATION MODULES
  # transformation <- c(transformation,rep("I",100)) # in case the user forgets to provide all the transformations
  if(length(trait) != length(transformation)){stop("The vector of transformations required need to be as many as the number of traits and viceversa.", call. = FALSE)}
  cbrt <- function(x){return(x^(1/3))}
  mpr <- function(x){ return(cgiarBase::replaceValues(x, c("[-]","[+]","?","[Het]","[Het]-[+/Bold]","[-]-Bold","[-]-Medium"),c(-1,1,NA,0,0,-1,-1)))}
  ctn <- function(x){return(as.numeric(as.character(x)))} # transform from character to number
  traitOrig <- trait
  if(is.null(phenoDTfile)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}

  ###################################
  # loading the dataset
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  mydata <- phenoDTfile

  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% colnames(mydata)){
      if(verbose){
        cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))
      }
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  trait <- setdiff(trait,traitToRemove)
  transformation <- transformation[which(traitOrig %in% trait)]
  #####################################
  # transformation
  counter <- 1
  for(iTrait in trait){ # iTrait=trait[1]
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    mydata[,paste(iTrait,transformation[counter], sep = "_")] <- do.call(transformation[counter],list(mydata[,iTrait]))
  }
  ##########################################

  return(mydata)
}
