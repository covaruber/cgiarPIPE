grm <- function(
    markerDTfile= NULL,#wd=NULL, 
    verbose=FALSE
){

  # if(is.null(wd)){wd <- getwd()}
  # md <- strsplit(wd,"/")[[1]]; md <- md[length(md)]
  # if(md != "DB"){stop("Please set your working directory to the DB folder", call. = FALSE)}

  id <- paste("grm",idGenerator(5,5),sep="")
  type <- "grm"

  ############################
  # loading the dataset
  if (is.null(markerDTfile)) stop("No input marker data file specified.")
  # ava.files <- dir(file.path(wd,"files_cleaned"))
  # if(!paste0(markerDTfile) %in% ava.files){stop("markerDTfile is not present in the files_cleaned folder. Please check the name of your file and its location",call. = FALSE)}
  M0 <- markerDTfile$cleaned #readRDS(file.path(wd,"files_cleaned",paste0(markerDTfile)))

  ############################
  # calculate the relationship matrix
  A <- sommer::A.mat(M0$M)

  #########################################
  ## update databases
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	NA,
    markerbooks	= NA,  markerDataFile =	NA,
    year = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = NA
  )
  # saveRDS(db.params, file = file.path(wd,"metadata",paste0(id,".rds")))
  # write the values used for cleaning to the modeling database
  # write predictions
  # write pipeline metrics

  # saveRDS(A, file = file.path(wd,"files_cleaned",paste0(id,".rds")))
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the files_cleaned folder under such id \n"))
  }
  result <- list(metrics=NA, predictions=NA, modeling=NA, metadata=db.params,
                 cleaned=A, outliers=NA, desire=NA, id=id)
  return(result)#paste("grm done:",id))
}
