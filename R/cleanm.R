cleanm <- function(
    markerDTfile= NULL, geno="geno", useColumns=NULL,
    verbose=FALSE,
    missingData=c("NN","")
    # wd=NULL
){
  # if(is.null(wd)){wd <- getwd()}
  # md <- strsplit(wd,"/")[[1]]; md <- md[length(md)]
  # if(md != "DB"){stop("Please set your working directory to the DB folder", call. = FALSE)}

  id <- paste("clm",idGenerator(5,5),sep="")
  type <- "cleaningm"


  ############################
  # loading the dataset
  if (is.null(markerDTfile)) stop("No input marker data file specified.")
  # modeling <- read.csv("modeling.csv")
  # ava.files <- dir(file.path(wd,"files_raw"))
  # if(!markerDTfile %in% ava.files){stop("markerDTfile is not present in the files_raw folder. Please check the name of your file and its location",call. = FALSE)}
  mydata <- markerDTfile # read.csv(file.path(wd,"files_raw",markerDTfile))

  if(!geno %in% colnames(mydata) ){
    stop("Please make sure that a column 'geno' with genotype names is included in the file and the rest of the columns is only markers",call. = FALSE)
    }else{
    colnames(mydata)[which(colnames(mydata) == geno)] <- "geno"
    }
  # print(useColumns)
  if(!is.null(useColumns)){
    mydata<- mydata[,useColumns]
  }
  # predictions <- read.csv("predictions.csv")
  # pipeline_metrics <- read.csv("pipeline_metrics.csv")


  ############################
  ## grm calculation

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

  M0 <- sommer::atcg1234(mydata)
  rownames(M0$M) <- genos
  # calculate the relationship matrix
  # A <- A.mat(M)

  #########################################
  ## update databases
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	NA,
    markerbooks	= NA,  markerDataFile =	markerDTfile$id,
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

  # saveRDS(M0, file = file.path(wd,"files_cleaned",paste0(id,".rds")))
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the files_cleaned folder under such id \n"))
  }
  result <- list(metrics=NA, predictions=NA, modeling=NA, metadata=db.params,
                 cleaned=M0, outliers=NA, desire=NA, id=id)
  return(result)#
  # return("cleanm done")
}
