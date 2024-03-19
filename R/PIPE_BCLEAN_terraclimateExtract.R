terraExtract <- function(phenoDTfile= NULL, verbose=FALSE){
  ## THIS FUNCTION EXTRATS WEATHER DATA FOR SOME TRIALS IN AN STA FILE THAT HAVE COORDINATE INFORMATION
  ## IS USED IN THE BCLEAN APP UNDER THE DATA EXTRACTION MODULES
  id <- paste( paste("wea",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "wea"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(length(grep("clp",phenoDTfile$id)) == 0){stop("This function can only be used in matched phenotypic files",call. = FALSE)}

  ############################
  # loading the dataset
  if(is.null(phenoDTfile$metadataFieldinst)){stop("There's no metadata for this file. Likely belongs to an older version of the cgiarPIPE package. Please match the columns of your data again.", call. = FALSE)}
  mydata <- phenoDTfile$metadataFieldinst # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  ############################
  ## index calculation
  badLats <- which(mydata$latitude == 1000)
  badLons <- which(mydata$longitude == 1000)
  if((length(badLons) > 0) | (length(badLats) > 0) ){
    stop("Some coordinates are not valid. Please correct and try again", call. = FALSE)
  }

  if( ( max(mydata$plantingDate) > "2020-12-31" ) | ( max(mydata$harvestingDate) > "2020-12-31") ){
    stop("No weather information was extracted. Dates provided are probably later than 2020-12-31. Terraclimate dataset can only extract data from 1958 to the end of 2020. Please try an earlier date. ", call. = FALSE)
  }else{
    wdata <- QBMS::get_terraclimate(lat=mydata$latitude, lon=mydata$longitude,
                           from=min(mydata$plantingDate), to=max(mydata$harvestingDate)
    )
  }
  wdata2 <- wdata
  toAdd <- list()
  for(k in 1:nrow(mydata)){
    ## subset to dates
    dateKp <- strsplit(as.character(mydata$plantingDate[k]), split="-")[[1]]
    yearKp <- as.numeric(dateKp[1])
    monthKp <- as.numeric(dateKp[2])
    dateKh <- strsplit(as.character(mydata$harvestingDate[k]), split="-")[[1]]
    yearKh <- as.numeric(dateKh[1])
    monthKh <- as.numeric(dateKh[2])
    # dates that match the needed dates
    wdata2$climate[[k]] <- wdata$climate[[k]][which( (wdata$climate[[k]]$year >= yearKp) & (wdata$climate[[k]]$month >= monthKp) & (wdata$climate[[k]]$year <= yearKh) & (wdata$climate[[k]]$month <= monthKh) ),]
    # growing degree days
    xmax <- ifelse(wdata2$climate[[k]]$tmax > 30, 30,wdata2$climate[[k]]$tmax)
    xmin <- ifelse(wdata2$climate[[k]]$tmin < 10, 10,wdata2$climate[[k]]$tmin)
    gdd <- ((xmax + xmin)/2) - 10
    wdata2$climate[[k]]$gdd <- gdd * 30 # average degree days times 30 days
    wdata2$climate[[k]]$cumgdd <- cumsum(wdata2$climate[[k]]$gdd)
    # summary for fieldinst
    res1 <- apply(wdata2$climate[[k]], 2, mean, na.rm=TRUE )
    res1 <- res1[-c(1:2)]
    res2 <- data.frame(matrix(res1, nrow=1))
    colnames(res2) <- names(res1)
    # res2$fieldinst <- mydata$fieldinst[k]
    toAdd[[k]] <- res2
    wdata2$climate[[k]]$fieldinst <- mydata$fieldinst[k]
    wdata2$climate[[k]]$latitude <- mydata$latitude[k]
    wdata2$climate[[k]]$longitude <- mydata$longitude[k]
  }
  wdataDf <- do.call(rbind, wdata2$climate)
  metadataFieldinst <- cbind(mydata, do.call(rbind,toAdd ))
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	phenoDTfile$id,
    markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = NA
  )
  ## write the values used for cleaning to the modeling database
  mod <- data.frame(
    trait = NA,
    traitLb = NA,
    traitUb = NA,
    outlierCoef = NA,
    analysisId = id,
    analysisType = type,
    fixedModel = NA,
    randomModel = NA,
    residualModel = NA,
    h2Threshold = NA
  )

  # write pipeline metrics
  pm <- data.frame(value=NA,
                   stdError=NA,
                   fieldinst="none",
                   trait="all",
                   analysisId=id, method= "terra",
                   traitUnits=NA,
                   stage = NA,
                   parameter= "many",
                   pipeline=NA
  )

  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  result <- list(metrics=pm, predictions=NA, modeling=mod, metadata=db.params,
                 cleaned=wdataDf, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst
  )
  return(result)

}
