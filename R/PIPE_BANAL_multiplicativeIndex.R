multIndex <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    trait= NULL,
    thresholds=NULL,
    verbose=TRUE
){
  ## THIS FUNCTION CALCULATED A MULTIPLICATIVE INDEX
  ## IS USED IN THE BANAL APP UNDER THE GENETIC EVALUATION MODULES
  id <- paste( paste("mid",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "mid"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(length(grep("met",phenoDTfile$id)) == 0){stop("Index can only be calculated on results from a MET analysis using across environment predictions",call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(length(trait) != length(thresholds)){stop("The number of traits and threshold values needs to be equal",call. = FALSE)}

  ############################
  # loading the dataset
  mydata <- phenoDTfile$predictions # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  mydata <- mydata[which(mydata$trait %in% trait),]
  names(thresholds) <- trait
  ############################
  ## index calculation

  for(iTrait in trait){ # iTrait <- trait[1]
    mydata[which(mydata$trait == iTrait),"predictedValue"] = apply(data.frame(mydata[which(mydata$trait == iTrait),"predictedValue"]),1,function(x){
      y <- ifelse(x > thresholds[iTrait], 1, 0)
      return(y)
    })
    mydata[which(mydata$trait == iTrait),"stdError"] = 0
  }

  baseToFill <-  mydata[which(mydata$trait == iTrait),]
  baseToFill$predictedValue <- 0
  genos <- unique(mydata$geno)
  for(iGeno in genos){ # iGeno = genos[1]
    baseToFill[which(baseToFill$geno == iGeno),"predictedValue"] <- sum(mydata[which(mydata$geno == iGeno),"predictedValue"])
  }
  predictionsBind <- rbind(mydata,baseToFill)

  # write predictions
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","mother","father","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")

  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	phenoDTfile$id,
    markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(mydata$stage)),collapse=", ")
  )
  ## write the values used for cleaning to the modeling database
  modeling <- data.frame(
    trait = trait,
    traitLb = thresholds,
    traitUb = thresholds,
    outlierCoef = NA,
    analysisId = id,
    analysisType = type,
    fixedModel = NA,
    randomModel = NA,
    residualModel = NA,
    h2Threshold = NA
  )

  # write pipeline metrics
  metrics <- data.frame(value=thresholds,
                   stdError=1e-6,
                   fieldinst="across",
                   trait=trait,
                   analysisId=id, method= "ifelse",
                   traitUnits=NA,
                   stage = paste(sort(unique(mydata$stage)),collapse=", "),
                   parameter= "threshold",
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", ")
  )

  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  if(is.null(phenoDTfile$metadataFieldinst)){
    metadataFieldinst=NA
  }else{
    metadataFieldinst=phenoDTfile$metadataFieldinst
  }
  result <- list(metrics=metrics, predictions=predictionsBind[,predcols], modeling=modeling, metadata=metadata,
                 cleaned=phenoDTfile$cleaned, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst
  )
  return(result)
}
