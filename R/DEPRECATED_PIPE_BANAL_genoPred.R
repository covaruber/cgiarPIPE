gpr <- function(
    markerEffectsDTfileTP= NULL, # analysis to be picked from predictions database
    markerDTfilePP= NULL,
    trait= NULL, # per trait
    fieldinst="across",
    verbose=FALSE
){
  ## THIS FUNCTION MAKES GENOMIC PREDICTION USING MARKER EFFECTS
  ## THIS IS NOT CURRENTLY IN USE IN ANY MODULE
  trait <- gsub("-","", trait)

  id <- paste( paste("gpr",idGenerator(5,5),sep=""), markerEffectsDTfileTP$idOriginal, sep = "_")
  type <- "gpr"
  if(is.null(markerEffectsDTfileTP)){stop("Please provide the marker effects from the training population", call. = FALSE)}
  if(is.null(markerDTfilePP)){stop("Please provide the markers from individuals to be predicted", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  # if new version make sure the ref alleles are the same
  if(!is.null(markerEffectsDTfileTP$ref.alleles)){
   co <- intersect(colnames(markerEffectsDTfileTP$ref.alleles), colnames(markerDTfilePP$ref.alleles))
   myTest <- markerEffectsDTfileTP$ref.alleles["Ref",co] ==  markerDTfilePP$ref.alleles["Ref",co]
   myTestNumber <- length(which(myTest == FALSE)) # number of different reference alleles
   if(myTestNumber > 0){
     stop("The reference alleles used in the TP and PP files are different. Analyses cannot proceed. Please match once
          again the marker file to be predicted using the reference alleles from the marker file used for the
          training population.", call=FALSE)
   }
  }else{
    cat("We couldn't check if the reference alleles used in the TP and PP files are the same. Please proceed carefully or match once
          again the marker file to be predicted using the reference alleles from the marker file used for the
          training population. This control is recent.")
  }
  ############################
  # loading the dataset
  mydataTP <- markerEffectsDTfileTP$predictions # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  mydataPP <- markerDTfilePP$cleaned$M

  mydataTP$trait <- gsub("-","", mydataTP$trait) # trait names shouldn't have strange characters
  utraits <- unique(mydataTP$trait)
  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% utraits){
      if(verbose){ cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))}
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  trait <- setdiff(trait,traitToRemove)

  if(is.null(fieldinst)){fieldinst <- unique(mydataTP$fieldinst)}
  missingData <- apply(mydataPP,2,function(x){length(which(is.na(x)))/length(x)}) # check percentage of missing data
  mydataPP <- mydataPP[,which(missingData < .50)] # keep markers with at least 50 % of information
  ############################
  ## genomic prediction analysis
  predictionsBindList <- list()
  for(iTrait in trait){ #  iTrait = trait[1]
    # print(iTrait)
    mydataTPprov <- mydataTP[which(mydataTP$trait %in% iTrait & mydataTP$fieldinst %in% fieldinst),]
    common <- unique(intersect(colnames(mydataPP), mydataTPprov[,"geno"])) # common markers
    if(length(common) == 0){
      print(paste("There was no intersection of markers in TP and markers in PP in trait",iTrait,". Skipping the trait. Please look at your input files. You may be picking the wrong file(s)."))
    }else{
      # markers in PP
      M <- mydataPP[,common] # PP subset with markers in common
      # intercept in TP
      # intercept <- mydataTPprov$predictedValue[which(mydataTPprov$genoType == "intercept")]
      # seIntercept <-  mydataTPprov$stdError[which(mydataTPprov$genoType == "intercept")]
      # marker effects from TP
      mydataTPprov <- mydataTPprov[which(mydataTPprov[,"geno"] %in% common),] # TP subset with markers in common
      rownames(mydataTPprov) <- mydataTPprov$geno


      predictedValue <- (M %*% as.matrix(mydataTPprov[colnames(M),"predictedValue"])) #+ intercept
      # stdError <- abs(M %*% as.matrix(mydataTPprov[colnames(M),"stdError"])) + seIntercept
      # standard error
      if(is.null(markerEffectsDTfileTP$PEV)){stop("PEV for marker effects not found. Please run the marker effects pipeline again. Some changes happened in past versions")}
      PEVm <- markerEffectsDTfileTP$PEV[[iTrait]]$genoF #$genoF#markerEffectsDTfileTP$PEV[[iTrait]]
      PEVm <- PEVm[common,common]
      PEVg <- M%*%PEVm%*%t(M)
      stdError <- sqrt(abs(diag(PEVg)))
      rel <- 1 - (diag(PEVg))/as.numeric(var(predictedValue, na.rm=TRUE))

      labelsGeno <- rep("predicted_TGV", length(rel)) # markerEffectsDTfileTP$genoMetaData$withPhenoOnly # tested_TGV
      tested_GEBV <- which(rownames(predictedValue) %in% markerEffectsDTfileTP$genoMetaData$withMarkandPheno)
      if(length(tested_GEBV) > 0){labelsGeno[tested_GEBV] <- "tested_GEBV"}
      predicted_GEBV <- which(rownames(predictedValue) %in% markerEffectsDTfileTP$genoMetaData$withMarkOnly)
      if(length(predicted_GEBV) > 0){labelsGeno[predicted_GEBV] <- "predicted_GEBV"}

      predictionsBindList[[iTrait]] <- data.frame(analysisId=id, pipeline= paste(sort(unique(mydataTPprov$pipeline)),collapse=", "),
                                                  trait=iTrait, genoCode=1:length(predictedValue), geno=rownames(predictedValue),
                                                  mother="unknown",father="unknown",
                                                  genoType=labelsGeno, genoYearOrigin=paste(sort(unique(mydataTPprov$genoYearOrigin)),collapse=", "),
                                                  genoYearTesting=paste(sort(unique(mydataTPprov$genoYearTesting)),collapse=", "),
                                                  fieldinst=paste(sort(unique(mydataTPprov$fieldinst)),collapse=", "),
                                                  predictedValue=predictedValue[,1], stdError=stdError, rel=rel,
                                                  stage=paste(sort(unique(mydataTPprov$stage)),collapse=", ")
      )
    }
  }
  if(length(predictionsBindList) == 0){stop("There was no intersection of markers in TP and markers in PP in any trait.",call. = FALSE)}
  predictionsBind <- do.call(rbind, predictionsBindList)
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","mother","father","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")

  #########################################
  ## update databases
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	markerEffectsDTfileTP$id,
    markerbooks	= NA,  markerDataFile =	markerDTfilePP$id,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  ## write the values used for cleaning to the modeling database
  mod <- data.frame(
    trait = trait,
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
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  if(is.null(markerEffectsDTfileTP$metadataFieldinst)){
    metadataFieldinst=NA
  }else{
    metadataFieldinst=markerEffectsDTfileTP$metadataFieldinst
  }
  result <- list(metrics=NA, predictions=predictionsBind[,predcols], modeling=mod, metadata=db.params,
                 cleaned=NA, outliers=NA, desire=NA, id=id, idOriginal=markerEffectsDTfileTP$idOriginal,
                 metadataFieldinst=metadataFieldinst
  )
  return(result)#paste("gs/gwas done:",id))
}
