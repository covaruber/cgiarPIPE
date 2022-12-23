gpr <- function(
    markerDTfileTP= NULL, # analysis to be picked from predictions database
    markerDTfilePP= NULL,
    trait= NULL, # per trait
    fieldinst="across",
    verbose=FALSE
){

  baseId <- idGenerator(5,5)
  id <- paste("gpr",baseId,sep="")
  type <- "gpr"
  if(is.null(markerDTfileTP)){stop("Please provide the predictions", call. = FALSE)}
  if(is.null(markerDTfilePP)){stop("Please provide the markers", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}

  ############################
  # loading the dataset
  mydataTP <- markerDTfileTP$predictions # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  mydataPP <- markerDTfilePP$cleaned$M

  utraits <- unique(mydataTP$trait)
  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% utraits){
      if(verbose){ cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))}
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  trait <- setdiff(trait,traitToRemove)

  ############################
  ## genomic prediction analysis
  predictionsBindList <- list()
  for(iTrait in trait){ #  iTrait = trait[1]
    print(iTrait)
    mydataTPprov <- mydataTP[which(mydataTP$trait %in% iTrait & mydataTP$fieldinst %in% fieldinst),]
    common <- unique(intersect(colnames(mydataPP), mydataTPprov[,"geno"])) # common markers
    M <- mydataPP[,common] # PP subset with markers in common
    mydataTPprov <- mydataTPprov[which(mydataTPprov[,"geno"] %in% common),] # TP subset with markers in common
    rownames(mydataTPprov) <- mydataTPprov$geno

    predictedValue <- M %*% as.matrix(mydataTPprov[colnames(M),"predictedValue"])
    stdError <- M %*% as.matrix(mydataTPprov[colnames(M),"stdError"])
    rel <- 1 - (stdError^2)/as.numeric(var(predictedValue, na.rm=TRUE))

    predictionsBindList[[iTrait]] <- data.frame(analysisId=id, pipeline= paste(sort(unique(mydataTPprov$pipeline)),collapse=", "),
                                                trait=iTrait, genoCode=1:length(predictedValue), geno=rownames(predictedValue),
                                                genoType="genomicPrediction", genoYearOrigin=paste(sort(unique(mydataTPprov$genoYearOrigin)),collapse=", "),
                                                genoYearTesting=paste(sort(unique(mydataTPprov$genoYearTesting)),collapse=", "),
                                                fieldinst=fieldinst, predictedValue=predictedValue[,1], stdError=stdError[,1], rel=rel[,1],
                                                stage=paste(sort(unique(mydataTPprov$stage)),collapse=", ")
    )
  }

  predictionsBind <- do.call(rbind, predictionsBindList)
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")

  #########################################
  ## update databases
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	markerDTfileTP$id,
    markerbooks	= NA,  markerDataFile =	markerDTfilePP$id,
    year = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  # saveRDS(db.params, file = file.path(wd,"metadata",paste0(id,".rds")))
  # db.params$analysisId <- idhits
  # db.params$analysisType <- typehits
  # saveRDS(db.params, file = file.path(wd,"metadata",paste0(idhits,".rds")))
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
  # saveRDS(mod, file = file.path(wd,"modeling",paste0(id,".rds")))
  # mod$analysisId <- idhits
  # mod$analysisType <- typehits
  # saveRDS(mod, file = file.path(wd,"modeling",paste0(idhits,".rds")))

  # write pipeline metrics
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  result <- list(metrics=NA, predictions=predictionsBind[,predcols], modeling=mod, metadata=db.params,
                 cleaned=NA, outliers=NA, desire=NA, id=id)
  return(result)#paste("gs/gwas done:",id))
}
