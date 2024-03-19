corLMM <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    fixedTerm= c("fieldinstF"),
    randomTerm=c("genoF"),
    residualBy=NULL,
    interactionsWithGeno=NULL,
    trait= NULL, # per trait
    traitFamily=NULL,
    fieldsToInclude=NULL,
    useWeights=TRUE,
    heritLB= 0.15,
    heritUB= 0.95,
    scaledDesire=TRUE,
    genoAmatrix=NULL,
    deregress=FALSE,
    maxit=50,
    nPC=0,
    markerEffects=FALSE,
    predictFullSet=FALSE,
    batchSizeToPredict=500,
    tolParInv=1e-4,
    verbose=TRUE
){
  ## THIS FUNCTION CALCULATES BLUPS OR GBLUPS FOR EACH FIELDINST LEVEL AND THEN THE CORRELATION BETWEEN THOSE (G)EBV ESTIMATES
  ## THIS FUNCTION I USED IN THE GENETIC EVALUATION MODULES FROM BANAL
  # phenoDTfile <<- phenoDTfile
  # genoAmatrix <<- genoAmatrix
  if(is.null(phenoDTfile)){stop("Please provide the phenotype file", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(is.null(traitFamily)){traitFamily <- rep("quasi(link = 'identity', variance = 'constant')", length(trait))}
  if(length(traitFamily) != length(trait)){stop("Trait distributions should have the same length than traits to be analyzed.", call. = FALSE)}
  
  traitFamily2 <<- traitFamily
  names(traitFamily2) <- trait
  id <- paste( paste("cor",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "cor"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  
  analysisTypeAmat <- substr(genoAmatrix$id,1,3)
  if(analysisTypeAmat %in% c("grm","nrm")){ # user provided a relationship matrix not markers
    if(predictFullSet){
    }else{
      uniqueGenos <- unique(phenoDTfile$predictions$geno)
      commonGenos <- intersect(rownames(genoAmatrix$cleaned), uniqueGenos)
      genoAmatrix$cleaned <- genoAmatrix$cleaned[commonGenos,commonGenos]
    }
  }
  
  listPhenoDTfile0 <- split(phenoDTfile$predictions, phenoDTfile$predictions$fieldinst)
  gList <- pList <- list() # to save genetic correlations for each trait
  for(jTrait in trait){ # jTrait <- trait[1]
    print(jTrait)
    resultField <- list() # to save predictions per field per trait
    listPhenoDTfile <- lapply(listPhenoDTfile0,function(x){x[which(x$trait == jTrait),]})
    for(iFieldinst in names(listPhenoDTfile)){ # iFieldinst <- names(listPhenoDTfile)[1]
      print(iFieldinst)
      phenoDTfileProv <- phenoDTfile
      phenoDTfileProv$predictions <- listPhenoDTfile[[iFieldinst]]
      provisionalResult <- try(cgiarPIPE::metLMM(
        phenoDTfile= phenoDTfileProv, fixedTerm=fixedTerm, randomTerm=randomTerm,
        residualBy=residualBy,interactionsWithGeno=interactionsWithGeno,
        trait= jTrait, 
        traitFamily=traitFamily2[jTrait],
        fieldsToInclude=fieldsToInclude,
        useWeights=useWeights,heritLB= heritLB,heritUB= heritUB,scaledDesire=scaledDesire,
        genoAmatrix=genoAmatrix,deregress=deregress,maxit=maxit,nPC=nPC,
        markerEffects=markerEffects,batchSizeToPredict=batchSizeToPredict,
        tolParInv=tolParInv, verbose=verbose
      ), silent = TRUE
      )
      if(!inherits(provisionalResult, "try-error")){
        provisionalResult$predictions$fieldinst <- iFieldinst
        resultField[[iFieldinst]] <- provisionalResult
      }
    }
    ## extract predictions
    pList[[jTrait]] <- do.call(rbind,lapply(resultField, function(x){unique(x$predictions)}) )
    ## to column bind predictions for each field
    allGenos <- unique(unlist(lapply(resultField, function(x){unique(x$predictions$geno)})))
    preds <- list()
    for(iFieldinst in names(resultField)){
      mm <- data.frame(matrix(NA,nrow=length(allGenos))); colnames(mm) <- iFieldinst
      rownames(mm) <- allGenos
      mm[resultField[[iFieldinst]]$predictions$geno,iFieldinst] <- resultField[[iFieldinst]]$predictions$predictedValue
      preds[[iFieldinst]] <- mm
    }
    predictionsWide <- do.call(cbind,preds) 
    
    Gfields <- cov(predictionsWide, use="pairwise.complete.obs")
    Cfields <- cor(predictionsWide, use="pairwise.complete.obs")
    extremes <- which(Cfields > 1, arr.ind = TRUE)
    if(nrow(extremes) > 0){ Cfields[extremes]=0.99}
    
    Cfields[lower.tri(Cfields)] <- NA
    Cdf <- as.data.frame(as.table(Cfields)) # converts a matrix in a data frame
    Cdf <- Cdf[which(!is.na(Cdf$Freq)),]
    ##
    Gfields[lower.tri(Gfields)] <- NA
    Gdf <- as.data.frame(as.table(Gfields)) # converts a matrix in a data frame
    Gdf <- Gdf[which(!is.na(Gdf$Freq)),]
    # write pipeline metrics
    # genetic covariance
    pm <- data.frame(value=Gdf$Freq,  stdError=NA,
                     fieldinst=paste(Gdf$Var1, Gdf$Var2, sep="###"),  trait=jTrait,
                     analysisId=id, method="CS",
                     traitUnits=NA, parameter="covG",
                     pipeline=paste(sort(unique(phenoDTfile$predictions$pipeline)),collapse=", "),
                     stage = paste(sort(unique(phenoDTfile$predictions$stage)),collapse=", ")
    )
    # genetic correlation
    pm2 <- data.frame(value=Cdf$Freq,  stdError=NA,
                      fieldinst=paste(Cdf$Var1, Cdf$Var2, sep="###"),  trait=jTrait,
                      analysisId=id, method="CS",
                      traitUnits=NA, parameter="corG",
                      pipeline=paste(sort(unique(phenoDTfile$predictions$pipeline)),collapse=", "),
                      stage = paste(sort(unique(phenoDTfile$predictions$stage)),collapse=", ")
    )
    gList[[jTrait]] <- rbind(pm,pm2);
  } # end of trait list
  
  predictions <- do.call(rbind,pList)
  pm2 <- do.call(rbind, gList)
  oldMetrics <- phenoDTfile$metrics
  oldMetrics <- oldMetrics[which(oldMetrics$trait %in% trait),]
  metrics <- rbind(pm2, oldMetrics)
  
  #########################################
  ## update databases: write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,  analysisType =	type,  fieldbooks	= NA,
    phenoDataFile =	phenoDTfile$id, markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(phenoDTfile$predictions$stage)),collapse=", ")
  )
  ## write the values used for cleaning to the modeling database
  modeling <- data.frame(
    trait = trait,
    traitLb = NA,
    traitUb = NA,
    outlierCoef = NA,
    analysisId = rep(id,length(trait)),
    analysisType =rep(type,length(trait)) ,
    fixedModel = rep("trait~1",length(trait)),
    randomModel = rep("geno",length(trait)), 
    residualModel = rep("units",length(trait)),
    h2Threshold = paste(heritLB,heritUB, sep=",")
  )
  ##
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  if(is.null(phenoDTfile$metadataFieldinst)){
    metadataFieldinst=NA
  }else{
    metadataFieldinst=phenoDTfile$metadataFieldinst
  }
  result <- list(metrics=metrics, predictions=predictions, modeling=modeling, metadata=metadata,
                 cleaned=phenoDTfile$cleaned, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst
  )
  return(result)
}
