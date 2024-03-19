ocs <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    relDTfile= NULL,
    trait= NULL, # per trait
    fieldinst="across",
    nCross=20,
    targetAngle=30, # in radians
    verbose=FALSE,
    maxRun=100,
    genoType=NULL,
    wd=NULL
){
  ## THIS FUNCTION CALCULATES THE OPTIMAL CROSS SELECTION USING A TRAIT AND A RELATIONSHIP MATRIX
  ## IS USED IN THE BANAL APP UNDER THE GENETIC EVALUATION MODULES
  id <- paste( paste("ocs",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "ocs"
  if(is.null(phenoDTfile)){stop("Please provide the predictions", call. = FALSE)}
  if(is.null(relDTfile)){stop("Please provide the markers", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(length(trait) > 1){
    stop(paste0(" Only one trait can be used for optimal contribution. We suggest using an index."), call. = FALSE)
  }
  if(length(fieldinst) > 1){
    stop(paste0(" Only one fieldinst can be used for optimal contribution. We suggest using an across environment value."), call. = FALSE)
  }

  otherTraits <- setdiff(unique(phenoDTfile$predictions$trait), trait)
  ############################
  # loading the dataset
  mydata <- phenoDTfile$predictions # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  if(!is.null(genoType)){
    mydata <- mydata[which(mydata$genoType %in% genoType),]
  }

  if(ncol(relDTfile$cleaned) != nrow(relDTfile$cleaned)){ # A matrix as data frame
    Adf <- relDTfile$cleaned
    A <- matrix(NA, nrow=max(Adf$Var1n), ncol = max(Adf$Var2n))
    A[as.matrix(Adf[,c("Var1n","Var2n")])] = Adf[,"Freq"]
    A[lower.tri(A)] <- t(A)[lower.tri(A)] # fill the lower triangular
    rownames(A) <- colnames(A) <- levels(Adf$Var1)
  }else{ # A matrix as an actual matrix
    A <- as.matrix(relDTfile$cleaned)
    A[lower.tri(A)] <- t(A)[lower.tri(A)] # fill the lower triangular
  }

  myrel <- A

  utraits <- unique(mydata$trait)
  if (!trait %in% utraits){
    stop(paste0("'", trait, "' is not present in the given dataset or the genoType doesn't correspond."), call. = FALSE)
  }
  mydata <- mydata[which(mydata$trait %in% trait),]
  mydata <- mydata[which(mydata$fieldinst %in% fieldinst),] # make sure only across
  if(nrow(mydata) == 0){stop("Please check the trait and fieldinst selected since there's no phenotypic data for that combination",call. = "FALSE")}
  # make sure you have same phenotypes and genotypes

  common <- intersect(rownames(myrel), mydata[,"geno"])
  if(length(common) == 0){
    stop("There was no intersection between the IDs in the relationship matrix and the IDs in the phenotypes provided. Please check your input files.",call. = FALSE)
  }
  myrel <- myrel[common,common]
  mydata <- mydata[which(mydata[,"geno"] %in% common),]

  ############################
  ## ocs analysis

  forLoop <- expand.grid(nCross, targetAngle)
  predictionsBindList <- list()
  meanCross <- meanFcross <- meanCrossSe <- meanFcrossSe <- numeric()

  for(iRow in 1:nrow(forLoop)){ # iRow=1

    ebv <- data.frame(mydata[,c("predictedValue")]); rownames(ebv) <- mydata[,"geno"]
    ebv <- data.frame(ebv[rownames(myrel),]); rownames(ebv) <- rownames(myrel)
    crossComb = t(combn(1:nrow(myrel), 2)) # all possible cross combintations
    eMP = (ebv[crossComb[,1],] +  ebv[crossComb[,2],])/2  # expected EBVs of all crosses based on # mean parent EBVs
    K <- as.matrix(myrel)
    # OCS: Determine a crossing plan
    plan = cgiarOcs::selectCrosses(nCross=forLoop[iRow,1], # number of crossed to be identified using OCS
                                   targetAngle=((forLoop[iRow,2])*pi)/180, # 30 degrees in radians
                                   u=eMP, # expected cross mean EBVs
                                   maxRun=maxRun,
                                   G=K)   # GRM
    dim(plan$crossPlan)

    crossPlan <- as.data.frame(plan$crossPlan) # list of crosses to be made already sorted by best
    crossPlan[ ,1] <- rownames(K)[crossPlan[ ,1]]
    crossPlan[ ,2] <- rownames(K)[crossPlan[ ,2]]
    colnames(crossPlan) <- c("Parent1", "Parent2", "OCS.merit")
    eMPsel = (ebv[crossPlan[ ,1],] +     # expected EBVs of selected crosses based on
                ebv[crossPlan[ ,2],])/2  # mean parent EBVs
    inbreeding = diag(K)
    inbreedingSel = (inbreeding[crossPlan[ ,1]] + inbreeding[crossPlan[ ,2]])/2

    predictionsBindList[[iRow]] <- data.frame(analysisId=id, pipeline= paste(sort(unique(mydata$pipeline)),collapse=", "),
                                  trait=trait, genoCode=1:nrow(crossPlan), geno=paste(crossPlan[,1],crossPlan[,2], sep=" x "),
                                  mother=crossPlan[,1],father=crossPlan[,2],
                                  genoType="predictedCross", genoYearOrigin=forLoop[iRow,1], genoYearTesting=forLoop[iRow,2],
                                  fieldinst=fieldinst, predictedValue=eMPsel, stdError=inbreedingSel, rel=crossPlan[,3], stage="unknown"
    )

    meanCross[iRow] <- mean(eMPsel)
    meanFcross[iRow] <- mean(inbreedingSel)
    meanCrossSe[iRow] <- sd(eMPsel)/sqrt(length(eMPsel))
    meanFcrossSe[iRow] <- sd(inbreedingSel)/sqrt(length(inbreedingSel))
    
    if(length(otherTraits) > 0){ # if there's more traits in the file, add the value of the crosses for those traits
      traitPredictions <- list()
      for(iTrait in otherTraits){ # iTrait <- otherTraits[1]
        
        provPredictions <- phenoDTfile$predictions
        provPredictions <- provPredictions[which(provPredictions$trait == iTrait), ]
        ebv2 <- data.frame(provPredictions[,c("predictedValue")]); rownames(ebv2) <- provPredictions[,"geno"]
        eMPtrait = (ebv2[crossPlan[ ,1],] +  ebv2[crossPlan[ ,2],])/2  #
        
        traitPredictions[[iTrait]] <- data.frame(  analysisId=id, pipeline= paste(sort(unique(mydata$pipeline)),collapse=", "),
                   trait=iTrait, genoCode=1:nrow(crossPlan), geno=paste(crossPlan[,1],crossPlan[,2], sep=" x "),
                   mother=crossPlan[,1],father=crossPlan[,2],
                   genoType="predictedCross", genoYearOrigin=forLoop[iRow,1], genoYearTesting=forLoop[iRow,2],
                   fieldinst=fieldinst, predictedValue=eMPtrait, stdError=inbreedingSel, rel=crossPlan[,3], stage="unknown"
        )
        
      }
      predictionsBindList[[iRow]] <- rbind(predictionsBindList[[iRow]], do.call(rbind, traitPredictions))
    }
  }

  predictionsBind <- do.call(rbind, predictionsBindList)

  #########################################
  ## write metrics
  metrics <- data.frame(value=c(meanCross,meanFcross),
                   stdError=c(meanCrossSe,meanFcrossSe),
                   fieldinst=fieldinst, trait=trait,
                   analysisId=id, method= "sum/n",
                   traitUnits=NA, stage = paste("nCross",forLoop[,1],"targetAngle",forLoop[,2], sep="_"),
                   parameter= c(rep("meanValue",length(meanCross)), rep("meanF",length(meanFcross)) ),
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", ")
  )

  #########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	phenoDTfile$id,
    markerbooks	= NA,  markerDataFile =	relDTfile$id,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  # saveRDS(db.params, file = file.path(wd,"metadata",paste0(id,".rds")))
  ## write the values used for cleaning to the modeling database
  modeling <- data.frame(
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

  # write predictions
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","mother","father","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")
  # write pipeline metrics
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the predictions database under such id \n"))
  }
  if(is.null(phenoDTfile$metadataFieldinst)){
    metadataFieldinst=NA
  }else{
    metadataFieldinst=phenoDTfile$metadataFieldinst
  }
  
  result <- list(metrics=metrics, predictions=predictionsBind[,predcols], modeling=modeling, metadata=metadata,
                 cleaned=phenoDTfile$cleaned, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst)
  return(result)
}
