gwas <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    markerDTfile= NULL,
    trait= NULL, # per trait
    fieldinst=NULL,
    verbose=FALSE
){
  ## THIS FUNCTION PERFORMS A GWAS ANALYSIS 
  ## IS USED IN THE BANAL APP UNDER THE MARKER MODULES
  baseId <- idGenerator(5,5)
  idmes <- paste( paste("mes",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")

  idhits <- paste("gwa",baseId,sep="")
  typemes <- "mes"
  typehits <- "gwa"
  if(is.null(phenoDTfile)){stop("Please provide the predictions", call. = FALSE)}
  if(is.null(markerDTfile)){stop("Please provide the markers", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}

  ############################
  # loading the dataset
  mydata <- phenoDTfile$predictions # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  myrel <- markerDTfile$cleaned # readRDS(file.path(wd,"files_cleaned",paste0(markerDTfile)))

  utraits <- unique(mydata$trait)
  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% utraits){
      if(verbose){ cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))}
      traitToRemove <- c(traitToRemove,trait[k])
      # stop(paste0("'", trait[k], "' is not a column in the given dataset"), call. = FALSE)
    }
  }
  trait <- setdiff(trait,traitToRemove)

  mydata <- mydata[which(mydata$trait %in% trait),]
  if(is.null(fieldinst)){
    fieldinst <- unique(mydata$fieldinst)
  }
  mydata <- mydata[which(mydata$fieldinst %in% fieldinst),] # make sure only across
  if(nrow(mydata) == 0){stop("Please check the trait and fieldinst selected since there's no phenotypic data for that combination",call. = "FALSE")}
  # make sure you have same phenotypes and genotypes

  common <- intersect(rownames(myrel$M), unique(mydata[,"geno"]))
  if(length(common) == 0){
    stop("There was no intersection of markers and phenotypes provided. Please check your input files.",call. = FALSE)
  }
  M <- myrel$M[common,]
  mydata <- mydata[which(mydata[,"geno"] %in% common),]
  # for files not properly cleaned by cleanm (TO BE REMOVED IN A MONTH WHEN THOSE FILES ARE GONE)
  missingData <- apply(M,2,function(x){length(which(is.na(x)))/length(x)}) # check percentage of missing data
  M <- M[,which(missingData < .50)] # keep markers with at least 50 % of information
  ############################
  ## ocs analysis
  predictionsBindMesList <- predictionsBindInterList <- predictionsBindHitsList <- listPev <- list()
  for(iTrait in trait){ #  iTrait = trait[1]
    # print(iTrait)
    mydataSub <- mydata[which(mydata$trait == iTrait),]
    n <- nrow(mydataSub) # to be used for degrees of freedom
    k <- 1 # to be used for degrees of freedom (number of levels in fixed effects)
    Mprov <- M[unique(mydataSub$geno),]
    MMT <- tcrossprod(Mprov) ## MM' = additive relationship matrix
    MMTinv<-solve(MMT + diag(1e-6, ncol(MMT), ncol(MMT))) ## inverse of MM'
    # MTMMTinv<-t(Mprov)%*%MMTinv # M' %*% (M'M)-
    mydataSub$w <- 1/(mydataSub$stdError^2)
    if(length(fieldinst) > 1){myFixedFormula <- as.formula(paste("predictedValue~fieldinst"))}else{myFixedFormula <- as.formula(paste("predictedValue~1"))}
    mc <- sommer::fixm(1)
    mi <- matrix(1/var(mydataSub$predictedValue, na.rm = TRUE))
    mixGBLUP <- try(sommer::mmer(myFixedFormula,
                                 random=~sommer::vsr(geno, Gu=MMT),
                                 rcov=~sommer::vsr(units, Gti=mi, Gtc=mc), # only if we were using STA
                                 nIters=30,
                                 weights = w,
                                 verbose = FALSE,
                                 data=mydataSub), silent = TRUE
    )
    if(!inherits(mixGBLUP,"try-error") ){ # if there's no error proceed
      myorder <- names(mixGBLUP$U$`u:geno`$predictedValue)
      Mprov2 <- Mprov[myorder,]
      MMT2 <- MMT[myorder,myorder]
      MMT2inv <- MMTinv[myorder,myorder]
      MTMMT2inv <-t(Mprov2)%*%MMT2inv
      a.from.g <-MTMMT2inv%*%matrix(mixGBLUP$U$`u:geno`$predictedValue,ncol=1)
      var.g <- kronecker(MMT2,mixGBLUP$sigma$`u:geno`) - mixGBLUP$PevU$`u:geno`$predictedValue
      # var.a.from.g <- t(Mprov2)%*%MMT2inv%*% (var.g) %*% t(MMT2inv)%*%Mprov2
      pev.a.from.g <- t(Mprov2)%*%MMT2inv%*% (mixGBLUP$PevU$`u:geno`$predictedValue) %*% t(MMT2inv)%*%Mprov2
      se.a.from.g <- sqrt(diag(pev.a.from.g) + min(diag(pev.a.from.g))+ 1e-7)
      t.stat.from.g <- a.from.g/se.a.from.g # t-statistic
      pvalGBLUP <- dt(t.stat.from.g,df=n-k-1) #
      minusLog10pvalGBLUP <- -log10(pvalGBLUP)

      rel <- 1 - (se.a.from.g^2)/as.numeric(mixGBLUP$sigma$`u:geno`)
      listPev[[iTrait]] <- pev.a.from.g
      # plot(minusLog10pvalGBLUP)


      predictionsBindInterList[[iTrait]] <- data.frame(analysisId=idmes, pipeline= paste(sort(unique(mydataSub$pipeline)),collapse=", "),
                                                     trait=iTrait, genoCode=1, geno="all", mother="unknown", father="unknown",
                                                     genoType="intercept", genoYearOrigin=paste(sort(unique(mydataSub$genoYearOrigin)),collapse=", "),
                                                     genoYearTesting=paste(sort(unique(mydataSub$genoYearTesting)),collapse=", "),
                                                     fieldinst="across",#paste(sort(unique(mydataSub$fieldinst)),collapse=", "),
                                                     predictedValue=mixGBLUP$Beta$Estimate[1], stdError=sqrt(diag(mixGBLUP$VarBeta))[1], rel=1,
                                                     stage=paste(sort(unique(mydataSub$stage)),collapse=", ")
      )
      predictionsBindMesList[[iTrait]] <- data.frame(analysisId=idmes, pipeline= paste(sort(unique(mydataSub$pipeline)),collapse=", "),
                                                     trait=iTrait, genoCode=1:nrow(a.from.g), geno=rownames(a.from.g), mother="unknown", father="unknown",
                                                     genoType="markerEffect", genoYearOrigin=paste(sort(unique(mydataSub$genoYearOrigin)),collapse=", "),
                                                     genoYearTesting=paste(sort(unique(mydataSub$genoYearTesting)),collapse=", "),
                                                     fieldinst="across",#paste(sort(unique(mydataSub$fieldinst)),collapse=", "),
                                                     predictedValue=a.from.g[,1], stdError=se.a.from.g, rel=rel, stage=paste(sort(unique(mydataSub$stage)),collapse=", ")
      )
      predictionsBindHitsList[[iTrait]] <- data.frame(analysisId=idhits, pipeline= paste(sort(unique(mydataSub$pipeline)),collapse=", "),
                                                      trait=iTrait, genoCode=1:nrow(pvalGBLUP), geno=rownames(a.from.g),mother="unknown", father="unknown",
                                                      genoType="gwasHits", genoYearOrigin=paste(sort(unique(mydataSub$genoYearOrigin)),collapse=", "),
                                                      genoYearTesting=paste(sort(unique(mydataSub$genoYearTesting)),collapse=", "),
                                                      fieldinst="across",#paste(sort(unique(mydataSub$fieldinst)),collapse=", "),
                                                      predictedValue=minusLog10pvalGBLUP[,1], stdError=se.a.from.g, rel=rel, stage=paste(sort(unique(mydataSub$stage)),collapse=", ")
      )
    } # end of if there was no error

  }

  if(length(predictionsBindMesList) == 0){stop("All models failed for your traits. Please look at your H2 values and check that marker data is correct for the predictions used.",call. = FALSE)}

  predictionsBindMes <- do.call(rbind, predictionsBindMesList)
  predictionsBindInt <- do.call(rbind, predictionsBindInterList)
  predictionsBindHits <- do.call(rbind, predictionsBindHitsList)
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","mother","father","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")
  finalPreds <- rbind(predictionsBindMes[,predcols],predictionsBindInt[,predcols])
  #########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= idmes,
    analysisType =	typemes,
    fieldbooks	= NA,
    phenoDataFile =	phenoDTfile$id,
    markerbooks	= NA,  markerDataFile =	markerDTfile$id,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(predictionsBindMes$stage)),collapse=", ")
  )
  ## write the values used for cleaning to the modeling database
  modeling <- data.frame(
    trait = trait,
    traitLb = NA,
    traitUb = NA,
    outlierCoef = NA,
    analysisId = idmes,
    analysisType = typemes,
    fixedModel = NA,
    randomModel = NA,
    residualModel = NA,
    h2Threshold = NA
  )
  # write pipeline metrics
  if(verbose){
    cat(paste("Your analysis id is:",idmes,"\n"))
    cat(paste("Your analysis id is:",idhits,"\n"))
    cat(paste("Your results will be available in the predictions database under such id \n"))
  }
  if(is.null(phenoDTfile$metadataFieldinst)){
    metadataFieldinst=NA
  }else{
    metadataFieldinst=phenoDTfile$metadataFieldinst
  }
  result <- list(metrics=NA, predictions=finalPreds, modeling=modeling, metadata=metadata,
                 cleaned=NA, outliers=NA, desire=NA, id=idmes, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst,
                 PEV=listPev
                 )
  return(result)#paste("gs/gwas done:",idmes))
}
