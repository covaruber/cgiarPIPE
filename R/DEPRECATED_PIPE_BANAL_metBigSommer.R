metBigSommer <- function(
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
    nPC=100,
    batchSizeToPredict=500,
    tolParInv=1e-4,
    verbose=TRUE
){
  ## THIS FUNCTION PERFORMS A MULTI TRIAL ANALYSIS FOR A BIG NUMBER OF INDIVIDUALS
  ## IS USED UNDER THE GENETIC EVALUATION MODULES
  type <- "met"
  id <- paste( paste(type,idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(is.null(traitFamily)){traitFamily <- rep("quasi(link = 'identity', variance = 'constant')", length(trait))}
  if(length(traitFamily) != length(trait)){stop("Trait distributions should have the same length than traits to be analyzed.", call. = FALSE)}
  names(traitFamily) <- trait
  heritLB <- rep(heritLB,100)
  heritUB <- rep(heritUB,100)
  traitOrig <- trait
  common <- intersect(fixedTerm,randomTerm)
  fixedTerm <- setdiff(fixedTerm,common)
  if("genoF" %in% randomTerm){}else{interactionsWithGeno <- NULL}
  ############################
  # loading the dataset
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  mydata <- phenoDTfile$predictions #readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  pipeline_metrics <- phenoDTfile$metrics #readRDS(file.path(wd,"metrics",paste0(phenoDTfile)))
  if(is.null(genoAmatrix)){
    surrogate <- "TGV"
    analysisTypeAmat <- "nnn"
  }else{
    analysisTypeAmat <- substr(genoAmatrix$id,1,3)
    surrogate <- ifelse(analysisTypeAmat %in% c("clm","grm"), "GEBV", "EBV") # if id has word grm then is a GEBV otherwise is EBV
  }

  utraits <- unique(mydata$trait)
  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% utraits){
      if(verbose){ cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))}
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  trait <- setdiff(trait,traitToRemove)
  traitToRemove <- character()
  if(!is.null(fieldsToInclude)){
    for(k in 1:length(trait)){
      if( length(which(fieldsToInclude[,trait[k]] > 0)) == 0 ){
        if(verbose){ cat(paste0("'", trait[k], "' has not fields selected. It will be removed from trait list \n"))}
        traitToRemove <- c(traitToRemove,trait[k])
        # stop(paste0("'", trait[k], "' is not a column in the given dataset"), call. = FALSE)
      }
    }
    trait <- setdiff(trait,traitToRemove)
  }
  heritLB <- heritLB[which(traitOrig %in% trait)]
  heritUB <- heritUB[which(traitOrig %in% trait)]
  ############################
  # remove outliers from each fieldinst
  fields <- unique(mydata$fieldinst)
  mydata$rowindex <- 1:nrow(mydata)
  # for(iTrait in trait){ # iTrait = trait[1]
  #   for(kField in 1:length(fields)){
  #     subData <- mydata[which((mydata$fieldinst == fields[kField]) & mydata$trait == iTrait ),]
  #     outlier <- boxplot.stats(x=subData$predictedValue,coef=2 )$out
  #     if(length(outlier) > 0){
  #       # print(outlier)
  #       mydata[subData[which(subData$predictedValue %in% outlier),"rowindex"],"predictedValue"] <- NA
  #     }
  #   }
  # }
  ##############################
  ## met analysis
  h2 <- vg <-vg.se <- h2.se <- mu <- numeric(); field <- trt <- vector()
  predictionsList <- list(); listPev <- list(); counter=counter2=1
  traitToRemove <- character()
  library(sommer)
  for(iTrait in trait){ # # iTrait = trait[1]  iTrait="value"
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    # subset data
    mydataSub <- droplevels(mydata[which(mydata$trait == iTrait),])
    if(!is.null(fieldsToInclude)){
      fieldsToIncludeTrait = rownames(fieldsToInclude)[which(fieldsToInclude[,iTrait] > 0)]
      mydataSub <- mydataSub[which(mydataSub$fieldinst %in% fieldsToIncludeTrait),]
    }
    # remove bad fieldinst
    pipeline_metricsSub <- pipeline_metrics[which(pipeline_metrics$trait == iTrait & pipeline_metrics$parameter %in% c("plotH2","H2","meanR2")),]
    goodFields <- unique(pipeline_metricsSub[which((pipeline_metricsSub$value > heritLB[counter2]) & (pipeline_metricsSub$value < heritUB[counter2])),"fieldinst"])
    mydataSub <- mydataSub[which(mydataSub$fieldinst %in% goodFields),]
    # remove bad records
    mydataSub <- mydataSub[which(mydataSub$predictedValue > 0),]

    if(nrow(mydataSub) == 0){
      print(paste("There was no predictions to work with in trait",iTrait,". Please look at your H2 boundaries. You may be discarding all fields."))
      traitToRemove <- c(traitToRemove,iTrait)
    }else{
      mydataSub$genoF <- as.factor(mydataSub$geno)
      mydataSub$fieldinstF <- as.factor(mydataSub$fieldinst)
      mydataSub$pipelineF <- as.factor(mydataSub$pipeline)
      mydataSub$stageF <- as.factor(mydataSub$stage)
      mydataSub$genoYearOriginF <- as.factor(mydataSub$genoYearOrigin)
      mydataSub$genoYearTestingF <- as.factor(mydataSub$genoYearTesting)
      mydataSub <- mydataSub[which(mydataSub$geno != ""),]
      ## add metadata from fieldinst
      if(!is.null(phenoDTfile$metadataFieldinst)){
        metas <- phenoDTfile$metadataFieldinst; #metas <- metas[,c("fieldinst","timePoint","latitude","longitude","altitude")]
        metasClass <- unlist(lapply(metas,class))
        numericMetas <- names(metasClass)[which(metasClass %in% c("integer","numeric"))]
        # center variables
        for(iMeta in numericMetas){
          metas[,iMeta] <- metas[,iMeta] - mean(metas[,iMeta])
        }
        mydataSub <- merge(mydataSub,metas, by="fieldinst", all.x = TRUE)
      }
      ## build the environmental index
      ei <- aggregate(predictedValue~fieldinstF, data=mydataSub,FUN=mean); colnames(ei)[2] <- "envIndex0"
      ei <- ei[with(ei, order(envIndex0)), ]
      ei$envIndex <- ei$envIndex0 - mean(ei$envIndex0)
      #
      colnames(ei) <- cgiarBase::replaceValues(colnames(ei), Search = "envIndex0", Replace = paste0("envIndex_",iTrait))
      phenoDTfile$metadataFieldinst <- merge(phenoDTfile$metadataFieldinst, ei[,c("fieldinstF",paste0("envIndex_",iTrait))], by.x="fieldinst", by.y="fieldinstF", all.x=TRUE)
      ## add the environmental index to the original dataset
      mydataSub <- droplevels(merge(mydataSub,ei[,c("fieldinstF","envIndex")], by="fieldinstF", all.x = TRUE))
      ## define the interactions to be used
      if(!is.null(interactionsWithGeno)){
        interactionsWithGenoTrait <- interactionsWithGeno
        interactionsWithGenoToRemove <- character()
        for(iInter in 1:length(interactionsWithGenoTrait)){
          if(interactionsWithGenoTrait[iInter] %in% colnames(mydataSub)){ # if trait is even there in dataset
            checkInters <- length(unique(mydataSub[,interactionsWithGenoTrait[iInter]]))
            if (checkInters < 2){ # there needs to be at least more than one level
              interactionsWithGenoToRemove <- c(interactionsWithGenoToRemove,interactionsWithGenoTrait[iInter])
            }
          }else{ # remove straight away
            interactionsWithGenoToRemove <- c(interactionsWithGenoToRemove,interactionsWithGenoTrait[iInter])
          }
        }
        interactionsWithGenoTrait <- setdiff(interactionsWithGenoTrait,interactionsWithGenoToRemove)
      }else{
        interactionsWithGenoTrait <- interactionsWithGeno
      }
      ##############
      # do analysis
      if(!is.na(var(mydataSub[,"predictedValue"],na.rm=TRUE))){ # if there's variance
        if( var(mydataSub[,"predictedValue"], na.rm = TRUE) > 0 ){
          # make sure the terms to be fitted have more than one level
          if(deregress){
            mydataSub$predictedValue <- mydataSub$predictedValue/mydataSub$rel
          }
          if(!is.null(randomTerm)){
            rTermsTrait <- randomTerm[which(apply(data.frame(randomTerm),1,function(x){length(table(mydataSub[,x]))}) > 1)]
            rTermsTrait <- setdiff(rTermsTrait,"genoF")
            rTermsTrait <- c("vsc(isc(xx$Z))",rTermsTrait)
          }else{rTermsTrait=NULL}

          if(!is.null(fixedTerm)){
            fixedTermTraitMinus <- setdiff(fixedTerm,"1") # remove intercept for a minute
            if(length(fixedTermTraitMinus) > 0){
              fixedTermTrait <- fixedTermTraitMinus[which(apply(data.frame(fixedTermTraitMinus),1,function(x){length(table(mydataSub[,x]))}) > 1)]
              fixedTermTrait <- c("1",fixedTermTrait)  # add it back
            }else{
              fixedTermTrait <- fixedTerm
            }
          }else{fixedTermTrait=NULL}

          if(!is.null(residualBy)){
            residualByTrait <- residualBy[which(apply(data.frame(residualBy),1,function(x){length(table(mydataSub[,x]))}) > 1)]
            if(length(residualByTrait) == 0){residualByTrait=NULL}
          }else{residualByTrait=NULL}
          #####
          if(length(interactionsWithGenoTrait) > 0){ # If there's interactions to be fitted build the formula terms
            interacs <- expand.grid(interactionsWithGenoTrait,"xx$Z")
            interacs <- as.data.frame(interacs[which(as.character(interacs[,1]) != as.character(interacs[,2])),])
            interacsUnlist <- apply(interacs,1,function(x){
              paste("vsc(",paste("dsc(",x[1],"),isc(",x[2],")"), ")")
            })
            rTermsTrait <- c(rTermsTrait,interacsUnlist)
          }else{
            rTermsTrait <- rTermsTrait
          }
          rTermsTrait <- setdiff(rTermsTrait, fixedTermTrait)

          if(length(rTermsTrait) == 0){
            ranran <- NULL#"~NULL"
          }else{
            ranran <- paste("~",paste(rTermsTrait, collapse=" + ") )
          }
          fix <- paste("predictedValueS ~",paste(fixedTermTrait, collapse=" + "))
          # Ai <- PED$ginv; attr(Ai, "INVERSE") <- TRUE
          if(!is.null(residualByTrait)){
            ranres <- as.formula(paste0("~",residualByTrait,""))
          }else{
            ranres <- NULL#"~units"
          }

          mydataSub=mydataSub[with(mydataSub, order(fieldinstF)), ]
          mydataSub$w <- 1/(mydataSub$stdError^2)
          if(verbose){
            cat(fix,"\n")
            cat(ranran,"\n")
          }

          if( is.null(genoAmatrix) ){
            stop("This form of genetic evaluation requires marker information",call. = FALSE)
          }else{
            genoFlevels <- as.character(unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"genoF"]))
            if(analysisTypeAmat %in% "clm"){ # user provided marker data
              commonBetweenMandP <- intersect(rownames(genoAmatrix$cleaned$M),genoFlevels)
              # print(genoFlevels)
              if(length(commonBetweenMandP) < 2){
                stop("Markers could not be matched with phenotypes. Please ensure that you have used the right marker file or check the rownames of your marker matrix and ids of your phenotypes.", call. = FALSE)
              }
              ## subset to data that can be used
              Md <- genoAmatrix$cleaned$M#[commonBetweenMandP,]
              mydataSub <- mydataSub[which(mydataSub$geno %in% commonBetweenMandP),]
              if(ncol(Md) > nrow(Md)){M <- tcrossprod(Md)}else{M <- Md}
              xx = with(mydataSub, cgiarPIPE::redmm(x=geno, M=M, nPC=nPC, returnLam = TRUE))
            }else{ # user provided a relationship matrix
              stop("This form of genetic evaluation cannot use a GRM. Please provide the marker calls instead. ",call. = FALSE)
            }
          }

          ##########
          ## stage 2
          ## use mmec for sparse equation
          ##########
          posib <- c(1,10,100,1000,10000,100000)
          closestFactor <- floor(mean((mydataSub[,"predictedValue"]))/posib)
          closestFactor <- which(closestFactor > 100)
          if(length(closestFactor) == 0){closestFactor <- 1}
          closestFactor <- closestFactor[max(closestFactor)]
          closestFactor <- posib[closestFactor]
          mydataSub[,"predictedValueS"] <- (mydataSub[,"predictedValue"])/closestFactor#min((mydataSub[,"predictedValue"]))
          mydataSub[,"stdErrorS"] <- (mydataSub[,"stdError"])/closestFactor#min((mydataSub[,"predictedValue"]))
          OM <- Matrix::Diagonal(x=1/(mydataSub[,"stdErrorS"]^2))
          m <- matrix(1/var(mydataSub$predictedValueS, na.rm = TRUE))
          ranFormulation=as.formula(ranran)

          if(verbose){
            print("Using weights in the analysis. Residual variance will be fixed to 1.")
          }
          
          mix <- try(
            mmec(as.formula(fix),
                 W=OM,
                 random = ranFormulation,
                 # emWeight = rep(1,100),
                 rcov=~vsc(isc(units,thetaC = matrix(3), theta = m)),
                 tolParConvLL = 1e-3,
                 tolParConvNorm = -1,
                 data=mydataSub,
                 nIters=maxit),
            silent = TRUE
          )
          # summary(mix)$varcomp
          # lapply(mix$theta, print)


          if(!inherits(mix,"try-error") ){ # if random model runs well try the fixed model
            ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # if model is random
            ######################################################
            ## do genetic evaluation only if there was genotype as random effect
            ######################################################
            Mfull <- genoAmatrix$cleaned$M
            if(ncol(Mfull) > nrow(Mfull)){M2 <- tcrossprod(Mfull)}else{M2 <- Mfull}
            xx2 = with(mydataSub, cgiarPIPE::redmm(x=geno, M=M2, nPC=nPC, returnLam = TRUE))

            predictedValue <- ( ( xx2$Lam %*% mix$uList$`vsc(isc(xx$Z))` ) + mix$b[1,1] ) * closestFactor
            start <- mix$partitions$`vsc(isc(xx$Z))`[1,1]
            end <- mix$partitions$`vsc(isc(xx$Z))`[1,2]
            pev <- xx2$Lam %*% mix$Ci[start:end,start:end] %*% t(xx2$Lam)  * closestFactor
            # print(dim(pev))
            stdError <- ( sqrt(Matrix::diag(pev)) )
            # print("all good 1")
            ######################################################
            ## end of genetic evaluation function
            ######################################################

            genoF <- rownames(predictedValue)
            pp <- data.frame(genoF=genoF,predictedValue=predictedValue[,1],stdError=stdError)
            ss = summary(mix)$varcomp #; rownames(ss) <- ss$VarComp
            Vg <- ss["xx:Z:isc:isc",1] * closestFactor; Vr <- ss["units:m:",1] * closestFactor
            pp$rel <- (Vg - Matrix::diag(pev))/Vg
            # print("all good 2")
            badRels <- which(pp$rel > 1); if(length(badRels) > 0){pp$rel[badRels] <- 0.9999}

            pp$trait <- iTrait # add trait
            ## heritabilities
            h2[counter] <- mean(pp$rel); h2.se[counter] <- 1e-6
            ## genetic variances
            vg[counter] <- Vg; vg.se[counter] <- 1e-6
            mu[counter] <- mean(pp$predictedValue, na.rm=TRUE)
            field[counter] <- "across"; trt[counter] <- iTrait

            # lpv <- sum(mix$EDdf$Model[1:which(mix$EDdf$Term == "genoF")])+1 # to be used as a starting point if random regression is requested
            # extract sensitivities if interaction is requested
            if(length(interactionsWithGenoTrait) > 0){ # if there's interactions
              if( length(intersect(interactionsWithGenoTrait, c("envIndex","timePoint","latitude","longitude","altitude","weather1","weather2"))) > 0 ){

                for(iInteractionTrait in c("envIndex","timePoint","latitude","longitude","altitude","weather1","weather2")){ # iInteractionTrait <- "envIndex"

                  if( (iInteractionTrait %in% interactionsWithGenoTrait) ){ # iInteractionTrait="envIndex"

                    counter <- counter+1

                    use <- paste0("vsc(",paste0("dsc(",iInteractionTrait,"), isc(","xx$Z",")"), ")")
                    predictedValue <- ( xx2$Lam %*% mix$uList[[use]] ) * closestFactor#+ mix$b[1,1]
                    start <- mix$partitions[[use]][1,1]
                    end <- mix$partitions[[use]][1,2]
                    pev <- ( xx2$Lam %*% mix$Ci[start:end,start:end] %*% t(xx2$Lam) ) * closestFactor
                    stdError <- sqrt(Matrix::diag(pev))

                    ######################################################
                    ## end of genetic evaluation function
                    ######################################################

                    genoF <- rownames(predictedValue)
                    pp2 <- data.frame(genoF=genoF,predictedValue=predictedValue[,1],stdError=stdError)
                    ss = summary(mix)$varcomp #; rownames(ss) <- ss$VarComp

                    vgname <- paste0(iInteractionTrait,":xx:Z:",iInteractionTrait,":",iInteractionTrait)
                    Vg <- ss[vgname,1]*closestFactor; Vr <- ss["units:m:",1]*closestFactor
                    pp2$rel <- (Vg - Matrix::diag(pev))/Vg
                    badRels <- which(pp2$rel > 1); if(length(badRels) > 0){pp2$rel[badRels] <- 0.9999}

                    pp2$trait <- paste(iTrait,iInteractionTrait,sep="-") # add trait
                    ## heritabilities
                    h2[counter] <- mean(pp2$rel); h2.se[counter] <- 1e-6
                    ## genetic variances
                    vg[counter] <- Vg; vg.se[counter] <- 1e-6
                    mu[counter] <- mean(pp2$predictedValue, na.rm=TRUE)
                    field[counter] <- "across"; trt[counter] <- paste(iTrait,iInteractionTrait,sep="-")
                    pp <- rbind(pp,pp2)

                  }

                }

              }
            } # if there's even interactions


            ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          }else{
            if(verbose){ cat(paste("Aggregating and assuming h2 = 0 \n"))}
            pp <- aggregate(predictedValue ~ genoF, FUN=mean, data=mydataSub)
            pp$stdError <- sd(pp$predictedValue)  # 1
            # pp$status <- "Aggregated"
            pp$rel <- 1e-6
            pp$trait <- iTrait
            ## heritabilities # h2Pred <- abs(round(mean(1 - (pp$stdError^2 / sum(summary(mix)$varcomp[1:2,1]))), 3))
            h2[counter] <- 1e-6; h2.se[counter] <- 1e-6
            ## genetic variances
            vg[counter] <- 1e-6; vg.se[counter] <- 1e-6
            mu[counter] <- mean(pp$predictedValue, na.rm=TRUE)
            field[counter] <- "across"; trt[counter] <- iTrait
          }

          pp$fieldinstF <- "across"
          mydataForGenoType <- droplevels(mydata[which(mydata$trait == iTrait),])
          pp$entryType <- apply(data.frame(pp$genoF),1,function(x){
            found <- which(mydataForGenoType$geno %in% x)
            if(length(found) > 0){
              x2 <- paste(sort(unique(toupper(trimws(mydataForGenoType[found,"genoType"])))), collapse = "#");
            }else{x2 <- "unknown"}
            return(x2)
          })
          mydataForGenoType <- NULL

          '%!in%' <- function(x,y)!('%in%'(x,y))
          pp[which(pp$genoF %in% mydataSub$geno),"entryType"] <- paste(pp[which(pp$genoF %in% mydataSub$geno),"entryType"],"tested",sep="_")
          pp[which(pp$genoF %!in% mydataSub$geno),"entryType"] <- paste(pp[which(pp$genoF %!in% mydataSub$geno),"entryType"],"predicted",sep="_")

          predictionsList[[counter2]] <- pp;
          counter=counter+1
        }
      }
    }
    counter2 = counter2+1
  }
  #
  trait <- setdiff(trait,traitToRemove)
  heritLB <- heritLB[which(traitOrig %in% trait)]
  heritUB <- heritUB[which(traitOrig %in% trait)]
  if(length(predictionsList) == 0){stop("There was no predictions to work with. Please look at your H2 boundaries. You may be discarding all fields.",call. = FALSE)}
  predictionsBind <- do.call(rbind, predictionsList)
  predictionsBind$analysisId <- id
  ###
  colnames(predictionsBind) <- cgiarBase::replaceValues(Source=colnames(predictionsBind), Search=c("genoF","fieldinstF","entryType"), Replace=c("geno","fieldinst","genoType"))
  predictionsBind$pipeline <- paste(sort(unique(mydata$pipeline)),collapse=", ")
  ##########################################
  ## add timePoint of origin, stage and geno code
  entries <- unique(mydata[,"geno"])
  vals <- apply(data.frame(entries),1,function(x){
    out <- (sort(mydata[which(mydata$geno %in% x),"genoYearOrigin"], decreasing = FALSE))[1]
    return(out)
  })
  vals2 <- apply(data.frame(entries),1,function(x){
    out <- paste(sort(unique(mydata[which(mydata$geno %in% x),"stage"], decreasing = FALSE)), collapse = ".")
    return(out)
  })
  vals3 <- apply(data.frame(entries),1,function(x){
    out <- (sort(mydata[which(mydata$geno %in% x),"genoCode"], decreasing = FALSE))[1]
    return(out)
  })
  vals4 <- apply(data.frame(entries),1,function(x){
    out <- (sort(mydata[which(mydata$geno %in% x),"genoYearTesting"], decreasing = FALSE))[1]
    return(out)
  })
  vals5 <- apply(data.frame(entries),1,function(x){
    out <- (sort(mydata[which(mydata$geno %in% x),"mother"], decreasing = FALSE))[1]
    return(out)
  })
  vals6 <- apply(data.frame(entries),1,function(x){
    out <- (sort(mydata[which(mydata$geno %in% x),"father"], decreasing = FALSE))[1]
    return(out)
  })
  baseOrigin <- data.frame(geno=entries,mother=vals5, father=vals6,genoYearOrigin=vals, genoYearTesting=vals4, stage=vals2, genoCode=vals3)
  predictionsBind <- merge(predictionsBind,baseOrigin, by="geno", all.x=TRUE)
  ## in case singleCross predictions are available
  if(!is.null(genoAmatrix$singleCrosses)){
    rownames(genoAmatrix$singleCrosses) <- genoAmatrix$singleCrosses$hybrid
    toGetParents <- which(predictionsBind$geno %in% rownames(genoAmatrix$singleCrosses))
    if(length(toGetParents) > 0){
      predictionsBind$mother[toGetParents] <- as.character(genoAmatrix$singleCrosses[predictionsBind$geno[toGetParents], "Var1"])
      predictionsBind$father[toGetParents] <- as.character(genoAmatrix$singleCrosses[predictionsBind$geno[toGetParents], "Var2"])
    }
  }
  # should we set to NA predictions lower than zero?
  # zeros <- intersect( which(predictionsBind$predictedValue < 0), setdiff(1:nrow(predictionsBind),grep("adaptability", predictionsBind$trait)))
  # if(length(zeros) > 0){predictionsBind[zeros,] <- NA}
  wide <- reshape(predictionsBind[,c("geno","trait","predictedValue")], direction = "wide", idvar = "geno",
                  timevar = "trait", v.names = "predictedValue", sep= "_")
  wide <- as.matrix(wide[,-1]); colnames(wide) <- unique(predictionsBind$trait)
  wide <- apply(wide,2,sommer::imputev)
  if(scaledDesire){
    if(verbose){cat(paste("scaledDesire has been set to",scaledDesire,"'desirev' values are expected to be the desired change in std. deviations \n"))}
    wide <- apply(wide,2,scale)
  }else{
    if(verbose){cat(paste("scaledDesire has been set to",scaledDesire,"'desirev' values are expected to be desired change in original units \n"))}
  }
  G <- cov(wide, use="pairwise.complete.obs")
  #########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	phenoDTfile$id,
    markerbooks	= NA,  markerDataFile =	ifelse(is.null(genoAmatrix),NA, genoAmatrix$id),
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  ## write the values used for cleaning to the modeling database
  modeling <- data.frame(
    trait = trait,
    traitLb = NA,
    traitUb = NA,
    outlierCoef = NA,
    analysisId = rep(id,length(trait)),
    analysisType =rep(type,length(trait)) ,
    fixedModel = rep(fix,length(trait)),
    randomModel = rep(ifelse(is.null(ranran),"~",ranran),length(trait)),
    residualModel = ifelse(is.null(ranres),"~units",rep(as.character(ranres),length(trait))),
    h2Threshold = paste(heritLB,heritUB, sep=",")
  )

  # write predictions
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","mother","father","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")

  # write pipeline metrics
  metrics <- data.frame(value=c(h2,vg,mu),  stdError=c(h2.se,vg.se,rep(1e-6,length(mu))),
                   fieldinst=c(field,field,field),  trait=c(trt,trt,trt),
                   analysisId=id, method=c(rep("(G-PEV)/G",length(h2)),rep("REML",length(vg)), rep("mean",length(mu)) ),
                   traitUnits=NA, parameter=c(rep("r2",length(h2)),rep("Vg",length(vg)), rep("mean",length(mu)) ),
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", "),
                   stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
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
  result <- list(metrics=metrics, predictions=predictionsBind[,predcols], modeling=modeling, metadata=metadata,
                 cleaned=NA, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst, PEV=listPev, genoMetaData=NA
  )
  return(result)#paste("met done:",id))
}
