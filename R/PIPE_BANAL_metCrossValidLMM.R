metCvLMM <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    fixedTerm= c("fieldinstF"),
    randomTerm=c("genoF"),
    residualBy=NULL,
    interactionsWithGeno=NULL,
    trait= NULL, # per trait
    fieldsToInclude=NULL,
    heritLB= 0.15,
    heritUB= 0.95,
    scaledDesire=TRUE,
    genoAmatrix=NULL,
    deregress=FALSE,
    maxit=50,
    fixErrorVar=TRUE,
    markerEffects=FALSE,
    batchSizeToPredict=500,
    tolParInv=1e-4,
    verbose=TRUE,
    strata="none", kfold=2, nCvSamples=10
){
  ## THIS FUNCTION PERFORMS CROSS VALIDATION FOR A MET MODEL TO DETERMINE THE ACCURACY OF A MODEL
  ## IS USED IN THE BANAL APP UNDER THE METRICS MODULES
  if(markerEffects){type <- "cvg"}else{type <- "cvg"}
  id <- paste( paste(type,idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  heritLB <- rep(heritLB,100)
  heritUB <- rep(heritUB,100)
  traitOrig <- trait
  common <- intersect(fixedTerm,randomTerm)
  fixedTerm <- setdiff(fixedTerm,common)
  if("genoF" %in% randomTerm){}else{
    stop("genotype has to be random for cross validations.", call. = FALSE)
    interactionsWithGeno <- NULL
  }
  if(is.null(genoAmatrix)){ stop("relationship matrix or markers have to be provided for cross validations.", call. = FALSE)}
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

  ###########################################
  ## define cross validation
  baseMeta <- phenoDTfile$metadataFieldinst
  for(iTrait in trait){ # # iTrait = trait[2]  iTrait="value"
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    # subset data
    mydataSub0 <- droplevels(mydata[which(mydata$trait == iTrait),])
    if(!is.null(fieldsToInclude)){
      fieldsToIncludeTrait = rownames(fieldsToInclude)[which(fieldsToInclude[,iTrait] > 0)]
      mydataSub0 <- mydataSub0[which(mydataSub0$fieldinst %in% fieldsToIncludeTrait),]
    }
    # remove bad fieldinst
    pipeline_metricsSub <- pipeline_metrics[which(pipeline_metrics$trait == iTrait & pipeline_metrics$parameter %in% c("plotH2","H2")),]
    goodFields <- pipeline_metricsSub[which((pipeline_metricsSub$value > heritLB[counter2]) & (pipeline_metricsSub$value < heritUB[counter2])),"fieldinst"]
    mydataSub0 <- mydataSub0[which(mydataSub0$fieldinst %in% goodFields),]
    mydataSub0$genoF <- as.factor(mydataSub0$geno)
    mydataSub0$fieldinstF <- as.factor(mydataSub0$fieldinst)
    mydataSub0$pipelineF <- as.factor(mydataSub0$pipeline)
    mydataSub0$stageF <- as.factor(mydataSub0$stage)
    mydataSub0$genoYearOriginF <- as.factor(mydataSub0$genoYearOrigin)
    mydataSub0$genoYearTestingF <- as.factor(mydataSub0$genoYearTesting)

    cvList <- cvListR <- list()
    if(strata != "none"){nCvSamples=1}
    for(iSample in 1:nCvSamples){ ## iSample =1 for loop for the ith cv sample
      if(strata == "none"){
        cvsdf <- paste("strata",sample(1:kfold, nrow(mydataSub0), replace = TRUE), sep="_")
      }else{cvsdf <- mydataSub0[,strata]}
      strataLevels <- unique(cvsdf)
      TPVP <- TPVPr <- list()
      for(iStrata in 1:length(strataLevels)){ # iStrata=1  ## for each strata level

        TPVP[[iStrata]] <- list(VP=which(cvsdf == strataLevels[iStrata]), TP=which(cvsdf != strataLevels[iStrata]) )
        mydataSub <- mydataSub0[TPVP[[iStrata]]$TP,]

        if(nrow(mydataSub) == 0){
          print(paste("There was no predictions to work with in trait",iTrait,". Please look at your H2 boundaries. You may be discarding all fields."))
          traitToRemove <- c(traitToRemove,iTrait)
        }else{

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
          phenoDTfile$metadataFieldinst <- merge(baseMeta, ei[,c("fieldinstF",paste0("envIndex_",iTrait))], by.x="fieldinst", by.y="fieldinstF", all.x=TRUE)
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
                interacs <- expand.grid("genoF",interactionsWithGenoTrait)
                interacs<- as.data.frame(interacs[which(as.character(interacs[,1]) != as.character(interacs[,2])),])
                interacsUnlist <- apply(interacs,1,function(x){paste(x,collapse = ":")})
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
              fix <- paste("predictedValue ~",paste(fixedTermTrait, collapse=" + "))
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
                # make sure the matrix only uses the leves for individuals with data
                genoFlevels <- unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"genoF"])
                Ainv <- A <- diag(length(genoFlevels))
                colnames(Ainv) <- rownames(Ainv) <- genoFlevels
                colnames(A) <- rownames(A) <- genoFlevels
                inter <- character()
                onlyInA <- character() # genotypes only present in A and not in dataset
                differ <- character()
                genoMetaData <- list(withMarkandPheno=inter, withPhenoOnly=genoFlevels, withMarkOnly=onlyInA)
              }else{
                genoFlevels <- as.character(unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"genoF"]))

                if(analysisTypeAmat %in% "clm"){ # user provided marker data

                  commonBetweenMandP <- intersect(rownames(genoAmatrix$cleaned$M),genoFlevels)
                  # print(genoFlevels)
                  if(length(commonBetweenMandP) < 2){
                    stop("Markers could not be matched with phenotypes. Please ensure that you have used the right marker file or check the rownames of your marker matrix and ids of your phenotypes.", call. = FALSE)
                  }
                  M <- genoAmatrix$cleaned$M[commonBetweenMandP,]
                  if(markerEffects){ # user wants to save marker effects instead of doing genetic evaluation
                    A <- tcrossprod(M) ## MM' = additive relationship matrix
                  }else{
                    if(ncol(M) > 5000){
                      A <- sommer::A.mat(M[,sample(1:ncol(M), 5000)])
                    }else{
                      A <- sommer::A.mat(M)
                    }
                  }
                  M <- NULL
                }else{ # user provided a relationship matrix
                  if(ncol(genoAmatrix$cleaned) != nrow(genoAmatrix$cleaned)){ # A matrix as data frame
                    Adf <- genoAmatrix$cleaned
                    A <- matrix(NA, nrow=max(Adf$Var1n), ncol = max(Adf$Var2n))
                    A[as.matrix(Adf[,c("Var1n","Var2n")])] = Adf[,"Freq"]
                    A[lower.tri(A)] <- t(A)[lower.tri(A)] # fill the lower triangular
                    rownames(A) <- colnames(A) <- levels(Adf$Var1)
                  }else{ # A matrix as an actual matrix
                    A <- as.matrix(genoAmatrix$cleaned)
                    A[lower.tri(A)] <- t(A)[lower.tri(A)] # fill the lower triangular
                  }
                }

                badGeno <- which(rownames(A) == "") # should have no covariance with anyone
                if(length(badGeno) > 0){A[badGeno,2:ncol(A)]=0; A[2:nrow(A),badGeno]=0} # make zero covariance with this genotype

                badBlankGenotype <- which(colnames(A)=="")
                if(length(badBlankGenotype) > 0){A <- A[-badBlankGenotype,-badBlankGenotype]}

                inter <- intersect(genoFlevels,colnames(A)) # go for sure
                onlyInA <- setdiff(colnames(A),genoFlevels) # genotypes only present in A and not in dataset
                differ <- setdiff(genoFlevels,inter) # are missing in A matrix
                genoMetaData <- list(withMarkandPheno=inter, withPhenoOnly=differ, withMarkOnly=onlyInA)
                # get inverse matrix
                if(length(inter) > 0){ #
                  A1inv <- solve(A[inter,inter] + diag(tolParInv,length(inter), length(inter)))
                }else{A1inv <- matrix(0,0,0)}
                if(length(differ) > 0){ # we have to add individuals without markers or not being part of the GRM?
                  if(length(differ) > 1){ # there's at least 2 inds to be added
                    A2inv <- diag(x=rep(mean(diag(A1inv)),length(differ)))
                  }else{ A2inv <- diag(1)*mean(diag(A1inv)) } # there's only one individual to be added
                  colnames(A2inv) <- rownames(A2inv) <- differ
                }else{A2inv <- matrix(0,0,0)}
                Ainv <- sommer::adiag1(A1inv,A2inv)
                Ainv[lower.tri(Ainv)] <- t(Ainv)[lower.tri(Ainv)] # fill the lower triangular
                colnames(Ainv) <- rownames(Ainv) <- c(colnames(A1inv), colnames(A2inv))
                A1inv <- NULL; A2inv <- NULL;
              }
              levelsInAinv <- colnames(Ainv)
              myGinverse <- ifelse("genoF" %in% rTermsTrait, list(genoF=Ainv), NA) # if genoF is random prepare the Ainv, otherwise just set to NA
              if(!is.na(myGinverse)){names(myGinverse) <- "genoF"} # if Ainv is used assign the name genoF to the relathionship matrix
              if(is.na(myGinverse)){myGinverse <- NULL} # make it NULL for the actual modeling
              if(length(ranran) == 0){ranFormulation=NULL}else{ranFormulation=as.formula(ranran)}
              if(fixErrorVar){myFamily = quasi(link = "identity", variance = "constant")}else{myFamily=gaussian()}
              mix <- try(
                LMMsolver::LMMsolve(fixed =as.formula(fix),
                                    random = ranFormulation,
                                    residual=ranres,
                                    weights = "w",
                                    ginverse = myGinverse,
                                    family = myFamily, # depending on whether we fix the error variance or not we fit use different family
                                    #trace = TRUE,
                                    data = mydataSub, maxit = maxit),
                silent = TRUE
              )
              myGinverse <- NULL      #
              # print(mix$VarDf)
              if(!inherits(mix,"try-error") ){ # if random model runs well try the fixed model
                if(is.null(phenoDTfile$metadataFieldinst)){numericMetas <- character()}
                for(iIndex in c(numericMetas,"envIndex")){
                  if( (iIndex %in% interactionsWithGenoTrait) ){names(mix$ndxCoefficients[[paste0("genoF:",iIndex)]]) <- names(mix$ndxCoefficients$genoF) } # copy same names than main geno effect
                }

                iGenoUnit <- "genoF" # in MET iGenoUnit is always "geno" only in STA we allow for different

                ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(markerEffects & (analysisTypeAmat %in% "clm")){ # user wants to predict marker effects and save those effects for posterior genomic prediction
                  ##
                  AGE <- A
                  mydataSub2 <- mydataSub[which(!is.na(mydataSub$predictedValue)),]
                  mydataSub2 <- droplevels(mydataSub2[which(mydataSub2$geno %in% colnames(A)),])
                  genEvaProv <- cgiarPIPE::geneticEvaluation(
                    fixed =as.formula(fix),
                    random = as.formula(ranran),
                    rcov=ranres,
                    weights = "w",
                    data = mydataSub2,
                    varComp=mix$VarDf,
                    # keep=ttGrouping[iGroup],
                    AFGE=AGE
                  )
                  ##
                  ss <- mix$VarDf; rownames(ss) <- ss$VarComp
                  nEffects <- mix$EDdf$Model; names(nEffects) <- mix$EDdf$Term; nEffects <- nEffects[-length(nEffects)]
                  end <- numeric()
                  for(i in 1:(length(nEffects))){end[i]=sum(nEffects[1:(i)])+1}
                  start <- end - nEffects
                  end <- start + nEffects - 1
                  # extract marker effects
                  pp <- data.frame(matrix(NA,ncol=5,nrow=1)); colnames(pp) <- c("genoF","predictedValue","stdError","rel","trait"); pp <- pp[-1,]
                  listPEVeffects <- list()
                  for(iIndex in c("genoF", paste("genoF",numericMetas,sep=":"), "genoF:envIndex") ){ # iIndex="genoF:envIndex"
                    if( (iIndex %in% names(mix$ndxCoefficients)) ){
                      # names(mix$coefficients[[iIndex]]) <- gsub("genoF_","",names(mix$coefficients[[iIndex]]))
                      names(genEvaProv[[iIndex]]$predictedValue) <- gsub("genoF_","",names(genEvaProv[[iIndex]]$predictedValue) )
                      # myorder <- intersect(names(mix$coefficients[[iIndex]]), rownames( genoAmatrix$cleaned$M))
                      myorder <- intersect( names(genEvaProv[[iIndex]]$predictedValue) , rownames( genoAmatrix$cleaned$M))
                      # toTakePEV <- (start[iIndex]:end[iIndex])[which(names(mix$coefficients[[iIndex]]) %in% myorder )]
                      # toTakePEV <- (start[iIndex]:end[iIndex])[which(names(genEvaProv[[iIndex]]$predictedValue) %in% myorder )]
                      M <- genoAmatrix$cleaned$M[myorder,]
                      MTMMTinv <-t(M)%*%Ainv[myorder,myorder]
                      # a.from.g <-MTMMTinv%*%matrix(mix$coefficients[[iIndex]][myorder],ncol=1)
                      a.from.g <-MTMMTinv%*%matrix(genEvaProv[[iIndex]]$predictedValue[myorder],ncol=1)
                      pev <- genEvaProv[[iIndex]]$pev
                      pev.a.from.g <- t(M)%*%Ainv[myorder,myorder]%*% pev %*% t(Ainv[myorder,myorder])%*%M # [toTakePEV,toTakePEV]
                      se.a.from.g <- sqrt(diag(pev.a.from.g))
                      rel <- 1 - (diag(pev.a.from.g)/ss[iIndex,"Variance"])
                      if(iIndex == "genoF"){
                        newTraitName <- iTrait
                      }else{
                        newTraitName <- paste(iTrait,gsub("genoF","",gsub("genoF:","",iIndex)),sep="-")
                      }
                      pp2 <- data.frame(genoF=rownames(a.from.g), predictedValue=a.from.g, stdError=se.a.from.g, rel=rel,
                                        trait=newTraitName )
                      h2[counter] <- mean(rel); h2.se[counter] <- 1e-6; vg[counter] <- ss[iIndex,"Variance"]; vg.se[counter] <- NA; field[counter] <- "across"; trt[counter] <- paste(iTrait,gsub("-geno","",gsub("genoF:","",iIndex)),sep="-")
                      mu[counter] <- mean(pp2$predictedValue, na.rm=TRUE)
                      pp <- rbind(pp,pp2)
                      listPEVeffects[[iIndex]] <- pev.a.from.g
                      counter <- counter+1
                    }
                  }
                  # pp2 <- data.frame(genoF="(Intercept)", predictedValue=mix$coefficients$`(Intercept)`, stdError=NA, rel=NA, trait=paste(iTrait,"(Intercept)",sep="-") )
                  # pp <- rbind(pp,pp2)
                  Ainv <- NULL
                  listPev[[iTrait]] <- listPEVeffects
                  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                }else{ # user wants to do genetic evaluation
                  if(length(grep("genoF", ranran)) > 0){
                    ######################################################
                    ## do genetic evaluation only if there was genotype as random effect
                    ######################################################
                    if(analysisTypeAmat %in% "clm"){ # if user provided markers
                      allIndsNames <- rownames(genoAmatrix$cleaned$M) # individuals in the marker matrix
                    }else{# if user provided relationship matrix
                      allIndsNames <- colnames(A) # individuals in the relationship matrix
                    }
                    indNamesInAinv <- levelsInAinv # reduced set of individuals in the met analysis
                    indNamesInAinvNoDiffer <- setdiff(indNamesInAinv,differ) # individuals with no markers data that were just expanded
                    indsToBePredicted <- c(indNamesInAinv, setdiff(allIndsNames, indNamesInAinv) )# additional individuals to be predicted in the A matrix
                    # Ainv <- NULL
                    if(length(indsToBePredicted) > 0 ){ # if there's additional inds to be predicted
                      groups0 <- sort(rep(paste0("g",1:1000),batchSizeToPredict)) # groups of 500
                      grouping <- data.frame(id=indsToBePredicted,grp=groups0[1:length(indsToBePredicted)])
                      grouping$grp[1:length(indNamesInAinv)] <- "g1"
                      groups <- unique(grouping$grp)
                    }else{ # no additional inds
                      grouping <- data.frame(id=indNamesInAinv,grp="g1")
                      groups <- unique(grouping$grp)
                    }
                    nEffectsToStore <- length(interactionsWithGenoTrait)+1
                    genEva <- lapply(vector(mode="list",nEffectsToStore), function(x,nn){vector(mode="list",nn)}, nn=4) #list of lists
                    names(genEva) <- names(mix$ndxCoefficients)[grep("genoF",names(mix$ndxCoefficients))]
                    genEva <- lapply(genEva, function(x){names(x) <- c("predictedValue","stdError","pev","R2"); return(x)})
                    ttGrouping <- table(grouping$grp)
                    for(iGroup in 1:length(groups)){ # iGroup=1  # for each group do the genetic evaluation
                      if(verbose){print(paste("Predicting batch",iGroup))}
                      iGrupNames <- grouping[which(grouping$grp == groups[iGroup]),"id"]
                      # make a new A for the ith group
                      present <- unique(c( intersect(colnames(A), indNamesInAinv), iGrupNames ))
                      mydataSub2 <- mydataSub[which(!is.na(mydataSub$predictedValue)),]
                      mydataSub2 <- droplevels(mydataSub2[which(mydataSub2$genoF %in% c(present)),])
                      if(analysisTypeAmat %in% "clm"){ # if user provided markers
                        AGE <- sommer::A.mat(genoAmatrix$cleaned$M[intersect(present,rownames(genoAmatrix$cleaned$M)),])
                      }else{
                        present2 <- intersect(colnames(A),present)
                        AGE <- A[present2,present2]
                      }
                      ## do the genetic evaluation
                      genEvaProv <- cgiarPIPE::geneticEvaluation(
                        fixed =as.formula(fix),
                        random = as.formula(ranran),
                        rcov=ranres,
                        weights = "w",
                        data = mydataSub2,
                        varComp=mix$VarDf,
                        keep=ttGrouping[iGroup],
                        AFGE=AGE
                      )
                      AGE <- NULL

                      # fill the list
                      for(uu in 1:length(genEvaProv)){ # for each genetic effect
                        genEva[[uu]]$predictedValue <- c(genEva[[uu]]$predictedValue,genEvaProv[[uu]]$predictedValue )
                        genEva[[uu]]$stdError <- c(genEva[[uu]]$stdError,genEvaProv[[uu]]$stdError )
                        if(iGroup == 1){
                          genEva[[uu]]$pev <- genEvaProv[[uu]]$pev
                          genEva[[uu]]$R2 <- genEvaProv[[uu]]$R2
                        }else{
                          genEva[[uu]]$pev <- Matrix::bdiag(genEva[[uu]]$pev,genEvaProv[[uu]]$pev )
                          genEva[[uu]]$R2 <- c(genEva[[uu]]$R2,genEvaProv[[uu]]$R2 ) #Matrix::bdiag(genEva[[uu]]$R2,genEvaProv[[uu]]$R2 )
                        }
                      }
                      genEvaProv <- NULL
                    }
                    ######################################################
                    ## end of genetic evaluation function
                    ######################################################
                  }
                  ## now extract the needed values
                  if(iGenoUnit %in% fixedTermTrait){ # user wants fixed effect predictions for genotype
                    shouldBeOne <- which(mix$ndxCoefficients[[iGenoUnit]] == 0)
                    if(length(shouldBeOne) > 0){mix$ndxCoefficients[[iGenoUnit]][shouldBeOne] = 1}
                    predictedValue <- mix$coefMME[mix$ndxCoefficients[[iGenoUnit]]] + mix$coefMME[mix$ndxCoefficients$`(Intercept)`]
                    names(predictedValue) <- names(mix$ndxCoefficients[[iGenoUnit]])
                    dims <- mix$EDdf
                    start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) # we don't add a one because we need the intercept
                    pev <- as.matrix(solve(mix$C))[start:(start+length(predictedValue)-1),start:(start+length(predictedValue)-1)]
                    stdError <- (sqrt(diag(pev)))
                  }else{ # user wants random effect predictions for genotype (main effect)
                    predictedValue <- genEva$genoF$predictedValue
                    stdError <- genEva$genoF$stdError
                    pev <- genEva$genoF$pev
                  }
                  genoF <- gsub("genoF_","", names(predictedValue))
                  pp <- data.frame(genoF,predictedValue,stdError)
                  ss = mix$VarDf; rownames(ss) <- ss$VarComp
                  Vg <- ss["genoF",2]; Vr <- ss["residual",2]
                  if(iGenoUnit %in% fixedTermTrait){ # add reliabilities to the data frame
                    pp$rel <- NA
                  }else{ # if random, reliability can be calculated for main effect
                    # print(names(genEva))
                    pp$rel <- genEva$genoF$R2
                    badRels <- which(pp$rel > 1); if(length(badRels) > 0){pp$rel[badRels] <- 0.9999}
                  }
                  pp$trait <- iTrait # add trait
                  ## heritabilities
                  h2[counter] <- mean(pp$rel); h2.se[counter] <- 1e-6
                  ## genetic variances
                  vg[counter] <- Vg; vg.se[counter] <- 1e-6
                  mu[counter] <- mean(pp$predictedValue, na.rm=TRUE)
                  field[counter] <- "across"; trt[counter] <- iTrait
                  lpv <- sum(mix$EDdf$Model[1:which(mix$EDdf$Term == "genoF")])+1 # to be used as a starting point if random regression is requested

                  ########################################
                  ########################################
                  ########################################
                  ## Validation
                  # if(length(TPVP[[iStrata]]$VP) > 0){}
                  vpData <- mydataSub0[TPVP[[iStrata]]$VP,]
                  vpData <- aggregate(predictedValue~geno, FUN=mean, data=vpData)
                  vpData <- merge(vpData, pp, by.x = "geno", by.y="genoF", all.x = TRUE )
                  plot(vpData$predictedValue.x, vpData$predictedValue.y)
                  TPVPr[[iStrata]] <- cor(vpData$predictedValue.x, vpData$predictedValue.y, use = "complete.obs")
                  ########################################
                  ########################################
                  ########################################
                  # extract sensitivities if interaction is requested
                  if(length(interactionsWithGenoTrait) > 0){ # if there's interactions
                    if( length(intersect(interactionsWithGenoTrait, c("envIndex","timePoint","latitude","longitude","altitude","weather1","weather2"))) > 0 ){

                      for(iInteractionTrait in c("envIndex","timePoint","latitude","longitude","altitude","weather1","weather2")){

                        if( (iInteractionTrait %in% interactionsWithGenoTrait) ){ # iInteractionTrait="envIndex"

                          counter <- counter+1
                          iGenoUnit <- paste0("genoF:",iInteractionTrait)#"genoF:envIndex"
                          if(iGenoUnit %in% fixedTermTrait){ # user wants fixed effect predictions for genotype:envIndex
                            shouldBeOne <- which(mix$ndxCoefficients[[iGenoUnit]] == 0)
                            if(length(shouldBeOne) > 0){mix$ndxCoefficients[[iGenoUnit]][shouldBeOne] = 1}
                            predictedValue <- mix$coefMME[mix$ndxCoefficients[[iGenoUnit]]] #+ mix$coefficients$`(Intercept)`
                            dims <- mix$EDdf
                            start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) # we don't add a one because we need the intercept
                            stdError <- (sqrt(diag(as.matrix(solve(mix$C)))))[start:(start+length(predictedValue)-1)]
                          }else{ # user wants random effect predictions for genotype:envIndex
                            predictedValue <- genEva[[iGenoUnit]]$predictedValue #+ mix$coefficients$`(Intercept)`
                            stdError <- genEva[[iGenoUnit]]$stdError
                          }
                          genoF <- gsub("genoF_","", names(predictedValue))
                          pp2 <- data.frame(genoF,predictedValue,stdError)
                          Vg <- ss[iGenoUnit,2];
                          if(iGenoUnit %in% fixedTermTrait){pp$rel <- 1e-6}else{pp2$rel <- genEva[[iGenoUnit]]$R2}
                          pp2$trait <- paste(iTrait,iInteractionTrait,sep="-")
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

                }
                ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

              }else{
                if(markerEffects){
                  # do nothing
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

              }

              pp$fieldinstF <- "across"
              mydataForGenoType <- droplevels(mydata[which(mydata$trait == iTrait),])
              pp$entryType <- apply(data.frame(pp$genoF),1,function(x){
                found <- which(mydataForGenoType$geno %in% x)
                if(length(found) > 0){
                  x2 <- paste(unique(mydataForGenoType[found,"genoType"]), collapse = "_");
                }else{x2 <- "unknown"}
                return(x2)
              })
              mydataForGenoType <- NULL
              if(!is.null(genoAmatrix)){
                if(markerEffects){
                  pp$entryType <- "markerEffect"
                }else{
                  pp$entryType <- paste(ifelse(pp$genoF %in% differ, "TGV", surrogate),
                                        pp$entryType,
                                        ifelse(pp$genoF %in% onlyInA, "predicted", "tested"),
                                        sep="_")
                }
              }
              ###
              # predictionsList[[counter2]] <- pp;
              counter=counter+1
            }
          }
        }
        


      } ### enf of loop for each strata level
      names(TPVP) <- strataLevels
      cvList[[iSample]] <- TPVP
      cvListR[[iSample]] <- data.frame(r=unlist(TPVPr), vp=names(TPVP), iSample=iSample, trait=iTrait )



    } ## end of loop for the ith cv sample
    names(cvList) <- names(cvListR) <- paste("sample",1:nCvSamples, sep="_")
    counter2 = counter2+1

  } # end of for loop for each trait
  parameters <- do.call(rbind, cvListR)
  #

  #########################################
  ## update databases
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
    stage = NA
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

  # write pipeline metrics
  metrics <- data.frame(value=c(h2,vg,mu, parameters$r),  stdError=c(h2.se,vg.se,rep(1e-6,length(mu)), rep(1e-6, length(parameters$r)) ),
                   fieldinst=c(field,field,field, parameters$vp),  trait=c(trt,trt,trt, parameters$trait),
                   analysisId=id, method=c(rep("(G-PEV)/G",length(h2)),rep("REML",length(vg)), rep("mean",length(mu)), rep("k-fold",length(parameters$trait)) ),
                   traitUnits=NA, parameter=c(rep("r2",length(h2)),rep("Vg",length(vg)), rep("mean",length(mu)), rep("CV",length(parameters$trait)) ),
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", "),
                   stage = NA
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
  result <- list(metrics=metrics, predictions=NA, modeling=modeling, metadata=metadata,
                 cleaned=NA, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst, PEV=NA, genoMetaData=NA
  )
  return(result)#paste("met done:",id))
}
