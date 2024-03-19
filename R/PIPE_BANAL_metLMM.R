metLMM <- function(
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
    # fixErrorVar=TRUE,
    markerEffects=FALSE,
    batchSizeToPredict=500,
    tolParInv=1e-4,
    verbose=TRUE
){
  ## THIS FUNCTION PERFORMS A MULT TRIAL ANALYSIS USING LMM SOLVER
  ## IS USED IN THE BANAL APP UNDER THE GENETIC EVALUATION MODULES
  type <- "met"
  id <- paste( paste(type,cgiarPIPE::idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  if(is.null(phenoDTfile)){stop("Please provide the phenotype file", call. = FALSE)}
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
  # cleaning <- phenoDTfile$outliers 
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
  ##############################
  ## met analysis
  h2 <- vg <-vg.se <- h2.se <- mu <- numeric(); field <- trt <- vector()
  predictionsList <- list(); listPev <- list(); counter=counter2=1
  traitToRemove <- character()
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
    if(verbose){print(paste("Fields included:",paste(goodFields,collapse = ",")))}
    ## remove records without marker data if marker effects
    if(markerEffects){
      if(analysisTypeAmat != "clm"){stop("The marker effects (rrBLUP) method is only possible providing a marker matrix", call. = FALSE)}
      mydataSub <- mydataSub[which(mydataSub$geno %in% rownames(genoAmatrix$cleaned$M)),]
      Mtrait <- genoAmatrix$cleaned$M[which(rownames(genoAmatrix$cleaned$M) %in% mydataSub$geno),]
    }
    LGrp <- list();   groupTrait <- NULL
    ## next step
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
            if(markerEffects){
              reduced <- with(mydataSub,cgiarPIPE::redmm(x=geno,M=Mtrait, nPC=nPC, returnLam = TRUE))
              LGrp[["QTL"]] <- c((ncol(mydataSub)+1):(ncol(mydataSub)+ncol(reduced$Z)))
              mydataSub <- cbind(mydataSub,reduced$Z)
              rTermsTrait <- randomTerm[which(apply(data.frame(randomTerm),1,function(x){length(table(mydataSub[,x]))}) > 1)]
              rTermsTrait <- setdiff(rTermsTrait,"genoF")
              rTermsTrait <- c("grp(QTL)",rTermsTrait)
            }else{
              rTermsTrait <- randomTerm[which(apply(data.frame(randomTerm),1,function(x){length(table(mydataSub[,x]))}) > 1)]
            }
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
            if(markerEffects){
              for(iInteraction in unique(interactionsWithGenoTrait)){ # iInteraction <- unique(interactionsWithGenoTrait)[1]
                LGrp[[paste0("QTL",iInteraction)]] <- c((ncol(mydataSub)+1):(ncol(mydataSub)+ncol(reduced$Z)))
                mydataSub <- cbind(mydataSub,reduced$Z*mydataSub[,iInteraction])
                rTermsTrait <- c(rTermsTrait,paste0("grp(QTL",iInteraction,")"))
              }
            }else{ 
              interacs <- expand.grid("genoF",interactionsWithGenoTrait)
              interacs<- as.data.frame(interacs[which(as.character(interacs[,1]) != as.character(interacs[,2])),])
              interacsUnlist <- apply(interacs,1,function(x){paste(x,collapse = ":")})
              rTermsTrait <- c(rTermsTrait,interacsUnlist)
            }
          }else{
            rTermsTrait <- rTermsTrait
          }
          rTermsTrait <- setdiff(rTermsTrait, fixedTermTrait)
          
          if(length(rTermsTrait) == 0){
            ranran <- NULL
            myGinverse <- NULL # if no random effects we don't have a g inverse
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
            myGinverse <- list(genoF=Ainv)
            levelsInAinv <- colnames(Ainv)
            genoMetaData <- list(withMarkandPheno=inter, withPhenoOnly=genoFlevels, withMarkOnly=onlyInA)
          }else{
            if(markerEffects){
              genoFlevels <- unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"genoF"])
              inter <- character()
              onlyInA <- character() # genotypes only present in A and not in dataset
              differ <- character()
              myGinverse <- NULL
              groupTrait <- LGrp
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
              if(analysisTypeAmat %in% "clm"){ # user provided marker data
                onlyInA <- setdiff(rownames(genoAmatrix$cleaned$M),genoFlevels) # genotypes only present in A and not in dataset (not tested)
              }else{
                onlyInA <- setdiff(rownames(A),genoFlevels) # genotypes only present in A and not in dataset (not tested)
              }
              differ <- setdiff(genoFlevels,inter) # are missing in A matrix
              genoMetaData <- list(withMarkandPheno=inter, withPhenoOnly=differ, withMarkOnly=onlyInA)
              # get inverse matrix
              if(length(inter) > 0){ #
                A1inv <- solve(A[inter,inter] + diag(tolParInv,length(inter), length(inter)))
              }else{A1inv <- matrix(0,0,0)}
              if(length(differ) > 0){ # we have to add individuals without markers or not being part of the GRM?
                if(length(differ) > 1){ # there's at least 2 inds to be added
                  A2inv <- diag(x=rep(median(diag(A1inv)),length(differ)))
                }else{ A2inv <- diag(1)*median(diag(A1inv)) } # there's only one individual to be added
                colnames(A2inv) <- rownames(A2inv) <- differ
              }else{A2inv <- matrix(0,0,0)}
              Ainv <- sommer::adiag1(A1inv,A2inv)
              Ainv[lower.tri(Ainv)] <- t(Ainv)[lower.tri(Ainv)] # fill the lower triangular
              colnames(Ainv) <- rownames(Ainv) <- c(colnames(A1inv), colnames(A2inv))
              A1inv <- NULL; A2inv <- NULL;
              levelsInAinv <- colnames(Ainv)
              myGinverse <- list(genoF=Ainv)
            }
          }
          
          if(length(ranran) == 0){ranFormulation=NULL}else{ranFormulation=as.formula(ranran)}
          if(useWeights){
            weightsFormulation="w"
            if(verbose){
              print("Using weights in the analysis. Residual variance will be fixed to 1.")
            }
          }else{
            weightsFormulation=NULL
            if(verbose){
              print("Ignoring weights in the analysis. Residual variance will be estimated.")
            }
          }
          mix <- try(
            LMMsolver::LMMsolve(fixed =as.formula(fix),
                                random = ranFormulation,
                                residual=ranres,
                                weights = weightsFormulation,
                                ginverse = myGinverse,
                                group = groupTrait,
                                family = eval(parse(text = traitFamily[iTrait])),
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
            if(markerEffects){ # user wants to predict using marker effects
              ##
              Mfull <- genoAmatrix$cleaned$M
              if((ncol(Mfull) < nrow(Mfull)) | nPC==0){M2 <- Mfull}else{  M2 <- tcrossprod(Mfull)}
              xx2 = with(mydataSub, cgiarPIPE::redmm(x=geno, M=M2, nPC=nPC, returnLam = TRUE)) # we need the new lambda for the fullt marker matrix
              ss <- mix$VarDf;  rownames(ss) <- ss$VarComp
              pp <- list()
              for(iGroup in names(LGrp)){ # iGroup <- names(LGrp)[1]  for each rrBLUP effect
                shouldBeOne <- which(mix$ndxCoefficients[[iGroup]] == 0)
                if(length(shouldBeOne) > 0){mix$ndxCoefficients[[iGroup]][shouldBeOne] = 1}
                blup <- mix$coefMME[mix$ndxCoefficients[[iGroup]]]
                names(blup) <- names(mix$ndxCoefficients[[iGroup]])
                names(blup) <- gsub(paste0(iGroup,"_"),"",names(blup))
                predictedValue <- (xx2$Lam %*% blup[colnames(xx2$Lam)]) + mix$coefMME[mix$ndxCoefficients$`(Intercept)`]
                if(length(shouldBeOne) > 0){predictedValue[1] = mix$coefMME[mix$ndxCoefficients$`(Intercept)`]}
                # names(predictedValue) <- names(mix$ndxCoefficients[[iGroup]])
                dims <- mix$EDdf
                start <- sum(dims[1:(which(dims$Term == iGroup) - 1),"Model"]) # we don't add a one because we need the intercept
                nEffects <- length(mix$coefMME[mix$ndxCoefficients[[iGroup]]])
                pev <- as.matrix(solve(mix$C))[start:(start+nEffects-1),start:(start+nEffects-1)]
                pev <- xx2$Lam %*% pev %*% t(xx2$Lam) 
                stdError <- (sqrt(Matrix::diag(pev)))
                rel <- 1 - (stdError/as.numeric(var(predictedValue)))
                pp[[iGroup]] <- data.frame(genoF=rownames(predictedValue), predictedValue=predictedValue, stdError=stdError, rel=rel,
                                           trait=paste0(iTrait,"_",iGroup) )
                h2[counter] <- mean(rel); h2.se[counter] <- 1e-6; vg[counter] <- ss[iGroup,"Variance"]; vg.se[counter] <- NA; field[counter] <- "across"; trt[counter] <- paste0(iTrait,"_",iGroup)
                mu[counter] <- mean(predictedValue, na.rm=TRUE)
                counter <- counter+1
              }
              pp <- do.call(rbind,pp)
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
                  grouping$grp[1:length(indNamesInAinv)] <- "g1" ## all the ones in the model go to group 1
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
                  AGE <- AGE + diag(1e-3,nrow(AGE),nrow(AGE)) 
                  ######################################
                  ## do the genetic evaluation
                  genEvaProv <- cgiarPIPE::geneticEvaluation2(
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
                  # print("finished")
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
                if(length(shouldBeOne) > 0){predictedValue[1] = mix$coefMME[mix$ndxCoefficients$`(Intercept)`]}
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
              h2[counter] <- median(pp$rel); h2.se[counter] <- 1e-6
              ## genetic variances
              vg[counter] <- Vg; vg.se[counter] <- 1e-6
              mu[counter] <- mean(pp$predictedValue, na.rm=TRUE)
              field[counter] <- "across"; trt[counter] <- iTrait
              lpv <- sum(mix$EDdf$Model[1:which(mix$EDdf$Term == "genoF")])+1 # to be used as a starting point if random regression is requested
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
                      h2[counter] <- median(pp2$rel); h2.se[counter] <- 1e-6
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
            
          }else{ # if model failed
            if(markerEffects){
              # do nothing
            }else{
              if(verbose){ cat(paste("Mixed model failed for this combination. Aggregating and assuming h2 = 0 \n"))}
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
              x2 <- paste(sort(unique(toupper(trimws(mydataForGenoType[found,"genoType"])))), collapse = "#");
            }else{x2 <- "unknown"}
            return(x2)
          })
          mydataForGenoType <- NULL
          if(!is.null(genoAmatrix)){
            if(markerEffects){
              pp$entryType <- "markerEffect"
            }else{
              pp$entryType <- paste(ifelse(as.character(pp$genoF) %in% differ, "TGV", surrogate),
                                    pp$entryType,
                                    ifelse(as.character(pp$genoF) %in% onlyInA, "predicted", "tested"),
                                    sep="_")
            }
          }
          ###
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
  # save desire file
  des <- desire(trait=trt,h2=h2, G=G
                # pathFile=file.path(wd,"desire",paste0("desire_",id,".txt"))
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
                 cleaned=phenoDTfile$cleaned, outliers=NA, desire=des, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst, PEV=listPev, genoMetaData=genoMetaData
  )
  return(result)#paste("met done:",id))
}
