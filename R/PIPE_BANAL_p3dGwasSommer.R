p3dSommer <- function(
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
    tolParInv=1e-4,
    verbose=TRUE
){
  ## THIS FUNCTION CALCULATES A P3D GWAS 
  ## IS USED IN THE BANAL APP IN THE MARKER MODULES
  type <- "gwa"
  id <- paste( paste(type,idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  heritLB <- rep(heritLB,100)
  heritUB <- rep(heritUB,100)
  traitOrig <- trait
  common <- intersect(fixedTerm,randomTerm)
  fixedTerm <- setdiff(fixedTerm,common)
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
  ############################
  ## gwa analysis
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
    pipeline_metricsSub <- pipeline_metrics[which(pipeline_metrics$trait == iTrait & pipeline_metrics$parameter %in% c("plotH2","H2")),]
    goodFields <- pipeline_metricsSub[which((pipeline_metricsSub$value > heritLB[counter2]) & (pipeline_metricsSub$value < heritUB[counter2])),"fieldinst"]
    mydataSub <- mydataSub[which(mydataSub$fieldinst %in% goodFields),]

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
      ei <- aggregate(predictedValue~fieldinstF, data=mydataSub,FUN=mean); colnames(ei)[2] <- "envIndex"
      ei <- ei[with(ei, order(envIndex)), ]
      ei$envIndex <- ei$envIndex - mean(ei$envIndex)
      ## add the environmental index to the original dataset
      mydataSub <- droplevels(merge(mydataSub,ei, by="fieldinstF", all.x = TRUE))
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
          rTermsTrait <- cgiarBase::replaceValues(rTermsTrait, Search = "genoF", Replace = "vsr(genoF, Gu=A)")

          if(length(rTermsTrait) == 0){
            ranran <- NULL#"~NULL"
          }else{
            ranran <- paste("~",paste(rTermsTrait, collapse=" + ") )
          }
          fix <- paste("predictedValue ~",paste(fixedTermTrait, collapse=" + "))
          # Ai <- PED$ginv; attr(Ai, "INVERSE") <- TRUE
          if(!is.null(residualByTrait)){
            ranres <- as.formula(paste0("~vsr(dsr(",residualByTrait,") , units)"))
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
            stop("A matrix is required")
          }else{
            genoFlevels <- as.character(unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"genoF"]))

            if(analysisTypeAmat %in% "clm"){ # user provided marker data
              commonBetweenMandP <- intersect(rownames(genoAmatrix$cleaned$M),genoFlevels)
              if(length(commonBetweenMandP) < 2){
                stop("Markers could not be matched with phenotypes. Please ensure that you have used the right marker file or check the rownames of your marker matrix and ids of your phenotypes.", call. = FALSE)
              }
              M <- genoAmatrix$cleaned$M[sort(commonBetweenMandP),]
              set.seed(1234)
              pickForA <- sample(colnames(M), min(c(10000,ncol(M))))
              A <<- sommer::A.mat(M[,pickForA])

            }else{ # user provided a relationship matrix
              stop("Gwas only takes a marker matrix")
            }

            badGeno <- which(rownames(A) == "") # should have no covariance with anyone
            if(length(badGeno) > 0){A[badGeno,2:ncol(A)]=0; A[2:nrow(A),badGeno]=0} # make zero covariance with this genotype

            badBlankGenotype <- which(colnames(A)=="")
            if(length(badBlankGenotype) > 0){A <- A[-badBlankGenotype,-badBlankGenotype]}

            ## ensure GWAS is safe (only inds with markers)
            mydataSub <- droplevels(mydataSub[which(mydataSub$geno %in% colnames(A)),])
            genoFlevels <- as.character(unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"genoF"]))


            inter <- intersect(genoFlevels,colnames(A)) # go for sure
            onlyInA <- setdiff(colnames(A),genoFlevels) # genotypes only present in A and not in dataset
            differ <- setdiff(genoFlevels,inter) # are missing in A matrix
            genoMetaData <- list(withMarkandPheno=inter, withPhenoOnly=differ, withMarkOnly=onlyInA)
          }
          levelsInAinv <- colnames(A)
          if(length(ranran) == 0){ranFormulation=NULL}else{ranFormulation=as.formula(ranran)}
          # myFamily=gaussian()
          W <- diagonal(1/(mydataSub$sdtError^2))
          if(!is.null(residualByTrait)){
            mix <- try(
              sommer::GWAS(fixed =as.formula(fix),
                           random = ranFormulation,
                           rcov=ranres,
                           # W = W,
                           gTerm = "u:genoF",
                           M=M,
                           data = mydataSub, nIters = maxit),
              silent = TRUE
            )
          }else{
            mix <- try(
              sommer::GWAS(fixed =as.formula(fix),
                           random = ranFormulation,
                           gTerm = "u:genoF",
                           # W = W,
                           M=M,
                           data = mydataSub, nIters = maxit),
              silent = TRUE
            )
          }
          if(!inherits(mix,"try-error") ){ # if random model runs well try the fixed model

            iGeno <- length(mix$U)
            blup <- mix$U[[iGeno]]
            mix$pheno <- data.frame(geno=colnames(A), pheno=blup)
            # remove some objects to make the object lighter
            mix$Vi <- NULL
            mix$P <- NULL
            mix$VarU <- NULL
            mix$PevU <- NULL
            predictionsList[[iTrait]] <- mix

          }else{
            cat("trait failed")
          }

        }
      }
    }
    counter2 = counter2+1
  }
  #

  #########################################
  ## update databases
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
  # pm <- data.frame(value=c(h2,vg,mu),  stdError=c(h2.se,vg.se,rep(1e-6,length(mu))),
  #                  fieldinst=c(field,field,field),  trait=c(trt,trt,trt),
  #                  analysisId=id, method=c(rep("(G-PEV)/G",length(h2)),rep("REML",length(vg)), rep("mean",length(mu)) ),
  #                  traitUnits=NA, parameter=c(rep("r2",length(h2)),rep("Vg",length(vg)), rep("mean",length(mu)) ),
  #                  pipeline=paste(sort(unique(mydata$pipeline)),collapse=", "),
  #                  stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  # )
  # save desire file
  # des <- desire(trait=trt,h2=h2, G=G
  #               # pathFile=file.path(wd,"desire",paste0("desire_",id,".txt"))
  # )
  ##
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  if(is.null(phenoDTfile$metadataFieldinst)){
    metadataFieldinst=NA
  }else{
    metadataFieldinst=phenoDTfile$metadataFieldinst
  }
  result <- list(metrics=NA, predictions=predictionsList, modeling=mod, metadata=db.params,
                 cleaned=NA, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst, genoAmatrix=genoAmatrix
  )
  return(result)#paste("met done:",id))
}
