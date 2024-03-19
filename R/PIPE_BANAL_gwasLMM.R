gwasLMM <- function(
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
    verbose=TRUE
){
  ## THIS FUNCTION PERFORMS A GWAS ANALYSIS 
  ## IS USED IN THE BANAL APP UNDER THE MARKER MODULES
  id <- paste( paste("mes",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "mes"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(is.null(genoAmatrix)){stop("Please provide the markers", call. = FALSE)}
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
  }else{
    surrogate <- "GEBV"#ifelse(length(grep("grm",genoAmatrix$id)) > 0, "GEBV", "EBV") # if id has word grm then is a GEBV otherwise is EBV
  }
  ####

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
  for(iTrait in trait){ # iTrait = trait[1]
    for(kField in 1:length(fields)){
      subData <- mydata[which((mydata$fieldinst == fields[kField]) & mydata$trait == iTrait ),]
      outlier <- boxplot.stats(x=subData$predictedValue,coef=2 )$out
      if(length(outlier) > 0){
        # print(outlier)
        mydata[subData[which(subData$predictedValue %in% outlier),"rowindex"],"predictedValue"] <- NA
      }
    }
  }
  ############################
  # ensure only inds with both stay in the datset
  myrel <- genoAmatrix$cleaned
  common <- intersect(rownames(myrel$M), unique(mydata[,"geno"]))
  if(length(common) == 0){
    stop("There was no intersection of markers and phenotypes provided. Please check your input files.",call. = FALSE)
  }
  M <- myrel$M[common,]
  mydata <- mydata[which(mydata[,"geno"] %in% common),]
  #########################
  ############################
  ## met analysis
  h2 <- vg <-vg.se <- h2.se <- mu <- numeric(); field <- trt <- vector()
  predictionsList <- predictionsBindInterList <- predictionsBindMesList <- listPev <- list(); counter=counter2=1
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

      ## build the environmental index
      ei <- aggregate(predictedValue~fieldinstF, data=mydataSub,FUN=mean); colnames(ei)[2] <- "envIndex"
      ei <- ei[with(ei, order(envIndex)), ]
      ei$envIndex <- ei$envIndex - mean(ei$envIndex)
      ## add the environmental index to the original dataset
      mydataSub <- droplevels(merge(mydataSub,ei, by="fieldinstF"))
      ## define the interactions to be used
      if(!is.null(interactionsWithGeno)){
        interactionsWithGenoTrait <- interactionsWithGeno
        interactionsWithGenoToRemove <- character()
        for(iInter in 1:length(interactionsWithGenoTrait)){
          checkInters <- length(unique(mydataSub[,interactionsWithGenoTrait[iInter]]))
          if (checkInters < 2){ # there needs to be at least more than one level
            interactionsWithGenoToRemove <- c(interactionsWithGenoToRemove,interactionsWithGenoTrait[iInter])
          }
        }
        interactionsWithGenoTrait <- setdiff(interactionsWithGenoTrait,interactionsWithGenoToRemove)
      }else{
        interactionsWithGenoTrait <- interactionsWithGeno
      }

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

          #
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
            stop("Markers are needed to estimate marker effects", call. = FALSE)
          }else{
            genoFlevels <- as.character(unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"genoF"]))
            Msub <- M[genoFlevels,]
            MMT <- tcrossprod(Msub) ## MM' = additive relationship matrix
            MMTinv<-solve(MMT + diag(1e-3, ncol(MMT), ncol(MMT))) ## inverse of MM'
            MMTinv[lower.tri(MMTinv)] <- t(MMTinv)[lower.tri(MMTinv)] # ensure is symmetric

          }

          myGinverse <- ifelse("genoF" %in% rTermsTrait, list(genoF=MMTinv), NA) # if genoF is random prepare the MMTinv, otherwise just set to NA
          if(!is.na(myGinverse)){names(myGinverse) <- "genoF"} # if MMTinv is used assign the name genoF to the relathionship matrix
          if(is.na(myGinverse)){myGinverse <- NULL} # make it NULL for the actual modeling
          if(length(ranran) == 0){ranFormulation=NULL}else{ranFormulation=as.formula(ranran)}
          if(fixErrorVar){myFamily = quasi(link = "identity", variance = "constant")}else{myFamily=gaussian()}
          mix <- try(
            LMMsolver::LMMsolve(fixed =as.formula(fix),
                                random = ranFormulation,
                                residual=ranres,
                                weights = "w",
                                ginverse = myGinverse,#list(genoF=MMTinv),
                                family = myFamily, # quasi(link = "identity", variance = "constant"),
                                #trace = TRUE,
                                data = mydataSub, maxit = maxit),
            silent = TRUE
          )
          # print(mix$VarDf)
          if(!inherits(mix,"try-error") ){ # if random model runs well try the fixed model
            if( ("envIndex" %in% interactionsWithGenoTrait) ){names(mix$coefficients$`genoF:envIndex`) <- names(mix$coefficients$genoF) }

            iGenoUnit <- "genoF" # in MET iGenoUnit is always "geno" only in STA we allow for different
            ## do genetic evaluation only if there was random effects
            if(!is.null(ranFormulation)){

              # reorder the matrices
              names(mix$coefficients$genoF) <- gsub("genoF_","",names(mix$coefficients$genoF))
              myorder <- names(mix$coefficients$genoF);
              Msub <- Msub[myorder,]
              MMT <- MMT[myorder,myorder]
              MMTinv <- MMTinv[myorder,myorder]
              MTMMTinv <-t(Msub)%*%MMTinv
              a.from.g <-MTMMTinv%*%matrix(mix$coefficients$genoF,ncol=1)
              pev <- solve(mix$C)

              # extract predicted values, standard errors and PEV
              vgeno <- which(mix$EDdf$Term == "genoF")
              vfixed <- 1:(vgeno-1)
              nEffects <- c(sum(mix$EDdf$Model[vfixed]), mix$EDdf$Model[vgeno])
              end <- numeric()
              for(i in 1:(length(nEffects))){end[i]=sum(nEffects[1:(i)])+1}
              start <- end - nEffects
              end <- start + nEffects - 1
              # remove fixed effects
              start2 <- start[-1]
              end2 <- end[-1]
              # var.g <- kronecker(MMT,mixGBLUP$sigma$`u:geno`) - mixGBLUP$PevU$`u:geno`$predictedValue
              pev.a.from.g <- t(Msub)%*%MMTinv%*% pev[start2:end2,start2:end2] %*% t(MMTinv)%*%Msub
              se.a.from.g <- sqrt(diag(pev.a.from.g))# + min(diag(pev.a.from.g))+ 1e-7)
              # t.stat.from.g <- a.from.g/se.a.from.g # t-statistic
              # pvalGBLUP <- dt(t.stat.from.g,df=n-k-1) #
              # minusLog10pvalGBLUP <- -log10(pvalGBLUP)
              rownames(mix$VarDf) <- mix$VarDf$VarComp
              freq = apply(Msub+1,2,mean)/2
              deno = (2 * mean(freq*(1-freq)) ) / ncol(Msub)
              rel <- 1 - (se.a.from.g^2)/(as.numeric(mix$VarDf["genoF","Variance"])/deno )
              listPev[[iTrait]] <- pev.a.from.g

            }
            if(iGenoUnit %in% fixedTermTrait){ # user wants fixed effect predictions for genotype
              stop("Not enabled", call. = FALSE)
            }else{ # user wants random effect predictions for genotype

              predictedValue <- a.from.g # genEva$genoF$predictedValue
              stdError <- se.a.from.g # genEva$genoF$stdError
              pev <- pev.a.from.g# genEva$genoF$pev
            }
            genoF <- gsub("genoF_","", rownames(predictedValue))
            # pp <- data.frame(genoF,predictedValue,stdError)
            ss = mix$VarDf; rownames(ss) <- ss$VarComp
            Vg <- ss["genoF",2]/deno; Vr <- ss["residual",2]
            if(iGenoUnit %in% fixedTermTrait){ # for purely fixed effect models
              stop("Not enabled", call. = FALSE)
            }else{ # if random reliability can be calculated
            }
            # pp$trait <- iTrait
            ## heritabilities
            h2[counter] <- mean(rel); h2.se[counter] <- 1e-6
            ## genetic variances
            vg[counter] <- Vg; vg.se[counter] <- 1e-6
            mu[counter] <- mean(predictedValue, na.rm=TRUE)
            field[counter] <- "across"; trt[counter] <- iTrait
            lpv <- length(predictedValue) + 1 # to be used as a starting point if random regression is requested
            # extract sensitivities if interaction is requested
            predictionsBindInterList[[iTrait]] <- data.frame(analysisId=id, pipeline= paste(sort(unique(mydataSub$pipeline)),collapse=", "),
                                                             trait=iTrait, genoCode=1, geno="all", mother="unknown", father="unknown",
                                                             genoType="intercept", genoYearOrigin=paste(sort(unique(mydataSub$genoYearOrigin)),collapse=", "),
                                                             genoYearTesting=paste(sort(unique(mydataSub$genoYearTesting)),collapse=", "),
                                                             fieldinst="across",#paste(sort(unique(mydataSub$fieldinst)),collapse=", "),
                                                             predictedValue=mix$coefficients$`(Intercept)`, stdError=sqrt(pev.a.from.g[1,1]), rel=1,
                                                             stage=paste(sort(unique(mydataSub$stage)),collapse=", ")
            )
            predictionsBindMesList[[iTrait]] <- data.frame(analysisId=id, pipeline= paste(sort(unique(mydataSub$pipeline)),collapse=", "),
                                                           trait=iTrait, genoCode=1:nrow(a.from.g), geno=rownames(a.from.g), mother="unknown", father="unknown",
                                                           genoType="markerEffect", genoYearOrigin=paste(sort(unique(mydataSub$genoYearOrigin)),collapse=", "),
                                                           genoYearTesting=paste(sort(unique(mydataSub$genoYearTesting)),collapse=", "),
                                                           fieldinst="across",#paste(sort(unique(mydataSub$fieldinst)),collapse=", "),
                                                           predictedValue=a.from.g[,1], stdError=se.a.from.g, rel=rel, stage=paste(sort(unique(mydataSub$stage)),collapse=", ")
            )

            if( ("envIndex" %in% interactionsWithGenoTrait)  ){
              counter <- counter+1
              iGenoUnit <- "genoF:envIndex"
              if(iGenoUnit %in% fixedTermTrait){ # user wants fixed effect predictions for genotype:envIndex
                predictedValue <- mix$coefficients[[iGenoUnit]] #+ mix$coefficients$`(Intercept)`
                dims <- mix$EDdf
                start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) # we don't add a one because we need the intercept
                stdError <- (sqrt(diag(as.matrix(solve(mix$C)))))[start:(start+length(predictedValue)-1)]
              }else{ # user wants random effect predictions for genotype:envIndex
                predictedValue <- genEva$`genoF:envIndex`$predictedValue #+ mix$coefficients$`(Intercept)`
                stdError <- genEva$`genoF:envIndex`$stdError
              }


              genoF <- gsub("genoF_","", names(predictedValue))
              pp2 <- data.frame(genoF,predictedValue,stdError)
              Vg <- ss["genoF:envIndex",2]; #Vr <- ss["residual",2]
              if(iGenoUnit %in% fixedTermTrait){pp$rel <- 1e-6}else{pp2$rel <- diag(genEva$`genoF:envIndex`$R2)}

              # pp2$rel <- abs(1 - (pp2$stdError^2)/(Vg))
              pp2$trait <- paste(iTrait,"adaptability",sep="-")
              ## heritabilities
              h2[counter] <- mean(pp2$rel); h2.se[counter] <- 1e-6
              ## genetic variances
              vg[counter] <- Vg; vg.se[counter] <- 1e-6
              mu[counter] <- mean(pp2$predictedValue, na.rm=TRUE)
              field[counter] <- "across"; trt[counter] <- paste(iTrait,"adaptability",sep="-")
              pp <- rbind(pp,pp2)
            }
          }else{
            cat("Model was singular for trait",iTrait,". Skipping this trait")
          }
          ###
          # predictionsList[[counter2]] <- pp;
          counter=counter+1
        }
      }
    }
    counter2 = counter2+1
  }
  # print(vg)
  trait <- setdiff(trait,traitToRemove)
  heritLB <- heritLB[which(traitOrig %in% trait)]
  heritUB <- heritUB[which(traitOrig %in% trait)]
  if(length(predictionsBindMesList) == 0){stop("There was no predictions to work with. Please look at your H2 boundaries. You may be discarding all fields.",call. = FALSE)}
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","mother","father","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")
  predictionsBindMes <- do.call(rbind, predictionsBindMesList)
  predictionsBindInt <- do.call(rbind, predictionsBindInterList)
  finalPreds <- rbind(predictionsBindMes[,predcols],predictionsBindInt[,predcols])

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
    stage = paste(sort(unique(predictionsBindMes$stage)),collapse=", ")
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
  metrics <- data.frame(value=c(h2,vg,mu),  stdError=c(h2.se,vg.se,rep(1e-6,length(mu))),
                   fieldinst=c(field,field,field),  trait=c(trt,trt,trt),
                   analysisId=id, method=c(rep("(G-PEV)/G",length(h2)),rep("REML",length(vg)), rep("mean",length(mu)) ),
                   traitUnits=NA, parameter=c(rep("r2",length(h2)),rep("Vg",length(vg)), rep("mean",length(mu)) ),
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", "),
                   stage = paste(sort(unique(predictionsBindMes$stage)),collapse=", ")
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
  result <- list(metrics=metrics, predictions=finalPreds, modeling=modeling, metadata=metadata,
                 cleaned=NA, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst, PEV=listPev,
                 ref.alleles=genoAmatrix$ref.alleles
  )
  return(result)
}
