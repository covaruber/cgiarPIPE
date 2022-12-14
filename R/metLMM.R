metLMM <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    fixedTerm= c("fieldinstF"),
    randomTerm=c("genoF"),
    residualBy=NULL,
    interactionsWithGeno=NULL,
    trait= NULL, # per trait
    heritLB= 0.15,
    heritUB= 0.95,
    scaledDesire=TRUE,# wd=NULL,
    genoAmatrix=NULL,
    maxit=50,
    verbose=TRUE
){

  id <- paste("met",idGenerator(5,5),sep="")
  type <- "met"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}

  ############################
  # loading the dataset
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  mydata <- phenoDTfile$predictions #readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  pipeline_metrics <- phenoDTfile$metrics #readRDS(file.path(wd,"metrics",paste0(phenoDTfile)))

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
  ## met analysis
  h2 <- vg <- cl <-vg.se <- h2.se <- cl.se <- mu <- numeric(); field <- trt <- vector()
  predictionsList <- list(); counter=1
  for(iTrait in trait){ # # iTrait = trait[1]  iTrait="GYKGPHA"
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    # subset data
    mydataSub <- droplevels(mydata[which(mydata$trait == iTrait),])
    # remove bad fieldinst
    pipeline_metricsSub <- pipeline_metrics[which(pipeline_metrics$trait == iTrait & pipeline_metrics$parameter == "H2"),]
    goodFields <- pipeline_metricsSub[which((pipeline_metricsSub$value > heritLB) & (pipeline_metricsSub$value < heritUB)),"fieldinst"]
    mydataSub <- mydataSub[which(mydataSub$fieldinst %in% goodFields),]

    if(nrow(mydataSub) == 0){stop(paste("There was no predictions to work with in trait",iTrait,". Please look at your H2 boundaries. You may be discarding all fields."),call. = FALSE)}

    mydataSub$genoF <- as.factor(mydataSub$geno)
    mydataSub$fieldinstF <- as.factor(mydataSub$fieldinst)
    mydataSub$pipelineF <- as.factor(mydataSub$pipeline)
    mydataSub$stageF <- as.factor(mydataSub$stage)
    mydataSub$genoYearOriginF <- as.factor(mydataSub$genoYearOrigin)
    mydataSub$genoYearTestingF <- as.factor(mydataSub$genoYearTesting)
    mydataSub <- mydataSub[which(mydataSub$geno != ""),]

    # do analysis
    if(!is.na(var(mydataSub[,"predictedValue"],na.rm=TRUE))){ # if there's variance
      if( var(mydataSub[,"predictedValue"], na.rm = TRUE) > 0 ){
        checks <- mydataSub[which(mydataSub[,"genoType"] == "check"),"geno"]
        # make sure the terms to be fitted have more than one level
        if(!is.null(randomTerm)){
          rTermsTrait <- randomTerm[which(apply(data.frame(randomTerm),1,function(x){length(table(mydataSub[,x]))}) > 1)]
        }else{rTermsTrait=NULL}
        if(!is.null(fixedTerm)){
          fixedTermTraitMinus <- setdiff(fixedTerm,"1")
          if(length(fixedTermTraitMinus) > 0){
            fixedTermTrait <- fixedTerm[which(apply(data.frame(fixedTermTraitMinus),1,function(x){length(table(mydataSub[,x]))}) > 1)]
            fixedTermTrait <- c("1",fixedTermTrait)
          }else{
            fixedTermTrait <- fixedTerm
          }
        }else{fixedTermTrait=NULL}
        if(!is.null(residualBy)){
          residualByTrait <- residualBy[which(apply(data.frame(residualBy),1,function(x){length(table(mydataSub[,x]))}) > 1)]
          if(length(residualByTrait) == 0){residualByTrait=NULL}
        }else{residualByTrait=NULL}


        if(!is.null(interactionsWithGeno)){
          interacs <- expand.grid("genoF",interactionsWithGeno)
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

        if(is.null(genoAmatrix)){
          # make sure the matrix only uses the leves for individuals with data
          genoFlevels <- unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"genoF"])
          Ainv <- diag(length(genoFlevels))
          colnames(Ainv) <- rownames(Ainv) <- genoFlevels
        }else{
          genoFlevels <- as.character(unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"genoF"]))
          A <- genoAmatrix$cleaned
          inter <- intersect(genoFlevels,colnames(A)) # go for sure
          differ <- setdiff(inter,genoFlevels) # are missing
          if(length(inter) > 0){
            A1 <- A[inter,inter]
            A1inv <- solve(A1 + diag(1e-5,ncol(A1), ncol(A1)))
          }else{A1inv <- matrix(0,0,0)}
          if(length(differ) > 0){
            A2inv <- diag(length(differ))
            colnames(A2inv) <- rownames(A2inv) <- differ
          }else{A2inv <- matrix(0,0,0)}
          Ainv <- sommer::adiag1(A1inv,A2inv)
          colnames(Ainv) <- rownames(Ainv) <- c(colnames(A1inv), colnames(A2inv))
        }

        myGinverse <- ifelse("genoF" %in% rTermsTrait, list(genoF=Ainv), NULL)
        if(!is.null(myGinverse)){names(myGinverse) <- "genoF"}
        # mydataSub2 <- mydataSub[which(mydataSub$predictedValue > 0),]
        mix <- try(
          LMMsolver::LMMsolve(fixed =as.formula(fix),
                              random = as.formula(ranran),
                              residual=ranres,
                              weights = "w",
                              ginverse = myGinverse,#list(genoF=Ainv),
                              #trace = TRUE,
                              data = mydataSub, maxit = maxit),
          silent = TRUE
        )
        if(!inherits(mix,"try-error") ){ # if random model runs well try the fixed model
          predictedValue <- mix$coefficients$genoF + mix$coefficients$`(Intercept)`
          stdError <- (sqrt(diag(as.matrix(solve(mix$C)))))[1:length(predictedValue)]
          genoF <- gsub("genoF_","", names(predictedValue))
          pp <- data.frame(genoF,predictedValue,stdError)
          ss = mix$VarDf; rownames(ss) <- ss$VarComp
          vg <- ss["genoF",2]; vr <- ss["residual",2]
          pp$rel <- abs(1 - (pp$stdError^2)/(vg))
          ## heritabilities
          h2[counter] <- mean(pp$rel); h2.se[counter] <- 1e-6
          ## genetic variances
          vg[counter] <- vg; vg.se[counter] <- 1e-6
        }else{
          if(verbose){ cat(paste("Aggregating and assuming h2 = 0 \n"))}
          pp <- aggregate(predictedValue ~ genoF, FUN=mean, data=mydataSub)
          pp$stdError <- 1
          pp$status <- "Aggregated"
          pp$rel <- 1e-6
          ## heritabilities # h2Pred <- abs(round(mean(1 - (pp$stdError^2 / sum(summary(mix)$varcomp[1:2,1]))), 3))
          h2[counter] <- 1e-6; h2.se[counter] <- 1e-6
          ## genetic variances
          vg[counter] <- 1e-6; vg.se[counter] <- 1e-6
        }
        pp$trait <- iTrait
        pp$fieldinstF <- "across"
        pp$entryType <- "test";  areChecks <- which(pp$genoF %in% checks)
        if(length(areChecks) > 0){pp$entryType[areChecks] <- "check"}
        predictionsList[[counter]] <- pp;
        ## cycle times
        genoCodes <- unique(mydataSub$genoCode)
        crossYears <- as.numeric(substr(as.character(as.numeric(gsub('[^[:digit:] ]', '', genoCodes) )),start=1, stop=2))
        #crossYears <- as.numeric(substr(as.character(as.numeric(gsub("GE","",genoCodes))),start=1, stop=2))
        datee <- Sys.Date()
        year.mo.day <- as.numeric(strsplit(as.character(datee),"-")[[1]])# <- as.numeric(strsplit(gsub("....-","",datee),"-")[[1]])
        your.year <- as.numeric(substr(year.mo.day[1],3,4))
        crossYears2 <- as.numeric(paste0(ifelse(crossYears < your.year, "20","19"), as.character(crossYears)))
        cl[counter] <- abs(mean(mydataSub$genoYearOrigin, na.rm=TRUE) - mean(crossYears2, na.rm=TRUE))
        cl.se[counter] <- 1e-6
        ## mean
        mu[counter] <- mean(pp$predictedValue, na.rm=TRUE)
        ## others
        field[counter] <- "across"; trt[counter] <- iTrait
        counter=counter+1
      }
    }
  }
  if(length(predictionsList) == 0){stop("There was no predictions to work with. Please look at your H2 boundaries. You may be discarding all fields.",call. = FALSE)}
  predictionsBind <- do.call(rbind, predictionsList)
  predictionsBind$analysisId <- id
  predictionsBind$id.geno <- NA
  colnames(predictionsBind) <- cgiarBase::replaceValues(Source=colnames(predictionsBind), Search=c("genoF","fieldinstF","entryType"), Replace=c("geno","fieldinst","genoType"))
  predictionsBind$pipeline <- paste(sort(unique(mydata$pipeline)),collapse=", ")
  ##########################################
  ## add year of origin, stage and geno code
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
  baseOrigin <- data.frame(geno=entries,genoYearOrigin=vals, genoYearTesting=vals4, stage=vals2, genoCode=vals3)
  predictionsBind <- merge(predictionsBind,baseOrigin, by="geno", all.x=TRUE)
  zeros <- which(predictionsBind$predictedValue < 0)
  if(length(zeros) > 0){predictionsBind[zeros,] <- NA}
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
  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	phenoDTfile$id,
    markerbooks	= NA,  markerDataFile =	NA,
    year = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  # saveRDS(db.params, file = file.path(wd,"metadata",paste0(id,".rds")))
  ## write the values used for cleaning to the modeling database
  mod <- data.frame(
    trait = trait,
    traitLb = NA,
    traitUb = NA,
    outlierCoef = NA,
    analysisId = rep(id,length(trait)),
    analysisType =rep(type,length(trait)) ,
    fixedModel = rep(fix,length(trait)),
    randomModel = rep(ranran,length(trait)),
    residualModel = ifelse(is.null(ranres),NA,rep(as.character(ranres),length(trait))),
    h2Threshold = rep(paste(c(heritLB,heritUB),collapse=" , "),length(trait))
  )
  # saveRDS(mod, file = file.path(wd,"modeling",paste0(id,".rds")))

  # write predictions
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")
  # saveRDS(predictionsBind[,predcols], file = file.path(wd,"predictions",paste0(id,".rds")))

  # write pipeline metrics
  pm <- data.frame(value=c(h2,vg,mu),  stdError=c(h2.se,vg.se,rep(1e-6,length(mu))),
                   fieldinst=c(field,field,field),  trait=c(trt,trt,trt),
                   analysisId=id, method=c(rep("cullis",length(h2)),rep("reml",length(vg)), rep("mean",length(mu)) ),
                   traitUnits=NA, parameter=c(rep("H2",length(h2)),rep("VG",length(vg)), rep("mean",length(mu)) ),
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", "),
                   stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  # saveRDS(pm, file = file.path(wd,"metrics",paste0(id,".rds")))
  # save desire file
  des <- desire(trait=trait,h2=h2, G=G
                # pathFile=file.path(wd,"desire",paste0("desire_",id,".txt"))
  )
  ##
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    # cat(paste("Your results will be available in the predictions database under such id \n"))
    # cat(paste("Your desire file will be available in the desire folder under such id \n"))
  }
  result <- list(metrics=pm, predictions=predictionsBind[,predcols], modeling=mod, metadata=db.params,
                 cleaned=NA, outliers=NA, desire=des, id=id)
  return(result)#paste("met done:",id))
}
