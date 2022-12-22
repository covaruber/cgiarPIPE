staLMM <- function(
    phenoDTfile= NULL,
    trait=NULL, # per trait
    fixedTerm=c("1","genoF"),
    maxit=50,
    verbose=FALSE
){

  id <- paste("sta",idGenerator(5,5),sep="")
  type <- "sta"
  if(is.null(phenoDTfile)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}

  ###################################
  # loading the dataset
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  mydata <- phenoDTfile$cleaned #readRDS(file.path(wd,"files_cleaned",paste0(phenoDTfile)))
  cleaning <- phenoDTfile$outliers #readRDS(file.path(wd,"outliers",paste0(phenoDTfile)))

  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% colnames(mydata)){
      if(verbose){
        cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))
      }
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  trait <- setdiff(trait,traitToRemove)
  #####################################
  # single trial analysis
  fields <- as.character(unique(mydata$fieldinstF))
  h2 <- se <- numeric(); field <- trt <- vector()
  predictionsList <- list(); counter=1
  for(iTrait in trait){ # iTrait=trait[1]
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    for(iField in fields){ # iField = fields[1]# "2019_WS_BINA_Regional_Station_Barishal"
      if(verbose){cat(paste("Analyzing field", iField,"\n"))}
      # subset data
      mydataSub <- droplevels(mydata[which(as.character(mydata$fieldinstF) %in% iField),])
      mydataSub$trait <- mydataSub[,iTrait]
      # remove outliers
      cleaningSub <- cleaning[which(cleaning$traitName %in% iTrait),]
      out <- which(mydataSub$rowindex %in% cleaningSub$indexRow )

      if(length(out) > 0){mydataSub[out,"trait"] <- NA}
      # do analysis
      if(!is.na(var(mydataSub[,"trait"],na.rm=TRUE))){ # if there's variance
        if( var(mydataSub[,"trait"], na.rm = TRUE) > 0 ){
          checks <- mydataSub[which(mydataSub[,"entryType"] == "check"),"geno"]
          # remove repeated row:col in each trial
          gridCheck <- with(mydataSub,table(rowcoord,colcoord))
          if(nrow(gridCheck) > 1){
            badRecords <- which(gridCheck > 1, arr.ind = TRUE)
            while(nrow(badRecords) > 0){
              rowsRemove <- list()
              for(iRow in 1:nrow(badRecords)){
                rowsRemove[[iRow]] <- which((mydataSub$rowcoord == badRecords[iRow,1]) & (mydataSub$colcoord == badRecords[iRow,2]) )[1]
              }
              if(verbose){cat("Removing replicated records assigned to same row and column in the same trial \n")}
              mydataSub <- droplevels(mydataSub[-unlist(rowsRemove),])
              gridCheck <- with(mydataSub,table(rowcoord,colcoord))
              badRecords <- which(gridCheck > 1, arr.ind = TRUE)
            }
          }

          # find best formula
          mde <- cgiarFTDA::asremlFormula(fixed=as.formula(paste("trait","~ 1")),
                                          random=~ at(fieldinstF):rowcoordF + at(fieldinstF):colcoordF + at(fieldinstF):trialF + at(fieldinstF):repF + at(fieldinstF):blockF,
                                          rcov=~at(fieldinstF):id(rowcoordF):id(colcoordF),
                                          dat=droplevels(mydataSub[which(!is.na(mydataSub[,"trait"])),]),

                                          minRandomLevels=list(rowcoordF= 3, colcoordF=3, trialF=2,repF=2, blockF=4),
                                          minResidualLevels=list(rowcoordF=5, colcoordF=5),

                                          exchangeRandomEffects=list(rowcoordF="colcoordF", colcoordF="rowcoordF"),

                                          exchangeResidualEffects=list(rowcoordF="colcoordF", colcoordF="rowcoordF"),

                                          customRandomLevels=NULL, customResidualLevels=NULL,

                                          xCoordinate= "rowcoordF",yCoordinate ="colcoordF",
                                          doubleConstraintRandom=c("rowcoordF","colcoordF"), verbose=verbose)

          factorsFitted <- unlist(lapply(mde$used$fieldinstF,length))
          factorsFittedGreater <- which(factorsFitted > 0)

          newRandom <- ifelse(length(factorsFittedGreater) > 0, names(factorsFitted)[factorsFittedGreater], NA  )
          if(is.na(newRandom)){newRandom=NULL}
          #
          if((length(mde$used$fieldinstF$rowcoordF) == 0) & (length(mde$used$fieldinstF$colcoordF) == 0)){
            newSpline <- NULL
          }else{
            newSpline = as.formula("~LMMsolver::spl2D(x1 = rowcoord, x2 = colcoord, nseg = c(10, 10))")
          }

          randomTermForRanModel <- c(newRandom,"genoF")
          fixedTermForRanModel <- setdiff(fixedTerm,randomTermForRanModel)
          fix <- paste("trait ~",paste(fixedTermForRanModel, collapse = " + "))
          #
          ranran <- paste(unique(c(randomTermForRanModel, newRandom)), collapse=" + ")
          ranran <- paste("~",ranran)
          # do we have replication in the fixed factors to be fitted?
          repFixedTerm <- apply(data.frame(setdiff(fixedTerm,"1")),1,function(x){length(which(table(mydataSub[,x]) > 1))})
          repFixedTermGreater <- which(repFixedTerm > 0)

          # at least one condition met
          if( (length(factorsFittedGreater) > 0) | (length(repFixedTermGreater) > 0) ){

            mixRandom <- try(
              LMMsolver::LMMsolve(fixed =as.formula(fix),
                                  random = as.formula(ranran),
                                  spline = newSpline, #trace = TRUE,
                                  # ginverse = list(genoF=Ainv),
                                  data = mydataSub, maxit = maxit),
              silent = TRUE
            );

            # only keep variance components that were greater than zero
            if(!inherits(mixRandom,"try-error") & ((length(factorsFittedGreater) > 0) | (length(repFixedTermGreater) > 0)) ){ # if random model runs well try the fixed model
              sm <- summary(mixRandom, which = "variances")
              newRanran <- setdiff((sm[,1])[which(sm[,2] >0.05)],c("residual","genoF"))
              ranran <- paste("~",paste(c(newRanran), collapse=" + "))
              if(ranran=="~ "){ranran=NULL}else{ranran <- as.formula(ranran)}
              rownames(sm) <- NULL
              if(verbose){
                print(sm)
                cat(paste(iTrait," ~",paste(fixedTerm, collapse = " + ")," \n"))
                cat(paste(ranran,"\n"))
              }
              fix <- paste("trait ~",paste(fixedTerm, collapse = " + "))
              mixFixed <- try(
                LMMsolver::LMMsolve(fixed =as.formula(fix),
                                    random = ranran,
                                    spline = newSpline, #trace = TRUE,
                                    data = mydataSub, maxit = maxit),
                silent = TRUE
              )
              if(!inherits(mixFixed,"try-error") ){ # if fixed model was not singular save all results

                predictedValue <- mixFixed$coefficients$genoF + mixFixed$coefficients$`(Intercept)`
                stdError <- (sqrt(diag(as.matrix(solve(mixFixed$C)))))[1:length(predictedValue)]
                genoF <- gsub("genoF_","", names(predictedValue))
                pp <- data.frame(genoF,predictedValue,stdError)
                pp$trait <- iTrait
                pp$fieldinstF <- iField
                pp$entryType <- "test";  areChecks <- which(pp$genoF %in% checks)
                if(length(areChecks) > 0){pp$entryType[areChecks] <- "check"}
                pp$genoYearTesting <- unique(mydataSub$year)[1]
                ## heritabilities
                ss = mixRandom$VarDf#summary(mixRandom, which = "variances")
                rownames(ss) <- ss$VarComp
                vg <- ss["genoF",2]; vr <- ss["residual",2]
                h2[counter] <-  vg / (vg+vr) ; se[counter] <- 0
                ## reliability
                pp$rel <- abs(1 - (pp$stdError^2)/(vg))
                predictionsList[[counter]] <- pp;
                field[counter] <- iField; trt[counter] <- iTrait
                counter=counter+1

              } # end of if fixed model run well
            }else{ # if there was singularities we just take means and assigna h2 of zero
              if(verbose){cat(paste("No design to fit, aggregating and assuming h2 = 0 \n"))}
              pp <- aggregate(trait ~ genoF, FUN=mean, data=mydataSub)
              colnames(pp)[2] <- "predictedValue"
              pp$stdError <- 1
              pp$status <- "averaged"
              pp$rel <- 1e-6
              pp$trait <- iTrait
              pp$fieldinstF <- iField
              pp$entryType <- "test";  areChecks <- which(pp$genoF %in% checks)
              if(length(areChecks) > 0){pp$entryType[areChecks] <- "check"}
              pp$genoYearTesting <- unique(mydataSub$year)[1]
              predictionsList[[counter]] <- pp;
              h2[counter] <- 1e-6; se[counter] <- 1e-6 # lm(rr$predictedValue~pp$predictedValue)$coefficients[2]
              field[counter] <- iField; trt[counter] <- iTrait
              counter=counter+1
            } # end of is mixed model run well

          }else{

            if(verbose){cat(paste("No design to fit, aggregating and assuming h2 = 0 \n"))}
            pp <- aggregate(trait ~ genoF, FUN=mean, data=mydataSub)
            colnames(pp)[2] <- "predictedValue"
            pp$stdError <- 1
            pp$status <- "averaged"
            pp$rel <- 1e-6
            pp$trait <- iTrait
            pp$fieldinstF <- iField
            pp$entryType <- "test";  areChecks <- which(pp$genoF %in% checks)
            if(length(areChecks) > 0){pp$entryType[areChecks] <- "check"}
            pp$genoYearTesting <- unique(mydataSub$year)[1]
            predictionsList[[counter]] <- pp;
            h2[counter] <- 1e-6; se[counter] <- 1e-6 # lm(rr$predictedValue~pp$predictedValue)$coefficients[2]
            field[counter] <- iField; trt[counter] <- iTrait
            counter=counter+1

          }

        }
      }
    }
  }
  predictionsBind <- do.call(rbind, predictionsList)
  predictionsBind$analysisId <- id
  colnames(predictionsBind) <- cgiarBase::replaceValues(Source=colnames(predictionsBind), Search=c("genoF","fieldinstF","entryType"), Replace=c("geno","fieldinst","genoType"))
  predictionsBind$pipeline <- paste(sort(unique(mydata$pipeline)),collapse=", ")
  datee <- Sys.Date()
  year.mo.day <- as.numeric(strsplit(as.character(datee),"-")[[1]])# <- as.numeric(strsplit(gsub("....-","",datee),"-")[[1]])
  your.year <- year.mo.day[1]
  predictionsBind <- predictionsBind[which(predictionsBind$genoYearTesting <= your.year),] # remove typos of individuals from future years
  ##########################################
  ## add year of origin and stage and geno code
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
  baseOrigin <- data.frame(geno=entries, genoCode=vals3, genoYearOrigin=vals, stage=vals2)
  predictionsBind <- merge(predictionsBind,baseOrigin, by="geno", all.x=TRUE)
  ##########################################
  ## update databases
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= id,analysisType =	type,fieldbooks	= NA,
    phenoDataFile =	phenoDTfile$id, markerbooks	= NA,  markerDataFile =	NA,
    year = NA,  season =	NA,  location =	NA,country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  ## write the values used for cleaning to the modeling database
  mod <- data.frame(
    trait = trait, traitLb = NA,traitUb = NA,outlierCoef = NA,
    analysisId = id,analysisType = type,fixedModel = fix,
    randomModel = ifelse(is.null(ranran),NA,setdiff(as.character(ranran),"~")),residualModel = "units",h2Threshold = NA
  )
  # write predictions
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")

  # write pipeline metrics
  pm <- data.frame(value=h2,stdError=se, fieldinst=field,trait=trt,
                   analysisId=id, method="ratio",traitUnits=NA,
                   parameter="H2",
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", "),
                   stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    # cat(paste("Your results will be available in the predictions database under such id \n"))
  }
  result <- list(metrics=pm, predictions=predictionsBind[,predcols], modeling=mod, metadata=db.params,
                 cleaned=NA, outliers=NA, desire=NA, id=id)
  return(result)#paste("sta done:", id))
}
