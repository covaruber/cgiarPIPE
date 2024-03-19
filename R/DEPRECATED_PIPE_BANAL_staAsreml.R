sta <- function(
    phenoDTfile= NULL,
    trait=NULL, # per trait
    fixedTerm=c("1","genoF"),
    workspace="360mb",
    pworkspace="360mb",#wd=NULL,
    verbose=FALSE
){

  id <- paste( paste("sta",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "sta"
  if(is.null(phenoDTfile)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}

  ###################################
  # loading the dataset
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  mydata <- phenoDTfile$cleaned #readRDS(file.path(wd,"files_cleaned",paste0(phenoDTfile)))
  cleaning <- phenoDTfile$outliers #readRDS(file.path(wd,"outliers",paste0(phenoDTfile)))
  # parameters <- read.csv(file.path(wd,"parameters.csv"))
  # modeling <- read.csv(file.path(wd,"modeling.csv"))
  # predictions <- read.csv(file.path(wd,"predictions.csv"))
  # pipeline_metrics <- read.csv(file.path(wd,"pipeline_metrics.csv"))

  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% colnames(mydata)){
      if(verbose){
        cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))
      }
      traitToRemove <- c(traitToRemove,trait[k])
      # stop(paste0("'", trait[k], "' is not a column in the given dataset"))
    }
  }
  trait <- setdiff(trait,traitToRemove)
  #####################################
  # single trial analysis
  tableFields <- table(mydata$fieldinstF)
  fields <- as.character(unique(names(tableFields[which(tableFields > 1)])))
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

          # Ai <- PED$ginv; attr(Ai, "INVERSE") <- TRUE
          if((length(mde$used$fieldinstF$rowcoordF) == 0) & (length(mde$used$fieldinstF$colcoordF) == 0)){
            newdat <- mydataSub
            glist <<- NULL
            funny <- NULL
          }else{
            Z <- sommer::spl2Db(x.coord=mydataSub$rowcoord,y.coord=mydataSub$colcoord, at.var=NULL,at.levels=NULL, nsegments = c(10,10),
                        degree = c(3,3), penaltyord = c(2,2),nestorder = c(1,1),
                        minbound=NULL, maxbound=NULL, method="Lee", what="base")
            newdat <- cbind(mydataSub,Z)
            glist <<- list(field=1:ncol(Z) + ncol(mydataSub)); #names(glist) <- "field"
            funny <- paste0("grp(","field",")")
          }
          fix <- paste("trait ~",paste(fixedTerm, collapse = " + "))
          ranran0 <-"~genoF" # simple formula, can be more complex

          ranran <- paste(c(ranran0, mde$random, funny), collapse=" + ")
          ranres <- "~dsum(~units | fieldinstF)"
          mixRandom <- try(
            asreml::asreml(as.formula(fix),
                   random= as.formula(ranran),
                   group=glist,
                   residual=as.formula(ranres),
                   na.action = na.method(x="include",y="include"),
                   trace=FALSE,
                   data=newdat, maxiter=50),
            silent = TRUE
          );  # mixRandom <- update(mixRandom, maxiter=5)


          # only keep variance components that were greater than zero
          if(!inherits(mixRandom,"try-error") ){ # if random model runs well try the fixed model
            sm <- summary(mixRandom)$varcomp;
            fix <- paste("trait ~",paste(fixedTerm, collapse = " + "))
            newRanran <- setdiff(rownames(sm)[which(sm[,1] >0.05)],c("fieldinstF!R","genoF"))
            newRanran <- gsub("):","' ):",gsub(", ",", '", newRanran))
            ranran <- paste("~",paste(c(newRanran), collapse=" + "))
            if(ranran=="~ "){ranran=~NULL; glist<<- NULL}else{ranran <- as.formula(ranran)}
            rownames(sm) <- NULL
            if(verbose){
              print(sm)
              cat(paste(iTrait," ~",paste(fixedTerm, collapse = " + ")," \n"))
              cat(paste(ranran,"\n"))
            }
            mixFixed <- try(
              asreml::asreml(as.formula(fix),
                     random= ranran,
                     group=glist,
                     residual=as.formula(ranres),
                     na.action = na.method(x="include",y="include"),
                     trace=FALSE,
                     data=newdat, maxiter=50, workspace=workspace),
              silent = TRUE
            )
            if(!inherits(mixFixed,"try-error") ){ # if fixed model was not singular save all results

              pp <- asreml::predict(mixFixed, classify = "genoF", pworkspace=pworkspace,data=newdat, trace=FALSE, maxit=1)$pvals#, aliased=TRUE)
              colnames(pp) <- cgiarBase::replaceValues(Source=colnames(pp), Search=c("predicted.value","std.error"), Replace=c("predictedValue","stdError"))
              pp$trait <- iTrait
              pp$fieldinstF <- iField
              pp$entryType <- "test";  areChecks <- which(pp$genoF %in% checks)
              if(length(areChecks) > 0){pp$entryType[areChecks] <- "check"}
              pp$genoYearTesting <- unique(newdat$year)[1]
              # pp$genoYearOrigin<- unique(newdat$genoYearOrigin)[1]
              ## heritabilities
              vr <- paste0("V",grep("fieldinstF!R",rownames(summary(mixRandom)$varcomp)))
              vg <- paste0("V",grep("genoF",rownames(summary(mixRandom)$varcomp)))
              h2Pred <- vpredict(mixRandom,  as.formula(paste0("~(",vg,")/(",vg,"+",vr,")")) )
              ## reliability
              pp$rel <- 1 - (pp$stdError^2)/(2*summary(mixRandom)$varcomp[as.numeric(gsub("V","",vg)),"component"])
              predictionsList[[counter]] <- pp;
              h2[counter] <- h2Pred$Estimate[1]; se[counter] <- h2Pred$SE[1] # lm(rr$predictedValue~pp$predictedValue)$coefficients[2]
              field[counter] <- iField; trt[counter] <- iTrait
              counter=counter+1

            } # end of if fixed model run well
          }else{ # if there was singularities we just take means and assigna h2 of zero
            if(verbose){cat(paste("No design to fit, aggregating and assuming h2 = 0 \n"))}
            pp <- aggregate(trait ~ genoF, FUN=mean, data=newdat)
            colnames(pp)[2] <- "predictedValue"
            pp$stdError <- 1
            pp$status <- "averaged"
            pp$rel <- 1e-6
            pp$trait <- iTrait
            pp$fieldinstF <- iField
            pp$entryType <- "test";  areChecks <- which(pp$genoF %in% checks)
            if(length(areChecks) > 0){pp$entryType[areChecks] <- "check"}
            pp$genoYearTesting <- unique(newdat$year)[1]
            # pp$genoYearOrigin<- unique(newdat$genoYearOrigin)[1]
            predictionsList[[counter]] <- pp;
            h2[counter] <- 1e-6; se[counter] <- 1e-6 # lm(rr$predictedValue~pp$predictedValue)$coefficients[2]
            field[counter] <- iField; trt[counter] <- iTrait
            counter=counter+1
          } # end of is mixed model run well

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
  # saveRDS(db.params, file = file.path(wd,"metadata",paste0(id,".rds")))
  ## write the values used for cleaning to the modeling database
  mod <- data.frame(
    trait = trait, traitLb = NA,traitUb = NA,outlierCoef = NA,
    analysisId = id,analysisType = type,fixedModel = fix,
    randomModel = setdiff(as.character(ranran),"~"),residualModel = ranres,h2Threshold = NA
  )
  # saveRDS(mod, file = file.path(wd,"modeling",paste0(id,".rds")))

  # write predictions
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")
  # saveRDS(predictionsBind[,predcols], file = file.path(wd,"predictions",paste0(id,".rds")))

  # write pipeline metrics
  pm <- data.frame(value=h2,stdError=se, fieldinst=field,trait=trt,
                   analysisId=id, method="ratio",traitUnits=NA,
                   parameter="H2",
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", "),
                   stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  # saveRDS(pm, file = file.path(wd,"metrics",paste0(id,".rds")))

  if(verbose){
  cat(paste("Your analysis id is:",id,"\n"))
  # cat(paste("Your results will be available in the predictions database under such id \n"))
  }
  result <- list(metrics=pm, predictions=predictionsBind[,predcols], modeling=mod, metadata=db.params,
                 cleaned=NA, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal)
  return(result)#paste("sta done:", id))
}
