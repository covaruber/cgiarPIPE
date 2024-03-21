staLMM <- function(
    phenoDTfile= NULL,
    trait=NULL, # per trait
    traitFamily=NULL,
    fixedTerm=c("1","genoF"),
    genoUnit = c("genoF"),
    maxit=50,
    returnFixedGeno=TRUE,
    verbose=TRUE
){
  ## THIS FUNCTION CALCULATES A SINGLE TRIAL ANALYSIS FOR MULTIPLE FIELDS AND TRAITS
  ## IS USED IN THE BANAL APP UNDER THE GENETIC EVALUATION MODULES
  id <- paste( paste("sta",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")

  fixedTerm <- unique(c(fixedTerm,"genoF"))
  type <- "sta"
  if(is.null(phenoDTfile)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(is.null(traitFamily)){traitFamily <- rep("gaussian(link = 'identity')", length(trait))}
  if(length(traitFamily) != length(trait)){stop("Trait distributions should have the same length than traits to be analyzed.", call. = FALSE)}
  names(traitFamily) <- trait
  ###################################
  # loading the dataset
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  mydata <- phenoDTfile$cleaned #readRDS(file.path(wd,"files_cleaned",paste0(phenoDTfile)))

  cleaning <- phenoDTfile$outliers #readRDS(file.path(wd,"outliers",paste0(phenoDTfile)))
  for(iName in c("geno","mother","father")){
    # also the ones that are already factors
    mydata[,paste0(iName,"F")] <- as.factor(mydata[,iName])
  }

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
  fixedFormulaForFixedModel = NULL
  randomFormulaForFixedModel = NULL
  fields <- as.character(unique(mydata$fieldinstF))

  if(length(fields) == nrow(mydata)){
    stop("The number of fieldinst levels is equal to the number of records/rows in the dataset.
          This means that probably you didn't select the right columns to define the fieldinst column.
          Please match again your raw file checking carefully the columns defining the fieldinst.", call.=FALSE )
  }

  h2 <- se <- cv <- r2 <- r2se <- numeric(); field <- trt <- vector(); diagnostic <- character()
  predictionsList <- list(); counter=1
  library(LMMsolver)
  for(iTrait in trait){ # iTrait=trait[1]
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    mydata[,paste(iTrait,"residual",sep="-")] <- NA
    for(iField in fields){ # iField = fields[1]# "2019_WS_BINA_Regional_Station_Barishal"
      if(verbose){cat(paste("Analyzing field", iField,"\n"))}
      # subset data
      mydataSub <- droplevels(mydata[which(as.character(mydata$fieldinstF) %in% iField),])
      mydataSub$trait <- mydataSub[,iTrait]
      # remove outliers
      if(!is.null(nrow(cleaning))){ # if there's outliers
        cleaningSub <- cleaning[which(cleaning$traitName %in% iTrait),]
        out <- which(mydataSub$rowindex %in% cleaningSub$indexRow )
      }else{out <- numeric()}
      # check the genetic units
      nLevelsGenounit <- apply(data.frame(genoUnit),1,function(x){length(table(mydataSub[,x])) }); names(nLevelsGenounit) <- genoUnit
      genoUnitTraitField <- names(nLevelsGenounit)[which(nLevelsGenounit > 1)]

      if(length(out) > 0){mydataSub[out,"trait"] <- NA}
      # do analysis
      if(!is.na(var(mydataSub[,"trait"],na.rm=TRUE))){ # if there's variance
        if( var(mydataSub[,"trait"], na.rm = TRUE) > 0 ){
          gridCheck <- with(mydataSub,table(rowcoord,colcoord))
          if(nrow(gridCheck) > 1){
            # try to fix assuming they jusy have the rows and cols by replicate
            badRecords <- which(gridCheck > 1, arr.ind = TRUE)
            if(nrow(badRecords) > 0){
              mydataSub <- cgiarBase::fixCoords( mydataSub=mydataSub, rowcoord="rowcoord", colcoord="colcoord", rep="rep" )
            }
            # make the check once again in case that didn't help
            gridCheck <- with(mydataSub,table(rowcoord,colcoord))
            badRecords <- which(gridCheck > 1, arr.ind = TRUE)
            if(nrow(badRecords) > 0){
              if(verbose){cat("Replicated records assigned to same row and column in the same trial. Ignoring row and column information for this trial. \n")}
              mydataSub$rowcoord <- NA
              mydataSub$rowcoordF <- NA
              mydataSub$colcoord <- NA
              mydataSub$colcoordF <- NA
            }
          }

          # find best experimental design formula
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

          newRandom <- ifelse(length(factorsFittedGreater) > 0, TRUE, FALSE  )
          if(newRandom){newRandom <- names(factorsFitted)[factorsFittedGreater]}else{newRandom<-NULL}
          #
          if((length(mde$used$fieldinstF$rowcoordF) == 0) & (length(mde$used$fieldinstF$colcoordF) == 0)){
            newSpline <- NULL
          }else{
            newSpline = as.formula("~spl2D(x1 = rowcoord, x2 = colcoord, nseg = c(10, 10))")
          }

          for(iGenoUnit in genoUnitTraitField){ # iGenoUnit <- genoUnitTraitField[1]
            myGeneticUnit <-  iGenoUnit
            randomTermForRanModel <- c(newRandom,genoUnit)
            fixedTermForRanModel <- setdiff(fixedTerm,randomTermForRanModel)
            fixedFormulaForRanModel <- paste("trait ~",paste(fixedTermForRanModel, collapse = " + "))
            #
            ranran <- paste(c(myGeneticUnit, unique(intersect(randomTermForRanModel,newRandom)) ), collapse = " + ")
            randomFormulaForRanModel <- paste("~",ranran)
            # do we have replication in the fixed factors to be fitted?
            toTabulate <- setdiff(fixedTermForRanModel,"1")
            if(length(toTabulate) > 0){
              repFixedTerm <- apply(data.frame(toTabulate),1,function(x){length(which(table(mydataSub[,x]) > 1))})
            }else{
              repFixedTerm <- character()
            }
            repFixedTermGreater <- which(repFixedTerm > 0)

            # impute fixed effect columns if they are numeric
            toImpute <- unlist(lapply(mydataSub[,setdiff(fixedTermForRanModel,"1"),drop=FALSE], class))
            keepToImpute <- which(toImpute %in% "numeric")
            if(length(keepToImpute) > 0){
              toImpute <- names(toImpute[keepToImpute])
              for(iImpute in toImpute){mydataSub[, iImpute] <- sommer::imputev(mydataSub[, iImpute])}
            }

            # at least one condition met: replicated random terms, or replicated fixed terms
            if( (length(factorsFittedGreater) > 0) | (length(repFixedTermGreater) > 0) ){

              mixRandom <- try(
                LMMsolver::LMMsolve(fixed =as.formula(fixedFormulaForRanModel),
                                    random = as.formula(randomFormulaForRanModel),
                                    spline = newSpline, #trace = TRUE,
                                    family = eval(parse(text = traitFamily[iTrait])),
                                    data = mydataSub, maxit = maxit),
                silent = TRUE
              );  # mixRandom$VarDf
              # if random model run
              # only keep variance components that were greater than zero
              if(!inherits(mixRandom,"try-error") & ((length(factorsFittedGreater) > 0) | (length(repFixedTermGreater) > 0)) ){ # if random model runs well try the fixed model
                ## save residuals
                whereResidualGoes <- mydataSub[which(!is.na(mydataSub$trait)),"rowindex"]
                mydata[whereResidualGoes,paste(iTrait,"residual",sep="-")] <- mixRandom$residuals[,1]
                ##
                sm <- summary(mixRandom, which = "variances")
                newRanran <- setdiff( (sm[,1])[which(sm[,2] >0.05)] , c("residual",genoUnitTraitField,"s(rowcoord, colcoord)"))

                ranran <- paste("~",paste(c(newRanran), collapse=" + "))
                if(ranran=="~ "){randomFormulaForFixedModel=NULL}else{randomFormulaForFixedModel <- as.formula(ranran)}
                rownames(sm) <- NULL
                if(verbose){
                  print(sm)
                  cat(paste(iTrait," ~",paste(iGenoUnit, collapse = " + ")," \n"))
                  cat(paste(randomFormulaForFixedModel,"\n"))
                }
                otherFixed <- setdiff(fixedTerm,genoUnitTraitField)
                fixedFormulaForFixedModel <- paste("trait ~",paste(c(iGenoUnit,otherFixed), collapse = " + "))
                mixFixed <- try(
                  LMMsolver::LMMsolve(fixed =as.formula(fixedFormulaForFixedModel),
                                      random = randomFormulaForFixedModel,
                                      spline = newSpline, #trace = TRUE,
                                      family = eval(parse(text = traitFamily[iTrait])),
                                      data = mydataSub, maxit = maxit),
                  silent = TRUE
                ) # mixFixed$VarDf
                ##################################################################
                ## if model fails try the simplest model, no random and no spatial
                if(inherits(mixFixed,"try-error") ){
                  mixFixed <- try(
                    LMMsolver::LMMsolve(fixed =as.formula(fixedFormulaForFixedModel),
                                        family = eval(parse(text = traitFamily[iTrait])),
                                        data = mydataSub, maxit = maxit),
                    silent = TRUE
                  )
                  if( inherits(mixFixed,"try-error") ){
                    if(verbose){cat(paste("Fixed effects models failed. Returning deregressed BLUPs \n"))}
                    predictedValue <-  mixRandom$coefMME[mixRandom$ndxCoefficients[[iGenoUnit]]] +  mixRandom$coefMME[mixRandom$ndxCoefficients$`(Intercept)`]
                    dims <- mixRandom$EDdf
                    start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) + 1 # we add the one when is random
                    pev <- as.matrix(solve(mixRandom$C))[start:(start+length(predictedValue)-1),start:(start+length(predictedValue)-1)]
                    stdError <- stdErrorRandom <- (sqrt(diag(pev)))

                    badSEs <- which( stdError < (sd(predictedValue, na.rm = TRUE)/100) )
                    if(length(badSEs) > 0){stdError[badSEs] <- sd(predictedValue, na.rm = TRUE)}

                    genoF <- gsub(paste0(iGenoUnit,"_"),"", names(mixRandom$ndxCoefficients[[iGenoUnit]]))
                    pp <- data.frame(genoF,predictedValue,stdError)
                    pp$trait <- iTrait
                    pp$fieldinstF <- iField
                    pp$entryType <- apply(data.frame(pp$genoF),1,function(x){found <-which(mydataSub$geno %in% x); x2 <- ifelse(length(found) > 0, paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "##"),"unknown_type"); return(x2)})
                    # print(head(mydataSub))
                    if(is.null(mydataSub$timePoint)){
                      pp$genoYearTesting <- unique(mydataSub$year)[1]
                    }else{
                      pp$genoYearTesting <- unique(mydataSub$timePoint)[1]
                    }
                    ## heritabilities
                    ss = mixRandom$VarDf#summary(mixRandom, which = "variances")
                    rownames(ss) <- ss$VarComp
                    vg <- ss[iGenoUnit,2]; vr <- ss["residual",2]
                    h2[counter] <-  vg / (vg+vr) ; se[counter] <- 0
                    cv[counter] <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
                    ## reliability
                    A <- diag(nrow(pev))
                    G <- A*vg # G matrix
                    R2 = (G - pev)/G
                    pp$rel <- diag(R2)
                    pp$predictedValue <- pp$predictedValue/pp$rel
                    r2[counter] <- mean(pp$rel)
                    r2se[counter] <- sd(pp$rel, na.rm = TRUE)/sqrt(length(pp$rel))

                    # pp$rel <- abs(1 - (stdErrorRandom^2)/(vg))
                    badRels <- which(pp$rel > 1); if(length(badRels) > 0){pp$rel[badRels] <- 0.9999}
                    predictionsList[[counter]] <- pp;
                    field[counter] <- iField; trt[counter] <- iTrait
                    counter=counter+1
                  }
                }
                ##################################################################
                if(!inherits(mixFixed,"try-error") ){ # if fixed model was not singular save all results

                  if(returnFixedGeno){ # user wants fixed effect predictions for genotype
                    shouldBeOne <- which(mixFixed$ndxCoefficients[[iGenoUnit]] == 0)
                    if(length(shouldBeOne) > 0){mixFixed$ndxCoefficients[[iGenoUnit]][shouldBeOne] = 1}
                    predictedValue <- mixFixed$coefMME[mixFixed$ndxCoefficients[[iGenoUnit]]] + mixFixed$coefMME[mixFixed$ndxCoefficients$`(Intercept)`]
                    if(length(shouldBeOne) > 0){predictedValue[1] = mixFixed$coefMME[mixFixed$ndxCoefficients$`(Intercept)`]}
                    dims <- mixFixed$EDdf
                    start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) # we don't add a one because we need the intercept
                    pev <- as.matrix(solve(mixFixed$C))[start:(start+length(predictedValue)-1),start:(start+length(predictedValue)-1)]
                    stdError <- sqrt(diag(pev))
                    stdError[1] <- mean(stdError[-1])
                    badSEs <- which( stdError < (sd(predictedValue, na.rm = TRUE)/100) )
                    if(length(badSEs) > 0){stdError[badSEs] <- sd(predictedValue, na.rm = TRUE)}
                    # just for reliability calculation
                    dims <- mixRandom$EDdf
                    start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) + 1 # we add the one when is random
                    pev <- as.matrix(solve(mixRandom$C))[start:(start+length(predictedValue)-1),start:(start+length(predictedValue)-1)]
                    # stdErrorRandom <- (sqrt(diag(pevRandom)))
                  }else{ # user wants random effect predictions for genotype
                    predictedValue <- mixRandom$coefMME[mixRandom$ndxCoefficients[[iGenoUnit]]] +  mixRandom$coefMME[mixRandom$ndxCoefficients$`(Intercept)`]
                    dims <- mixRandom$EDdf
                    start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) + 1 # we add the one when is random
                    pev <- as.matrix(solve(mixRandom$C))[start:(start+length(predictedValue)-1),start:(start+length(predictedValue)-1)]
                    stdError <- stdErrorRandom <- (sqrt(diag(pev)))
                    # move to std deviation if model is wrong and stdError is close to zero
                    badSEs <- which( stdError < (sd(predictedValue, na.rm = TRUE)/100) )
                    if(length(badSEs) > 0){stdError[badSEs] <- sd(predictedValue, na.rm = TRUE)}
                  }
                  genoF <- gsub(paste0(iGenoUnit,"_"),"", names(mixRandom$ndxCoefficients[[iGenoUnit]]))
                  pp <- data.frame(genoF,predictedValue,stdError)
                  pp$trait <- iTrait
                  pp$fieldinstF <- iField
                  pp$entryType <- apply(data.frame(pp$genoF),1,function(x){found <-which(mydataSub$geno %in% x); x2 <- ifelse(length(found) > 0, paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#"),"unknown_type"); return(x2)})
                  # print(head(mydataSub))
                  if(is.null(mydataSub$timePoint)){
                    pp$genoYearTesting <- unique(mydataSub$year)[1]
                  }else{
                    pp$genoYearTesting <- unique(mydataSub$timePoint)[1]
                  }
                  ## heritabilities
                  ss = mixRandom$VarDf#summary(mixRandom, which = "variances")
                  rownames(ss) <- ss$VarComp
                  vg <- ss[iGenoUnit,2]; vr <- ss["residual",2]
                  h2[counter] <-  vg / (vg+vr) ; se[counter] <- 0
                  cv[counter] <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100

                  ## reliability
                  A <- diag(nrow(pev))
                  G <- A*vg # G matrix
                  R2 = (G - pev)/G
                  pp$rel <- diag(R2)
                  r2[counter] <- mean(pp$rel)
                  r2se[counter] <- sd(pp$rel, na.rm = TRUE)/sqrt(length(pp$rel))

                  badRels <- which(pp$rel > 1); if(length(badRels) > 0){pp$rel[badRels] <- 0.9999}
                  predictionsList[[counter]] <- pp;
                  field[counter] <- iField; trt[counter] <- iTrait
                  counter=counter+1

                } # end of if fixed model run well
              }else{ # if there was singularities we just take means and assigna h2 of zero
                if(verbose){cat(paste("No design to fit or singularities encountered in the random model, aggregating and assuming h2 = 0 \n"))}
                pp <- aggregate(as.formula(paste("trait ~", iGenoUnit)), FUN=mean, data=mydataSub)
                colnames(pp)[2] <- "predictedValue"
                pp$stdError <- sd(pp$predictedValue) #1
                # pp$status <- "averaged"
                pp$rel <- 1e-6
                pp$trait <- iTrait
                pp$fieldinstF <- iField
                pp$entryType <- apply(data.frame(pp$genoF),1,function(x){found <-which(mydataSub$geno %in% x); x2 <- ifelse(length(found) > 0, paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#"),"unknown_type"); return(x2)})
                if(is.null(mydataSub$timePoint)){
                  pp$genoYearTesting <- unique(mydataSub$year)[1]
                }else{
                  pp$genoYearTesting <- unique(mydataSub$timePoint)[1]
                }
                predictionsList[[counter]] <- pp;
                h2[counter] <- 1e-6; se[counter] <- 1e-6
                cv[counter] <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
                field[counter] <- iField; trt[counter] <- iTrait
                r2[counter] <- 1e-6
                r2se[counter] <- 1e-6
                counter=counter+1
              } # end of is mixed model run well

            }else{

              if(verbose){cat(paste("No design to fit, aggregating and assuming h2 = 0 \n"))}
              pp <- aggregate(as.formula(paste("trait ~", iGenoUnit)), FUN=mean, data=mydataSub)
              colnames(pp)[2] <- "predictedValue"
              pp$stdError <- sd(pp$predictedValue)  # 1
              # pp$status <- "averaged"
              pp$rel <- 1e-6
              pp$trait <- iTrait
              pp$fieldinstF <- iField
              pp$entryType <- apply(data.frame(pp$genoF),1,function(x){found <-which(mydataSub$geno %in% x); x2 <- ifelse(length(found) > 0, paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#"),"unknown_type"); return(x2)})
              if(is.null(mydataSub$timePoint)){
                pp$genoYearTesting <- unique(mydataSub$year)[1]
              }else{
                pp$genoYearTesting <- unique(mydataSub$timePoint)[1]
              }
              predictionsList[[counter]] <- pp;
              h2[counter] <- 1e-6; se[counter] <- 1e-6 # lm(rr$predictedValue~pp$predictedValue)$coefficients[2]
              cv[counter] <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
              field[counter] <- iField; trt[counter] <- iTrait
              r2[counter] <- 1e-6
              r2se[counter] <- 1e-6
              counter=counter+1

            }

          }


        }else{

          for(iGenoUnit in genoUnitTraitField){

            if(verbose){
              cat(paste("No design to fit, aggregating for predicted values, std. errors assumed equal to std. deviation of the trial. In addition assuming h2 = 0 for the trial \n"))
            }
            pp <- aggregate(as.formula(paste("trait ~", iGenoUnit)), FUN=mean, data=mydataSub)
            colnames(pp)[2] <- "predictedValue"
            pp$stdError <- sd(pp$predictedValue)  # 1
            pp$trait <- iTrait
            pp$fieldinstF <- iField
            pp$entryType <- apply(data.frame(pp$genoF),1,function(x){found <-which(mydataSub$geno %in% x); x2 <- ifelse(length(found) > 0, paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#"),"unknown_type"); return(x2)})
            if(is.null(mydataSub$timePoint)){
              pp$genoYearTesting <- unique(mydataSub$year)[1]
            }else{
              pp$genoYearTesting <- unique(mydataSub$timePoint)[1]
            }
            pp$rel <- 1e-6
            predictionsList[[counter]] <- pp;
            h2[counter] <- 1e-6; se[counter] <- 1e-6 #
            cv[counter] <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
            field[counter] <- iField; trt[counter] <- iTrait
            r2[counter] <- 1e-6
            r2se[counter] <- 1e-6
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
  timePoint.mo.day <- as.numeric(strsplit(as.character(datee),"-")[[1]])# <- as.numeric(strsplit(gsub("....-","",datee),"-")[[1]])
  your.timePoint <- timePoint.mo.day[1]
  predictionsBind <- predictionsBind[which(predictionsBind$genoYearTesting <= your.timePoint),] # remove typos of individuals from future years
  ##########################################
  ## add timePoint of origin and stage and geno code
  entries <- unique(predictionsBind[,"geno"])
  vals <- apply(data.frame(entries),1,function(x){
    identified <- which(mydata$geno %in% x)
    out <- ifelse(length(identified) > 0,as.character((sort(mydata[identified,"genoYearOrigin"], decreasing = FALSE))[1]),"unknown")
    return(out)
  })
  vals2 <- apply(data.frame(entries),1,function(x){
    identified <- which(mydata$geno %in% x)
    out <- ifelse(length(identified) > 0,as.character((sort(mydata[identified,"stage"], decreasing = FALSE))[1]),"unknown")
    return(out)
  })
  vals3 <- apply(data.frame(entries),1,function(x){
    identified <- which(mydata$geno %in% x)
    out <- ifelse(length(identified) > 0,as.character((sort(mydata[identified,"genoCode"], decreasing = FALSE))[1]),"unknown")
    return(out)
  })
  vals4 <- apply(data.frame(entries),1,function(x){
    identified <- which(mydata$geno %in% x)
    out <- ifelse(length(identified) > 0,as.character((sort(mydata[identified,"mother"], decreasing = FALSE))[1]),"unknown")
    return(out)
  })
  vals5 <- apply(data.frame(entries),1,function(x){
    identified <- which(mydata$geno %in% x)
    out <- ifelse(length(identified) > 0,as.character((sort(mydata[identified,"father"], decreasing = FALSE))[1]),"unknown")
    return(out)
  })
  baseOrigin <- data.frame(geno=entries, mother=vals4, father=vals5, genoCode=vals3, genoYearOrigin=vals, stage=vals2)
  predictionsBind <- merge(predictionsBind,baseOrigin, by="geno", all.x=TRUE)
  if(length(which(is.na(predictionsBind$stage))) < length(predictionsBind$stage)){predictionsBind$stage <- sommer::imputev(predictionsBind$stage)}
  if(length(which(is.na(predictionsBind$genoCode))) < length(predictionsBind$genoCode)){predictionsBind$genoCode <- sommer::imputev(predictionsBind$genoCode)}
  if(length(which(is.na(predictionsBind$genoYearOrigin))) < length(predictionsBind$genoYearOrigin)){predictionsBind$genoYearOrigin <- sommer::imputev(predictionsBind$genoYearOrigin)}

  ##########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,analysisType =	type,fieldbooks	= NA,
    phenoDataFile =	phenoDTfile$id, markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  ## write the values used for cleaning to the modeling database
  modeling <- data.frame(
    trait = trait, traitLb = NA,traitUb = NA,outlierCoef = NA,
    analysisId = id,analysisType = type,fixedModel = ifelse(is.null(fixedFormulaForFixedModel),NA,setdiff(as.character(fixedFormulaForFixedModel),"~")) ,
    randomModel = ifelse(is.null(randomFormulaForFixedModel),NA,setdiff(as.character(randomFormulaForFixedModel),"~")),
    residualModel = "units",h2Threshold = NA
  )
  # write predictions
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","mother","father","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")

  # write pipeline metrics
  metrics <- data.frame(value=c(h2,cv, r2),stdError=c(se,rep(NA,length(cv)),r2se), fieldinst=c(field,field,field),trait=c(trt,trt,trt),
                   analysisId=id, method=c(rep("vg/(vg+ve)",length(h2)), rep("sd/mu", length(cv)), rep("(G-PEV)/G", length(r2)) ),
                   traitUnits=c(rep("unitless",length(h2)), rep("percentage",length(cv)), rep("unitless",length(r2)) ),
                   parameter=c(rep("plotH2",length(h2)), rep("CV",length(cv)), rep("meanR2",length(r2)) ),
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", "),
                   stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }

  if(is.null(phenoDTfile$metadataFieldinst)){
    metadataFieldinst=NA
  }else{
    metadataFieldinst=phenoDTfile$metadataFieldinst
  }
  result <- list(metrics=metrics, predictions=predictionsBind[,predcols], modeling=modeling, metadata=metadata,
                 cleaned=mydata, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst
  )
  return(result)
}
