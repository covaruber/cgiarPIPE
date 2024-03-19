rgg <- function(
    phenoDTfile= NULL,
    trait=NULL, # per trait
    fixedTerm="genoYearOrigin",
    deregressWeight=1,
    deregress=FALSE,
    partition=FALSE,
    yearsToUse=NULL,
    wd=NULL, verbose=TRUE
){
  ## THIS FUNCTION CALCULATES THE REALIZED GENETIC GAIN FOR SOME TRAITS
  ## IS USED IN THE BANAL APP UNDER THE METRICS MODULES
  id <- paste( paste("rgg",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "rgg"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}

  ############################
  # loading the dataset
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  parameters <- phenoDTfile$metadata # readRDS(file.path(wd,"metadata",paste0(phenoDTfile)))
  mydata <- phenoDTfile$predictions #readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  if(!is.null(yearsToUse)){
    yearsToUse <- as.character(yearsToUse)
    mydata <- mydata[which(as.character(mydata$genoYearOrigin) %in% yearsToUse),]
  } else{
    yearsToUse <- as.character(unique(mydata$genoYearOrigin))
  }
  if(nrow(mydata) == 0){stop("No data to work with with the specified parameters. You may want to check the yearsToUse parameter.",call. = FALSE)}
  if(length(unique(na.omit(mydata[,fixedTerm]))) <= 1){stop("Only one year of data. Realized genetic gain analysis cannot proceed.", call. = FALSE)}
  ############################
  ## gg analysis
  gg <- ggp <- inter <- gg.y1 <- gg.yn <- segp <- seb1 <- seb0 <- r2 <- pv <- ntrial <- ntrial.se <- numeric();
  field <- trt <- vector()
  counter=1
  for(iTrait in trait){ # iTrait="GY"
    if(verbose){
      cat(paste("Analyzing trait", iTrait,"\n"))
    }
    # subset data
    mydataSub <- droplevels(mydata[which(mydata$trait == iTrait),])
    mydataSub$fieldinstF <- as.factor(mydataSub$fieldinst)
    mydataSub$genoF <- as.factor(mydataSub$geno)
    mydataSub$predictedValue.d <- mydataSub$predictedValue/mydataSub$rel
    mydataSub[,fixedTerm] <- as.numeric(mydataSub[,fixedTerm])
    # do analysis
    if(!is.na(var(mydataSub[,"predictedValue"],na.rm=TRUE))){ # if there's variance
      if( var(mydataSub[,"predictedValue"], na.rm = TRUE) > 0 ){
        checks <- mydataSub[which(mydataSub[,"genoType"] == "check"),"geno"]
        ranran <- "~NULL"
        if(deregress){
          mydataSub$predictedValue <- mydataSub$predictedValue.d
          if(verbose){
            print("Deregressing predicted values using the reliability. Assuming you are providing BLUPs.")
          }
        }else{
          if(verbose){
            print("Using predicted values directly. Assuming you are providing BLUEs.")
          }
        }
        fix <- paste("predictedValue ~",paste(fixedTerm, collapse=" + "))
        # Ai <- PED$ginv; attr(Ai, "INVERSE") <- TRUE
        ranres <- "~units"#"~dsum(~units | fieldinstF)"
        mydataSub=mydataSub[with(mydataSub, order(fieldinstF)), ]
        mydataSub$w <- 1/(mydataSub$stdError)

        hh<-split(mydataSub,mydataSub[,fixedTerm])
        hh <- lapply(hh,function(x){
          outlier <- boxplot.stats(x=x[, "predictedValue"],coef=1.5 )$out
          bad <- which(x$predictedValue %in% outlier)
          if(length(bad) >0){out <- x[-bad,]}else{out<-x}
          return(out)
        })

        if(partition){
          p1 <- p2 <- p3 <- p4 <- p5 <- p6 <- p7 <- p8 <- numeric();cc <- 1
          for(u in 1:(length(hh))){
            for(w in 1:u){
              if(u != w){
                mydataSub2 <- do.call(rbind,hh[c(u,w)])
                mix <- lm(as.formula(fix), data=mydataSub2)
                sm <- summary(mix)
                p1[cc] <- sm$coefficients[2,1]*ifelse(deregress,deregressWeight,1)
                baseline <- mix$coefficients[1] + ( (mix$coefficients[2]*ifelse(deregress,deregressWeight,1))*min(as.numeric(mydataSub[which(mydataSub$trait == iTrait),fixedTerm]) ))
                p2[cc] <- round(( (mix$coefficients[2]*ifelse(deregress,deregressWeight,1)) /baseline) * 100,3)
                p3[cc] <- round((sm$coefficients[2,2]/baseline) * 100,3)
                p4[cc] <- sm$coefficients[1,1]
                p5[cc] <- sm$coefficients[2,2]
                p6[cc] <- sm$coefficients[1,2]
                p7[cc] <- sm$r.squared
                p8[cc] <- 1 - pf(sm$fstatistic[1], df1=sm$fstatistic[2], df2=sm$fstatistic[3])
                cc <- cc+1
              }
            }
          }
          gg[counter] <- mean(p1, na.rm=TRUE); ggp[counter] <- mean(p2, na.rm=TRUE); segp[counter] <- mean(p3, na.rm=TRUE)
          inter[counter] <- mean(p4, na.rm=TRUE); seb1[counter] <- mean(p5, na.rm=TRUE); seb0[counter] <- mean(p6, na.rm=TRUE)
          r2[counter] <- mean(p7, na.rm=TRUE); pv[counter] <- mean(p8, na.rm=TRUE)
        }else{
          mix <- lm(as.formula(fix), data=mydataSub)
          sm <- summary(mix)
          gg[counter] <- sm$coefficients[2,1]*ifelse(deregress,deregressWeight,1)
          baseline <- mix$coefficients[1] + ( (mix$coefficients[2]*ifelse(deregress,deregressWeight,1))*min(as.numeric(mydataSub[which(mydataSub$trait == iTrait),fixedTerm]) , na.rm=TRUE ))
          ggp[counter] <- round(( (mix$coefficients[2]*ifelse(deregress,deregressWeight,1)) /baseline) * 100,3)
          segp[counter] <- round((sm$coefficients[2,2]/baseline) * 100,3)
          inter[counter] <- sm$coefficients[1,1]
          seb1[counter] <- sm$coefficients[2,2]
          seb0[counter] <- sm$coefficients[1,2]
          r2[counter] <- sm$r.squared
          pv[counter] <- 1 - pf(sm$fstatistic[1], df1=sm$fstatistic[2], df2=sm$fstatistic[3])
        }

        # mix <- asreml(as.formula(fix),
        #                     random= as.formula(ranran),
        #                     # group=prov$glist,
        #                     residual=as.formula(ranres),
        #                     # weights = w,
        #                     workspace="900mb",
        #                     trace=FALSE,
        #                     # family = asreml::asr_gaussian(dispersion = 1),
        #                     na.action = na.method(x="exclude",y="include"),
        #                     data=mydataSub, maxiter=30)
        # summary(mix)$varcomp
        # gg[counter] <- mix$coefficients$fixed[1,1]
        # inter[counter] <- mix$coefficients$fixed[2,1]
        # seb1[counter] <- mix$vcoeff$fixed[1]
        # seb0[counter] <- mix$vcoeff$fixed[2]
        ntrial[counter] <- NA#mean(apply(dtwfs/dtwfs,2,sum, na.rm=TRUE))
        ntrial.se[counter] <- NA#sd(apply(dtwfs/dtwfs,2,sum, na.rm=TRUE))
        #
        field[counter] <- "across"; trt[counter] <- iTrait
        gg.y1[counter] <- sort(unique(mydataSub[,fixedTerm]), decreasing = FALSE)[1]
        gg.yn[counter] <- sort(unique(mydataSub[,fixedTerm]), decreasing = TRUE)[1]
        counter=counter+1
      }
    }
  }
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
    stage = paste(sort(unique(mydataSub$stage)),collapse=", ")
  )
  ## write the values used for cleaning to the modeling database
  modeling <- data.frame(
    trait = trait,
    traitLb = NA,
    traitUb = NA,
    outlierCoef = NA,
    analysisId = id,
    analysisType = type,
    fixedModel = fix,
    randomModel = ranran,
    residualModel = ranres,
    h2Threshold = NA
  )
  # write pipeline metrics
  metrics <- data.frame(value=c(gg,inter,ggp, gg.y1,r2,pv, ntrial),
                   stdError=c(seb1,seb0,segp,gg.yn,rep(1e-6,length(r2)),rep(1e-6,length(pv)), ntrial.se),
                   fieldinst=c(field,field,field,field,field,field,field),
                   trait=c(trt,trt,trt,trt,trt,trt,trt),
                   analysisId=id, method= ifelse(deregress,"blup+dereg","mackay"),traitUnits=NA,
                   stage = paste(sort(unique(mydataSub$stage)),collapse=", "),
                   parameter= c(rep("ggSlope",length(gg)), rep("ggInter",length(inter)), rep("gg%",length(gg)), rep("ggYear",length(gg.y1)), rep("r2",length(r2)), rep("pVal",length(pv)), rep("nTrialPerYear", length(ntrial)) ),
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", "),
                   yearOriginAnalyzed = rep(paste(sort(unique(yearsToUse)),collapse=", "),7)
                   )
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
                 metadataFieldinst=metadataFieldinst)
  return(result)#paste("rgg done:",id))
}
