ammi <- function(phenoDTfile, trait, verbose=TRUE){

  id <- paste("amm",idGenerator(5,5),sep="")
  type <- "ammi"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}

  mydata <- phenoDTfile$predictions

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

  percentExplained <- traitNames <- list()
  predictionsList <- list()
  for(iTrait in trait){ # iTrait=trait[1]
    mydataSub <- mydata[which(mydata$trait == iTrait),]

    if(verbose){print(paste("Only one environment for trait:",iTrait,". Trait will be skipped."))}

    if(length(table(mydataSub$fieldinst)) > 1){ # there's enough fields
      GEmedie <-  reshape(mydataSub[,c("geno","fieldinst","predictedValue")],
                          idvar = "geno", timevar = "fieldinst", v.names = "predictedValue",
                          direction = "wide", sep="_")
      colnames(GEmedie) <- gsub("predictedValue_","",colnames(GEmedie))
      GEmedie[,-c(1)] <- apply(as.data.frame(GEmedie[,-c(1)]),2,sommer::imputev)
      rownamesGEmedie <- GEmedie[,1]
      GEmedie <- as.matrix(GEmedie[,-1])
      rownames(GEmedie) <- rownamesGEmedie

      GE <- as.data.frame(t(scale( t(scale(GEmedie, center=T,scale=F)), center=T, scale=F)))
      sum(GE^2)

      U <- svd(GE)$u
      V <- svd(GE)$v
      D <- diag(svd(GE)$d)
      Sg <- U %*% sqrt(D)
      Se <- V %*% sqrt(D)
      rownames(Sg) <- rownames(GE)#levels(dataset$Genotype)
      rownames(Se) <- colnames(GE)#levels(dataset$Environment)
      colnames(Sg) <- colnames(Se) <- paste("PC", 1:ncol(Se), sep ="")

      mySum <- numeric()
      for(PC in 1:ncol(Se)){
        Sg2 <- as.matrix(Sg[,1:PC])
        Se2 <- as.matrix(Se[,1:PC])
        GE2 <- Sg2 %*% t(Se2)
        mySum[PC] <- sum(GE2^2)
      }
      percentExplained[[iTrait]] <- mySum/mySum[length(mySum)]
      genoType <- c(rep("geno",nrow(Sg)), rep("fieldinst",nrow(Se)))
      biplot <- rbind(Sg,Se)
      biplotList <- list(); colnamesBiplot <- colnames(biplot); counter1 <- 1
      for(iCol in 1:ncol(biplot)){
        biplotList[[counter1]] <- data.frame(analysisId=id, pipeline=paste(unique(mydata$pipeline), collapse = ", "),
                                             trait=paste(iTrait,colnamesBiplot[iCol], sep="-"), geno=rownames(biplot),
                                             genoType=genoType,fieldinst="across", predictedValue=biplot[,iCol],
                                             stdError=1, rel=1
        ); counter1 <- counter1+1
      }
      predictionsList[[iTrait]] <- do.call(rbind,biplotList)
      traitNames[[iTrait]] <- unique(predictionsList[[iTrait]]$trait)

    }
  }

  if(length(predictionsList) == 0){stop("No traits with more than 1 fieldinst. Analysis stopped.", call. = FALSE)}

  predictionsBind <- do.call(rbind, predictionsList)
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
  ## write the values used for cleaning to the modeling database
  mod <- data.frame(
    trait = trait,
    traitLb = NA,
    traitUb = NA,
    outlierCoef = NA,
    analysisId = rep(id,length(trait)),
    analysisType =rep(type,length(trait)) ,
    fixedModel = NA,
    randomModel = NA,
    residualModel = NA,
    h2Threshold = NA
  )

  # write predictions
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")

  # write pipeline metrics

  pm <- data.frame(value=unlist(percentExplained),  stdError=0,
                   fieldinst="across",  trait=unlist(traitNames),
                   analysisId=id, method="ammi-svd",
                   traitUnits=NA, parameter="%expl",
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", "),
                   stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  ##
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    # cat(paste("Your results will be available in the predictions database under such id \n"))
    # cat(paste("Your desire file will be available in the desire folder under such id \n"))
  }
  result <- list(metrics=pm, predictions=predictionsBind[,predcols], modeling=mod, metadata=db.params,
                 cleaned=NA, outliers=NA, desire=NA, id=id)

}
