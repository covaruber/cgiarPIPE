ammi <- function(phenoDTfile, genoAmatrix=NULL, 
                 trait, verbose=TRUE, nPC=5){
  ## FUNCTION TO PERFORM AMMI USING PHENOTYPES (POSSIBLE IMPUTED BY A GRM)
  ## IT IS USED IN THE AMMI MODULE IN BANAL
  if(is.null(genoAmatrix)){ # PCA will not use this information
    analysisTypeAmat <- "nnn"
  }else{ # we need to know if is a grm or numerator
    analysisTypeAmat <- substr(genoAmatrix$id,1,3)
  }
  id <- paste( paste("amm",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "ammi"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}

  mydata <- phenoDTfile$predictions # extract phenotypes
  utraits <- unique(mydata$trait)
  traitToRemove <- character()
  for(k in 1:length(trait)){ # check which traits are really available and subset
    if (!trait[k] %in% utraits){
      if(verbose){ cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))}
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  trait <- setdiff(trait,traitToRemove)
  ## start the loop for each trait do PCA
  percentExplained <- traitNames <- list()
  predictionsList <- list()
  for(iTrait in trait){ # iTrait=trait[1]
    mydataSub <- mydata[which(mydata$trait == iTrait),] # subset for the ith trait
    if(length(table(mydataSub$fieldinst)) > 1){ # there's enough fields
      GEmedie <-  reshape(mydataSub[,c("geno","fieldinst","predictedValue")],
                          idvar = "geno", timevar = "fieldinst", v.names = "predictedValue",
                          direction = "wide", sep="_") # go to wide format, traits in columns
      colnames(GEmedie) <- gsub("predictedValue_","",colnames(GEmedie)) # remove added label by reshape
      rownamesGEmedie <- GEmedie[,1] # make it a clean matrix
      GEmedie <- as.matrix(GEmedie[,-1])
      rownames(GEmedie) <- rownamesGEmedie
      ## if user provides a relationship matrix impute missing phenotypes using that information
      if(!is.null(genoAmatrix)){
        genoFlevels <- rownamesGEmedie 
        if(analysisTypeAmat %in% "clm"){ # user provided marker data
          commonBetweenMandP <- intersect(rownames(genoAmatrix$cleaned$M),genoFlevels)
          if(length(commonBetweenMandP) < 2){
            stop("Markers could not be matched with phenotypes. Please ensure that you have used the right marker file or check the rownames of your marker matrix and ids of your phenotypes.", call. = FALSE)
          }
          M <- genoAmatrix$cleaned$M[commonBetweenMandP,]
          A <- sommer::A.mat(M)
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
        if(length(differ) > 0){ # we have to add individuals without markers or not being part of the GRM?
          if(length(differ) > 1){ # there's at least 2 inds to be added
            A2 <- diag(x=rep(mean(diag(A)),length(differ)))
          }else{ A2 <- diag(1)*mean(diag(A)) } # there's only one individual to be added
          colnames(A2) <- rownames(A2) <- differ
        }else{A2 <- matrix(0,0,0)}
        A3 <- sommer::adiag1(A,A2)
        A3[lower.tri(A3)] <- t(A3)[lower.tri(A3)] # fill the lower triangular
        colnames(A3) <- rownames(A3) <- c(colnames(A), colnames(A2))
        nearest0 <- max(c(round(nrow(GEmedie)*.05),5))
        Gu <- A3[rownames(GEmedie), rownames(GEmedie)]
        hh <- corImputation(wide=GEmedie, Gu=Gu, nearest=nearest0, roundR=FALSE)
        GEmedie <- hh$imputed

      }else{ # user didn't have relationship matrix to impute missing phenotypes, do correlation-based imputation
        nearest0 <- max(c(round(nrow(GEmedie)*.05),5))
        hh <- corImputation(wide=GEmedie, Gu=NULL, nearest=nearest0, roundR=FALSE)
        GEmedie <- hh$corImputed
      }
      ## now do the actual AMMI with means
      GE <- as.data.frame(t(scale( t(scale(GEmedie, center=T,scale=F)), center=T, scale=F)));   sum(GE^2)
      U <- svd(GE)$u
      V <- svd(GE)$v
      D <- diag(svd(GE)$d)
      Sg <- U %*% sqrt(D)
      Se <- V %*% sqrt(D)
      # only keep desired #of PCs
      nPC1 <- min(c(nPC,ncol(Sg)))
      Sg <- Sg[,1:nPC1]
      nPC2 <- min(c(nPC,ncol(Se)))
      Se <- Se[,1:nPC2]
      ##
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
      percentExplained[[iTrait]] <- mySum/mySum[length(mySum)] # calculate %explained by each PC
      genoType <- c(rep("geno",nrow(Sg)), rep("fieldinst",nrow(Se)))
      biplot <- rbind(Sg,Se) # save information for biplot
      biplotList <- list(); colnamesBiplot <- colnames(biplot); counter1 <- 1
      for(iCol in 1:ncol(biplot)){ # we need to save each PC as a different trait to stick to the predictions format
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
  ## add timePoint of origin, stage and geno code to the predictions table
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
  #########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,nanalysisType =	type,
    fieldbooks	= NA,  phenoDataFile =	phenoDTfile$id,
    markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
  )
  ## write the values used for cleaning to the modeling database
  modeling <- data.frame(
    trait = trait, traitLb = NA, traitUb = NA,
    outlierCoef = NA,  analysisId = rep(id,length(trait)),
    analysisType =rep(type,length(trait)) ,
    fixedModel = NA, randomModel = NA,residualModel = NA, h2Threshold = NA
  )
  # write pipeline metrics
  metrics <- data.frame(value=unlist(percentExplained),  stdError=0,
                   fieldinst="across",  trait=unlist(traitNames),
                   analysisId=id, method="ammi-svd",
                   traitUnits=NA, parameter="%expl",
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", "),
                   stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
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
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","mother","father","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")
  result <- list(metrics=metrics, predictions=predictionsBind[,predcols], modeling=modeling, metadata=metadata,
                 cleaned=phenoDTfile$cleaned, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst
  )
}
