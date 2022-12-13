gwas <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    markerDTfile= NULL,
    trait= NULL, # per trait
    fieldinst="across",
    verbose=FALSE
    # wd=NULL
){

  # if(is.null(wd)){wd <- getwd()}
  # md <- strsplit(wd,"/")[[1]]; md <- md[length(md)]
  # if(md != "DB"){stop("Please set your working directory to the DB folder", call. = FALSE)}

  baseId <- idGenerator(5,5)
  idmes <- paste("mes",baseId,sep="")
  idhits <- paste("gwa",baseId,sep="")
  typemes <- "mes"
  typehits <- "gwa"
  if(is.null(phenoDTfile)){stop("Please provide the predictions", call. = FALSE)}
  if(is.null(markerDTfile)){stop("Please provide the markers", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  # if(length(trait) > 1){
  #   stop(paste0(" Only one trait can be used for optimal contribution. We suggest using an index"), call. = FALSE)
  # }

  ############################
  # loading the dataset
  mydata <- phenoDTfile$predictions # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))

  # ava.files <- dir(file.path(wd,"files_cleaned"))
  # if(!paste0(markerDTfile) %in% ava.files){stop("markerDTfile is not present in the files_cleaned folder. Please check the name of your file and its location",call. = FALSE)}
  myrel <- markerDTfile$cleaned # readRDS(file.path(wd,"files_cleaned",paste0(markerDTfile)))

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

  mydata <- mydata[which(mydata$trait %in% trait),]
  mydata <- mydata[which(mydata$fieldinst %in% fieldinst),] # make sure only across
  if(nrow(mydata) == 0){stop("Please check the trait and fieldinst selected since there's no phenotypic data for that combination",call. = "FALSE")}
  # make sure you have same phenotypes and genotypes

  common <- intersect(rownames(myrel$M), mydata[,"geno"])
  M <- myrel$M[common,]
  mydata <- mydata[which(mydata[,"geno"] %in% common),]

  ############################
  ## ocs analysis
  predictionsBindMesList <- predictionsBindHitsList <- list()
  for(iTrait in trait){ #  iTrait = trait[3]
    print(iTrait)
    mydataSub <- mydata[which(mydata$trait == iTrait),]
    n <- nrow(mydataSub) # to be used for degrees of freedom
    k <- 1 # to be used for degrees of freedom (number of levels in fixed effects)
    MMT <-tcrossprod(M) ## MM' = additive relationship matrix
    MMTinv<-solve(MMT + diag(1e-6, ncol(MMT), ncol(MMT))) ## inverse of MM'
    MTMMTinv<-t(M)%*%MMTinv # M' %*% (M'M)-
    mixGBLUP <- sommer::mmer(predictedValue~1,
                     random=~sommer::vsr(geno, Gu=MMT), rcov=~units, nIters=30,
                     verbose = FALSE,
                     data=mydataSub)
    myorder <- names(mixGBLUP$U$`u:geno`$predictedValue)
    M <- M[myorder,]
    MMT <- MMT[myorder,myorder]
    MMTinv <- MMTinv[myorder,myorder]
    MTMMTinv <- MTMMTinv[,myorder]
    a.from.g <-MTMMTinv%*%matrix(mixGBLUP$U$`u:geno`$predictedValue,ncol=1)
    var.g <- kronecker(MMT,mixGBLUP$sigma$`u:geno`) - mixGBLUP$PevU$`u:geno`$predictedValue
    var.a.from.g <- t(M)%*%MMTinv%*% (var.g) %*% t(MMTinv)%*%M
    se.a.from.g <- sqrt(diag(var.a.from.g) + min(diag(var.a.from.g))+ 1e-7)
    t.stat.from.g <- a.from.g/se.a.from.g # t-statistic
    pvalGBLUP <- dt(t.stat.from.g,df=n-k-1) #
    minusLog10pvalGBLUP <- -log10(pvalGBLUP)

    rel <- 1 - (se.a.from.g^2)/as.numeric(mixGBLUP$sigma$`u:geno`)
    # plot(minusLog10pvalGBLUP)
    # write predictions

    predictionsBindMesList[[iTrait]] <- data.frame(analysisId=idmes, pipeline= paste(sort(unique(mydataSub$pipeline)),collapse=", "),
                                  trait=trait, genoCode=1:nrow(a.from.g), geno=rownames(a.from.g),
                                  genoType="markerEffect", genoYearOrigin=paste(sort(unique(mydataSub$genoYearOrigin)),collapse=", "),
                                  genoYearTesting=paste(sort(unique(mydataSub$genoYearTesting)),collapse=", "),
                                  fieldinst=fieldinst, predictedValue=a.from.g[,1], stdError=se.a.from.g, rel=rel, stage=paste(sort(unique(mydataSub$stage)),collapse=", ")
    )
    predictionsBindHitsList[[iTrait]] <- data.frame(analysisId=idhits, pipeline= paste(sort(unique(mydataSub$pipeline)),collapse=", "),
                                                trait=trait, genoCode=1:nrow(pvalGBLUP), geno=rownames(a.from.g),
                                                genoType="gwasHits", genoYearOrigin=paste(sort(unique(mydataSub$genoYearOrigin)),collapse=", "),
                                                genoYearTesting=paste(sort(unique(mydataSub$genoYearTesting)),collapse=", "),
                                                fieldinst=fieldinst, predictedValue=minusLog10pvalGBLUP[,1], stdError=se.a.from.g, rel=rel, stage=paste(sort(unique(mydataSub$stage)),collapse=", ")
    )
  }

  predictionsBindMes <- do.call(rbind, predictionsBindMesList)
  predictionsBindHits <- do.call(rbind, predictionsBindHitsList)
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")
  # saveRDS(predictionsBindMes[,predcols], file = file.path(wd,"predictions",paste0(idmes,".rds")))
  # saveRDS(predictionsBindHits[,predcols], file = file.path(wd,"predictions",paste0(idhits,".rds")))

  #########################################
  ## update databases
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= idmes,
    analysisType =	typemes,
    fieldbooks	= NA,
    phenoDataFile =	NA,
    markerbooks	= NA,  markerDataFile =	markerDTfile,
    year = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(predictionsBindMes$stage)),collapse=", ")
  )
  # saveRDS(db.params, file = file.path(wd,"metadata",paste0(idmes,".rds")))
  # db.params$analysisId <- idhits
  # db.params$analysisType <- typehits
  # saveRDS(db.params, file = file.path(wd,"metadata",paste0(idhits,".rds")))
  ## write the values used for cleaning to the modeling database
  mod <- data.frame(
    trait = trait,
    traitLb = NA,
    traitUb = NA,
    outlierCoef = NA,
    analysisId = idmes,
    analysisType = typemes,
    fixedModel = NA,
    randomModel = NA,
    residualModel = NA,
    h2Threshold = NA
  )
  # saveRDS(mod, file = file.path(wd,"modeling",paste0(idmes,".rds")))
  # mod$analysisId <- idhits
  # mod$analysisType <- typehits
  # saveRDS(mod, file = file.path(wd,"modeling",paste0(idhits,".rds")))

  # write pipeline metrics
  if(verbose){
    cat(paste("Your analysis id is:",idmes,"\n"))
    cat(paste("Your analysis id is:",idhits,"\n"))
    cat(paste("Your results will be available in the predictions database under such id \n"))
  }
  result <- list(metrics=NA, predictions=predictionsBindMes[,predcols], modeling=mod, metadata=db.params,
                 cleaned=NA, outliers=NA, desire=NA, id=idmes)
  return(result)#paste("gs/gwas done:",idmes))
}
