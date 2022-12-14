index <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    trait= NULL,
    desirev = NULL,
    scaled=TRUE,
    # wd=NULL,
    verbose=TRUE
){

  # if(is.null(wd)){wd <- getwd()}
  # md <- strsplit(wd,"/")[[1]]; md <- md[length(md)]
  # if(md != "DB"){stop("Please set your working directory to the DB folder", call. = FALSE)}

  id <- paste("idx",idGenerator(5,5),sep="")
  type <- "idx"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(length(grep("met",phenoDTfile)) == 0){stop("Index can only be calculated on results from a MET analysis using across environment predictions",call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(length(trait) != length(desirev)){stop("The number of traits and desirev values needs to be equal",call. = FALSE)}

  ############################
  # loading the dataset
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  mydata <- phenoDTfile$predictions # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  mydata <- mydata[which(mydata$trait %in% trait),]

  ############################
  ## index calculation

  wide0 <- reshape(mydata[,c("geno","trait","predictedValue")], direction = "wide", idvar = "geno",
                   timevar = "trait", v.names = "predictedValue", sep= "_")
  wide <- as.matrix(wide0[,-1]); colnames(wide) <- gsub("predictedValue_","", colnames(wide0)[-1])#unique(mydata$trait)
  wide <- as.matrix(wide[,trait]); colnames(wide) <- trait # ensure order of the user so weights also match
  wide <- apply(wide,2,sommer::imputev)

  if(scaled){
    if(verbose){
      cat(paste("scaled has been set to",scaled,"'desirev' values are expected to be the desired change in std. deviations \n"))
    }
    wide <- apply(wide,2,scale)
    wide[which(is.na(wide), arr.ind = TRUE)] <- 0
  }else{
    if(verbose){
      cat(paste("scaled has been set to",scaled,"'desirev' values are expected to be desired change in original units \n"))
    }
  }
  G <- cov(wide, use="pairwise.complete.obs")
  G[which(is.na(G), arr.ind = TRUE)] <- 0
  b <- MASS::ginv(G)%*%desirev # desired weights Ginv*d, equivalent to knowing w (economic weights)
  merit <- wide %*% b

  #########################################
  ## update databases
  newped <- data.frame(analysisId=id,pipeline=unique(mydata$pipeline)[1], trait="index",
                       geno=wide0[,1], predictedValue=merit,stdError=1e-6,rel=1e-6, fieldinst="across")
  ## zdd genoType
  checks <- mydata[which(mydata[,"genoType"] == "check"),"geno"]
  newped$genoType <- "test";  areChecks <- which(newped$geno %in% checks)
  if(length(areChecks) > 0){pp$genoType[areChecks] <- "check"}
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
  predictionsBind <- merge(newped,baseOrigin, by="geno", all.x=TRUE)

  # write predictions
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")
  # saveRDS(predictionsBind[,predcols], file = file.path(wd,"predictions",paste0(id,".rds")))

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
    stage = paste(sort(unique(mydata$stage)),collapse=", ")
  )
  # saveRDS(db.params, file = file.path(wd,"metadata",paste0(id,".rds")))
  ## write the values used for cleaning to the modeling database
  mod <- data.frame(
    trait = colnames(wide),
    traitLb = desirev,
    traitUb = b[,1],
    outlierCoef = NA,
    analysisId = id,
    analysisType = type,
    fixedModel = NA,
    randomModel = NA,
    residualModel = NA,
    h2Threshold = NA
  )
  # saveRDS(mod, file = file.path(wd,"modeling",paste0(id,".rds")))

  # write pipeline metrics
  pm <- data.frame(value=c(b[,1],desirev),
                   stdError=1e-6,
                   fieldinst="across",
                   trait=c( colnames(wide), colnames(wide)),
                   analysisId=id, method= ifelse(scaled,"desireScaled","desireOriginal"),
                   traitUnits=NA,
                   stage = paste(sort(unique(mydata$stage)),collapse=", "),
                   parameter= c(rep("weight",length(colnames(wide))),rep("desire", length(colnames(wide)))),
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", ")
  )
  # saveRDS(pm, file = file.path(wd,"metrics",paste0(id,".rds")))

  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    # cat(paste("Your results will be available in the predictions database under such id \n"))
  }
  result <- list(metrics=pm, predictions=predictionsBind[,predcols], modeling=mod, metadata=db.params,
                 cleaned=NA, outliers=NA, desire=NA, id=id)
  return(result)#paste("index done:",id))
}
