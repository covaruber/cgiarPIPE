tpeDef <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    trait= NULL, # per trait
    fieldsToInclude=NULL,
    weightByAccuracy=TRUE,
    verbose=TRUE
){
  ## THIS FUNCTION CALCULATES THE TRADE OFF FOR THE DIFFERENT ENVIRONMENTS WHEN USING AN ACROSS ENV ESTIMATE
  ## IS USED IN THE BANAL APP UNDER THE WEATHER ANALYSIS MODULES
  id <- paste( paste("tpe",cgiarPIPE::idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "tpe"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(length(trait) > 1){stop("Only one trait at the time can be processed", call. = FALSE)}
  
  rgToDegree <- function(rg){
    # 180 minus formula to convert rg=c(-1, 1) into degree = c(0, 180)
    y <- 180 - ( (179*rg) + 1 ) 
    return(y)
  }
  
  trade <- function(degrees){
    # formula to convert degrees into radians
    radians = degrees * ( pi / 180.0 )
    # find the adjacent side of a triangle using adjacent = cos(theta) * hypothenuse
    # we assume always hypothenuse = 1
    a = cos(radians/2) # over 2 because is the middle between the 2 locations, forms the triangle
    return(a)
  }
  
  response <- function(tradeMatrix, pick){
    tradeMatrixSub <- tradeMatrix[pick,pick, drop=FALSE]
    notSelected <- setdiff(1:ncol(tradeMatrix),pick)
    result <- vector(mode="numeric", length = ncol(tradeMatrixSub))
    for(i in 1:ncol(tradeMatrixSub)){ # i=1
      col <- tradeMatrixSub[-i,i]
      row <- tradeMatrixSub[i,-i]
      result[i] <- ifelse(length(c(col,row)) > 0, mean(c(col,row)), 1)
    }
    fullResult <- vector(mode="numeric", length = ncol(tradeMatrix))
    fullResult[pick] <- result
    fullResult[notSelected] <- 0
    average <- mean(result)
    result2 <- list(tradeOffs=fullResult,averageGain=average,selectedFields=pick)
    return(result2)
  }
  
  ############################
  # loading the dataset
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  mydata <- phenoDTfile$metrics #readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  ############################
  # subset to a specific trait and parameter
  mydata = mydata[which(mydata$trait %in% trait ),]
  mydata = mydata[which(mydata$parameter %in%  "corG"),]
  # add field columns as character and numbers
  mydata4 <- apply(data.frame(mydata[,"fieldinst"]),1,function(x){strsplit(x,"###")[[1]]})
  mydata4 <- data.frame(t(mydata4))
  mydata4$Freq <- round(mydata[,"value"],3)
  mydata5 <- mydata4
  mydata5$X1 <- mydata4$X2
  mydata5$X2 <- mydata4$X1
  mydata6 <- rbind(mydata4,mydata5)
  mydata6 <- unique(mydata6)
  mydata6$Var1 <- as.numeric(as.factor(mydata6$X1))
  mydata6$Var2 <- as.numeric(as.factor(mydata6$X2))
  mydata6$mytext <- paste(mydata6$X1,mydata6$X2, sep="###")
  # make it a matrix
  A <- matrix(NA, nrow=max(mydata6$Var1), ncol = max(mydata6$Var2))
  A[as.matrix(mydata6[,c("Var1","Var2")])] = mydata6[,"Freq"]
  A[lower.tri(A)] <- t(A)[lower.tri(A)] # fill the lower triangular
  rownames(A) <- colnames(A) <- levels(as.factor(mydata6$X1))
  nas <- which(is.na(A), arr.ind = TRUE)
  if(nrow(nas) > 0){A[nas]=stats::median(A[lower.tri(A)], na.rm = TRUE)}
  # calculate rgToDegree and tradeMatrix
  tradeMatrix <- trade(rgToDegree(A));tradeMatrix
  # calculate response for selected fields
  pick <- which(cgiarBase::cleanChar(colnames(A)) %in% cgiarBase::cleanChar(fieldsToInclude) )
  result <- response(tradeMatrix,pick)
  # weight by reliability
  myR2 <- phenoDTfile$metrics 
  myR2 = myR2[which(myR2$trait %in% trait ),]
  myR2 = myR2[which(myR2$parameter %in%  "meanR2"),]
  rownames(myR2) <- myR2$fieldinst
  r2ToWeight <- myR2[fieldsToInclude,"value"]
  
  if(weightByAccuracy){
    result$tradeOffs <- result$tradeOffs * sqrt(r2ToWeight)
  }
  
  # save the new metrics of this analysis
  metrics <- data.frame(value=c(result$tradeOffs, result$averageGain, result$selectedFields),  stdError=NA,
                   fieldinst=c(cgiarBase::cleanChar(colnames(A)), "across", cgiarBase::cleanChar(colnames(A))[pick]),  
                   trait=trait,
                   analysisId=id, method=c(rep("cosine", length(result$tradeOffs)),
                                           rep("mean", length(result$averageGain)),
                                           rep("none", length(result$selectedFields))
                   ),
                   traitUnits=NA, parameter=c(rep("tradeOff", length(result$tradeOffs)),
                                              rep("averageGain", length(result$averageGain)),
                                              rep("selectedField", length(result$selectedFields))
                                              ),
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", "),
                   stage = paste(sort(unique(mydata$stage)),collapse=", ")
  )
  
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
    stage = paste(sort(unique(mydata$stage)),collapse=", ")
  )
  ## write the values used for cleaning to the modeling database
  modeling <- data.frame(
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
  ##
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    # cat(paste("Your results will be available in the predictions database under such id \n"))
    # cat(paste("Your desire file will be available in the desire folder under such id \n"))
  }
  if(is.null(phenoDTfile$metadataFieldinst)){
    metadataFieldinst=NA
  }else{
    metadataFieldinst=phenoDTfile$metadataFieldinst
  }
  result <- list(metrics=metrics, predictions=NA, modeling=modeling, metadata=metadata,
                 cleaned=NA, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst
  )
  return(result)
}
