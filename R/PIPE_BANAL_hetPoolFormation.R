hetf <- function(
    phenoDTfile=NULL,
    trait =NULL, 
    nPoolSearch=3,
    method="kmeans",
    verbose=FALSE
){
  ## THIS FUNCTION PERFORMS AN HETEROTIC POOL FORMATION
  ## IS USED IN THE BANAL APP UNDER STRUCTURE MODULES
  ## we still have to figure out how to add male and female information to pure predictions
  
  id <- paste( paste("hpf",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "hpf"
  
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  if (is.null(trait)) stop("No trait specified.")
  
  predictionsList <- list()
  for(iTrait in trait){ # iTrait <- trait[1]
    
    mydata=phenoDTfile$predictions
    mydata <- mydata[which(mydata$trait %in% iTrait),]
    
    mydata2 <- mydata
    mydata2$mother <- mydata$father
    mydata2$father <- mydata$mother
    mydata <- rbind(mydata,mydata2)
    mydata$cross <- paste(mydata$mother, mydata$father, sep="-")
    mydata <- mydata[which(!duplicated(mydata$cross)),]
    mydata$motherN <- as.numeric(as.factor(mydata$mother))
    mydata$fatherN <- as.numeric(as.factor(mydata$father))
    
    A <- matrix(NA, nrow=max(mydata$motherN), ncol = max(mydata$fatherN)) 
    A[as.matrix(mydata[,c("motherN","fatherN")])] = mydata[,"predictedValue"]
    A[lower.tri(A)] <- t(A)[lower.tri(A)] # fill the lower triangular
    rownames(A) <- colnames(A) <- levels(as.factor(mydata$mother))
    A <- corImputation(wide=A, Gu=NULL, nearest=10, roundR=FALSE)$imputed
    
    if(method == "kmeans"){
      cl <- kmeans(x=A, centers=nPoolSearch, nstart = 25)
    }else if(method == "hclust"){
      cl <- hclust(as.dist(A))
      nPerPool <- round(length(cl$order)/nPoolSearch)
      start <- 1
      cl$cluster <- rep(NA, length(cl$labels)); names(cl$cluster) <- cl$labels
      for(i in 1:nPoolSearch){
        cl$cluster[which(cl$order < (nPerPool*i) & cl$order > start )] <- i
        start <- (nPerPool*i)
      }
      cl$cluster[which(is.na(cl$cluster))] <- i
    }else{
      stop("Method not enabled yet", call. = FALSE)
    }
    
    predictionsList[[iTrait]] <- data.frame(analysisId=id, pipeline="none", trait=iTrait, genoCode=NA, 
                              geno=names(cl$cluster), mother=NA, father=NA, genoType=NA, genoYearOrigin=1,
                              genoYearTesting=1,fieldinst="across", 
                              predictedValue=cl$cluster, stdError=NA, rel=NA, stage=NA)
    
    
  }
  predictions <- do.call(rbind,predictionsList)
  
  
  #########################################
  metrics <- data.frame(value=nPoolSearch,  stdError=1e-3,
                   fieldinst="across",  trait=trait,
                   analysisId=id, method=method,
                   traitUnits=NA, parameter="pools",
                   pipeline="unknown",
                   stage = "unknown"
  )
  ########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	NA,
    markerbooks	= NA,  markerDataFile =	phenoDTfile$id,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = NA
  )
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the files_cleaned folder under such id \n"))
  }
  result <- list(metrics=metrics, predictions=predictions, modeling=NA, metadata=metadata,
                 cleaned=NA, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=NA)
  return(result)
}
