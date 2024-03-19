pcaPheno <- function(
    phenoDTfile=NULL,
    fieldinst= NULL, # fieldinst= "-"
    traits =NULL, # traits <- unique(mydata$trait);traits
    wideBy="trait", # column to be used for reshaping
    idvar="geno",
    nPC=3,
    verbose=FALSE
){
  ## THIS FUNCTION PERFORMS A PCA ON A PHENOTYPIC FILE IN WIDE FORMAT (FEATURES)
  ## IS USED IN THE BANAL APP UNDER THE STRUCTURE MODULE
  id <- paste( paste("pca",idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "pca"

  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  if (is.null(traits)) stop("No traits specified.")
  if (is.null(fieldinst)) stop("No fieldinst specified.")
  if(length(fieldinst) > 1 & length(traits) >1){
    stop("We can only handle multiple traits for one field or one trait for multiple fields.")
  }
  if(("trait"  %in% wideBy) & length(traits) < 2){
    stop("To cluster by traits we need more than one trait.")
  }
  if((fieldinst %in% wideBy) & length(fieldinst) < 2){
    stop("To cluster by fieldinst we need more than one fieldinst.")
  }

  mydata=phenoDTfile$predictions
  mydata <- mydata[which(mydata$fieldinst %in% fieldinst),]
  mydata <- mydata[which(mydata$trait %in% traits),]
  
  mydataWide <- reshape(mydata[,c(idvar,wideBy,"predictedValue")], direction = "wide", idvar = idvar,
                        timevar = wideBy, v.names = "predictedValue", sep= "_")
  colnames(mydataWide) <- gsub("predictedValue_","",colnames(mydataWide))
  # impute missing values
  mydataWide <- mydataWide[,which(colnames(mydataWide) != "NA")]
  mydataWide <- mydataWide[which(!is.na(mydataWide[,1])),]
  
  levelsOfWideBy <- unique(na.omit(mydata[,wideBy]))
  # rownamesOfWideBy <- unique(na.omit(mydata[,idvar]))

  
  ############################
  # do PCA
  mark6 <- mydataWide[,levelsOfWideBy]; rownames(mark6) <- mydataWide[,1] # mydataWide[,geno]
  mark6 <- as.matrix(mark6)
  # mark6 <- corImputation(mark6, nearest = 10)$imputed
  mark6 <- apply(mark6, 2, sommer::imputev)
  #PCA
  mark6.centered <- scale(mark6,center=T,scale=F)
  mark.svd <- svd(mark6.centered)
  percentExplained <- (mark.svd$d^2/sum(mark.svd$d^2))*100
  # plot(mark.svd$d^2/sum(mark.svd$d^2),ylab="Fraction Variation Explained")
  ix <- which(mark.svd$d>1e-9)
  #Project markers onto the PCs (M language to V language, although changes V grid to M grid
  # ... since the linear transformation M is applied to V)
  pca <- mark6.centered %*% mark.svd$v[,ix]
  colnames(pca) <- paste0("PC",1:ncol(pca))

  if(nPC > ncol(pca)){nPC <- ncol(pca)}

  pcaList <- list(); colnamespca <- colnames(pca); counter1 <- 1
  for(iCol in 1:nPC){
    pcaList[[counter1]] <- data.frame(analysisId=id, pipeline="unknown",
                                      trait=paste("marker",colnamespca[iCol], sep="-"), geno=rownames(pca),
                                      mother="unknown",father="unknown",
                                      genoType="test",fieldinst="across", predictedValue=pca[,iCol],
                                      stdError=1, rel=1, genoYearOrigin=1, genoYearTesting=1, stage="unknown",
                                      genoCode= 1:nrow(pca)
    ); counter1 <- counter1+1
  }
  predictionsBind <- do.call(rbind,pcaList)
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","mother","father","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")
  #########################################
  metrics <- data.frame(value=percentExplained,  stdError=1e-3,
                   fieldinst="across",  trait="marker",
                   analysisId=id, method="pca-svd",
                   traitUnits=NA, parameter="%expl",
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
  result <- list(metrics=metrics, predictions=predictionsBind[,predcols], modeling=NA, metadata=metadata,
                 cleaned=NA, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=NA)
  return(result)
}
