pca <- function(
    markerDTfile= NULL,
    minMAF = 0.01,
    nPC=3,
    verbose=FALSE
){
  ## THIS FUNCTION CALCULATES A PCA USING A MARKER MATRIX
  ## IS USED IN THE BANAL APP UNDER THE STRUCTURE MODULES
  id <- paste( paste("pca",idGenerator(5,5),sep=""), markerDTfile$idOriginal, sep = "_")
  type <- "pca"

  ############################
  # loading the dataset
  if (is.null(markerDTfile)) stop("No input marker data file specified.")
  M0 <- markerDTfile$cleaned #readRDS(file.path(wd,"files_cleaned",paste0(markerDTfile)))

  ############################
  # do PCA
  mark5 <- M0$M + 1
  MAF <- apply(mark5,2,function(x){y <- mean(x)/2; return(min(y,1-y))})
  mark6 <- mark5[,which(MAF>minMAF)]
  gid <- rownames(mark6)
  n.gid <- length(gid)

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
    markerbooks	= NA,  markerDataFile =	markerDTfile$id,
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
                 cleaned=NA, outliers=NA, desire=NA, id=id, idOriginal=markerDTfile$idOriginal,
                 metadataFieldinst=NA)
  return(result)
}
