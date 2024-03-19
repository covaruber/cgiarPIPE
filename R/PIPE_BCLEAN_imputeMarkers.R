imputeM <- function(
    markerDTfile = NULL, nearest=10, propNaForG=.9,  # proportion completeness
    verbose=TRUE, method="median"
){
  ## THIS FUNCTION IMPUTES A MARKER MATRIX
  ## IS USED IN THE BCLEAN APP UNDER THE QA/QC MODULES
  id <- paste( paste("clm",idGenerator(5,5),sep=""), markerDTfile$idOriginal, sep = "_")
  type <- "clm"
  
  ############################
  # loading the dataset
  if(method == "median"){
    propNAs <- apply(markerDTfile$cleaned$M,2,function(x){1-(length(which(is.na(x)))/length(x))})
    M <- markerDTfile$cleaned$M
    M <- M[,which(propNAs >= propNaForG)]
    Mi <- apply(M,2, sommer::imputev)
    rownames(Mi) <- rownames(M)
    markerDTfile$cleaned$M <- Mi
    markerDTfile$cleaned$ref.alleles <- markerDTfile$cleaned$ref.alleles[,which(propNAs > propNaForG)]
    
  }else if(method == "correlation"){
    propNAs <- apply(markerDTfile$cleaned$M,2,function(x){1-(length(which(is.na(x)))/length(x))})
    altValue <- quantile(propNAs,.8)
    X <- markerDTfile$cleaned$M[,which(propNAs >= min(c(propNaForG,altValue)) )]
    X <- apply(X, 2, sommer::imputev)
    Gu <- cor(t(X))
    Mi <- sommer::corImputation(markerDTfile$cleaned$M, nearest = nearest, Gu=Gu, roundR = TRUE)
    markerDTfile$cleaned$M <- Mi$imputed
    markerDTfile$cleaned$ref.alleles <- markerDTfile$cleaned$ref.alleles[,which(propNAs > propNaForG)]
  }else{
    stop("Method not implemented yet.",call. = FALSE)
  }
  
  #########################################
  markerDTfile$id <- gsub("rbm","clm", markerDTfile$id)
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  return(markerDTfile)
}
