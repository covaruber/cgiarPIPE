scm <- function( # single cross matrix function
    markerDTfile= NULL,#wd=NULL,
    markerDTfile2=NULL,
    hybridBatch=1000,
    separator=":",
    additive=TRUE,
    verbose=TRUE
){
  ## THIS FUNCTION CALCULATES A SINGLE-CROSS-HYBRID MARKER MATRIX USING MARKER FROM PARENTAL LINES
  ## IS USED IN THE BANAL APP UNDER THE STRUCTURE MODULES
  id <- paste( paste("clm",cgiarPIPE::idGenerator(5,5),sep=""), markerDTfile$idOriginal,markerDTfile2$idOriginal, sep = "_")
  type <- "clm"
  
  ############################
  # loading the dataset
  if (is.null(markerDTfile)) stop("No input marker data file specified.")
  if (is.null(markerDTfile2)) stop("No input marker data file 2 specified.")
  
  ## extract marker matrices and reference alleles
  M1 <- markerDTfile$cleaned$M 
  M2 <- markerDTfile2$cleaned$M 
  refM1 <- markerDTfile$cleaned$ref.alleles
  refM2 <- markerDTfile2$cleaned$ref.alleles
  # fill the missing references
  if(is.na(refM1) | is.na(refM2) ){
    if(verbose){cat("Reference alleles are not available for one of the marker files (maybe you started with numeric calls). \nCareful with results")}
  }else{
    badRefM1 <- which(is.na(refM1["Ref",]))
    if(length(badRefM1) > 0){ refM1["Ref",badRefM1] <- refM1["Alt",badRefM1]  }
    badRefM2 <- which(is.na(refM2["Ref",]))
    if(length(badRefM2) > 0){ refM2["Ref",badRefM2] <- refM2["Alt",badRefM2]  }
  }
  # reduce to common snps
  commonSNPs <- intersect(colnames(M1),colnames(M2))
  if (length(commonSNPs) == 0) stop("No SNPs in common. Please provide marker matrices that have SNPs in common")
  M1 <- M1[,commonSNPs]
  M2 <- M2[,commonSNPs]
  if(!is.na(refM1) & !is.na(refM2) ){
    refM1 <- refM1[,commonSNPs]
    refM2 <- refM2[,commonSNPs]
    # check that both matrices used the same reference alleles (at least 95% of them)
    lengthEqualRefs <- length(which(refM1["Ref",] == refM2["Ref",]))
    if (lengthEqualRefs < (length(commonSNPs)*.95) ) stop("The marker matrices provided have different reference alleles. Please call the marker matrices using the same reference alleles")
  }
  
  ############################
  ## first get all possible hybrids
  res1 <- build.HMM(M1, M2, 
                    return.combos.only = TRUE, separator=separator)
  
  ## build the marker matrix for batches of 1000 hybrids
  batches <- sort(rep(1:1000,min(c(nrow(res1$data.used),hybridBatch))))
  res1$data.used$batch <- batches[1:nrow(res1$data.used)]
  data.usedBatches <- split(res1$data.used, res1$data.used$batch)
  # res2 <- vector(mode="list", length = length(data.usedBatches))
  # start the loop
  for(i in 1:length(data.usedBatches)){
    prov <- build.HMM(M1, M2,
                      custom.hyb = data.usedBatches[[i]], 
                      separator=separator
    )
    if(i == 1){
      if(additive){
        M <- prov$HMM.add
      }else{
        M <- prov$HMM.dom
      }
    }else{
      if(additive){
        M <- rbind(M,prov$HMM.add)
      }else{
        M <- rbind(M,prov$HMM.dom)
      }
    }
    
  }
  # M <- do.call(rbind,res2)
  final <- list(M=as.matrix(M), ref.alleles=refM1)
  #########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	markerDTfile$id,
    markerbooks	= NA,  markerDataFile =	markerDTfile2$id,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = NA
  )
  
  # saveRDS(A, file = file.path(wd,"files_cleaned",paste0(id,".rds")))
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the files_cleaned folder under such id \n"))
  }
  result <- list(metrics=NA, predictions=NA, modeling=NA, metadata=metadata,
                 cleaned=final, outliers=NA, desire=NA, id=id, idOriginal=markerDTfile$idOriginal,
                 metadataFieldinst=NA, singleCrosses=res1$data.used)
  return(result)
}
