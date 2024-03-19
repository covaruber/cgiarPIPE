rbindm <- function(
    markerDTfile= NULL,
    markerDTfile2= NULL,
    verbose=FALSE
){
  ## THIS FUNCTION ROW BINDS TWO MARKER MATRICES AND ASSOCIATED TABLES
  ## IS USED IN THE BCLEAN APP UNDER THE TRANSFORMATIONS MODULES
  if(is.null(markerDTfile)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(markerDTfile2)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  
  start=substr(markerDTfile$id,1,3)
  start2=substr(markerDTfile2$id,1,3)
  id <- paste( paste("rbm",cgiarPIPE::idGenerator(5,5),sep=""), markerDTfile$idOriginal,"extended", sep="_")
  type <- "rbm"
  
  ###################################
  # loading the dataset
  res <- list()
  elementsToMerge <- intersect(names(markerDTfile), names(markerDTfile2))
  for(i in elementsToMerge){ # i = "cleaned"
    if(i == "cleaned"){
      
      ## extract marker matrices and reference alleles
      M1 <- markerDTfile$cleaned$M 
      colnames(M1) <- gsub("[[:punct:]]", "_", colnames(M1))
      M2 <- markerDTfile2$cleaned$M 
      colnames(M2) <- gsub("[[:punct:]]", "_", colnames(M2))
      refM1 <- markerDTfile$cleaned$ref.alleles
      colnames(refM1) <- gsub("[[:punct:]]", "_", colnames(refM1))
      refM2 <- markerDTfile2$cleaned$ref.alleles
      colnames(refM2) <- gsub("[[:punct:]]", "_", colnames(refM2))
      # fill the missing reference alleles
      if(length(refM1) > 1){ 
        badRefM1 <- which(is.na(refM1["Ref",]))
        refM1["Ref",badRefM1] <- refM1["Alt",badRefM1]
      }else{
        stop("Error in marker file 1. Binding marker files requires reference alleles to be available. Probably you provided numeric calls from the start. Those cannot be merged.")
      }
      if(length(refM2) > 1){ 
        badRefM2 <- which(is.na(refM2["Ref",]))
        refM2["Ref",badRefM2] <- refM2["Alt",badRefM2]  
      }else{
        stop("Error in marker file 2. Binding marker files requires reference alleles to be available. Probably you provided numeric calls from the start. Those cannot be merged.")
      }
      # reduce to common snps
      commonSNPs <- intersect(colnames(M1),colnames(M2))
      if (length(commonSNPs) == 0) stop("No SNPs in common. Please provide marker matrices that have SNPs in common or use similar SNP names")
      # which additional individuals exist in marker file 2
      keepM2 <- setdiff(rownames(M2), rownames(M1))
      if(length(keepM2) > 1){ # if the marker file 2 has some additional individuals compared to M1
        M2 <- M2[keepM2,,drop=FALSE] # only keep samples from M2 that are different from M1
        # check that both matrices used the same reference alleles (at least 95% of them)
        if(length(refM1) > 1){ # if there's reference alleles
          # take the numeric matrices and bring them back to the original letter format
          M2b <- cgiarBase::markerBackTransform(marks=M2+1, refs=refM2)
          M1b <- cgiarBase::markerBackTransform(marks=M1+1, refs=refM1)
          # form different sets depending on common or differential SNPs
          Mc <- rbind(M1b[,commonSNPs], M2b[, commonSNPs])
          Md1 <- M1b[,setdiff(colnames(M1b), colnames(M2b)), drop=FALSE]
          Md2 <- M2b[,setdiff(colnames(M2b), colnames(M1b)), drop=FALSE]
          # now make them numeric again
          if(ncol(Mc) > 0){McN <- cgiarBase::atcg1234(Mc, maf=-.1, silent = TRUE)}else{McN <- list(M=Mc, ref.alleles=NULL)}
          if(ncol(Md1) > 0){Md1N <- cgiarBase::atcg1234(Md1, maf=-.1, silent = TRUE)}else{Md1N <- list(M=Md1, ref.alleles=NULL)}
          if(ncol(Md2) > 0){Md2N <- cgiarBase::atcg1234(Md2, maf=-.1, silent = TRUE)}else{Md2N <- list(M=Md2, ref.alleles=NULL)}
          # form the final matrix using the numeric calls
          Mx <- matrix(NA, nrow=nrow(Md1N$M)+nrow(Md2N$M), ncol=ncol(Md1N$M)+ncol(Md2N$M))
          rownames(Mx) <- c(rownames(Md1N$M), rownames(Md2N$M))
          colnames(Mx) <- c(colnames(Md1N$M), colnames(Md2N$M))
          if(ncol(Md1N$M) > 0){Mx[rownames(Md1N$M),colnames(Md1N$M)] <- Md1N$M}
          if(ncol(Md2N$M) > 0){Mx[rownames(Md2N$M),colnames(Md2N$M)] <- Md2N$M}
          # bind the common markers and the differentials
          Mz <- cbind(McN$M, Mx[rownames(McN$M),])
        }
        res[[i]] <- list(M=Mz, ref.alleles= cbind(McN$ref.alleles, Md1N$ref.alleles, Md2N$ref.alleles))
      }else{
        res[[i]] <- markerDTfile$cleaned
      }
      
    }else{ # if is not the marker matrix but any other table
      dataList <- list(markerDTfile[[i]], markerDTfile2[[i]])
      dataListColNames <- lapply(dataList,colnames)
      commonNames <- Reduce(intersect, dataListColNames) #use intersect in a list of vectors
      if(is.null(commonNames)){
        res[[i]] <- NA
      }else{
        dataListCommon <- lapply(dataList,function(x){return(x[,commonNames])})
        res[[i]]<- do.call(rbind, dataListCommon)
      }
    }
  }
  names(res) <- elementsToMerge
  # identify which tables where pure NAs and shouldn't be merged
  badMerging <- which(unlist(lapply(res, function(y){length(which(is.na(y)))/(nrow(y)*ncol(y))})) == 1)
  if(length(badMerging) > 0){
    for(i in 1:length(badMerging)){
      res[names(badMerging)[i]] <- NA # transform a table of NAs to a single NA
    }
  }
  ###################################
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  res$id <- id
  res$idOriginal = markerDTfile$idOriginal
  
  return(res)
}
