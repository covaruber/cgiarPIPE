rbindp <- function(
    phenoDTfile= NULL,
    phenoDTfile2= NULL,
    verbose=FALSE
){
  ## THIS FUNCTION ROW BINDS PHENOTYPE FILES AND ASSOCIATED TABLES
  ## IS USED IN THE BCLEAN APP UNDER THE TRANSFORMATION MODULES
  if(is.null(phenoDTfile)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(phenoDTfile2)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  
  start=substr(phenoDTfile$id,1,3)
  start2=substr(phenoDTfile2$id,1,3)
  id <- paste( paste(start,idGenerator(5,5),sep=""), phenoDTfile$idOriginal,"extended", sep="_")
  type <- start
  
  ###################################
  # loading the dataset
  res <- list()
  elementsToMerge <- intersect(names(phenoDTfile), names(phenoDTfile2))
  for(i in elementsToMerge){ # i="cleaned"
    dataList <- list(phenoDTfile[[i]], phenoDTfile2[[i]])
    dataListColNames <- lapply(dataList,colnames)
    commonNames <- Reduce(intersect, dataListColNames) #use intersect in a list of vectors
    if(is.null(commonNames)){
      res[[i]] <- NA
    }else{
      # if((i == "cleaned") & (type == "clp")){ # for phenotype files in the cleaned tab make sure we don't loose traits
      differentNames <- c( setdiff(colnames(phenoDTfile[[i]]), colnames(phenoDTfile2[[i]])) , setdiff(colnames(phenoDTfile2[[i]]), colnames(phenoDTfile[[i]])) )
      allNames <- unique(unlist(dataListColNames))
      dfg1 <- data.frame(matrix(ncol = length(allNames), nrow = nrow(phenoDTfile[[i]]))); 
      colnames(dfg1) <- allNames
      dfg1[,colnames(phenoDTfile[[i]])] <- phenoDTfile[[i]]
      dfg2 <- data.frame(matrix(ncol = length(allNames), nrow = nrow(phenoDTfile2[[i]]))); 
      colnames(dfg2) <- allNames
      dfg2[,colnames(phenoDTfile2[[i]])] <- phenoDTfile2[[i]]
      newFile <- rbind(dfg1[,c(commonNames,differentNames)],dfg2[,c(commonNames,differentNames)])
      rownames(newFile) <- NULL
      isCleanedPhenoTable <- which(colnames(newFile) == "rowindex")
      if(length(isCleanedPhenoTable) > 0){newFile[,isCleanedPhenoTable] <- 1:nrow(newFile)}
      res[[i]]<- newFile
      # }else{ # only common names should be joined for any other table
      #   dataListCommon <- lapply(dataList,function(x){return(x[,commonNames])})
      #   res[[i]]<- do.call(rbind, dataListCommon)
      # }
    }
  }
  names(res) <- elementsToMerge
  # res=mapply(rbind, phenoDTfile, phenoDTfile2, SIMPLIFY=FALSE)
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
  res$idOriginal = phenoDTfile$idOriginal
  return(res)
}
