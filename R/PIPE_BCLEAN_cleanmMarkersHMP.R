cleanmHMP <- function(
    markerDTfile= NULL,
    linkFile=NULL, linkFileGenoCol="ABSCENT",
    linkFilesSampleCol="ABSCENT",
    markerColumn=NULL, # where marker ids are
    indsColumn=NULL, # where individuals start
    metadataColumns=NULL,
    idOriginal=NULL,
    refAlleles=NULL,
    singleLetter=FALSE,
    cleanCharacters=TRUE
){
  ## FUNCTION TO MATCH MARKER FILES TO DESIRED COLUMNS (GENO & USECOLUMNS) AND CONVERT TO NUMERIC CODE (-1,0,1) AND DO SOME BASIC CLEANING
  ## THIS ONE STARTS FROM AN .HMP FILE
  ## IT IS USED IN THE MATCH COLUMNS MARKERS HMP MODULE IN BCLEAN
  if(is.null(markerColumn)){stop("Marker column indices are needed")}
  if(is.null(indsColumn)){stop("Individual column index is needed")}
  
  newColnames <- markerDTfile[,markerColumn] # newColnames
  markerDTfile <- t(markerDTfile[,-markerColumn])
  if(!is.null(metadataColumns)){
    mymetadata <- markerDTfile[,c(markerColumn,metadataColumns)]
  }else{mymetadata <- character()}
  colnames(markerDTfile) <- newColnames
  take <- setdiff(1:nrow(markerDTfile), c(1:(indsColumn-2))) # -2 because we already reduced the column with column names
  markerDTfile <- markerDTfile[take,]
  replaceXs <- function(x){ifelse( startsWith(x, "X"), substr(x, 2, nchar(x)), x) }
  rownames(markerDTfile) <- replaceXs(rownames(markerDTfile))  # gsub("X","", rownames(markerDTfile))
  rownames(markerDTfile) <- stringi::stri_trans_general(rownames(markerDTfile), "Latin-ASCII")
  # rownames(markerDTfile) <- gsub("X","", rownames(markerDTfile))
  rownames(markerDTfile) <- gsub("[.]","-", rownames(markerDTfile))
  if(ncol(markerDTfile) < 10000){
    rownames(markerDTfile) <- cgiarBase::replaceValues(rownames(markerDTfile), Search=c("1-Jan","1-Feb","1-Mar","1-Apr","1-May","1-Jun","1-Jul","1-Aug","1-Sep","1-Oct","1-Nov","1-Dec"), Replace=paste(1:12,1,sep="-") )
  }
  
  markerDTfile <- cbind(rownames(markerDTfile), markerDTfile); colnames(markerDTfile)[1] <- "GID"
  if(!is.null(linkFile)){
    if (!linkFileGenoCol %in% colnames(linkFile)){
      stop(paste0("'", linkFileGenoCol, "' is not a column in the link dataset.Please check"))
    }
    if (!linkFilesSampleCol %in% colnames(linkFile)){
      stop(paste0("'", linkFilesSampleCol, "' is not a column in the link dataset.Please check"))
    }
    markerDTfile[,"GID"] <- cgiarBase::replaceValues(Source=markerDTfile[,"GID"],Search = linkFile[,linkFilesSampleCol],Replace = linkFile[,linkFileGenoCol] )
  }
  rownames(markerDTfile) <- NULL
  
  geno="GID"; useColumns=1:ncol(markerDTfile);
  verbose=FALSE;
  missingData=c("NN","FAIL","FAILED","Uncallable","Unused","NA","")
  specialToRemove=c("[:]","[|]")
  idOriginal=idOriginal;
  refAlleles=refAlleles
  
  id <- paste( paste("clm",cgiarPIPE::idGenerator(5,5),sep=""), gsub(".csv","",gsub(" ","",idOriginal)), sep = "_")
  type <- "cleaningm"
  ############################
  # loading the dataset
  if (is.null(markerDTfile)) stop("No input marker data file specified.")
  if(!geno %in% colnames(markerDTfile) ){
    stop("Please make sure that a column 'geno' with genotype names is included in the file and the rest of the columns is only markers",call. = FALSE)
  }else{
    colnames(markerDTfile)[which(colnames(markerDTfile) == geno)] <- "geno"
  }
  ############################
  # transform to numbers
  genos <- markerDTfile[,"geno"]
  markerDTfile <- markerDTfile[,-which(colnames(markerDTfile) %in% "geno")]
  # replace missing data foraneous call
  if(!is.null(missingData)){
    print("setting missing data")
    for(u in 1:length(missingData)){
      miss <- which(markerDTfile == missingData[u], arr.ind = TRUE)
      if(nrow(miss) > 0){
        markerDTfile[miss] = NA
      }
    }
  }
  # replace special characters with empty space
  if(!is.null(specialToRemove)){
    for(u in 1:length(specialToRemove)){
      markerDTfile <- gsub(specialToRemove[u],"",markerDTfile)
    }
  }
  shouldBeTransformed=FALSE # is input letters or numbers
  isLetter <- which(unlist(lapply(markerDTfile[,1:min(c(10,ncol(markerDTfile)))], class)) == "character")
  if(length(isLetter) > 1){shouldBeTransformed=TRUE}  # if we found these are nor numbers but letters we transform
  if(shouldBeTransformed){ # is letter coded
    ## remove markers with too much missing indivisuals
    print("removing markers with too much missing data")
    xxx <- which(is.na(markerDTfile), arr.ind = TRUE)
    badCols <- as.numeric(names(which( table(xxx[,"col"]) > (nrow(markerDTfile)*.5))))
    badRows <- as.numeric(names(which( table(xxx[,"row"]) > (ncol(markerDTfile)*.5))))
    if(length(badCols) > 0){
      markerDTfile <- markerDTfile[,-badCols]
    }
    if(length(badRows) > 0){
      markerDTfile <- markerDTfile[-badRows,]
      genos <- genos[-badRows]
    }
    if(median(nchar(markerDTfile[,1]), na.rm=TRUE) <= 1){
      stop("This function only allows double-letter coded markers. Please look at other function to transform it first", call. = FALSE)
    }
    if(singleLetter){
      rowNamesM <- rownames(markerDTfile)
      markerDTfile <- apply(markerDTfile,2,function(x){
        nDigits <- nchar(x)
        toTransform <- which(nDigits == 1)
        if(length(toTransform) > 0){
          x[toTransform] <- paste0(x[toTransform],x[toTransform])
        }
        toDelete <- which(nDigits > 2)
        if(length(toDelete) > 0){
          x[toDelete] <- NA
        }
        return(x)
      })
      rownames(markerDTfile) <- rowNamesM
    }
    
    M0 <- cgiarBase::atcg1234(markerDTfile, ref.alleles = refAlleles, silent = TRUE)
    
    rownames(M0$M) <- genos
  }else{ # is coded numerically
    
    ## remove markers with too much missing indivisuals
    xxx <- which(is.na(markerDTfile), arr.ind = TRUE)
    badCols <- as.numeric(names(which( table(xxx[,"col"]) > (nrow(markerDTfile)*.5))))
    badRows <- as.numeric(names(which( table(xxx[,"row"]) > (ncol(markerDTfile)*.5))))
    if(length(badCols) > 0){
      markerDTfile <- markerDTfile[,-badCols]
    }
    if(length(badRows) > 0){
      markerDTfile[-badRows,]
      genos <- genos[-badRows]
    }
    #  impute the rest
    markerDTfile <- apply(markerDTfile,2,sommer::imputev) # impute
    ## center
    markerDTfile <- apply(markerDTfile,2,function(x){y <- x - (max(x)+min(x))/2; return(y)}) # make sure is centered (if monomorphic it will be zero)
    M0 <- list(M=markerDTfile, ref.alleles=NA)
    
    rownames(M0$M) <- genos
  }
  # ensure no weird names in individuals
  if(cleanCharacters){
    rownames(M0$M) <- cgiarBase::cleanChar(as.character(rownames(M0$M))) # remove special characters
  }
  #########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,  analysisType =	type,
    fieldbooks	= NA,  phenoDataFile =	NA,
    markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,  stage = NA
  )
  #
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the files_cleaned folder under such id \n"))
  }
  result <- list(metrics=NA, predictions=NA, modeling=NA, metadata=metadata,
                 cleaned=M0, outliers=NA, desire=NA, id=id, idOriginal=idOriginal,
                 metadataFieldinst=NA)
  result$hmpMetadata <- mymetadata
  return(result)
}
