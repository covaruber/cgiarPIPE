cleanFields <- function(
    phenoDTfile= NULL,
    fieldsToClean=NULL,
    verbose=FALSE
){
  ## FUNCTION TO SET TO MISSING DATA SPECIFIC LEVELS FOR A VARIABLE (E.G. BLOCK) WITHIN AN SPECIFIC FIENDINST (E.G., IBADAN_2010)
  ## IT IS USED IN THE CLEAN FIELDS MODULE IN BCLEAN
  id <- phenoDTfile$id
  id <- paste(id,"FC",sep="_")
  if(is.null(phenoDTfile)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(fieldsToClean)){stop("Please provide table of fields to be cleaned", call. = FALSE)}

  #####################################
  # cleaning
  mydata <- phenoDTfile$cleaned
  for(i in 1:ncol(fieldsToClean)){ # for each fieldinstance
    for(j in 1:nrow(fieldsToClean)){ # for each experimental design factor
      if(fieldsToClean[j,i] == 0){
        fieldI <- which(mydata[,"fieldinstF"] == colnames(fieldsToClean)[i])
        mydata[fieldI,rownames(fieldsToClean)[j]] = NA
      }
    }
  }
  # move results to the factor column
  for(j in 1:nrow(fieldsToClean)){ # for each experimental design factor
    mydata[,paste0(rownames(fieldsToClean)[j],"F")] <- as.factor(mydata[,rownames(fieldsToClean)[j]])
  }
  phenoDTfile$cleaned <- mydata
  ##########################################
  ## update databases
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  phenoDTfile$id <- id
  return(phenoDTfile)
}
