pcheck <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    checksAndTraitsToInclude=NULL,
    verbose=TRUE
){
  ## THIS FUNCTION CALCULATES A PERCENTAGE ADVANTAGE OVER CHECKS 
  ## IS USED IN THE BANAL APP UNDER THE TRANSFORMATION MODULES
  id <- paste0(phenoDTfile$id,"_pcheck")
  type <- "pch"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(length(grep("met",phenoDTfile)) == 0){stop("Percentage check comparison can only be calculated on results from a MET analysis using across environment predictions",call. = FALSE)}
  if(is.null(checksAndTraitsToInclude)){stop("Please provide traits to be analyzed", call. = FALSE)}

  ############################
  # loading the dataset
  mydata <- phenoDTfile$predictions # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))

  ############################
  ## percentage check calculation

  myNewData <- list(); counter <- 1
  for(i in 1:nrow(checksAndTraitsToInclude)){
    for(j in 1:ncol(checksAndTraitsToInclude)){
      if(!is.na(checksAndTraitsToInclude[i,j])){
        if(checksAndTraitsToInclude[i,j] > 0){
          checkName <- rownames(checksAndTraitsToInclude)[i]
          traitName <- colnames(checksAndTraitsToInclude)[j]
          vTrait <- which(mydata$trait == traitName)
          if(length(vTrait) > 0){# there's information for that trait
            provData <- mydata[vTrait,]
            vCheck <- which(provData$geno == checkName)
            if(length(vCheck) > 0){# there's information for the check
              provDataCheck <- provData[vCheck,]
              meanCheck <- mean(provDataCheck$predictedValue)
              if(!is.na(meanCheck)){# there's an actual value
                provData0 <- provData
                provData[,"predictedValue"] <- ((provData0[,"predictedValue"] - meanCheck)/meanCheck) * 100
                provData[,"stdError"] <- ((provData0[,"stdError"])/meanCheck) * provData0[,"predictedValue"]
                provData$trait <- paste(provData$trait,"percOver",checkName,sep="-")
                myNewData[[counter]] <- provData
                counter <- counter+1
              }
            }
          }
        }
      }
    }
  }

  myNewDataRbind <- do.call(rbind, myNewData)
  ###########
  phenoDTfile$predictions <- rbind(phenoDTfile$predictions,myNewDataRbind)
  phenoDTfile$id <- id
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  return(phenoDTfile)#paste("index done:",id))
}
