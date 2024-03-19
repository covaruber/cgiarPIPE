cleanp <- function(
  phenoDTfile= NULL, pipeline="ABSENT",
  stage= "ABSENT", timePoint= "ABSENT", season= "ABSENT",
  location= "ABSENT", country= "ABSENT",
  trial= "ABSENT", geno= "ABSENT", # design= "ABSENT",
  mother="ABSENT", father="ABSENT",
  genoCode="ABSENT", genoYearOrigin="ABSENT",
  rep= "ABSENT", block= "ABSENT",rowcoord= "ABSENT",
  colcoord= "ABSENT",entryType= "ABSENT",
  # trait parameters
  trait= NULL, fieldinst= c("timePoint", "season", "location"),
  # outlier cleaning parameters
  outlierCoef= NA,
  traitLB = 0.0001, traitUB= Inf,# wd=NULL,
  verbose=FALSE,  idOriginal=NULL,
  # fieldinst metadata parameters
  latitude="ABSENT", longitude="ABSENT", altitude="ABSENT",
  plantingDate = "ABSENT", harvestingDate="ABSENT",
  cleanCharacters=TRUE
){
  ## THIS FUNCTION TAKES A RAW PHENOTYPIC FILE AND MATCHES THE COLUMN NAMES TO CREATE A STANDARD FORMAT
  ## IT IS USED IN THE MATCH COLUMNS PHENOTYPE MODULE IN BCLEAN
  idOriginal <- gsub(".csv","",cgiarBase::cleanChar(as.character(gsub(" ","_",idOriginal))))
  id <- paste( paste("clp",idGenerator(5,5),sep=""), idOriginal, sep = "_")
  if(is.null(trait)){stop("Please provide traits to be carried through the pipeline.", call. = FALSE)}
  if(geno=="ABSENT"){stop("Please provide the geno column. This cannot be ABSENT", call. = FALSE)}
  outlierCoef <- rep(outlierCoef,100); outlierCoef <- outlierCoef[1:length(trait)]
  traitLB <- rep(traitLB,100); traitLB <- traitLB[1:length(trait)]
  traitUB <- rep(traitUB,100); traitUB <- traitUB[1:length(trait)]
  type <- "cleaningp"
  ###################################
  ###################################
  ## shape the data
  if (is.null(phenoDTfile)){stop("No input phenotypic data file specified.")}
  mydata <- as.data.frame(phenoDTfile)
  datacolparams <- c("pipeline","timePoint", "season", "location", "country", "trial", #"design",
                     "geno","mother","father", "genoCode","genoYearOrigin","rep", "block", "rowcoord", "colcoord", "entryType", "stage",
                     "latitude","longitude","altitude","plantingDate","harvestingDate"
                     )
  colsFieldsNotAllowed <- setdiff(fieldinst,datacolparams)
  if(length(colsFieldsNotAllowed) > 0){
    stop("Columns to form the fieldinst levels cannot be different to the new standardized names.")
  }
  provided <- c(pipeline,timePoint, season, location, country, trial, geno, mother, father, genoCode, genoYearOrigin, rep, block, rowcoord, colcoord, entryType, stage,
                latitude,longitude,plantingDate,harvestingDate,altitude) # design
  counter=1
  columnClasses <- unlist(lapply(mydata,class))
  for(j in 1:length(datacolparams)){ # we check parameters belonging to data columns
    if(provided[j] == "ABSENT"){ # if column doesn't exist
      # if ABSENT column should be numeric
      if(datacolparams[j] %in% c("rep", "block", "rowcoord", "colcoord","timePoint","genoYearOrigin","altitude")){
        mydata[,datacolparams[counter]] <- 1
      }else if(datacolparams[j] %in% c("plantingDate","harvestingDate")){
        if(datacolparams[j] %in% c("plantingDate")){
          mydata[,datacolparams[counter]] <- Sys.Date() #"2011-01-01" #Sys.Date()-180 # 6 months ago
        }else{
          mydata[,datacolparams[counter]] <- Sys.Date() # "2020-12-31"  #mydata[,"plantingDate"]
        }
      }else if(datacolparams[j] %in% c("longitude","latitude")){
        mydata[,datacolparams[counter]] <- 1000
      }else{# if ABSENT column should be character
        mydata[,datacolparams[counter]] <- paste("base",datacolparams[counter],sep="_")
      }
    }else{ # if it does exist
      if (!provided[j] %in% colnames(mydata)){
        stop(paste0("'", provided[j], "' is not a column in the given dataset or you have repeated a column name"))
      }else{
        mydata[,datacolparams[j]] <- mydata[,provided[j]]
        # ensure the followin columns are numeric
        if(datacolparams[j] %in% c("rep", "block", "rowcoord", "colcoord","timePoint","genoYearOrigin","latitude","longitude","altitude")){
          if(columnClasses[provided[j]] %in% c("character","factor") ){
            mydata[,datacolparams[j]] <- as.numeric(as.factor(mydata[,datacolparams[j]]))
          }
        }else if(datacolparams[j] %in% c("plantingDate","harvestingDate")){
          mydata[,datacolparams[j]] <- gsub('[^[:digit:] ]', '-', mydata[,datacolparams[j]]) # leave only digits
          mydata[,datacolparams[j]] <- as.Date(mydata[,datacolparams[j]])
        }
      }
    }; counter=counter+1
  }
  # check if traits are present
  for(k in 1:length(trait)){
    if (!trait[k] %in% colnames(mydata)){
      stop(paste0("'", trait[k], "' is not a column in the given dataset or you have repeated a column name"))
    }else{
      mydata[,trait[k]] <- as.numeric(mydata[,trait[k]])
    }
  }
  # subset to relevant columns
  datacolparamsData <- setdiff(datacolparams, c("latitude","longitude","altitude"))
  pick <- which(colnames(mydata) %in% c(datacolparamsData,trait))
  mydatan <- mydata[,pick]
  # make sure the geno column is always uppercase and no spaces
  if(cleanCharacters){
    for(iName in c("geno","mother","father")){ # iName = c("geno","mother","father")[1]
      mydatan[,iName] <- cgiarBase::cleanChar(as.character(mydatan[,iName])) # change to upper case
      mydatan[,iName] <- stringi::stri_trans_general( mydatan[,iName] , "Latin-ASCII")
    }
    # make sure all the columns the user wants to use to form fieldinst have good names
    for(iPar in setdiff(fieldinst,"timePoint")){ # except timepoint which is numeric in nature
      mydatan[,iPar] <- cgiarBase::cleanCharField(as.character(mydatan[,iPar]))
    }
  }
  # entryType
  et <- which(is.na(mydatan$entryType))
  if(length(et) > 0){mydatan$entryType[et]= "test"}

  # convert to factors
  for(iVar in datacolparamsData){
    mydatan[,paste0(iVar,"F")] <- as.factor(mydatan[,iVar])
  }
  # add field instance or environment column where the model will be run on
  if(length(fieldinst) == 0){stop("fieldinst forming columns are a required argument", call. = FALSE)}
  if(length(fieldinst) > 1){
    mydatan[,"fieldinstF"] <- apply(mydatan[,fieldinst],1,function(x){paste(x,collapse = "_")})
  }else{
    mydatan[,"fieldinstF"] <- mydatan[,fieldinst]
  }
  if(length(unique(mydatan[,"fieldinstF"])) == nrow(mydatan)){
    stop("The number of fieldinst levels is equal to the number of records/rows in the dataset. 
          This means that probably you didn't select the right columns to define the fieldinst column. 
          Please match again your raw file checking carefully the columns defining the fieldinst.", call.=FALSE )
  }
  mydatan[,"rowindex"] <- 1:nrow(mydatan)
  # metadata
  metadata <-  unique(mydata[,unique(c(fieldinst,"timePoint","location","country","latitude","longitude","altitude","plantingDate","harvestingDate"))])
  if(is.null(dim(metadata[,fieldinst]))){
    metadata$fieldinst <- metadata[,fieldinst]
  }else{
    metadata$fieldinst <- as.factor(apply(metadata[,fieldinst],1,function(x){paste(x,collapse = "_")}))
  }
  metadata <- metadata[which(!duplicated(metadata[,"fieldinst"])),]
  ######################################
  ## clean the data: removing outliers by fieldinst
  mydata_cleaned <- mydatan
  if(!is.null(trait)){
    outList <- list(); counter=1
    for (j in 1:length(trait)) {
      for (i in 1:nlevels(mydatan[, "fieldinstF"])) {
        sampleDT <- mydatan[mydatan[, "fieldinstF"] == levels(mydatan[, "fieldinstF"])[i], ]
        if(length(which(!is.na(sampleDT[, trait[j]]))) > 10){ # if enough data then clean
          if(!is.na(outlierCoef[j])){
            outlier <- boxplot.stats(x=sampleDT[, trait[j]],coef=outlierCoef[j] )$out
            toSilence <- sampleDT[which(sampleDT[,trait[j]] %in% outlier),"rowindex"]
            typeOut <- rep("outlier",length(toSilence))
          }else{
            toSilence <- numeric()
            typeOut <- character()
          }
          outOfBounds <- which((sampleDT[, trait[j]] < traitLB[j]) | (sampleDT[, trait[j]] > traitUB[j]) )
          if(length(outOfBounds) > 0){toSilence <- c(toSilence, sampleDT[outOfBounds,"rowIndex"]); typeOut <- c(typeOut, rep("bound",length(outOfBounds))) }
          if(length(toSilence) > 0){
            outList[[counter]] <- data.frame(indexRow=toSilence, traitName=trait[j],type="outlier",discardDecision=TRUE);
            counter=counter+1
          }
        }# end of if enough data
      }# end for each trial
    } # end for each trait
    outBind <- do.call(rbind, outList)
    if(!is.null(outBind)){outBind$analysisId <- id}
  }
  ####################################################
  ## update database: write the parameters to the parameter database
  metadata0 <- data.frame(
    analysisId	= id, analysisType =	type, fieldbooks	= NA,
    phenoDataFile =	NA,markerbooks	= NA,
    markerDataFile =	NA,timePoint = timePoint,season =	season,
    location =	location,country	= country,trial	= trial,
    geno = geno,rep	= rep,block =	block,
    rowcoord =	rowcoord,colcoord = colcoord,stage = stage
  )
  ## write the values used for cleaning to the modeling database
  modeling <- data.frame(
    trait = trait,traitLb = NA, traitUb = NA,
    outlierCoef = outlierCoef, analysisId = id,
    analysisType = type, fixedModel = NA,randomModel = NA,
    residualModel = NA,h2Threshold = NA
  )
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  result <- list(metrics=NA, predictions=NA, modeling=modeling, metadata=metadata0,
                 cleaned=mydata_cleaned, outliers=outBind, desire=NA, id=id, idOriginal=idOriginal,
                 metadataFieldinst=metadata)
  return(result)
}
