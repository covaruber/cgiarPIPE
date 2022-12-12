cleanp <- function(
  phenoDTfile= NULL, pipeline="ABSCENT",
  stage= "ABSCENT", year= "ABSCENT", season= "ABSCENT",
  location= "ABSCENT", country= "ABSCENT",
  trial= "ABSCENT", design= "ABSCENT",geno= "ABSCENT",
  genoCode="ABSCENT", genoYearOrigin="ABSCENT",
  rep= "ABSCENT", block= "ABSCENT",rowcoord= "ABSCENT",
  colcoord= "ABSCENT",entryType= "ABSCENT",
  # trait parameters
  trait= NULL, fieldinst= c("year", "season", "location"),
  # outlier cleaning parameters
  outlierCoef= 1.5,
  traitLB = 0.01, traitUB= Inf,# wd=NULL, 
  verbose=FALSE
){
  #
  # if(is.null(wd)){wd <- getwd()}
  # md <- strsplit(wd,"/")[[1]]; md <- md[length(md)]
  # if(md != "DB"){stop("Please set your working directory to the DB folder", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be cleaned", call. = FALSE)}
  # traitClasses <- unlist(lapply(xx,class))
  outlierCoef <- rep(outlierCoef,100); outlierCoef <- outlierCoef[1:length(trait)]
  traitLB <- rep(traitLB,100); traitLB <- traitLB[1:length(trait)]
  traitUB <- rep(traitUB,100); traitUB <- traitUB[1:length(trait)]
  id <- paste("clp",idGenerator(5,5),sep="")
  type <- "cleaningp"
  ###################################

  ###################################
  ## shape the data
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  mydata <- phenoDTfile # read.csv(file.path(wd,"files_raw",phenoDTfile))
  # cleaning <- read.csv(file.path(wd,"cleaning.csv"))
  # parameters <- read.csv(file.path(wd,"parameters.csv"))
  # modeling <- read.csv(file.path(wd,"modeling.csv"))
  # defining the columns to be present - checking of params
  datacolparams <- c("pipeline","year", "season", "location", "country", "trial", "design",
                     "geno", "genoCode","genoYearOrigin","rep", "block", "rowcoord", "colcoord", "entryType", "stage")
  provided <- c(pipeline,year, season, location, country, trial, design, geno, genoCode, genoYearOrigin, rep, block, rowcoord, colcoord, entryType, stage)
  counter=1
  for(j in 1:16){ # we check parameters belonging to data columns
    if(provided[j] == "ABSCENT"){ # if column doesn't exist
      if(datacolparams[j] %in% c("rep", "block", "rowcoord", "colcoord","year","genoYearOrigin")){ # if abscent column should be numeric
        mydata[,datacolparams[counter]] <- 1
      }else{# if abscent column should be character
        mydata[,datacolparams[counter]] <- paste("dummy",datacolparams[counter],sep="_")
      }
    }else{ # if it does exist
      if (!provided[j] %in% colnames(mydata)){
        stop(paste0("'", provided[j], "' is not a column in the given dataset or you have repeated a column name"))
      }else{
        mydata[,datacolparams[j]] <- mydata[,provided[j]]
        # colnames(mydata) <- replaceValues(colnames(mydata),provided[j], datacolparams[counter])
        if(datacolparams[j] %in% c("rep", "block", "rowcoord", "colcoord")){
          mydata[,datacolparams[j]] <- gsub('[^[:digit:] ]', '', mydata[,datacolparams[j]]) # leave only digits
          mydata[,datacolparams[j]] <- as.numeric(mydata[,datacolparams[j]])
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
  mydatan <- mydata[,c(datacolparams,trait)]
  # make sure fieldinst or fieldinst columns don't have messy identifiers
  for(iPar in setdiff(fieldinst,"year")){
    mydatan[,iPar] <- gsub(" ","_",mydatan[,iPar])
  }
  # entryType
  mydatan$entryType2 <- "test"
  areChecks <- grep("check",mydatan$entryType, ignore.case = TRUE)
  if(length(areChecks) > 0){mydatan$entryType2[areChecks] <- "check"}
  # convert to factors
  for(iVar in datacolparams){
    mydatan[,paste0(iVar,"F")] <- as.factor(mydatan[,iVar])
  }
  # add field instance or environment column where the model will be run on
  mydatan[,"fieldinstF"] <- apply(mydatan[,fieldinst],1,function(x){paste(x,collapse = "_")})
  mydatan[,"fieldinstF"] <- gsub(",","_",mydatan[,"fieldinstF"]) # can cause issues in model fitting
  mydatan[,"fieldinstF"] <- gsub("!","_",mydatan[,"fieldinstF"])
  mydatan[,"fieldinstF"] <- gsub(" ",".",mydatan[,"fieldinstF"])
  mydatan[,"fieldinstF"] <- as.factor(mydatan[,"fieldinstF"])
  mydatan[,"rowindex"] <- 1:nrow(mydatan)

  # identifying checks for further analysis
  checks <- unique(mydatan[which(mydatan[, "entryType2"] != "test"), "geno"])

  ######################################
  ## clean the data
  # removing outliers by fieldinst
  mydata_cleaned <- mydatan

  if(!is.null(trait)){
    outList <- list(); counter=1
    for (j in 1:length(trait)) {
      for (i in 1:nlevels(mydatan[, "fieldinstF"])) {
        sampleDT <- mydatan[mydatan[, "fieldinstF"] == levels(mydatan[, "fieldinstF"])[i], ]
        if(length(which(!is.na(sampleDT[, trait[j]]))) > 10){ # if enough data then clean
          outlier <- boxplot.stats(x=sampleDT[, trait[j]],coef=outlierCoef[j] )$out
          toSilence <- sampleDT[which(sampleDT[,trait[j]] %in% outlier),"rowindex"]
          typeOut <- rep("outlier",length(toSilence))
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
    # print(str(outBind))
    if(!is.null(outBind)){outBind$analysisId <- id}
  }
  ####################################################
  ## update database
  # saveRDS(mydata_cleaned, file = file.path(wd,"files_cleaned",paste0(id,".rds")))
  ## write the outliers to the cleaning database
  # if(!is.null(trait)){saveRDS(outBind, file = file.path(wd,"outliers",paste0(id,".rds")))}
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= id, analysisType =	type, fieldbooks	= NA,
    phenoDataFile =	NA,markerbooks	= NA,
    markerDataFile =	NA,year = year,season =	season,
    location =	location,country	= country,trial	= trial,
    design =	design,geno = geno,rep	= rep,block =	block,
    rowcoord =	rowcoord,colcoord = colcoord,stage = stage
  )
  # saveRDS(db.params, file = file.path(wd,"metadata",paste0(id,".rds")))
  # parameters <- unique(rbind(parameters, db.params[,colnames(parameters)]))
  # write.csv(parameters, file = file.path(wd,"parameters.csv"),row.names = FALSE)
  ## write the values used for cleaning to the modeling database
  mod <- data.frame(
    trait = trait,traitLb = NA, traitUb = NA,
    outlierCoef = outlierCoef, analysisId = id,
    analysisType = type, fixedModel = NA,randomModel = NA,
    residualModel = NA,h2Threshold = NA
  )
  # saveRDS(mod, file = file.path(wd,"modeling",paste0(id,".rds")))
  # modeling <- unique(rbind(modeling, mod[,colnames(modeling)]))
  # write.csv(modeling, file = file.path(wd,"modeling.csv"),row.names = FALSE)

  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    # cat(paste("Your results will be available in the 'files_cleaned' folder under such id \n"))
  }
  result <- list(metrics=NA, predictions=NA, modeling=mod, metadata=db.params,
                 cleaned=mydata_cleaned, outliers=outBind, desire=NA, id=id)
  return(result)#paste("cleaning done:",id))
}
