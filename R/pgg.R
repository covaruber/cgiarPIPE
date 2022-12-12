pgg <- function(
    phenoDTfile= NULL,
    trait=NULL, # per trait
    fieldinst=NULL,
    year=NULL,
    percentage=10,
    lIdeal=NA,# wd=NULL, 
    verbose=TRUE
){

  # if(is.null(wd)){wd <- getwd()}
  # md <- strsplit(wd,"/")[[1]]; md <- md[length(md)]
  # if(md != "DB"){stop("Please set your working directory to the DB folder", call. = FALSE)}

  id <- paste("pgg",idGenerator(5,5),sep="")
  type <- "pgg"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}

  ############################
  # loading the dataset
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  mydata <- phenoDTfile$predictions #readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))

  ############################
  ## gg analysis
  p <- percentage/100
  i <- dnorm(qnorm(1 - p))/p
  R <- ggAge <- ggIdeal <- age <- sigma <- r <- numeric();
  counter=1
  for(iTrait in trait){ # iTrait="YLD_TON"
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    # subset data
    mydataSub <- droplevels(mydata[which((mydata$trait == iTrait) & (mydata$fieldinst %in% fieldinst) & (mydata$genoYearTesting == year) ),])
    # calculate parameters
    rels <- mydataSub$rel;
    badrels <- which(rels < 0); if(length(badrels) > 0){ rels[badrels] <- NA}
    r[iTrait] <- ifelse(length(na.omit(rels)) > 0, mean(sqrt(na.omit(rels)), na.rm=TRUE), 0)

    sigma[iTrait] <- sd(mydataSub$predictedValue, na.rm = TRUE)

    mydataSubSorted <- mydataSub[with(mydataSub, order(-predictedValue)), ]
    mydataSubSortedSel <- mydataSubSorted[1:round(nrow(mydataSubSorted) * p),]
    age[iTrait] <- mean(mydataSubSortedSel$genoYearTesting, na.rm=TRUE) - mean(mydataSubSortedSel$genoYearOrigin, na.rm=TRUE)

    R[iTrait] <- r[iTrait] * sigma[iTrait] * i

    ggAge[iTrait] =  R[iTrait]/age[iTrait]

    ggIdeal[iTrait] =  R[iTrait]/lIdeal
  }
  #########################################
  ## update databases
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	NA,
    markerbooks	= NA,  markerDataFile =	NA,
    year = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(mydataSub$stage)),collapse=", ")
  )
  # saveRDS(db.params, file = file.path(wd,"metadata",paste0(id,".rds")))
  ## write the values used for cleaning to the modeling database
  mod <- data.frame(
    trait = trait,
    traitLb = NA,
    traitUb = NA,
    outlierCoef = NA,
    analysisId = id,
    analysisType = type,
    fixedModel = NA,
    randomModel = NA,
    residualModel = NA,
    h2Threshold = NA
  )
  # saveRDS(mod, file = file.path(wd,"modeling",paste0(id,".rds")))

  # write predictions

  # write pipeline metrics
  pm <- data.frame(value=c(rep(i,length(trait)),r,sigma, age, lIdeal, R, ggAge, ggIdeal),
                   stdError=NA,
                   fieldinst=fieldinst,
                   trait=rep(trait,8), analysisId=id,
                   method= c(rep("qnorm",length(r)), rep("reliability",length(r)), rep("sd",length(r)), rep("age",length(r)), rep("ideal",length(r)), rep("product",length(r)), rep("age", length(r)), rep("ideal", length(r)) ),
                   traitUnits=NA, stage = paste(sort(unique(mydataSub$stage)),collapse=", "),
                   parameter= c(rep("i",length(r)), rep("r",length(r)), rep("sigma",length(r)), rep("La",length(r)), rep("Li",length(r)), rep("R",length(r)), rep("gga", length(r)), rep("ggi", length(r)) ),
                   pipeline=paste(sort(unique(mydata$pipeline)),collapse=", ")
  )
  # saveRDS(pm, file = file.path(wd,"metrics",paste0(id,".rds")))
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the pipeline_metrics database under such id \n"))
  }
  result <- list(metrics=pm, predictions=NA, modeling=mod, metadata=db.params,
                 cleaned=NA, outliers=NA, desire=NA, id=id)
  return(result)#paste("pgg done:",id))
}
