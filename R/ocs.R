ocs <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    relDTfile= NULL,
    trait= NULL, # per trait
    fieldinst="across",
    nCross=20,
    targetAngle=30, # in radians
    verbose=FALSE,
    wd=NULL
){

  # if(is.null(wd)){wd <- getwd()}
  # md <- strsplit(wd,"/")[[1]]; md <- md[length(md)]
  # if(md != "DB"){stop("Please set your working directory to the DB folder", call. = FALSE)}

  id <- paste("ocs",idGenerator(5,5),sep="")
  type <- "ocs"
  if(is.null(phenoDTfile)){stop("Please provide the predictions", call. = FALSE)}
  if(is.null(relDTfile)){stop("Please provide the markers", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(length(trait) > 1){
    stop(paste0(" Only one trait can be used for optimal contribution. We suggest using an index"), call. = FALSE)
  }

  ############################
  # loading the dataset
  mydata <- phenoDTfile$predictions # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))

  # ava.files <- dir(file.path(wd,"files_cleaned"))
  # if(!paste0(relDTfile) %in% ava.files){stop("relDTfile is not present in the files_cleaned folder. Please check the name of your file and its location",call. = FALSE)}
  myrel <- relDTfile$cleaned #readRDS(file.path(wd,"files_cleaned",paste0(relDTfile)))

  utraits <- unique(mydata$trait)
  if (!trait %in% utraits){
    stop(paste0("'", trait[k], "' is not present in the given dataset"), call. = FALSE)
  }
  mydata <- mydata[which(mydata$trait %in% trait),]
  mydata <- mydata[which(mydata$fieldinst %in% fieldinst),] # make sure only across
  if(nrow(mydata) == 0){stop("Please check the trait and fieldinst selected since there's no phenotypic data for that combination",call. = "FALSE")}
  # make sure you have same phenotypes and genotypes

  common <- intersect(rownames(myrel), mydata[,"geno"])
  myrel <- myrel[common,common]
  mydata <- mydata[which(mydata[,"geno"] %in% common),]

  ############################
  ## ocs analysis
  ebv <- data.frame(mydata[,c("predictedValue")]); rownames(ebv) <- mydata[,"geno"]
  ebv <- data.frame(ebv[rownames(myrel),]);
  crossComb = t(combn(1:nrow(myrel), 2)) # all possible cross combintations
  eMP = (ebv[crossComb[,1],] +     # expected EBVs of all crosses based on
           ebv[crossComb[,2],])/2  # mean parent EBVs
  K <- as.matrix(myrel)
  # OCS: Determine a crossing plan
  plan = selectCrosses(nCross=nCross, # number of crossed to be identified using OCS
                       targetAngle=targetAngle*pi/180, # 30 degrees in radians
                       u=eMP, # expected cross mean EBVs
                       G=K)   # GRM
  dim(plan$crossPlan)

  crossPlan <- as.data.frame(plan$crossPlan) # list of crosses to be made already sorted by best
  crossPlan[ ,1] <- rownames(K)[crossPlan[ ,1]]
  crossPlan[ ,2] <- rownames(K)[crossPlan[ ,2]]
  colnames(crossPlan) <- c("Parent1", "Parent2", "OCS.merit")
  crossPlan

  predictionsBind <- data.frame(analysisId=id, pipeline= paste(sort(unique(mydata$pipeline)),collapse=", "),
                                trait=trait, genoCode=1:nrow(crossPlan), geno=paste(crossPlan[,1],crossPlan[,2], sep=" x "),
                                genoType="predictedCross", genoYearOrigin=1, genoYearTesting=1,
                                fieldinst=fieldinst, predictedValue=crossPlan[,3], stdError=1e-6, rel=1e-6, stage="unknown"
                                )

  #########################################
  ## update databases
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	phenoDTfile,
    markerbooks	= NA,  markerDataFile =	relDTfile,
    year = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = paste(sort(unique(predictionsBind$stage)),collapse=", ")
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
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")
  # saveRDS(predictionsBind[,predcols], file = file.path(wd,"predictions",paste0(id,".rds")))

  # write pipeline metrics
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the predictions database under such id \n"))
  }
  result <- list(metrics=NA, predictions=predictionsBind[,predcols], modeling=mod, metadata=db.params,
                 cleaned=NA, outliers=NA, desire=NA, id=id)
  return(result)#paste("ocs done:",id))
}
