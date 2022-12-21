cleanr <- function(
    pedigreeDTfile= NULL, indiv=NULL, dam=NULL, sire=NULL,
    verbose=FALSE
){

  id <- paste("clr",idGenerator(5,5),sep="")
  type <- "cleaningr"


  ############################
  # loading the dataset
  if (is.null(pedigreeDTfile)) stop("No input pedigree data file specified.")
  if (is.null(indiv)) stop("No column specifying the name of the individual specified.")
  if (is.null(sire)) stop("No column specifying the name of the sire/male specified.")
  if (is.null(dam)) stop("No column specifying the name of the dam/female specified.")

  mydata <- pedigreeDTfile # read.csv(file.path(wd,"files_raw",pedigreeDTfile))

  if(!indiv %in% colnames(mydata) ){
    stop("Please make sure that a column 'geno' with genotype names is included in the file and the rest of the columns is only markers",call. = FALSE)
  }
  if(!sire %in% colnames(mydata) ){
    stop("Please make sure that a column 'geno' with genotype names is included in the file and the rest of the columns is only markers",call. = FALSE)
  }
  if(!dam %in% colnames(mydata) ){
    stop("Please make sure that a column 'geno' with genotype names is included in the file and the rest of the columns is only markers",call. = FALSE)
  }

  ############################
  ## nrm calculation

  ids <- unique(c(mydata[,dam], mydata[,sire], mydata[,indiv]))
  idsDf <- data.frame(idsn = 1:length(ids)); rownames(idsDf) <- ids
  orPed <- unique(mydata[,c(indiv,dam,sire)])

  orPedN <- apply(orPed,2,function(x){ # not sure why is working for the  "" genotypes but is working :)
    idsDf[as.character(x),]
  })
  # which dams and sires are not part of the df
  damsToAdd <- setdiff(orPedN[,dam],orPedN[,indiv])
  siresToAdd <- setdiff(orPedN[,sire],orPedN[,indiv])

  damsToAddDf <- data.frame(damsToAdd,NA,NA); colnames(damsToAddDf) <- colnames(orPedN)
  siresToAddDf <- data.frame(siresToAdd,NA,NA); colnames(siresToAddDf) <- colnames(orPedN)
  orPedN2 <- rbind(damsToAddDf,siresToAddDf,orPedN)

  ped <- pedigreemm::pedigree(sire = orPedN2[,sire],
                  dam  = orPedN2[,dam], label= orPedN2[,indiv])

  A <- as.matrix(pedigreemm::getA(ped))
  # Ai <-pedigreemm::getAInv(ped)

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
    stage = NA
  )
  # write the values used for cleaning to the modeling database
  # write predictions
  # write pipeline metrics

  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  result <- list(metrics=NA, predictions=NA, modeling=NA, metadata=db.params,
                 cleaned=A, outliers=NA, desire=NA, id=id)
  return(result)#
}
