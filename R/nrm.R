nrm <- function(
    pedigreeDTfile= NULL,
    verbose=FALSE
){

  # if(is.null(wd)){wd <- getwd()}
  # md <- strsplit(wd,"/")[[1]]; md <- md[length(md)]
  # if(md != "DB"){stop("Please set your working directory to the DB folder", call. = FALSE)}

  id <- paste("nrm",idGenerator(5,5),sep="")
  type <- "nrm"

  ############################
  # loading the dataset
  if (is.null(pedigreeDTfile)) stop("No input marker data file specified.")
  mydata <- pedigreeDTfile$cleaned #readRDS(file.path(wd,"files_cleaned",paste0(pedigreeDTfile)))

  ids <- unique(c(mydata[,"dam"], mydata[,"sire"], mydata[,"indiv"]))
  idsDf <- data.frame(idsn = 1:length(ids)); rownames(idsDf) <- ids
  idsDfInverse <- data.frame(ids = ids); rownames(idsDf) <- 1:length(ids)
  orPedN <- apply(mydata,2,function(x){ # not sure why is working for the  "" genotypes but is working :)
    idsDf[as.character(x),]
  })
  # which dams and sires are not part of the df
  damsToAdd <- setdiff(orPedN[,"dam"],orPedN[,"indiv"])
  siresToAdd <- setdiff(orPedN[,"sire"],orPedN[,"indiv"])

  damsToAddDf <- data.frame(damsToAdd,NA,NA); colnames(damsToAddDf) <- colnames(orPedN)
  siresToAddDf <- data.frame(siresToAdd,NA,NA); colnames(siresToAddDf) <- colnames(orPedN)
  orPedN2 <- rbind(damsToAddDf,siresToAddDf,orPedN)

  ped <- pedigreemm::pedigree(sire = orPedN2[,"sire"],
                              dam  = orPedN2[,"dam"], label= orPedN2[,"indiv"])

  A <- as.matrix(pedigreemm::getA(ped))

  rownames(A) <- colnames(A) <- idsDfInverse[rownames(A),]
  # Ai <-pedigreemm::getAInv(ped)

  #########################################
  ## update databases
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	pedigreeDTfile$id,
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
    cat(paste("Your results will be available in the files_cleaned folder under such id \n"))
  }
  result <- list(metrics=NA, predictions=NA, modeling=NA, metadata=db.params,
                 cleaned=A, outliers=NA, desire=NA, id=id)
  return(result)#paste("grm done:",id))
}
