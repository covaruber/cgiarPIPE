nrm <- function(
    pedigreeDTfile= NULL,
    verbose=FALSE
){
  ## THIS FUNCTION CALCULATES A NUMERATOR RELATIONSHIP MATRIX
  ## IS USED IN THE BANAL APP UNDER THE STRUCTURE MODULES
  id <- paste( paste("nrm",idGenerator(5,5),sep=""), pedigreeDTfile$idOriginal, sep = "_")
  type <- "nrm"

  ############################
  # loading the dataset
  if (is.null(pedigreeDTfile)) stop("No input marker data file specified.")
  mydata <- pedigreeDTfile$cleaned #readRDS(file.path(wd,"files_cleaned",paste0(pedigreeDTfile)))
  mydata <- unique(mydata)
  for(j in c("sire","dam","indiv")){
    mydata[,j] <- cgiarBase::replaceValues(mydata[,j], Search = "", Replace = NA)
    mydata[,j] <- cgiarBase::replaceValues(mydata[,j], Search = " ", Replace = NA)
    mydata[,j] <- cgiarBase::replaceValues(mydata[,j], Search = "  ", Replace = NA)
  }
  mydata <- unique(mydata)
  # "IR78761-B-SATB2-17-1"
  # the ones that should have NA in dam and sire
  badSire <-  which(is.na(mydata[,2]) & (mydata[,1] == mydata[,3]) )
  if(length(badSire) > 0){ mydata[badSire,3]=NA}
  badDam <- which(is.na(mydata[,3]) & (mydata[,1] == mydata[,2]) )
  if(length(badDam) > 0){ mydata[badDam,2]=NA}
  # the ones who cannot be true
  unrealisticDam <- which(!is.na(mydata[,3]) & (mydata[,1] == mydata[,2]) )
  if(length(unrealisticDam) > 0){ mydata[unrealisticDam,1]=paste0(mydata[unrealisticDam,1],"-2")}
  unrealisticSire <- which(!is.na(mydata[,2]) & (mydata[,1] == mydata[,3]) )
  if(length(unrealisticSire) > 0){ mydata[unrealisticSire,1]=paste0(mydata[unrealisticSire,1],"-2")}
  # the ones that are identical id to mother and father
  unrealisticInd <- which((mydata[,1] == mydata[,2]) & (mydata[,1] == mydata[,3]) )
  if(length(unrealisticInd) > 0){ mydata[unrealisticInd,2]=NA;mydata[unrealisticInd,3]=NA}
  # unique
  mydata <- unique(mydata)
  first <- which(is.na(mydata[,2]) & is.na(mydata[,3])) # no dams and sire
  # second <- which(is.na(mydata[,3]) & !is.na(mydata[,2])) # no sire
  # third <- which(is.na(mydata[,2]) & !is.na(mydata[,3])) # no dam 9313
  foundersName <- na.omit(unique(mydata[first,"indiv"]))
  damsToAddName <- na.omit(setdiff(unique(mydata[,"dam"]),c(unique(mydata[,"indiv"]),foundersName ) )) # dams not present in inviduals or founders
  siresToAddName <- na.omit(setdiff(unique(mydata[,"sire"]),c(damsToAddName, unique(mydata[,"indiv"]), foundersName) ))
  progenyToAddName <- na.omit(setdiff(unique(mydata[,"indiv"]),c(damsToAddName, siresToAddName, foundersName) ))
  ids <- na.omit(unique(c(foundersName,damsToAddName, siresToAddName,progenyToAddName)))
  #
  idsDf <- data.frame(idsn = 1:length(ids)); rownames(idsDf) <- as.character(ids)
  #
  idsDfInverse <- data.frame(ids = ids); rownames(idsDfInverse) <- 1:length(ids)
  # matrix of numerical ids for individuals with information for father, mother or both
  mydataWithPedigreeOnly <- mydata[which(mydata[,"indiv"] %in% progenyToAddName),]
  mydataWithPedigreeOnly <- mydataWithPedigreeOnly[which(!duplicated(mydataWithPedigreeOnly[,"indiv"])),]
  orPedN <- apply(mydataWithPedigreeOnly,2,function(x){ # not sure why is working for the  "" genotypes but is working :)
    idsDf[as.character(x),]
  })
  orPedN <- orPedN[ order(orPedN[,3],orPedN[,2]), ]
  # matrix of numerical ids for individuals without father and mother
  orPedNfounders <- matrix(NA, ncol=3, nrow=length(c(foundersName,damsToAddName, siresToAddName)))
  orPedNfounders[,1] <- 1:length(c(foundersName,damsToAddName, siresToAddName))
  colnames(orPedNfounders) <- colnames(orPedN)
  # bind both
  orPedN2 <- rbind(orPedNfounders,orPedN)
  pede <- data.frame(sire = as.character(orPedN2[,3]),
                     dam  = as.character(orPedN2[,2]),
                     label= as.character(orPedN2[,1])
  )
  pede2<- pedigreemm::editPed(sire=pede$sire, dam= pede$dam, label=pede$label, verbose = FALSE)
  ped<- with(pede2, pedigreemm::pedigree(label=label, sire=sire, dam=dam))

  A <- as.matrix(pedigreemm::getA(ped))
  rownames(A) <- colnames(A) <- as.character(idsDfInverse[rownames(A),])
  A[lower.tri(A)] <- NA

  #########################################
  ## update databases
  ## write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,  analysisType =	type,
    fieldbooks	= NA, phenoDataFile =	pedigreeDTfile$id,
    markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA, stage = NA
  )

  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the files_cleaned folder under such id \n"))
  }
  result <- list(metrics=NA, predictions=NA, modeling=NA, metadata=metadata,
                 cleaned=A, outliers=NA, desire=NA, id=id, idOriginal=pedigreeDTfile$idOriginal,
                 metadataFieldinst=NA)
  return(result)
}
