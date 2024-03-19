mabc <- function(
    markerDTfile= NULL,
    linkFile=NULL,
    idGenotypeLink = NULL, # columns specifying the geotype IDs matching the marker file
    nameGenotypeLink = NULL, # column specifying the genotype name
    crossGroupLink = NULL,
    parentsLink = NULL, # column specifying the recurrent parent
    verbose=FALSE
){
  ## THIS FUNCTION PERFORMS F1 VERIFICATION USING MARKER INFORMATION
  ## THIS IS USED IN THE BCLEAN APP UNDER QA/QC MODULES FOR MARKER DATA
  id <- paste( paste("mbc", cgiarPIPE::idGenerator(5,5),sep=""), markerDTfile$idOriginal, sep = "_")
  type <- "mbc"
  ############################
  # loading the dataset
  if (is.null(markerDTfile)) stop("No input marker data file specified.")
  if (is.null(linkFile)) stop("No link file specified.")
  if (is.null(idGenotypeLink)) stop("No genotype ID column specified. Please do so.")
  if (is.null(nameGenotypeLink)) stop("No genotype name column specified. Please do so.")
  if (is.null(crossGroupLink)) stop("No cross column specified. Please do so.")
  if (is.null(parentsLink)) stop("No parent indicator column specified. Please do so.")
  M0 <- markerDTfile$cleaned$M + 1 
  RA <- markerDTfile$cleaned$ref.alleles
  data <- linkFile$cleaned[,c(idGenotypeLink,nameGenotypeLink,parentsLink,crossGroupLink)]
  data[,idGenotypeLink] <- toupper(data[,idGenotypeLink])
  ############################
  # calculate metrics for MABC or parent recovery
  '%!in%' <- function(x,y)!('%in%'(x,y)) 
  groups <- unique(data[,crossGroupLink])
  predictionsList <- vector(mode="list",length = length(groups))
  for(iGroup in groups){ # iGroup <- groups[1]
    dataG <- data[which(data[,crossGroupLink] == iGroup),] # subset of data
    parentsId <- dataG[which(!is.na(dataG[,parentsLink])),idGenotypeLink ] # id of parents
    parentsName <- dataG[which(!is.na(dataG[,parentsLink])),nameGenotypeLink ]; # name of parents
    progenyId <- dataG[which(is.na(dataG[,parentsLink])),idGenotypeLink ]
    progenyName <- dataG[which(is.na(dataG[,parentsLink])),nameGenotypeLink ]; 
    names(parentsId) <- parentsName ; names(progenyId) <- progenyName 
    unmatchedParents <- which(parentsId %!in% rownames(M0))
    if(length(unmatchedParents) > 0){ # this group cannot be processed
      print(paste("In group",iGroup,"parents", paste(parentsId[unmatchedParents], collapse = ","), "were not found" ))
    }else{ # this group can be processed
      Mpar <- M0[parentsId,, drop=FALSE]
      Mf1 <- apply(Mpar,2,mean)
      # check progeny
      unmatchedProgeny <- which(progenyId %!in% rownames(M0))
      if(length(unmatchedProgeny) > 0){
        print(paste("In group",iGroup,"progeny", paste(progenyId[unmatchedProgeny], collapse = ","), "were not found and will be missing from results" ), call. = FALSE)
      }
      # update progeny IDs and names
      progenyId0 <- progenyId
      progenyId <- intersect(rownames(M0), progenyId )
      progenyName <- names(progenyId0)[which(progenyId0 %in% progenyId)]
      # if Aa in the recurrent parent 
      # hetParentsN <- apply(Mpar,2,function(x){length(which(x == 1))})
      # bothParHet <- which(hetParentsN == 1)
      # if(length(bothParHet) > 0){Mf1[bothParHet] = 3}
      # now calculate match and other parameters
      if(length(progenyId) > 0){ # there was progeny to work with
        Mpro <- M0[progenyId,]
        resultMatch <- vector(mode="list", length = ncol(Mpro))
        for(iMark in 1:ncol(Mpro)){ # iMark=1
          # if(Mf1[iMark] == 3){ # if both parents were heterozygotes
          #   matches <- rep(1,nrow(Mpro)); names(matches) <- rownames(Mpro) # everyone is a match, since all genotypes are possible
          # }else{ # if at least one parent was homozygoues
          match1 <- Mpro[,iMark] - Mf1[iMark] 
          match2 <- Mpro[,iMark] - floor(Mf1[iMark]) # to detect heterozygotes
          match3 <- Mpro[,iMark] - ceiling(Mf1[iMark])  # to detect heterozygotes
          matches <- apply(rbind(match1,match2,match3),2,function(x){ifelse(length(which(x == 0)) > 0,1,0)}) # 0 is no match, 1 is match
          # }
          resultMatch[[iMark]] <- matches
        }
        resultMatchAll <- do.call(cbind,resultMatch)
        # parameters
        percMatch <- apply(resultMatchAll,1,mean)
        hetero <- apply(Mpro,1,function(x){length(which(x == 1))/length(x)})
        heteroF1 <- length(which(Mf1 == 1))/length(Mf1)
        heteroDeviation <- hetero - heteroF1
        percComplete <- apply(Mpro,1,function(x){length(which(!is.na(x)))/length(x)})
        nMarkers <- apply(Mpro,1,function(x){length(which(!is.na(x)))})
        # decision <- ifelse(percMatch > threshold, 1, 0)
        predictionsList[[iGroup]] <- data.frame(analysisId=id, pipeline= paste(sort(unique(linkFile$cleaned$pipeline)),collapse=", "),
                                                trait=c(rep("percentMatch", nrow(Mpro)), rep("nMarkers", nrow(Mpro)), rep("percentComplete", nrow(Mpro)), rep("heterozygocity", nrow(Mpro)) ), 
                                                genoCode=c(progenyId,progenyId,progenyId,progenyId), 
                                                geno=c(progenyName,progenyName,progenyName,progenyName),
                                                mother=c( rep(parentsName[1], nrow(Mpro)), rep(parentsName[1], nrow(Mpro)), 
                                                          rep(parentsName[1], nrow(Mpro)), rep(parentsName[1], nrow(Mpro)) ),
                                                father=c( rep(parentsName[2], nrow(Mpro)), rep(parentsName[2], nrow(Mpro)), 
                                                          rep(parentsName[2], nrow(Mpro)), rep(parentsName[2], nrow(Mpro)) ),
                                                genoType=iGroup, genoYearOrigin=unique(linkFile$cleaned$genoYearOrigin)[1], genoYearTesting=unique(linkFile$cleaned$timePoint)[1],
                                                fieldinst=paste(sort(unique(linkFile$cleaned$fieldinstF)),collapse=", "), 
                                                predictedValue=c(percMatch,nMarkers,percComplete,hetero), 
                                                stdError=NA, 
                                                rel=1, stage= paste(parentsId, collapse=" x ")
        )
      }else{print(paste("No progeny found for cross",iGroup))}
    }
  } # of of iGroup
  predictions <- do.call(rbind, predictionsList); rownames(predictions) <- NULL
  head(predictions)
  #########################################
  ## update databases: write the parameters to the parameter database
  metadata <- data.frame(
    analysisId	= id,  analysisType =	type,
    fieldbooks	= NA,  phenoDataFile =	linkFile$id,
    markerbooks	= NA,  markerDataFile =	markerDTfile$id,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,stage = NA
  )
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the files_cleaned folder under such id \n"))
  }
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","mother","father","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")
  result <- list(metrics=NA, predictions=predictions[,predcols], modeling=NA, metadata=metadata,
                 cleaned=markerDTfile$cleaned, outliers=NA, desire=NA, id=id, idOriginal=linkFile$idOriginal,
                 metadataFieldinst=NA)
  return(result)
}
