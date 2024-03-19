cleanOutl <- function(
    phenoDTfile= NULL,
    newOutliers=NULL,
    verbose=FALSE
){
  ## THIS FUNCTION TAKES OUTLIERS IDENTIFIED IN THE VISUALIZATION APP AND BINDS THEM TO THE ORIGINAL FILE IN A CLEAN WAY
  ## IT IS USED IN THE OUTLIER CLEANING MODULE IN BCLEAN
  if(is.null(phenoDTfile)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(newOutliers)){stop("Please provide new outliers to be added", call. = FALSE)}
  id <- paste(phenoDTfile$id,unique(newOutliers$traitName), "I",sep="")
  #####################################
  # transformation
  phenoDTfile$outliers <- unique(rbind(phenoDTfile$outliers,newOutliers))
  ##########################################
  ## update databases
  ## write the parameters to the parameter database
  phenoDTfile$id <- id
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  return(phenoDTfile)
}
