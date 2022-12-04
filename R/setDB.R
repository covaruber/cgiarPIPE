setDB <- function(path=NULL, refreshDB=FALSE, refreshDirCleaned=FALSE, refreshDirRaw=FALSE,
                  refreshDirPredictions=FALSE, refreshDirMetrics=FALSE, refreshDirModeling=FALSE,
                  refreshDirOutliers=FALSE, refreshDirMetadata=FALSE){
  
  if(refreshDB){
    if(is.null(path)){stop("Please provide the path where the DB folder can be found", call. = FALSE)}
    refreshDirRaw=TRUE
    refreshDirCleaned=TRUE
    refreshDirPredictions=TRUE
    refreshDirMetrics=TRUE
    refreshDirModeling=TRUE
    refreshDirOutliers=TRUE
    refreshDirMetadata=TRUE
    cat("Databases have been reset. \n")
  }
  if(refreshDirRaw){
    if(is.null(path)){stop("Please provide the path where the DB folder can be found", call. = FALSE)}
    # save initial analysis DBs
    files <- dir(file.path(path,"DB","files_raw"))
    files <- setdiff(files,"example1.csv")
    unlink(file.path(path,"DB","files_raw", files))
    cat("Directory with raw files has been cleaned \n")
  }
  if(refreshDirCleaned){
    if(is.null(path)){stop("Please provide the path where the DB folder can be found", call. = FALSE)}
    # save initial analysis DBs
    files <- dir(file.path(path,"DB","files_cleaned"))
    files <- setdiff(files,"example1.rds")
    unlink(file.path(path,"DB","files_cleaned", files))
    cat("Directory with cleaned files has been cleaned \n")
  }
  if(refreshDirPredictions){
    if(is.null(path)){stop("Please provide the path where the DB folder can be found", call. = FALSE)}
    # save initial analysis DBs
    files <- dir(file.path(path,"DB","predictions"))
    files <- setdiff(files,"example1.rds")
    unlink(file.path(path,"DB","predictions", files))
    cat("Directory with prediction files has been cleaned \n")
  }
  if(refreshDirMetrics){
    if(is.null(path)){stop("Please provide the path where the DB folder can be found", call. = FALSE)}
    # save initial analysis DBs
    files <- dir(file.path(path,"DB","metrics"))
    files <- setdiff(files,"example1.rds")
    unlink(file.path(path,"DB","metrics", files))
    cat("Directory with metrics files has been cleaned \n")
  }
  if(refreshDirModeling){
    if(is.null(path)){stop("Please provide the path where the DB folder can be found", call. = FALSE)}
    # save initial analysis DBs
    files <- dir(file.path(path,"DB","modeling"))
    files <- setdiff(files,"example1.rds")
    unlink(file.path(path,"DB","modeling", files))
    cat("Directory with modeling files has been cleaned \n")
  }
  if(refreshDirOutliers){
    if(is.null(path)){stop("Please provide the path where the DB folder can be found", call. = FALSE)}
    # save initial analysis DBs
    files <- dir(file.path(path,"DB","outliers"))
    files <- setdiff(files,"example1.rds")
    unlink(file.path(path,"DB","outliers", files))
    cat("Directory with outlier files has been cleaned \n")
  }
  if(refreshDirMetadata){
    if(is.null(path)){stop("Please provide the path where the DB folder can be found", call. = FALSE)}
    # save initial analysis DBs
    files <- dir(file.path(path,"DB","metadata"))
    files <- setdiff(files,"example1.rds")
    unlink(file.path(path,"DB","metadata", files))
    cat("Directory with metadata files has been cleaned \n")
  }
  if(!is.null(path)){ # if user provides a path then create the DB
    if(dir.exists(file.path(path,"DB"))){
      cat("Database folder 'DB' already exists please delete it first if you want to create it again.")
    }else{
      dir.create(file.path(path,"DB")) # create the root database
      dir.create(file.path(path,"DB","desire")) #
      dir.create(file.path(path,"DB","files_raw")) # create raw files folder
      dir.create(file.path(path,"DB","files_cleaned")) # create cleaned files folder
      dir.create(file.path(path,"DB","predictions")) # create cleaned files folder
      dir.create(file.path(path,"DB","metrics")) # create cleaned files folder
      dir.create(file.path(path,"DB","modeling")) # create cleaned files folder
      dir.create(file.path(path,"DB","outliers")) # create cleaned files folder
      dir.create(file.path(path,"DB","metadata")) # create cleaned files folder
      
      data("example1")
      write.csv(example1,file=file.path(path,"DB","files_raw","example1.csv"), row.names = FALSE)
      saveRDS(example1,file=file.path(path,"DB","files_cleaned","example1.rds"))
      
      cat("Database and all subdirectories have been created. \n")
    }
  }
  
}

