traitNames <- function(id){
  md <- strsplit(getwd(),"/")[[1]]; md <- md[length(md)]
  if(md != "DB"){stop("Please set your working directory to the DB folder", call. = FALSE)}
  mydata <- readRDS(file.path(wd,"predictions",paste0(id)))
  res <- unique(mydata$trait)
  return(res)
}

preds <- function(id){
  md <- strsplit(getwd(),"/")[[1]]; md <- md[length(md)]
  if(md != "DB"){stop("Please set your working directory to the DB folder", call. = FALSE)}
  mydata <- readRDS(file.path(wd,"predictions",paste0(id)))
  return(mydata)
}

metrics <- function (id){
  md <- strsplit(getwd(), "/")[[1]]
  md <- md[length(md)]
  if (md != "DB") {
    stop("Please set your working directory to the DB folder", 
         call. = FALSE)
  }
  mydata <- readRDS(file.path(wd,"predictions",paste0(id)))
  return(mydata)
}

traitMeans <- function(id){
  mydata <- readRDS(file.path(wd,"predictions",paste0(id)))
  trait <- unique(mydata$trait)
  mm <- vv <- numeric()
  for(iTrait in 1:length(trait)){
    sub <- mydata[which(mydata$trait == trait[iTrait]),]
    mm[iTrait] <- mean(sub$predictedValue, na.rm=TRUE)
  }
  names(mm) <- trait
  return(mm)
}