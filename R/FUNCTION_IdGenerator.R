idGenerator <- function(nn=5, nl=5){
  rm(.Random.seed, envir=globalenv())
  
  # out1 <- apply(data.frame(1:nn),1,function(x){sample(0:9,1)})
  out2 <- apply(data.frame(1:nl),1,function(x){sample(letters,1)})
  # out1 <- paste(out1,collapse = "")
  out2 <- paste(out2,collapse = "")
  # out <- paste(out2,out1,sep="")
  out <- paste( gsub(":",".",gsub("-",".",gsub(" ","-",Sys.time()))) , out2, sep=".")
  return(out)
}
