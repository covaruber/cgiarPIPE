
extractSommer <- function(object,object2,rTermsTrait){
  ## THIS FUNCTION EXTRACTS PREDICTED VALUES, STD ERROR AND RELIABILITIES FROM A SOMMER OBJECT]
  ## THIS FUNCTION IS NOT CURRENTLY IN USE IN ANY MODULE
  result <- list()
  for(ire in 1:length(rTermsTrait)){ # for each random effect
    predictedValue <- object$uList[[ire]]
    stdError <- sqrt(object$uPevList[[ire]])
    vs = object$partitions[[ire]]
    vs1=vs[1,1];vs2=vs[1,2]
    pev=object$Ci[vs1:vs2,vs1:vs2]
    G = solve(object2$Ai[[ire]])*as.numeric(object$theta[[ire]])
    R2 = (G - pev)/G
    result[[ire]] <- list(predictedValue=predictedValue,stdError=stdError,
                           pev=pev, R2=diag(R2))
  }
  names(result) <- rTermsTrait
  return(result)
}
