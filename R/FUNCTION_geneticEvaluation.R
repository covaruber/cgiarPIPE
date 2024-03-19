geneticEvaluation <- function(fixed=NULL,random=NULL,rcov=NULL,weights=NULL, 
                              AFGE=NULL, varComp=NULL, data=NULL, keep=NULL, tolParInv=1e-5){

  # AFGE <<- AFGE
  # rebuild random part
  yuyu <- strsplit(as.character(random[2]), split = "[+]")[[1]] # random parts
  rtermss <- apply(data.frame(yuyu),1,function(x){ # split random terms
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  rtermssCopy <- rtermss
  for(u in 1:length(rtermss)){
    if(rtermss[u] == "genoF"){ # main effect
      rtermss[u] <- paste("sommer::vsr( ",rtermss[u],", Gu=Ac )")
    }else{ # environmental index interaction for FW
      if(rtermss[u] == "genoF:envIndex" | rtermss[u] == "envIndex:genoF"){
        rtermss[u] <- paste("sommer::vsr( sommer::dsr(envIndex), genoF, Gu=Ac )")
      }else{
        rtermss[u] <- rtermss[u] #paste("vsr( ",rtermss[u]," )")
      }
    }
  }
  random2 <- as.formula(paste("~",paste(rtermss,collapse = "+")))
  # print(random2)
  # rebuild residual part
  if(!is.null(rcov)){
    # print("hey here")
    yuyu <- strsplit(as.character(rcov[2]), split = "[+]")[[1]] # random parts
    rtermss <- apply(data.frame(yuyu),1,function(x){ # split random terms
      strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    for(u in 1:length(rtermss)){
      rtermss[u] <- paste("sommer::vsr( sommer::dsr(",rtermss[u],"), units )")
    }
    rcov2 <- as.formula(paste("~",paste(rtermss,collapse = "+")))
  }

  # ensure AFGE is complete
  genoFlevels <- as.character(unique(data[,"genoF"]))
  inter <- intersect(genoFlevels,colnames(AFGE)) # go for sure
  differ <- setdiff(genoFlevels,inter) # are missing in AFGE matrix
  if(length(differ) > 0){
    A2 <- diag(length(differ))
    colnames(A2) <- rownames(A2) <- differ
    Ac <- sommer::adiag1(AFGE,A2)
  }else{Ac <- AFGE}

  # print(setdiff(levels(data$genoF),colnames(Ac)))
  # print(colnames(Ac))
  # get matrices from sommer
  if(is.null(rcov)){
    mixSommer <- sommer::mmer(fixed=fixed, random=random2, data=data, returnParam = TRUE, dateWarning = FALSE)
  }else{
    mixSommer <- sommer::mmer(fixed=fixed, random=random2, rcov=rcov2, data=data, returnParam = TRUE, dateWarning = FALSE)
  }
  # form H
  # print(str(mixSommer))
  H <- Matrix::Diagonal(x=data[,weights])
  Hs <- chol(H)
  # print(str(Hs))
  # form R
  if(is.null(rcov)){
    ve <- as.list(varComp[grep("residual",varComp$VarComp),"Variance"])
  }else{
    ve <- as.list(varComp[grep("!R",varComp$VarComp),"Variance"])
  }
  veInv <- lapply(ve, function(x){return(1/x)})
  listRinv <- mapply("*", mixSommer$R, veInv, SIMPLIFY=FALSE)
  Rinv <- Reduce("+",listRinv)
  Ri = Hs %*%  Rinv %*% t(Hs);
  # form W
  X <- do.call(cbind,mixSommer$X)
  Z <- do.call(cbind, mixSommer$Z)
  W <- cbind(X,Z)
  # print(str(W))
  # form G
  Gx <- Matrix::Matrix(0,nrow=ncol(X),ncol=ncol(X))
  nonResidual <- setdiff(1:nrow(varComp), c(grep("residual",varComp$VarComp),grep("!R",varComp$VarComp) ) )
  vz <- as.list(varComp[nonResidual,"Variance"])
  listG <- mapply("*", mixSommer$K, vz, SIMPLIFY=FALSE)
  listGinv <-  lapply(listG, function(x){solve(x + diag(tolParInv,ncol(x), ncol(x)))})
  Ginvz <- do.call(Matrix::bdiag, listGinv); listGinv <- NULL
  Gi <- Matrix::bdiag(Gx,Ginvz); Ginvz <- NULL
  # print(str(Gi))
  # form M or LHS
  RiW <- Ri %*% W
  M = Matrix::crossprod(W,RiW)
  M = M + Gi
  # get coefficient matrix
  Ci = solve(M); M <- NULL
  # print("here1")
  # form RHS
  y <- Matrix::Matrix(mixSommer$yvar)
  Riy <- Ri %*% y
  RHS = Matrix::crossprod(W,Riy)
  # print("here2")

  # get solutions
  bu <- Ci %*% RHS
  # Ci <- NULL
  # extract predicted values, standard errors and PEV
  nEffects <- c(ncol(X), unlist(lapply( mixSommer$Z,ncol)))
  end <- numeric()
  for(i in 1:(length(nEffects))){end[i]=sum(nEffects[1:(i)])+1}
  start <- end - nEffects
  end <- start + nEffects - 1
  # remove fixed effects
  start2 <- start[-1]
  end2 <- end[-1]
  resultList <- list()
  geneticTerms <- intersect(1:nrow(varComp), grep("genoF",varComp$VarComp) )
  for(i in geneticTerms){ # i=1
    wgi <- (start2[i]):(end2[i])
    if(rtermssCopy[i] == "genoF"){
      predictedValue <- bu[wgi,] + bu[1,1]
    }else{
      predictedValue <- bu[wgi,]
    }
    pev <- Ci[wgi,wgi]
    stdError <- (sqrt(diag(pev)))
    R2 <- (listG[[i]] - pev)/listG[[i]]
    if(!is.null(keep) & (length(grep("genoF",rtermssCopy[i])) > 0) ){
      keep2 <- (length(predictedValue)-keep+1):length(predictedValue)
      predictedValue <- predictedValue[keep2]
      pev <- pev[keep2,keep2, drop=FALSE]
      stdError <- (sqrt(diag(pev)))
      R2 <- R2[keep2,keep2, drop=FALSE]
    }
    resultList[[rtermssCopy[i]]] <- list(predictedValue=predictedValue, stdError=stdError, pev=pev, R2=diag(R2))
  }
  R2 <- NULL; pev <- NULL; predictedValue <- NULL; stdError <- NULL; Ci <- NULL

  return(resultList)
}
