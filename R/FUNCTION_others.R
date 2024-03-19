corImputation <- function(wide, Gu=NULL, nearest=2, roundR=FALSE){
  if(is.null(Gu)){
    X <- apply(wide, 2, sommer::imputev)
    Gu <- cor(t(X))
  }
  wide2 <- wide
  rowNamesWide <-  rownames(wide)
  # for each feature
  for(iEnv in 1:ncol(wide)){ # iEnv=10
    withData <- which(!is.na(wide[,iEnv]))
    withoutData <- which(is.na(wide[,iEnv]))
    toPredict <- 1:nrow(wide)
    # if(length(toPredict) > 1){
    likelihood=Gu[as.character(rowNamesWide)[toPredict],as.character(rowNamesWide)[withData], drop=FALSE]
    # }else{
    #   likelihood=matrix(Gu[as.character(rowNamesWide)[toPredict],as.character(rowNamesWide)[withData]], nrow=1)
    #   rownames(likelihood) <- as.character(rowNamesWide)[toPredict]
    #   colnames(likelihood) <- as.character(rowNamesWide)[withData]
    # }
    replacement <- numeric()
    for(iInd in 1:nrow(likelihood)){ # iInd=1
      # averaging only the positively correlated to avoid decrease   # wide[iInd, iEnv]
      indLik <- sort(abs(likelihood[iInd,]), decreasing = TRUE)
      toAverage <- indLik[1:min(c(nearest,length(withData)))]
      indLikToAverage <- likelihood[iInd,names(toAverage)]
      replacement[iInd] <- mean(wide[names(which(indLikToAverage > 0)),iEnv]) 
    }
    names(replacement) <- rownames(likelihood)
    # time to replace the missing data
    dd=data.frame(replacement=replacement, index=1:length(replacement),
                  imputed=1, id=rownames(likelihood),orVal=wide[,iEnv])
    dd$imputed[which(dd$id %in% names(withoutData))]=0
    head(dd)
    if(roundR){
      wide2[toPredict,iEnv] <- round(replacement)
    }else{
      wide2[toPredict,iEnv] <- replacement
    }
  }
  # for each individual
  for(jRow in 1:nrow(wide)){
    miss <- which(is.na(wide[jRow,]))
    if(length(miss) > 0){
      dd=data.frame(full=as.vector(unlist(wide2[jRow,])), partial=as.vector(unlist(wide[jRow,])))
      model <- RcppArmadillo::fastLm(partial~full,data=dd[which(!is.na(dd$partial)),])
      # model <- lm(partial~full,data=dd[which(!is.na(dd$partial)),])
      pp=model$coefficients[1]+(dd[which(is.na(dd$partial)),"full"]*model$coefficients[2])
      if(roundR){
        wide[jRow,miss] <- round(pp)#round(predict(model,newdata = dd[which(is.na(dd$partial)),]))
      }else{
        wide[jRow,miss] <- pp#predict(model,newdata = dd[which(is.na(dd$partial)),])
      }
    }
  }
  stillEmpty <- which(is.na(wide), arr.ind = TRUE)
  if(nrow(stillEmpty) > 0){wide[stillEmpty] <- mean(wide, na.rm=TRUE)}
  stillEmpty <- which(is.na(wide2), arr.ind = TRUE)
  if(nrow(stillEmpty) > 0){wide2[stillEmpty] <- mean(wide2, na.rm=TRUE)}
  
  return(list(imputed=wide, corImputed=wide2))
}

neMarker <- function(M, maxNe=100, maxMarker=1000, nSamples=5){
  # function
  v <- sample(1:ncol(M), min(c(maxMarker, ncol(M))))
  M <- M[,v]
  nAllelesPop <- apply(M,2,function(x){round((length(table(x))/3) + 1)})
  # plot(nAllelesPop)
  nAllelesPopTotal <- sum(nAllelesPop)
  maxNe <- min(c(maxNe, nrow(M)))
  
  allelesCovered <- allelesCoveredSe <- vector(mode="numeric", length = maxNe)
  for(i in 2:maxNe){
    print(i)
    allelesCoveredSample <- vector(mode="numeric", length = nSamples)
    for(j in 1:nSamples){
      ii <- sample(1:nrow(M),i)
      nAllelesPopI <- apply(M[ii,],2,function(x){round((length(table(x))/3) + 1)})
      allelesCoveredSample[j] <- sum(nAllelesPopI)
    }
    allelesCovered[i] <- mean(allelesCoveredSample)/nAllelesPopTotal
    allelesCoveredSe[i] <- ( sd(allelesCoveredSample)/nSamples ) /nAllelesPopTotal
  }
  result <- data.frame(allelesCovered=allelesCovered, allelesCoveredSe=allelesCoveredSe, Ne=1:maxNe)
  return(result)
}

##
build.HMM <- function(M1,M2, custom.hyb=NULL, return.combos.only=FALSE, separator=":"){
  # build hybrid marker matrix
  
  if(!is.null(custom.hyb)){
    pheno <- custom.hyb
    found <- length(which(colnames(pheno) %in% c("Var1","Var2","hybrid")))
    if(found != 3){
      stop("Column names Var1, Var2, hybrid need to be present when you provide \n       a data table to customize the hybrid genotypes to be build.\n", call. = FALSE)
    }
    return.combos.only=FALSE
  }else{
    a <- rownames(M1)
    b <- rownames(M2)
    pheno <- expand.grid(a,b)
    pheno <- pheno[!duplicated(t(apply(pheno, 1, sort))),]
    pheno$hybrid <- paste(pheno$Var1, pheno$Var2, sep=separator)
  }
  
  if(!return.combos.only){
    # check that marker matrices are in -1,0,1 format
    checkM1 <- c(length(which(M1 == -1)),length(which(M1 == 1)),length(which(M1 == 2)))
    checkM2 <- c(length(which(M2 == -1)),length(which(M2 == 1)),length(which(M2 == 2)))
    
    checkM1[which(checkM1 > 0)] <- 1
    checkM2[which(checkM2 > 0)] <- 1
    
    if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
    }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
      cat("Either -1 or 1 alleles not detected in M1, we assume you have coded homozygotes \n       as 0 and 1 instead of -1 and 1. We'll fix it.\n")
    }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
      cat("Either -1 or 1 alleles not detected in M1, we assume you have coded homozygotes \n       as 0 and 2 instead of -1 and 1. We'll fix it.\n")
    }
    
    if(all(checkM2 == c(1,1,0))){ # homo markers were coded correctly as -1,1
    }else if(all(checkM2 == c(0,1,0)) | all(checkM2 == c(1,0,0))){ # homo markers were coded as 0 1
      cat("Either -1 or 1 alleles not detected in M2, we assume you have coded homozygotes \n       as 0 and 1 instead of -1 and 1. We'll fix it.\n")
    }else if(all(checkM2 == c(0,0,1))){ # homo markers were coded as 0 2
      cat("Either -1 or 1 alleles not detected in M2, we assume you have coded homozygotes \n       as 0 and 2 instead of -1 and 1. We'll fix it.\n")
    }
    
    
    
    ## add markers coming from parents M1
    Z1 <- Matrix::sparse.model.matrix(~Var1-1,pheno);dim(Z1); 
    colnames(Z1) <- gsub("Var1","",colnames(Z1))
    M1 <- M1[colnames(Z1),]
    #M1[1:4,1:4]; Z1[1:4,1:4]; 
    ## add markers coming from parents M2
    Z2 <- Matrix::sparse.model.matrix(~Var2-1,pheno);dim(Z2); 
    colnames(Z2) <- gsub("Var2","",colnames(Z2))
    M2 <- M2[colnames(Z2),]
    #M2[1:4,1:4]; Z2[1:4,1:4];  
    
    ## create the 
    # Z3 <- model.matrix(~hybrid-1,pheno);dim(Z3);
    # colnames(Z3) <- gsub("hybrid","",colnames(Z3))
    # hyb.names <- colnames(Z3)[as.vector(apply(Z3,1,function(x){which(x==1)}))] # names of hybrids
    hyb.names <- pheno$hybrid
    ## marker matrix for hybrids one for each parent
    cat(paste("Building hybrid marker matrix for",nrow(Z1),"hybrids\n"))
    
    # M1 <- as(M1, Class="dgCMatrix")
    # M2 <- as(M2, Class="dgCMatrix")
    # Z1 <- as(Z1, Class="dgCMatrix")
    # Z2 <- as(Z2, Class="dgCMatrix")
    
    cat("Extracting M1 contribution\n")
    if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
      Md <- Z1 %*% M1;  # was already converted to -1,1
    }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
      Md <- 2*Z1 %*% M1 - 1;  # 2*Z.dent %*% M.dent - 1   # convert to -1,1
    }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
      Md <- Z1 %*% M1 - 1;  # Z.dent %*% M.dent - 1   # convert to -1,1
    }
    
    cat("Extracting M2 contribution\n")
    if(all(checkM2 == c(1,1,0))){ # homo markers were coded correctly as -1,1
      Mf <- Z2 %*% M2;  # was already converted to -1,1
    }else if(all(checkM2 == c(0,1,0)) | all(checkM2 == c(1,0,0))){ # homo markers were coded as 0 1
      Mf <- 2*Z2 %*% M2 - 1;  # 2*Z.dent %*% M.dent - 1   # convert to -1,1
    }else if(all(checkM2 == c(0,0,1))){ # homo markers were coded as 0 2
      Mf <- Z2 %*% M2 - 1;  # Z.dent %*% M.dent - 1   # convert to -1,1
    }
    
    ## marker matrix coded as additive -1,0,1
    Mdf <- (Md + Mf)*(1/2) # normal marker matrix for the hybrids
    rownames(Mdf) <- hyb.names
    #hist(Mdf)
    
    ## dominance matrix for hybrids (0,1 coded)
    Delta <- 1/2*(1 - Md * Mf) #performs element wise multiplication = Hadamard product
    rownames(Delta) <- hyb.names
    #hist(Delta)
    cat("Done!!\n")
    return(list(HMM.add=Mdf, HMM.dom=Delta, data.used=pheno))
    
  }else{
    return(list(HMM.add=NA, HMM.dom=NA, data.used=pheno))
  }
}

redmm <- function (x, M = NULL, Lam=NULL, nPC=50, cholD=FALSE, returnLam=FALSE) {
  
  if(system.file(package = "RSpectra") == ""){
    stop("Please install the RSpectra package to use the redmm() function.",call. = FALSE)
  }else{
    requireNamespace("RSpectra",quietly=TRUE)
  }
  
  if(is.null(M)){
    stop("M cannot be NULL. We need a matrix of features that defines the levels of x")
  }else{
    
    if (inherits(x, "dgCMatrix") | inherits(x, "matrix")) {
      notPresentInM <- setdiff(colnames(Z),rownames(M))
      notPresentInZ <- setdiff(rownames(M),colnames(x))
    }else{
      notPresentInM <- setdiff(unique(x),rownames(M))
      notPresentInZ <- setdiff(rownames(M),unique(x))
    }
    if(is.null(Lam)){ # user didn't provide a Lambda matrix
      if(nPC == 0){ # user wants to use the full marker matrix
        Lam <- Lam0 <- M
      }else{ # user wants to use the PCA method
        nPC <- min(c(nPC, ncol(M)))
        if(cholD){
          smd <- try(chol(M) , silent = TRUE)
          if(inherits(smd, "try-error")){smd <- try(chol((M+diag(1e-5,nrow(M),nrow(M))) ) , silent = TRUE)}
          Lam0 = t(smd)
        }else{
          smd <- RSpectra::svds(M, k=nPC, which = "LM")
          Lam0 <- smd$u
        }
        Lam = Lam0[,1:min(c(nPC,ncol(M))), drop=FALSE]
        rownames(Lam) <- rownames(M)
        colnames(Lam) <- paste0("nPC",1:nPC)
      }
    }else{ # user provided it's own Lambda matrix
      Lam0=Lam
      Lam = Lam0[,1:min(c(nPC,ncol(M))), drop=FALSE]
      rownames(Lam) <- rownames(M)
      colnames(Lam) <- paste0("nPC",1:nPC)
    }
  }
  if (inherits(x, "dgCMatrix") | inherits(x, "matrix")) {
    Z <- x
  }else{
    if (!is.character(x) & !is.factor(x)) {
      namess <- as.character(substitute(list(x)))[-1L]
      Z <- Matrix(x, ncol = 1)
      colnames(Z) <- namess
    }else {
      dummy <- x
      levs <- na.omit(unique(dummy))
      if (length(levs) > 1) {
        Z <- Matrix::sparse.model.matrix(~dummy - 1, na.action = na.pass)
        colnames(Z) <- gsub("dummy", "", colnames(Z))
      } else {
        vv <- which(!is.na(dummy))
        Z <- Matrix(0, nrow = length(dummy))
        Z[vv, ] <- 1
        colnames(Z) <- levs
      }
    }
  }
  
  Zstar <- as.matrix(Z %*% Lam[colnames(Z),])
  if(returnLam){
    return(list(Z = Zstar, Lam=Lam, Lam0=Lam0)) 
  }else{return(Zstar)}
  
}