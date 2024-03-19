
desire <- function(trait,h2, G){ 
  ## THIS FUNCTION CREATES THE INPUT FILE FOR THE DESIRE SOFTWARE
  ## THIS IS USED IN THE GENETIC EVALUATION MODULES IN BANAL UNDER THE MULTI TRIAL ANALYSIS MODULE
  if(missing(trait)){stop("trait missing.", call. = FALSE)}
  if(missing(h2)){stop("h2 missing.", call. = FALSE)}
  if(missing(G)){stop("G missing.", call. = FALSE)}
  s1 <- data.frame("  Input file for program Desire."); colnames(s1)<-"m"
  s2 <- data.frame("  Number of traits ...  "); colnames(s2)<-"m"
  s3 <- data.frame(length(trait)); colnames(s3) <- "m"
  s4 <- data.frame("  Trait names ..."); colnames(s4)<-"m"
  s5 <- data.frame(trait); colnames(s5)<-"m"
  s6 <- data.frame("  Starting economic weights ..."); colnames(s6)<-"m"
  s7 <- data.frame(rep(1,length(trait))); colnames(s7)<-"m" # rep(1,length(traitNamesSTGx))
  s8 <- data.frame("  Heritabilities ..."); colnames(s8)<-"m"
  s9 <- data.frame(round(h2,3)); colnames(s9)<-"m"
  s10 <- data.frame("  Standard Deviations ..."); colnames(s10)<-"m"
  s11 <- data.frame(rep(1,length(trait))); colnames(s11)<-"m"
  s12 <- data.frame("  Correlation matrix (rp on upper diagonal, rg on lower)."); colnames(s12)<-"m"
  s13 <- round(G,3)
  
  s.1.12 <- rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
  s.1.12plus <- cbind(s.1.12,matrix("",ncol=ncol(s13) - 1, nrow=nrow(s.1.12)))
  colnames(s.1.12plus) <- colnames(s13)
  desireResult <- rbind(s.1.12plus,s13)
  
  return(desireResult)
}
