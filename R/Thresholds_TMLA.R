Thresholds_TMLA <- function(Eval,
                            Eval_JS,
                            sensV){
  
  #Calculates Model Threhsolds
  res <- NULL
  
  #Dismo Thresholds
  ThrDis <- c("kappa","spec_sens","no_omission","sensitivity")
  nom <- c("MAX_KAPPA","MAX_TSS","LPT","SENSITIVITY")
  if(is.null(sensV)){
    ThrDis <- unlist(sapply(ThrDis,function(x) threshold(Eval)[x]))
  }else{
    ThrDis <- unlist(sapply(ThrDis,function(x) threshold(Eval,sensitivity=sensV)[x]))
  }
  names(ThrDis) <- NULL
  PosI <- lapply(as.list(ThrDis), function(x) which(Eval@t%in%x))
  TPR <- sapply(PosI, function(x) Eval@TPR[x])
  names(TPR) <- NULL
  TNR <- sapply(PosI, function(x) Eval@TNR[x])
  names(TNR) <- NULL
  resD <- cbind.data.frame(THR=nom,THR_VALUE=ThrDis,TPR,TNR)
  res <- rbind(res,resD)
  
  #Jaccard Thresholds
  PosII <- which(Eval@t == Eval_JS$JaccardTHR)
  TP_TN <- c(Eval_JS$TPR[PosII],Eval_JS$TNR[PosII])
  resJ <- cbind.data.frame(THR="JACCARD",THR_VALUE=Eval_JS$JaccardTHR,TPR=TP_TN[1],TNR=TP_TN[2])
  res <- rbind(res,resJ)
  
  #Sorensen Thresholds
  PosIII <- which(Eval@t == Eval_JS$SorensenTHR)
  TP_TN <- c(Eval_JS$TPR[PosIII],Eval_JS$TNR[PosIII])
  resS <- cbind.data.frame(THR="SORENSEN",THR_VALUE=Eval_JS$SorensenTHR,TPR=TP_TN[1],TNR=TP_TN[2])
  res <- rbind(res,resS)
  
  return(res)
}