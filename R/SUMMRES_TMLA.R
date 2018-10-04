#written by Santiago Velazco

# Function for summary of models performance
SUMMRES<-function(Eval, N, Thr){
  TSS <- sapply(Eval, function(x) max(x@TPR + x@TNR) - 1)
  summres <- data.frame(matrix(0, N, 4))
  for (i in 1:N) {
    Pos <- which(Eval[[i]]@t == Thr[i])
    summres[i, ] <- cbind(Eval[[i]]@auc,
                          Eval[[i]]@TPR[Pos],
                          Eval[[i]]@TNR[Pos],
                          Eval[[i]]@kappa[Pos])
  }
  colnames(summres) <- c("AUC", "TPR", "TNR", "Kappa")
  summres <- cbind(THR=as.numeric(Thr), TSS=TSS, summres)
  res <- data.frame(matrix(round(colMeans(summres), 4), 1, 6))
  colnames(res) <- colnames(summres)
  return(res)
}