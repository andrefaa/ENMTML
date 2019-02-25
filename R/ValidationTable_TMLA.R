#written by Santiago Velazco & Andre Andrade

# Function for summary of models performance
Validation_Table_TMLA<-function(Eval,
                                Eval_JS,
                                N){
  
  #Creates model validation table
  
  summres <- data.frame(matrix(0, N, 6))
  for(i in 1:N){
    #Calculate Evaluation Metrics
    summres[i, ] <- cbind(AUC=Eval[[i]]@auc,
                          Kappa=max(Eval[[i]]@kappa),
                          TSS=sapply(Eval, function(x) max(x@TPR + x@TNR) - 1),
                          Jaccard=max(Eval_JS[[i]]@Jaccard),
                          Sorensen=max(Eval_JS[[i]]@Sorensen),
                          Fpb=max(Eval_JS[[i]]@Fpb))
    colnames(summres) <- c( "AUC","Kappa", "TSS","Jaccard","Sorensen","Fpb")
  }
  res <- data.frame(matrix(round(colMeans(summres), 3),nrow=1,ncol=6))
  colnames(res) <- colnames(summres)
  return(res)
}
