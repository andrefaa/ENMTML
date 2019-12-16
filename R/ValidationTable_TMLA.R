#written by Santiago Velazco & Andre Andrade

# Function for summary of models performance
Validation_Table_TMLA <- function(Eval,
                                  Eval_JS,
                                  N,
                                  Thr) {
  #Creates model validation table
  
  # summres <- data.frame(matrix(0, (N*length(Thr)), 6))
  summres <- list()
  for (i in 1:N) {
    #Calculate Evaluation Metrics
    # summres[i,] <- cbind(
    #   AUC = Eval[[i]]@auc,
    #   Kappa = max(Eval[[i]]@kappa),
    #   TSS = sapply(Eval, function(x)
    #     max(x@TPR + x@TNR) - 1),
    #   Jaccard = max(Eval_JS[[i]]$Jaccard),
    #   Sorensen = max(Eval_JS[[i]]$Sorensen),
    #   Fpb = max(Eval_JS[[i]]$Fpb)
    # )-->OLD
    summres[[i]] <- cbind(
      AUC = Eval[[i]]@auc,
      Kappa = Eval[[i]]@kappa,
      TSS = (Eval[[i]]@TPR + Eval[[i]]@TNR) - 1,
      Jaccard = Eval_JS[[i]]$Jaccard,
      Sorensen = Eval_JS[[i]]$Sorensen,
      Fpb = Eval_JS[[i]]$Fpb,
      OR = (1-Eval[[i]]@TPR)
    )
    colnames(summres[[i]]) <-
      c("AUC", "Kappa", "TSS", "Jaccard", "Sorensen", "Fpb","OR")
  }
  summres <- plyr::ldply(summres,data.frame,.id=NULL)
  summres <- cbind(Threshold=Thr,summres)
  if (N != 1) {
    resSD <- aggregate(.~Threshold, data=summres, function(x) stats::sd(x))
    resSD <- resSD[,-1]
    res <- aggregate(.~Threshold, data=summres, mean)
    # resSD <-
    #   data.frame(matrix(round(apply(summres, 2, stats::sd), 3), nrow = 1, ncol = 6))
    # res <-
    #   data.frame(matrix(round(colMeans(summres), 3), nrow = 1, ncol = 6))
    # colnames(res) <- colnames(summres)
    colnames(resSD) <- paste0(colnames(resSD), "_SD")
    res <- cbind(res, resSD)
  } else{
    res <- aggregate(.~Threshold, data=summres, mean)
    # res <- data.frame(matrix(round(colMeans(summres), 3), nrow = 1, ncol = 6))
    # colnames(res) <- colnames(summres)
  }
  return(res)
}
