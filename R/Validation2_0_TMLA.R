Validation2_0 <- function(Eval,
                          Thr,
                          PredPoint,
                          RastPart,
                          N) {
  #Evaluation Metrics 2.0
  Jac <- list()
  OPR <- list()
  Fcpb <- list()
  for (i in 1:N) {
    Pos <- which(Eval[[i]]@t == Thr[i])
    Conf <- Eval[[i]]@confusion[Pos, ]
    Nste <- length(PredPoint$PresAbse)
    P <- sum(PredPoint$PresAbse == 1)
    A <- sum(PredPoint$PresAbse == 0)
    Prev <-
      sum(na.omit((RastPart[[i]] > Thr[i])[])) / length(na.omit(RastPart[[i]][]))
    UPR <- Conf["fn"] / (Conf["tp"] + Conf["fn"])
    Xis <- P / A * (1 - Prev) / Prev
    c <- P / (Prev * A)
    
    #Corrected Jaccard, OPR & Boyce
    FFP <- (1 - UPR) * (Prev * Nste - P)
    FP <- A * ((Conf["fp"] / A) * (N - P) - FFP) / ((1 - Prev) * Nste)
    Jac[[i]] <- Conf["tp"] / (Conf["tp"] + Xis * FP + Conf["fn"])
    OPR[[i]] <- (Xis * FP) / (Conf["tp"] + Xis * FP)
    Fcpb <- 2 * Conf["tp"] / (Conf["fn"] + Conf["tp"] + c * Conf["fp"])
    #Fcpb correction [Scale between 0-1]
    Fcpb[[i]] <- Fcpb / (2 * c / (1 + c))
  }
  
  res <- c(Jac, OPR, Fcpb)
  names(res) <- rep(c("JAC", "OPR", "FCPB"), N)
  return(res)
}
