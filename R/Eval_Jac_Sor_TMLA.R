#written by Santiago Velazco & Andre Andrade

# Function for ENM evaluation by Jaccard, Sorensen & Fpb
Eval_Jac_Sor_TMLA <- function(p,a){
  #Evaluate ENMs using Jaccard and Sorensen Indexes
  
  #Parameters:
    #p:presence points suitability
    #a:absence points suitability
  #Initialisation
  np <- length(p)
  na <- length(a)
  if (na == 0 | np == 0) {
    stop('cannot evaluate a model without absence and presence data that are not NA')
  }
  #Threshold breaks:
  if (length(p) > 1000) {
    tr <- as.vector(quantile(p, 0:1000/1000))
  } else {
    tr <- p
  }
  if (length(a) > 1000) {
    tr <- c(tr, as.vector(quantile(a, 0:1000/1000)))
  } else {
    tr <- c(tr, a)
  }
  tr <- sort(unique( round(tr, 8)))
  tr <- c( tr - 0.0001, tr[length(tr)] + c(0, 0.0001))
  
  res <- matrix(ncol=4, nrow=length(tr))
  colnames(res) <- c('tp', 'fp', 'fn', 'tn')
  #Confusion Matrix
  for (i in 1:length(tr)) {
    res[i,1] <- length(p[p>=tr[i]])  # a  true positives
    res[i,2] <- length(a[a>=tr[i]])  # b  false positives
    res[i,3] <- length(p[p<tr[i]])    # c  false negatives
    res[i,4] <- length(a[a<tr[i]])    # d  true negatives
  }
  a = res[,1]
  b = res[,2]
  c = res[,3]
  d = res[,4]

  #Sorensen Index
  SOR <- 2*a/(c + 2*a + b)
  SorTHR <- tr[which(SOR==max(SOR))][1]
  
  #Jaccard Index
  JAC <- a/(c + a + b)
  JacTHR <- tr[which(JAC==max(JAC))][1]
  
  #Fpb
  Fpb <- 2*JAC
  FpbTHR <- tr[which(Fpb==max(Fpb))][1]
  
  #Final Result Object
  xc <- list()
  xc$presences <- np
  xc$absences <- na
  xc$Sorensen <- SOR
  xc$SorensenTHR <- SorTHR
  xc$Jaccard <- JAC
  xc$JaccardTHR <- JacTHR
  xc$Fpb <- Fpb
  xc$t <- tr
  xc$TPR <- a / (a + c)
  xc$TNR <- d / (b + d)
  
  #Return final result
  return(xc)
}
