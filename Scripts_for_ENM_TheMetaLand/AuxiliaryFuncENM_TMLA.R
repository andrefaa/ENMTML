# Standarize raster values between 0 and 1------
STANDAR <- function(x){
  result <- (x-cellStats(x, min))/(cellStats(x, max)-cellStats(x, min))
  return(result)
}

# Function for summary of models performance------
SUMMRES<-function(Eval, N, Thr){
  res <- data.frame(matrix(0, N, 6))
    for(h in 1:length(Thr)){
    summres <- data.frame(matrix(0, N, 5))
      for (i in 1:N) {
        summres[i, ] <- cbind(Eval[[i]]@auc,
                            Eval[[i]]@TPR[which(Eval[[i]]@t == Thr[h])],
                            Eval[[i]]@TNR[which(Eval[[i]]@t == Thr[h])],
                            Eval[[i]]@kappa[which(Eval[[i]]@t == Thr[h])],
                            ((Eval[[i]]@TPR[which(Eval[[i]]@t == Thr[h])] + Eval[[i]]@TNR[which(Eval[[i]]@t == Thr[h])]) - 1))
      }
    summres <- cbind(Thr[h], summres)
    summres <- data.frame(matrix(round(colMeans(summres), 5), 1, 6))
    res[h,] <- summres
    }
  Nom <- NULL
  if("no_omission"%in%names(Thr)){
    Nom <- c(Nom,"LPT")
  }
  if("spec_sens"%in%names(Thr)){
    Nom <- c(Nom,"MAX")
  }
  res <- cbind(Nom,res)
  colnames(res) <- c("TYPE","THR","AUC", "TPR", "TNR", "Kappa","TSS")
  return(res)
}

SUMMRES_ENS<-function(Eval, Thr){
  res <- data.frame(matrix(0, length(Thr), 5))
  for(h in 1:length(Thr)){
    res[h,] <- cbind(Eval@auc,
                        Eval@TPR[which(Eval@t == Thr[[h]])],
                        Eval@TNR[which(Eval@t == Thr[[h]])],
                        Eval@kappa[which(Eval@t == Thr[[h]])],
                        ((Eval@TPR[which(Eval@t == Thr[[h]])] + Eval@TNR[which(Eval@t == Thr[[h]])]) - 1))
    }
    colnames(res) <- c("AUC", "TPR", "TNR", "Kappa","TSS")
    TYPE <- NULL
    if("no_omission"%in%names(Thr)){
      TYPE <- c(TYPE,"LPT")
    }
    if("spec_sens"%in%names(Thr)){
      TYPE <- c(TYPE,"MAX")
    }
    THR <- unlist(Thr)
    names(THR) <- NULL
    res <- cbind(TYPE, THR, res)
  return(res)
}



# Prediction of different algorithm------
PREDICT <- function(Variables, Models_List){
ListRaster <- as.list(names(Models_List))
names(ListRaster) <- names(Models_List)
Algorithm <- names(Models_List)

if(any(Algorithm=="BIOCLIM")==TRUE){
  ListRaster[["BIOCLIM"]] <- round(predict(Variables, Models_List[["BIO"]]),4)
}
if(any(Algorithm=="MAXENTD")==TRUE){
  ListRaster[["MAXENTD"]] <- round(predict(Variables, Models_List[["MAXENTD"]], args='outputformat=cloglog'), 4)
}
if(any(Algorithm=="MAXENTS")==TRUE){
  ListRaster[["MAXENTS"]] <- round(predict(Variables, Models_List[["MAXENTS"]], args='outputformat=cloglog'), 4)
}
if(any(Algorithm=="MAXENTD_NEW")==TRUE){
  ListRaster[["MAXENTD_NEW"]] <- round(predict(Variables, Models_List[["MAXENTD_NEW"]], clamp=F, type="cloglog"), 4)
}
if(any(Algorithm=="MAXENTS_NEW")==TRUE){
  ListRaster[["MAXENTS_NEW"]] <- round(predict(Variables, Models_List[["MAXENTS_NEW"]], clamp=F, type="cloglog"), 4)
}
if(any(Algorithm=="MAXLIKE")==TRUE){
  ListRaster[["MAXLIKE"]] <- round(predict(Variables, Models_List[["MAXLIKE"]]), 4)
}
if(any(Algorithm=="SVM")==TRUE){
  ListRaster[["SVM"]] <- round(predict(Variables, Models_List[["SVM"]]), 4)
}
if(any(Algorithm=="RF")==TRUE){
  ListRaster[["RF"]] <- round(predict(Variables, Models_List[["RF"]]), 4)
}
if(any(Algorithm=="GLMC")==TRUE){
  ListRaster[["GLMC"]] <- round(predict(Variables, Models_List[["GLMC"]]), 4)
}
if(any(Algorithm=="GAMC")==TRUE){
  ListRaster[["GAMC"]] <- round(predict(Variables, Models_List[["GAMC"]]), 4)
}
if(any(Algorithm=="GAUSS")==TRUE){
  Pred <- predict.graf.raster(Models_List[["GAUSS"]], Variables, type = "response", 
                              CI = 0.95, maxn = NULL)
  ListRaster[["GAUSS"]] <- round(Pred$posterior.mode, 4)
}
return(ListRaster)
}

# Predict Raster wiht GRaF package GAUSS algorithm
predict.graf.raster <- function (object, x, type, CI, maxn, ...) {
    # the following adapted from dismo:::predict
    variables <- colnames(object$x)
    if (is.null(CI)) {
      # if CIs aren't required
      out <- raster(x)
    } else {
      if (CI == 'std') {
        # if the SD is needed (latent case)
        out <- brick(x,
                     nl = 2)
      } else {
        # if proper CIs are required
        out <- brick(x,
                     nl = 3)
      }
    }
    ncols <- ncol(out)
    
    firstrow <- 1
    firstcol <- 1
    
    if (!canProcessInMemory(out, 3)) {
      
      filename <- rasterTmpFile()
      out <- writeStart(out, filename = filename, ...)
      inMemory <- FALSE
    } else {
      v <- array(NA, dim = c(nrow(out), ncol(out), nlayers(out)))
      inMemory <- TRUE
    }
    
    tr <- blockSize(out, n = nlayers(x) + 2)
    pb <- pbCreate(tr$n)
    
    for (i in 1:tr$n) {
      rr <- firstrow + tr$row[i] - 1
      
      rowvals <- getValuesBlock(x,
                                row = rr,
                                nrows = tr$nrows[i], 
                                firstcol,
                                ncols)
      rowvals <- rowvals[, variables, drop = FALSE]
      
      res <- matrix(NA,
                    nrow = nrow(rowvals),
                    ncol = nlayers(out))
      rowvals <- na.omit(rowvals)
      if (length(rowvals) > 0) {
        newdata <- data.frame(rowvals)
        # loop through and convert factors to factors
        for (col in ncol(newdata)) {
          if (is.factor(object$obsx[, col])) {
            newdata[, col] <- factor(newdata[, col])
          }
        }
        
        # predict to this batch of pixels
        p <- predict.graf(object,
                          newdata = newdata,
                          type = type,
                          CI = CI,
                          maxn = maxn,
                          ...)
        naind <- as.vector(attr(rowvals, "na.action"))
        if (!is.null(naind)) {
          res[-naind, ] <- p
        } else {
          res <- p
        }
      }
      
      if (inMemory) {
        res <- array(res, dim = c(ncol(out), tr$nrows[i], nlayers(out)))
        res <- aperm(res, c(2, 1, 3))
        rows <- tr$row[i]:(tr$row[i] + tr$nrows[i] - 1)
        v[rows, , ] <- res
      } else {
        out <- writeValues(out, res, tr$row[i])
      }
      pbStep(pb, i)
    }
    pbClose(pb)
    
    if (inMemory) {
      # permute v so that it gets vectorized in the right direction
      v <- aperm(v, c(2, 1, 3))
      for (i in 1:nlayers(out)) {
        out <- setValues(out,
                         as.vector(v[, , i]),
                         layer = i)
      }
    } else {
      out <- writeStop(out)
    }
    names(out) <- colnames(p)
    return (out)
}

hingeval <-
  
  function(x, min, max)
    
  {
    
    pmin(1, pmax(0, (x-min)/(max-min)))
    
  }
