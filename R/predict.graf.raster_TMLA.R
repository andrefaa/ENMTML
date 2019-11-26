#source https://github.com/goldingn/GRaF/blob/master/R/predict.graf.raster.R

# Predict Raster wiht GRaF package GAUSS algorithm
predict.graf.raster <- function (object, x, type, CI, maxn, ...) {
  # the following adapted from dismo:::predict
  variables <- colnames(object$x)
  if (is.null(CI)) {
    # if CIs aren't required
    out <- raster::raster(x)
  } else {
    if (CI == 'std') {
      # if the SD is needed (latent case)
      out <- raster::brick(x,
                           nl = 2)
    } else {
      # if proper CIs are required
      out <- raster::brick(x,
                           nl = 3)
    }
  }
  ncols <- ncol(out)
  
  firstrow <- 1
  firstcol <- 1
  
  if (!canProcessInMemory(out, 3)) {
    filename <- raster::rasterTmpFile()
    out <- raster::writeStart(out, filename = filename, ...)
    inMemory <- FALSE
  } else {
    v <- array(NA, dim = c(nrow(out), ncol(out), nlayers(out)))
    inMemory <- TRUE
  }
  
  tr <- raster::blockSize(out, n = raster::nlayers(x) + 2)
  pb <- raster::pbCreate(tr$n)
  
  for (i in 1:tr$n) {
    rr <- firstrow + tr$row[i] - 1
    
    rowvals <- raster::getValuesBlock(x,
                                      row = rr,
                                      nrows = tr$nrows[i],
                                      firstcol,
                                      ncols)
    rowvals <- rowvals[, variables, drop = FALSE]
    
    res <- matrix(NA,
                  nrow = nrow(rowvals),
                  ncol = raster::nlayers(out))
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
      p <- predict.graf(
        object,
        newdata = newdata,
        type = type,
        CI = CI,
        maxn = maxn,
        ...
      )
      naind <- as.vector(attr(rowvals, "na.action"))
      if (!is.null(naind)) {
        res[-naind,] <- p
      } else {
        res <- p
      }
    }
    
    if (inMemory) {
      res <- array(res, dim = c(ncol(out), tr$nrows[i], nlayers(out)))
      res <- aperm(res, c(2, 1, 3))
      rows <- tr$row[i]:(tr$row[i] + tr$nrows[i] - 1)
      v[rows, ,] <- res
    } else {
      out <- raster::writeValues(out, res, tr$row[i])
    }
    pbraster::Step(pb, i)
  }
  raster::pbClose(pb)
  
  if (inMemory) {
    # permute v so that it gets vectorized in the right direction
    v <- aperm(v, c(2, 1, 3))
    for (i in 1:raster::nlayers(out)) {
      out <- raster::setValues(out,
                               as.vector(v[, , i]),
                               layer = i)
    }
  } else {
    out <- raster::writeStop(out)
  }
  names(out) <- colnames(p)
  return(out)
}
