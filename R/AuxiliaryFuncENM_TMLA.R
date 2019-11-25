# Standarize raster values between 0 and 1------
STANDAR <- function(x) {
  result <-
    (x - raster::cellStats(x, min)) / (raster::cellStats(x, max) - raster::cellStats(x, min))
  return(result)
}

STANDAR_FUT <- function(ModelFut, ModelPre) {
  result <-
    (ModelFut - raster::cellStats(ModelPre, min)) / (raster::cellStats(ModelPre, max) -
                                                       raster::cellStats(ModelPre, min))
  return(result)
}


# Function to remove outliers from MAH and DOMAIN preiction
rem_out <- function(r) {
  ss <- quantile(r[], na.rm = T)
  me <- median(r[], na.rm = T)
  out <- boxplot.stats(r[])[[4]]
  r[which(r[] %in% out[out <= me])] <- ss[2]
  return(r)
}

# Prediction for Mahalanobis and Domaint------
PREDICT_DomainMahal <- function(mod, variables) {
  df <- na.omit(rasterToPoints(variables))
  pred <- dismo::predict(mod, df[, -c(1:2)])
  result <- variables[[1]]
  result[which(!is.na(result[]))] <- pred
  return(result)
}


# Prediction of different algorithm------
PREDICT <- function(Variables, Models_List) {
  ListRaster <- as.list(names(Models_List))
  names(ListRaster) <- names(Models_List)
  Algorithm <- names(Models_List)

  if (any(Algorithm == "BIOCLIM") == TRUE) {
    ListRaster[["BIOCLIM"]] <-
      round(predict(Variables, Models_List[["BIO"]]), 4)
  }
  if (any(Algorithm == "MAXENTD") == TRUE) {
    ListRaster[["MAXENTD"]] <-
      round(predict(Variables, Models_List[["MAXENTD"]], args = 'outputformat=cloglog'),
            4)
  }
  if (any(Algorithm == "MAXENTS") == TRUE) {
    ListRaster[["MAXENTS"]] <-
      round(predict(Variables, Models_List[["MAXENTS"]], args = 'outputformat=cloglog'),
            4)
  }
  if (any(Algorithm == "MAXENTD_NEW") == TRUE) {
    ListRaster[["MAXENTD_NEW"]] <-
      round(predict(
        Variables,
        Models_List[["MAXENTD_NEW"]],
        clamp = F,
        type = "cloglog"
      ),
      4)
  }
  if (any(Algorithm == "MAXENTS_NEW") == TRUE) {
    ListRaster[["MAXENTS_NEW"]] <-
      round(predict(
        Variables,
        Models_List[["MAXENTS_NEW"]],
        clamp = F,
        type = "cloglog"
      ),
      4)
  }
  if (any(Algorithm == "MAXLIKE") == TRUE) {
    ListRaster[["MAXLIKE"]] <-
      round(predict(Variables, Models_List[["MAXLIKE"]]), 4)
  }
  if (any(Algorithm == "SVM") == TRUE) {
    ListRaster[["SVM"]] <-
      round(predict(Variables, Models_List[["SVM"]]), 4)
  }
  if (any(Algorithm == "RF") == TRUE) {
    ListRaster[["RF"]] <-
      round(predict(Variables, Models_List[["RF"]]), 4)
  }
  if (any(Algorithm == "GLMC") == TRUE) {
    ListRaster[["GLMC"]] <-
      round(predict(Variables, Models_List[["GLMC"]]), 4)
  }
  if (any(Algorithm == "GAMC") == TRUE) {
    ListRaster[["GAMC"]] <-
      round(predict(Variables, Models_List[["GAMC"]]), 4)
  }
  if (any(Algorithm == "GAUSS") == TRUE) {
    Pred <-
      predict.graf.raster(
        Models_List[["GAUSS"]],
        Variables,
        type = "response",
        CI = 0.95,
        maxn = NULL
      )
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
    pmin(1, pmax(0, (x - min) / (max - min)))

  }

############################################################
#                                                          #
#             function used for data partition             ####
#    in  BlockPartition_TMLA() & BandsPartition_TMLA()     #
#                                                          #
############################################################

KM <- function(cell_coord, variable, NumAbsence) {
  # cell_env = cell environmental data
  # variable = a stack raster with variables
  # NumAbsence = number of centroids sellected as absences
  # This function will be return a list whith the centroids sellected
  # for the clusters
  var <- raster::extract(variable, cell_coord)
  Km <- kmeans(cbind(cell_coord, var), centers = NumAbsence)
  ppa2 <- (list(
    Centroids = Km$centers[, 1:2],
    Clusters = Km$cluster
  ))

  abse <- ppa2$Centroids
  colnames(abse) <- c("lon", "lat")

  ppaNA <- raster::extract(variable[[1]], abse)

  if (sum(is.na(ppaNA)) != 0) {
    ppaNA1 <- which(is.na(ppaNA))
    # Wich pseudo absence is a environmental NA
    ppaNA2 <-
      as.data.frame(cbind(ppa2$Clusters, raster::rasterToPoints(variable[[1]])[, -3]))
    ppa.clustres <- split(ppaNA2[, 2:3], ppaNA2[, 1])
    ppa.clustres <- ppa.clustres[ppaNA1]
    nearest.point <- list()

    if (length(ppaNA1) < 2) {
      error <- matrix(abse[ppaNA1, ], 1)
      colnames(error) <- colnames(abse)
    } else{
      error <- abse[ppaNA1, ]
    }

    nearest.point <- list()
    for (eee in 1:length(ppaNA1)) {
      x <-
        flexclust::dist2(ppa.clustres[[eee]], error, method = "euclidean", p = 2)
      x <- as.data.frame(x)
      nearest.point[[eee]] <-
        as.data.frame(ppa.clustres[[eee]][which.min(x[, eee]), ])
    }
    nearest.point <- t(matrix(unlist(nearest.point), 2))
    abse[ppaNA1, ] <- nearest.point
  }
  return(abse)
}

# Inverse bioclim
inv_bio <- function(e, p) {
  Model <- dismo::bioclim(e, p)
  r <- dismo::predict(Model, e)
  names(r) <- "Group"
  r <- round(r, 5)
  r <- (r - minValue(r)) /
    (maxValue(r) - minValue(r))
  r <- (1 - r) >= 0.99 #environmental constrain
  r[which(r[, ] == FALSE)] <- NA
  return(r)
}

# Inverse geo
inv_geo <- function(e, p, d) {
  Model <- dismo::circles(p, lonlat = T, d = d)
  r <- mask(e[[1]], Model@polygons, inverse = T)
  names(r) <- "Group"
  r[is.na(r) == F] <- 1
  r[which(r[, ] == FALSE)] <- NA
  return(r)
}
