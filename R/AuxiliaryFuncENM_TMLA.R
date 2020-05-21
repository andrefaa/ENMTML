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
  me <- stats::median(r[], na.rm = T)
  out <- grDevices::boxplot.stats(r[])[[4]]
  r[which(r[] %in% out[out <= me])] <- ss[2]
  return(r)
}

# Prediction for Mahalanobis and Domaint------
PREDICT_DomainMahal <- function(mod, variables) {
  df <- stats::na.omit(raster::rasterToPoints(variables))
  pred <- dismo::predict(mod, df[, -c(1:2)])
  result <- variables[[1]]
  result[which(!is.na(result[]))] <- pred
  return(result)
}

#Prediction for ENFA----
PREDICT_ENFA <- function(mod,prediction_dataset,train_dataset=NULL){
  if (class(prediction_dataset)=="data.frame"){
    Zli <- as.matrix(prediction_dataset) %*% as.matrix(mod$co)
    ZER <- 1:nrow(prediction_dataset)
    f1 <- function(x) rep(x, ZER)
    Sli <- apply(Zli, 2, f1)
    m <- apply(Sli, 2, mean)
    cov <- t(as.matrix(Zli)) %*% as.matrix(Zli)/nrow(Zli)
    res <- (data.frame(MD = stats::mahalanobis(Zli, center = m, cov = cov,inverted = F)))*-1
    return(res)
  }else if(class(prediction_dataset)=="RasterBrick"){
    PredRas <- raster::values(prediction_dataset)
    POS <- which(is.na(PredRas[,1]))
    Zli <- as.matrix(stats::na.omit(raster::values(prediction_dataset)) %*% as.matrix(mod$co))
    POSPRE <- raster::cellFromXY(prediction_dataset[[1]],train_dataset[train_dataset$PresAbse==1,c("x","y")])
    ZER <- rep(0,nrow(PredRas))
    ZER[POSPRE] <- 1
    if(length(POS) > 0) {
      ZER <- ZER[-POS]
    }
    f1 <- function(x) rep(x, ZER)
    Sli <- apply(Zli, 2, f1)
    m <- apply(Sli, 2, mean)
    cov <- t(as.matrix(Zli)) %*% as.matrix(Zli)/nrow(Zli)
    PredRas <- (data.frame(MD = stats::mahalanobis(Zli, center = m, cov = cov,inverted = F)))*-1
    XY <- raster::xyFromCell(prediction_dataset[[1]],1:raster::ncell(prediction_dataset[[1]]))
    PredRAS <- data.frame(cbind(XY,ENF=NA))
    if(length(POS) > 0) {
      PredRAS[-POS, "ENF"] <- PredRas
    } else{
      PredRAS[, "ENF"] <- PredRas
    }
    sp::gridded(PredRAS) <- ~x+y
    FinalModelT <- rem_out(raster::raster(PredRAS))
    FinalModel <- FinalModelT
    return(FinalModel)
  }
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
  Km <- stats::kmeans(cbind(cell_coord, var), centers = NumAbsence)
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
      error <- matrix(abse[ppaNA1,], 1)
      colnames(error) <- colnames(abse)
    } else{
      error <- abse[ppaNA1,]
    }
    
    nearest.point <- list()
    for (eee in 1:length(ppaNA1)) {
      x <-
        flexclust::dist2(ppa.clustres[[eee]], error, method = "euclidean", p = 2)
      x <- as.data.frame(x)
      nearest.point[[eee]] <-
        as.data.frame(ppa.clustres[[eee]][which.min(x[, eee]),])
    }
    nearest.point <- t(matrix(unlist(nearest.point), 2))
    abse[ppaNA1,] <- nearest.point
  }
  return(abse)
}

# Inverse bioclim
inv_bio <- function(e, p) {
  Model <- dismo::bioclim(e, p)
  r <- dismo::predict(Model, e)
  names(r) <- "Group"
  r <- round(r, 5)
  r <- (r - raster::minValue(r)) /
    (raster::maxValue(r) - raster::minValue(r))
  r <- (1 - r) >= 0.99 #environmental constrain
  r[which(r[,] == FALSE)] <- NA
  return(r)
}


# Inverse geo
inv_geo <- function(e, p, d) {
  Model <- dismo::circles(p, lonlat = T, d = d)
  r <- raster::mask(e[[1]], Model@polygons, inverse = T)
  names(r) <- "Group"
  r[is.na(r) == F] <- 1
  r[which(r[,] == FALSE)] <- NA
  return(r)
}


# MESS function modified from modEvA package 
MESS <- function (V, P, id.col = NULL) 
{
  index.V <- 1:nrow(V)
  index.P <- 1:nrow(P)
  if (is.null(id.col)) {
    if (ncol(V) != ncol(P)) 
      stop("The number of variables in V and P does not match.")
  }
  else {
    if (ncol(V) != ncol(P) - 1) 
      stop("'id.col' should refer to a column in P; P should therefore have one more column than V.")
    P.input <- P
    P <- P[, -id.col]
  }
  n.vars <- ncol(V)
  n.varP <- ncol(P)
  nrow.P <- nrow(P)
  results <- matrix(nrow = nrow.P, ncol = n.vars, dimnames = list(NULL, 
                                                                  colnames(P)))
  for (i in 1:n.vars) {
    min.Vi <- min(V[, i], na.rm = TRUE)
    max.Vi <- max(V[, i], na.rm = TRUE)
    SIM <- vector("numeric", nrow.P)
    for (j in 1:nrow.P) {
      VV <- V[, i]
      VV[VV < P[j, i]] <- 1
      VV[VV >= P[j, i]] <- 0
      Fj <- sum(VV, na.rm = TRUE) * 100/length(VV)
      if (Fj == 0) 
        SIM[j] <- (P[j, i] - min.Vi)/(max.Vi - min.Vi) * 
        100
      else if (Fj > 0 && Fj <= 50) 
        SIM[j] <- 2 * Fj
      else if (Fj > 50 && Fj < 100) 
        SIM[j] <- 2 * (100 - Fj)
      else if (Fj == 100) 
        SIM[j] <- (max.Vi - P[j, i])/(max.Vi - min.Vi) * 
        100
    }
    results[, i] <- SIM
  }
  results <- data.frame(results)
  results$TOTAL <- apply(results[, 1:n.vars], 1, min)
  results$MoD <- as.factor(colnames(results)[apply(results[, 
                                                           1:n.vars], 1, which.min)])
  if (!is.null(id.col)) {
    results <- data.frame(P.input[, id.col], results)
    colnames(results)[1] <- colnames(P.input)[id.col]
  }
  return(results)
}

