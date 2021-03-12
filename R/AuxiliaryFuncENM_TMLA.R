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


# Function to remove outliers from MAH and DOMAIN prediction
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
  }else if(class(prediction_dataset)=="RasterBrick" || class(prediction_dataset)=="RasterStack"){
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

madifa2 <- function (dudi, pr, scannf = TRUE, nf = 2) 
{
  if (!inherits(dudi, "dudi")) 
    stop("object of class dudi expected")
  call <- match.call()
  if (any(is.na(dudi$tab))) 
    stop("na entries in table")
  if (!is.vector(pr)) 
    stop("pr should be a vector")
  prb <- pr
  pr <- pr/sum(pr)
  row.w <- dudi$lw
  col.w <- dudi$cw
  Z <- as.matrix(dudi$tab)
  n <- nrow(Z)
  f1 <- function(v) sum(v * pr)
  center <- apply(Z, 2, f1)
  Z <- sweep(Z, 2, center)
  f2 <- function(v) sum((v^2) * pr)
  Ze <- sweep(Z, 2, sqrt(col.w), "*")
  DpZ <- apply(Ze, 2, function(x) x * pr)
  Se <- crossprod(Ze, DpZ)
  Ge <- crossprod(Ze, apply(Ze, 2, function(x) x * row.w))
  eS <- eigen(Se)
  b <- ifelse(is.nan(eS$values^(-0.5)),0,eS$values)
  S12 <- eS$vectors %*% diag(b) %*% t(eS$vectors)
  W <- S12 %*% Ge %*% S12
  s <- eigen(W)$values
  if (scannf) {
    barplot(s)
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  if (nf <= 0 | nf > ncol(Ze)) 
    nf <- 1
  tt <- as.data.frame((S12 %*% eigen(W)$vectors))
  ww <- apply(tt, 2, function(x) x/sqrt(col.w))
  norw <- sqrt(diag(t(as.matrix(tt)) %*% as.matrix(tt)))
  co <- sweep(ww, 2, norw, "/")
  li <- Z %*% apply(co, 2, function(x) x * col.w)
  co <- as.data.frame(co)
  li <- as.data.frame(li)
  varus <- apply(li, 2, function(x) sum((x^2) * pr))
  l1 <- sweep(li, 2, sqrt(unlist(varus)), "/")
  mahasu <- apply(l1, 1, function(x) sqrt(sum(x^2)))
  co <- data.frame(co[, 1:nf])
  li <- data.frame(li[, 1:nf])
  l1 <- data.frame(l1[, 1:nf])
  names(co) <- paste("Axis", (1:nf), sep = "")
  row.names(co) <- dimnames(dudi$tab)[[2]]
  names(li) <- paste("Comp.", (1:nf), sep = "")
  names(l1) <- paste("Comp.", (1:nf), sep = "")
  row.names(li) <- dimnames(dudi$tab)[[1]]
  row.names(l1) <- dimnames(dudi$tab)[[1]]
  corav <- cor(dudi$tab, li)
  madifa <- list(call = call, tab = data.frame(Z), pr = prb, 
                 cw = col.w, nf = nf, eig = s, lw = row.w, li = li, l1 = l1, 
                 co = co, mahasu = mahasu, cor = corav)
  class(madifa) <- "madifa"
  return(madifa)
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



KM_BLOCK <- function(cell_coord, variable, NumAbsence) {
  # cell_env = cell environmental data
  # variable = a stack raster with variables
  # NumAbsence = number of centroids sellected as absences
  # This function will be return a list whith the centroids sellected
  # for the clusters
  var <- raster::extract(variable, cell_coord)
  Km <- stats::kmeans(cbind(cell_coord, var), centers = NumAbsence)
  return(list(
    Centroids = Km$centers[, 1:2],
    Clusters = Km$cluster
  ))
  
  
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


##%######################################################%##
#                                                          #
####     MESS function modified from modEvA package     ####
#                                                          #
##%######################################################%##

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

##%######################################################%##
#                                                          #
####        RasterLayer NA CHECK By Boris Leroy         ####
#                                                          #
##%######################################################%##
synchroniseNA <- function(x)
{
  if(raster::canProcessInMemory(x, n = 2))
  {
    val <- raster::getValues(x)
    NA.pos <- unique(which(is.na(val), arr.ind = T)[, 1])
    val[NA.pos, ] <- NA
    x <- raster::setValues(x, val)
    return(x)
  } else
  {
    x <- raster::mask(x, calc(x, fun = sum))
    return(x)
  }
}

##%######################################################%##
#                                                          #
####                      I moran                       ####
#                                                          #
##%######################################################%##
# I'moran from ape
Moran.I <- function (x, weight, scaled = FALSE, na.rm = FALSE, alternative = "two.sided")
{
  if (dim(weight)[1] != dim(weight)[2])
    stop("'weight' must be a square matrix")
  n <- length(x)
  if (dim(weight)[1] != n)
    stop("'weight' must have as many rows as observations in 'x'")
  ei <- -1/(n - 1)
  nas <- is.na(x)
  if (any(nas)) {
    if (na.rm) {
      x <- x[!nas]
      n <- length(x)
      weight <- weight[!nas, !nas]
    }
    else {
      warning("'x' has missing values: maybe you wanted to set na.rm = TRUE?")
      return(list(observed = NA, expected = ei, sd = NA,
                  p.value = NA))
    }
  }
  ROWSUM <- rowSums(weight)
  ROWSUM[ROWSUM == 0] <- 1
  weight <- weight/ROWSUM
  s <- sum(weight)
  m <- mean(x)
  y <- x - m
  cv <- sum(weight * y %o% y)
  v <- sum(y^2)
  obs <- (n/s) * (cv/v)
  if (scaled) {
    i.max <- (n/s) * (sd(rowSums(weight) * y)/sqrt(v/(n -
                                                        1)))
    obs <- obs/i.max
  }
  S1 <- 0.5 * sum((weight + t(weight))^2)
  S2 <- sum((apply(weight, 1, sum) + apply(weight, 2, sum))^2)
  s.sq <- s^2
  k <- (sum(y^4)/n)/(v/n)^2
  sdi <- sqrt((n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * s.sq) -
                 k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s.sq))/((n -
                                                                     1) * (n - 2) * (n - 3) * s.sq) - 1/((n - 1)^2))
  alternative <- match.arg(alternative, c("two.sided",
                                          "less", "greater"))
  pv <- stats::pnorm(obs, mean = ei, sd = sdi)
  if (alternative == "two.sided")
    pv <- if (obs <= ei)
      2 * pv
  else 2 * (1 - pv)
  if (alternative == "greater")
    pv <- 1 - pv
  list(observed = obs, expected = ei, sd = sdi, p.value = pv)
}

##%######################################################%##
#                                                          #
####                    LatLong2grid                    ####
#                                                          #
##%######################################################%##
lat2grd <- function(input){
  toradians <- atan(1)/45
  radiusearth <- 0.5*(6378.2+6356.7)
  sine51 <- sin( 51.5*toradians)
  output <- data.frame(cbind(
    x=(input[,1]*toradians)*radiusearth*sine51,
    y=(input[,2]*toradians)*radiusearth
  ))
}

##%######################################################%##
#                                                          #
####                    RandomPoints                    ####
##%######################################################%##
OptimRandomPoints <- function(r, n, p) {
  v <- raster::getValues(r)
  v <- which(!is.na(v))
  v <- v[!v%in%raster::cellFromXY(r,p)]
  v <- sample(v, n)
  v <- raster::xyFromCell(r, v)
  return(v)
}
