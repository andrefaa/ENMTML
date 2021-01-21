PCAFuturo <- function(Env,
                      Dir,
                      DirP,
                      Save = '') {
  #0.Create PCA Folder
  DirPCA <- file.path(Dir, "PCA")
  dir.create(DirPCA)
  #PCA Tables Folder
  DirPCATab <- file.path(DirPCA, "Tables")
  dir.create(DirPCATab)
  
  #Projection Folders
  
  DirP_PCA <- file.path(dirname(Dir), 'Projection_PCA')
  dir.create(DirP_PCA)
  FoldersProj <- list.files(dirname(DirP[[1]]))
  FoldersProj <- as.list(file.path(DirP_PCA, FoldersProj))
  lapply(FoldersProj, function(x)
    dir.create(x))
  
  #1.Present PCA
  DF <- raster::rasterToPoints(Env)
  DF <- stats::na.omit(DF)
  PcaR <- DF[, -c(1:2)]
  
  means <- colMeans(PcaR)
  stds <- apply(PcaR, 2, stats::sd)
  
  #Scale transform
  DScale <- data.frame(apply(PcaR, 2, scale))
  
  # PCA
  DPca <- stats::prcomp(DScale,
                 retx = TRUE,
                 center = T,
                 scale = T)
  
  #Coefficients
  Coef <- DPca$rotation
  Coef2 <- data.frame(cbind(Variable = names(Env), Coef))
  utils::write.table(
    Coef2,
    file.path(DirPCATab, "Coeficient.txt"),
    sep = "\t",
    row.names = F
  )
  
  #Cummulative Variance
  NEixos <- length(summary(DPca)$importance[3, ])
  CumVar <- summary(DPca)$importance[3, ]
  VarEx <- data.frame(CumVar)
  utils::write.table(
    VarEx,
    file.path(DirPCATab, "CumulativeVariance.txt"),
    sep = "\t",
    row.names = F
  )
  
  #Save PCsthat account for 95% of total variability
  Eix <- as.data.frame(DPca$x)
  EixXY <- cbind(DF[, (1:2)], Eix)
  sp::gridded(EixXY) <- ~ x + y
  PCAPr <- raster::stack(EixXY)
  raster::crs(PCAPr) <- raster::crs(Env)
  PCA.95 <- PCAPr[[1:(sum(VarEx <= 0.95) + 1)]]
  if (Save == "Y") {
    raster::writeRaster(
      PCA.95,
      paste(DirPCA, names(PCA.95), sep = "/"),
      bylayer = T,
      format = "GTiff",
      overwrite = T
    )
  }
  
  #2.Project PCA
  
  DirP <- as.list(DirP)
  ProjEX <-
    lapply(DirP, function(x)
      unique(tools::file_ext(list.files(x))))
  form <- c('bil', 'asc', 'txt', 'tif')
  ProjEX <- unique(unlist(ProjEX)[unlist(ProjEX) %in% form])
  
  if (any(ProjEX %in% c('asc', 'bil', 'tif'))) {
    ProjT <-
      lapply(DirP, function(x)
        raster::brick(raster::stack(list.files(x, pattern=paste0('\\.', ProjEX, '$'),
                                       full.names=T))))
  }
  
  if (any(ProjEX %in% 'txt')) {
    ProjT <- list()
    for (j in DirP) {
      ProjT[[j]] <-
        utils::read.table(file.path(DirP[[j]], list.files(DirP[[j]], pattern = '.txt')), h = T)
      if(any(colnames(ProjT[[j]])%in%c("long","lat"))){
        colnames(ProjT[[j]])[colnames(ProjT[[j]])%in%c("long","lat")] <- c("x","y")
      }
      sp::gridded(ProjT[[j]]) <- ~ x + y
      ProjT[[j]] <- raster::brick(raster::stack(ProjT[[j]]))
    }
  }
  
  ProjE <- lapply(ProjT, function(x)
    raster::rasterToPoints(x))
  ProjE <- lapply(ProjE, function(x)
    stats::na.omit(x))
  ProjXY <- lapply(ProjE, "[", , c("x","y"))
  ProjER <-
    lapply(ProjE, function(z)
      z[, !(colnames(z) %in% c("x", "y"))])
  
  scale <- lapply(ProjER, function(x)
    sweep(x, 2, means))
  scale <- lapply(scale, function(x)
    x %*% diag(1 / stds))
  PCAFut <- lapply(scale, function(x)
    x %*% Coef)
  PCAFut <- Map(cbind,ProjXY, PCAFut)
  PCAFut <- lapply(PCAFut, function(x) as.data.frame(x))
  # PCAFut <-
  #   lapply(PCAFut, function(x)
  #     data.frame(cbind(ProjE[[1]][, (1:2)], x)))
  PCAFut.95 <- list()
  for (j in 1:length(PCAFut)) {
    sp::gridded(PCAFut[[j]]) <- ~ x + y
    PCAFut[[j]] <- raster::stack(PCAFut[[j]])
    raster::crs(PCAFut[[j]]) <- raster::crs(Env)
    names(PCAFut[[j]]) <- names(PCAPr)
    PCAFut.95[[j]] <- PCAFut[[j]][[1:raster::nlayers(PCA.95)]]
    if (Save == "Y") {
      raster::writeRaster(
        PCAFut.95[[j]],
        paste(FoldersProj[[j]], names(PCAFut.95[[j]]), sep = "/"),
        bylayer = T,
        format = "GTiff",
        overwrite = T
      )
    }
  }
  return(PCAFut.95)
}
