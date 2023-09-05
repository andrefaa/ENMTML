pca_raster <- function(env_layer,
                      Dir,
                      DirP = NULL) {
  #0.Create PCA Folder
  DirPCA <- file.path(Dir, "PCA")
  dir.create(DirPCA)
  #PCA Tables Folder
  DirPCATab <- file.path(DirPCA, "Tables")
  dir.create(DirPCATab)
  
  #Projection Folders
  
  DirP_PCA <- file.path(dirname(Dir), 'Projection_PCA')
  dir.create(DirP_PCA)
  
  if(!is.null(DirP)){
    FoldersProj <- list.files(dirname(DirP[[1]]))
    FoldersProj <- as.list(file.path(DirP_PCA, FoldersProj))
    lapply(FoldersProj, function(x)
      dir.create(x))
  }
  
  # mean
  means <- raster::cellStats(env_layer, mean)
  # SD
  stds <- raster::cellStats(env_layer, sd)
  # Standardize values
  env_layer <- raster::scale(env_layer, center = means, scale = stds)
  p0 <- raster::as.data.frame(env_layer, xy = FALSE, na.rm = TRUE)
  pca <- stats::prcomp(p0,
                       retx = TRUE,
                       scale. = FALSE,
                       center = FALSE
  )
  cof <- pca$rotation
  
  cvar <- summary(pca)$importance["Cumulative Proportion", ]
  naxis <- Position(function(x) {
    x >= 0.95
  }, cvar)
  cvar <- data.frame(cvar)
  
  pca <- stats::prcomp(p0, retx = TRUE, scale. = FALSE, center = FALSE, rank. = naxis)
  rm(p0)
  env_layer <- terra::predict(terra::rast(env_layer), pca) %>% 
    raster::stack() %>% raster::brick()
  
  result <- list(
    env_layer = env_layer,
    coefficients = data.frame(cof) %>% dplyr::tibble(variable = rownames(.), .),
    cumulative_variance = dplyr::tibble(PC = 1:nrow(cvar), cvar)
  )
  
  utils::write.table(
    result$coefficients,
    file.path(DirPCATab, "Coeficient.txt"),
    sep = "\t",
    row.names = F
  )
  
  #Cummulative Variance
  utils::write.table(
    result$cumulative_variance,
    file.path(DirPCATab, "CumulativeVariance.txt"),
    sep = "\t",
    row.names = F
  )
  
  #Save PCsthat account for 95% of total variability
    raster::writeRaster(
      result$env_layer,
      paste(DirPCA, names(result$env_layer), sep = "/"),
      bylayer = T,
      format = "GTiff",
      overwrite = T
    )
  
  if(!is.null(DirP)){
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
      for (j in 1:length(DirP)) {
        ProjT[[j]] <-
          utils::read.table(list.files(DirP[[j]], pattern = '.txt',full.names = T), h = T)
        if(any(colnames(ProjT[[j]])%in%c("long","lat"))){
          colnames(ProjT[[j]])[colnames(ProjT[[j]])%in%c("long","lat")] <- c("x","y")
        }
        sp::gridded(ProjT[[j]]) <- ~ x + y
        ProjT[[j]] <- raster::brick(raster::stack(ProjT[[j]]))
      }
    }
    
    for (j in 1:length(ProjT)) {
      scen <- raster::scale(ProjT[[j]], center = means, scale = stds) %>% 
        terra::rast()
      scen <- terra::predict(scen, pca) %>% raster::stack()
      raster::writeRaster(
        scen,
        paste(FoldersProj[[j]], paste0(names(scen), ".tif"), sep = "/"),
        bylayer = TRUE,
        format = "GTiff",
        overwrite = TRUE
      )
    }
  }
  
  return(result)
}
