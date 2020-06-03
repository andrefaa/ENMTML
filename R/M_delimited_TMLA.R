M_delimited <- function(var,
                        occ_xy,
                        method = c('BUFFER', 'MASK'),
                        BufferDistanceKm,
                        EcoregionsFile,
                        Dir,
                        spN,
                        Buffer_Opt,
                        SaveM = TRUE,
                        Mfolder = NULL) {
  var <- var[[1]]
  var <- !is.na(var)
  var[var[] == 0] <- NA
  
  Dir_M <- "Extent_Masks"
  if (file.exists(file.path(Dir, Dir_M))) {
    Dir_M <- file.path(Dir, Dir_M)
  } else {
    dir.create(file.path(Dir, Dir_M))
    Dir_M <- file.path(Dir, Dir_M)
  }
  
  #Check if GeoMasks already exist----
  if (any(paste0(spN, ".tif") %in% list.files(Dir_M, pattern = ".tif"))) {
    if (all(paste0(spN, ".tif") %in% list.files(Dir_M, pattern = ".tif"))) {
      print("GeoMasks already exist for all species! Using already-created restrictions!")
      return(Dir_M)
    } else{
      print("GeoMasks already exist for some species! Creating GeoMasks for the rest of species")
      spN <-
        spN[!paste0(spN, ".tif") %in% list.files(Dir_M, pattern = ".tif")]
    }
  }
  
  
  #Extent Restriction----
  if (method == 'BUFFER') {
    if (Buffer_Opt == 1) {
      occ_km <- list()
      for (i in 1:length(spN)) {
        occ_km[[i]] <- SpatialEpi::latlong2grid(occ_xy[[i]])
      }
      names(occ_km) <- names(occ_xy)
      occ_km <-
        lapply(occ_km, function(x)
          as.matrix(stats::dist(x, diag = F)))
      for (i in 1:length(occ_km)) {
        diag(occ_km[[i]]) <- Inf
      }
      
      BufferDistanceKm <-
        sapply(occ_km, function(x)
          max(apply(x, 2, min)))
      BufferDistanceKm <- BufferDistanceKm * 1000
    }
    
    if (Buffer_Opt == 2) {
      BufferDistanceKm <- BufferDistanceKm * 1000
    }
    # if (Buffer_Opt == 3) {
    #   cat("Please select the txt file with species specific distances in km")
    #   BufferDistanceKm <- utils::read.table(file.choose(),sep="\t",h=T)
    #   BufferDistanceKm <- BufferDistanceKm$Dist * 1000
    # }
    
    varCord <- data.frame(sp::coordinates(var))
    varCord <- varCord[which(!is.na(var[[1]][])), ]
    varCord$Cell <- raster::cellFromXY(var, xy = varCord)
    sp::coordinates(varCord) <- ~ x + y
    raster::crs(varCord) <- crs(var)
    
    M_list <- lapply(occ_xy, function(x) {
      sp::coordinates(x) <- ~ x + y
      raster::crs(x) <- crs(var)
      return(x)
    })
    
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    foreach(i = 1:length(M_list),
            .packages = c("raster")) %dopar% {
              if (Buffer_Opt == 2) {
                M <- dismo::circles(M_list[[i]], d = BufferDistanceKm, lonlat = T)
              } else{
                M <-
                  dismo::circles(M_list[[i]], d = BufferDistanceKm[i], lonlat = T)
              }
              M <- rgeos::gUnaryUnion(M@polygons)
              Filter <- sp::over(varCord, M)
              M <- varCord$Cell[which(is.na(Filter))]
              M2 <- var
              M2[M] <- NA
              M2
              # M2 <- crop(M2,extent(varCord@coords[which(!is.na(Filter)), ]))
              raster::writeRaster(M2,
                                  paste(Dir_M, paste(names(M_list[i]), ".tif", sep = ""), sep = "/"),
                                  format = "GTiff",
                                  overwrite = T)
            }
    parallel::stopCluster(cl)
  }
  
  if (method == 'MASK') {
    Dir <- EcoregionsFile
    EcoregionsFileExt <- unique(tools::file_ext(Dir))
    if (any(EcoregionsFileExt %in% c('bil', 'asc', 'tif'))) {
      EcoregionsFile <- raster::brick(stack(Dir))
      EcoregionsFile <- raster::crop(EcoregionsFile, var)
    }
    if (EcoregionsFileExt == 'shp') {
      EcoregionsFile <- raster::shapefile(Dir)
      EcoregionsFile <- raster::crop(EcoregionsFile, var)
      EcoregionsFile <- raster::rasterize(EcoregionsFile, var)
    }
    if (EcoregionsFileExt == 'txt') {
      EcoregionsFile <- utils::read.table(Dir, h = T)
      sp::gridded(EcoregionsFile) <- ~ x + y
      EcoregionsFile <-
        raster::brick(raster::stack(EcoregionsFile))
      EcoregionsFile <- raster::crop(EcoregionsFile, var)
    }
    
    sp.Ecoregions <-
      lapply(occ_xy, function(x)
        raster::extract(EcoregionsFile, x))
    sp.Ecoregions <- lapply(sp.Ecoregions, function(x) {
      x <- unique(stats::na.omit(x))
      x <- x[x != 0]
    })
    
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    foreach(i = 1:length(sp.Ecoregions),
            .packages = c("raster")) %dopar% {
              EcoregionSp <- EcoregionsFile
              EcoregionSp[!EcoregionSp[] %in% sp.Ecoregions[[i]]] <- NA
              EcoregionSp <- EcoregionSp > 0
              raster::writeRaster(
                EcoregionSp,
                paste(Dir_M, paste(
                  names(sp.Ecoregions[i]), ".tif", sep = ""
                ), sep = "/"),
                format = "GTiff",
                overwrite = T
              )
            }
    parallel::stopCluster(cl)
  }
  
  if (method == 'USER-DEFINED') {
    Dir <- EcoregionsFile
    EcoregionsFileExt <- unique(tools::file_ext(list.files(Dir)))
    if (any(EcoregionsFileExt %in% c('bil', 'asc', 'tif'))) {
      EcoregionsFile <-
        raster::brick(raster::stack(
          list.files(Dir, full.names = T, pattern = EcoregionsFileExt)
        ))
      if (is.na(raster::crs(EcoregionsFile))) {
        raster::crs(EcoregionsFile) <-
          "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
      }
      if (!raster::compareRaster(EcoregionsFile, var)) {
        if (raster::crs(EcoregionsFile) != raster::crs(var)) {
          EcoregionsFile <- raster::projectRaster(EcoregionsFile, var)
        }
        EcoregionsFile <- raster::resample(EcoregionsFile, var)
        EcoregionsFile <- raster::crop(EcoregionsFile, var)
      }
    }
    if (EcoregionsFileExt == 'shp') {
      EcoregionsFile <- raster::shapefile(Dir)
      EcoregionsFile <- raster::crop(EcoregionsFile, var)
      EcoregionsFile <- raster::rasterize(EcoregionsFile, var)
    }
    if (EcoregionsFileExt == 'txt') {
      EcoregionsFile <- utils::read.table(Dir, h = T)
      sp::gridded(EcoregionsFile) <- ~ x + y
      EcoregionsFile <- raster::brick(raster::stack(EcoregionsFile))
      if (is.na(raster::crs(EcoregionsFile))) {
        raster::crs(EcoregionsFile) <-
          "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
      }
      if (!raster::compareRaster(EcoregionsFile, var)) {
        if (raster::crs(EcoregionsFile) != raster::crs(var)) {
          EcoregionsFile <- raster::projectRaster(EcoregionsFile, var)
        }
        EcoregionsFile <- raster::resample(EcoregionsFile, var)
        EcoregionsFile <- raster::crop(EcoregionsFile, var)
      }
    }
    
    EcoregionsFile[EcoregionsFile != 1] <- NA
    raster::writeRaster(
      EcoregionsFile,
      paste(Dir_M, paste(names(sp.Ecoregions), ".tif", sep = ""), sep = "/"),
      format = "GTiff",
      overwrite = T,
      bylayer = T
    )
  }
  return(Dir_M)
}
