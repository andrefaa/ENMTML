M_delimited <- function(var,
                        occ_xy,
                        method = c('buffer', 'ecoregions'),
                        BufferDistanceKm,
                        EcoregionsFile,
                        Dir,
                        spN,
                        SaveM = TRUE) {
  sapply(c('raster', 'dismo', 'rgeos', 'foreach', 'doParallel', 'tools'),
         require,
         character.only = TRUE)
  
  var <- var[[1]]
  var <- !is.na(var)
  var[var[] == 0] <- NA
  
  Dir_M<-"Extent_Masks"
  if (file.exists(file.path(Dir,Dir_M))){
   Dir_M<-file.path(Dir,Dir_M)
  } else {
   dir.create(file.path(Dir,Dir_M))
   Dir_M<-file.path(Dir,Dir_M)
  }
  
  #Check if GeoMasks already exist----
  if(all(paste0(spN,".tif")%in%list.files(Dir_M,pattern=".tif"))){
    print("GeoMasks Already Exist! Using already-created restrictions! ")
    return(Dir_M)
  }
  
  #Extent Restriction----
  if (method == 'buffer') {
      if (is.null(BufferDistanceKm) == T) {
        print('Please define a source of buffer in km:')
        cat(("1-Default based on presences distances\n2-A single buffer in m for all species\n3-A txt tab delimited file"))
        Buffer_Opt <- as.integer(readLines(n = 1))
        
        if (Buffer_Opt == 1) {
          occ_km <- list()
          for(i in 1:length(spN)){
            occ_km[[i]] <- latlong2grid(occ_xy[[i]])
          }
          names(occ_km) <- names(occ_xy)
          occ_km <- lapply(occ_km, function(x) as.matrix(dist(x,diag=F)))
          for(i in 1:length(occ_km)){
            diag(occ_km[[i]]) <- Inf
          }
          
          BufferDistanceKm <- sapply(occ_km, function(x) max(apply(x,2,min)))
          BufferDistanceKm <- BufferDistanceKm*1000
        }

        if (Buffer_Opt == 2) {
          cat("Please difine a distance in km")
          BufferDistanceKm <- as.integer(readLines(n = 1))
          BufferDistanceKm <- BufferDistanceKm * 1000
        }
        if (Buffer_Opt == 3) {
          cat("Please select the txt file with species specific distances in km")
          BufferDistanceKm <- read.table(file.choose(),sep="\t",h=T)
          BufferDistanceKm <- BufferDistanceKm$Dist * 1000
        }
      }
      
      varCord <- data.frame(coordinates(var))
      varCord <- varCord[which(!is.na(var[[1]][])),]
      varCord$Cell <- cellFromXY(var, xy = varCord)
      head(varCord)
      coordinates(varCord) <- ~ x + y
      crs(varCord) <-
        "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      
      M_list <- lapply(occ_xy, function(x) {
        coordinates(x) <- ~ x + y
        crs(x) <-
          "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        return(x)
      })
      
      cl <- makeCluster(detectCores()-1)
      registerDoParallel(cl)
      foreach(i = 1:length(M_list), .packages = c("raster")) %dopar% {
        if(Buffer_Opt==2){
          M <- dismo::circles(M_list[[i]], d = BufferDistanceKm, lonlat = T)
        }else{
          M <- dismo::circles(M_list[[i]], d = BufferDistanceKm[i], lonlat = T)
        }
          M <- rgeos::gUnaryUnion(M@polygons)
          Filter <- sp::over(varCord, M)
          M <- varCord$Cell[which(is.na(Filter))]
          M2 <- var
          M2[M] <- NA
          M2
          # M2 <- crop(M2,extent(varCord@coords[which(!is.na(Filter)), ]))
          writeRaster(M2,paste(Dir_M,paste(names(M_list[i]),".tif",sep=""),sep="/"),format="GTiff",overwrite=T)  
      }
      stopCluster(cl)
    }
  
  if (method == 'ecoregions') {
      if(is.null(EcoregionsFile)==T){
        print("Please select the file to be used as Ecoregions mask")
        Dir <- file.choose()
        EcoregionsFileExt <- unique(file_ext(Dir))
        if (any(EcoregionsFileExt %in% c('bil', 'asc', 'tif'))) {
          EcoregionsFile <- brick(stack(Dir))
          EcoregionsFile <- crop(EcoregionsFile, var)
        }
        if(EcoregionsFileExt=='shp'){
          EcoregionsFile <- shapefile(Dir)
          EcoregionsFile <- crop(EcoregionsFile, var)
          EcoregionsFile <- rasterize(EcoregionsFile, var)
          }
        if (EcoregionsFileExt == 'txt') {
          EcoregionsFile <- read.table(Dir, h = T)
          gridded(EcoregionsFile) <- ~ x + y
          EcoregionsFile <- brick(stack(EcoregionsFile))
          EcoregionsFile <- crop(EcoregionsFile, var)
        }
        
        sp.Ecoregions <- lapply(occ_xy,function(x) extract(EcoregionsFile,x))
        sp.Ecoregions <- lapply(sp.Ecoregions,function(x) {
          x <- unique(na.omit(x))
          x <- x[x!=0]
          })
        
        cl <- makeCluster(detectCores()-1)
        registerDoParallel(cl)
        foreach(i = 1:length(sp.Ecoregions), .packages = c("raster")) %dopar% {
          EcoregionSp <- EcoregionsFile%in%sp.Ecoregions[[i]]
          EcoregionSp[EcoregionSp[]%in%0] <- NA
          # EcoregionSp <- crop(EcoregionSp,extent(rasterToPoints(EcoregionSp)[,-3]))
          writeRaster(EcoregionSp,paste(Dir_M,paste(names(sp.Ecoregions[i]),".tif",sep=""),sep="/"),format="GTiff",overwrite=T)  
        }
        stopCluster(cl)
      }
  }
  return(Dir_M)
}
  