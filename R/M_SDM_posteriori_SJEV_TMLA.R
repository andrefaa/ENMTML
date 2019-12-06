#written by Santiago Velazco

MSDM_Posterior <- function(RecordsData,
                           Threshold = thr,
                           cutoff = c('OBR', 'PRES', 'LR', 'MCP', 'MCP-B'),
                           CUT_Buf = CUT_Buf,
                           DirSave = NULL,
                           DirRaster = NULL) {
  # Arguments-----
  # RecordsData: data.frame. A data.frame with species,
  #              occurrences coordinates, and a column with
  #              presences and absences (1 and 0).
  # PresAbse: character. Names of the column with presences and absences (1 and 0)
  # Species: character. Names of the column species names
  # x: character. Names of the column with longitude
  # y: character. Names of the column with latitude
  # FromTo: numeric. A numeric vector to select which species have to be processed
  # Threshold: character. Type of threshold to be use, default is  equal sensitivity and specificity
  # DirRaster: character. Directory that contains species adequability.
  # DirSave: character. Directory to save the raster (in .tif format)
  
  #Create Binary folder
  foldCat <- paste(DirSave, "BIN", sep = "/")
  for (i in 1:length(foldCat)) {
    dir.create(foldCat[i])
  }
  
  # Vector with species names
  SpNames <-
    substr(list.files(DirRaster, pattern = '.tif$'),
           1,
           nchar(list.files(DirRaster, pattern = '.tif$')) - 4)
  
  # Data.frame wiht two columns 1-names of the species
  # 2-the directory of raster of each species
  if (is.null(DirRaster) == TRUE) {
    stop("Give a directory in the DirRaster argument")
  } else{
    RasterList <- list.files(DirRaster, pattern = '.tif$')
    sp <- gsub('.tif$', '', RasterList)
    RasterList <-
      list.files(DirRaster, pattern = '.tif$', full.names = T)
    RasterList <- data.frame(sp, RasterList, stringsAsFactors = F)
    colnames(RasterList) <- c("sp", 'RasterList')
    #Bianries
    RasterListBin <-
      list.files(file.path(DirRaster, Threshold),
                 pattern = '.tif$',
                 full.names = T)
    RasterListBin <-
      data.frame(sp, RasterListBin, stringsAsFactors = F)
    colnames(RasterListBin) <- c("sp", 'RasterList')
  }
  
  # loop to process each species
  for (s in 1:length(SpNames)) {
    print(paste(s, "from", length(SpNames), ":", SpNames[s]))
    # Read the raster of the species
    Adeq <-
      raster::raster(RasterList[RasterList[, "sp"] == SpNames[s], 'RasterList'])
    Bin <-
      raster::raster(RasterListBin[RasterListBin[, "sp"] == SpNames[s], 'RasterList'])
    if (is.na(crs(Adeq))) {
      crs(Adeq) <-
        "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    }
    crs(Bin) <- crs(Adeq)
    
    # # Extract values for one species and calculate the threshold
    SpData <-
      data.frame(RecordsData[RecordsData[, "sp"] == SpNames[s],])
    
    #### Cutoff MCP----
    if (cutoff == "MCP") {
      hull <-
        dismo::convHull(SpData[SpData[, "PresAbse"] == 1, c("x", "y")], lonlat = TRUE)
      hull <- dismo::predict(Adeq, hull, mask = TRUE)
      Adeq[(hull[] == 0)] <- 0
      Bin[(hull[] == 0)] <- 0
      raster::writeRaster(
        Adeq,
        paste(DirSave, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
        format = "GTiff",
        NAflag = -9999,
        overwrite = TRUE
      )
      raster::writeRaster(
        Bin,
        paste(foldCat, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
        format = "GTiff",
        NAflag = -9999,
        overwrite = TRUE
      )
    }
    if (cutoff == "MCP-B") {
      hull <-
        dismo::convHull(SpData[SpData[, "PresAbse"] == 1, c("x", "y")], lonlat = TRUE)
      hull <- dismo::predict(Adeq, hull, mask = TRUE)
      
      pts1 <- SpData[SpData[, "PresAbse"] == 1, c("x", "y")]
      spraster <- raster::rasterize(pts1, Adeq, field = 1)
      sps <- methods::as(spraster, 'SpatialPixels')@coords
      dist <- flexclust::dist2(sps, sps, method = 'euclidean', p = 2)
      dist[dist == 0] <- NA
      distmin <- apply(dist, 1, function(x)
        min(x, na.rm = T))
      
      hull2 <- hull
      hull2[hull2[] == 0] <- NA
      hull2 <- raster::boundaries(hull2)
      hull2[hull2[] == 0] <- NA
      df <- raster::rasterToPoints(hull2)
      df <- df[df[, 3] == 1, -3]
      buf <- dismo::circles(df, lonlat = T, d = CUT_Buf)
      buf <- dismo::predict(Adeq, buf,  mask = TRUE)
      buf[(hull[] == 1)] <- 1
      buf[(!is.na(Adeq[]) & is.na(buf[]))] <- 0
      Adeq[which(buf[] != 1)] <- 0
      Bin[which(buf[] != 1)] <- 0
      
      raster::writeRaster(
        Adeq,
        paste(DirSave, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
        format = "GTiff",
        NAflag = -9999,
        overwrite = TRUE
      )
      raster::writeRaster(
        Bin,
        paste(foldCat, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
        format = "GTiff",
        NAflag = -9999,
        overwrite = TRUE
      )
    }
    
    if (cutoff %in% c("OBR", "LR", "PRES")) {
      # Transform coordinate in a SpatialPoints object
      pts1 <- SpData[SpData[, "PresAbse"] == 1, c("x", "y")]
      sp::coordinates(pts1) <- ~ x + y
      raster::crs(pts1) <- raster::crs(Adeq)
      
      # Raster with areas equal or grater than the threshold
      AdeqBin <- Adeq * Bin
      AdeqBin[AdeqBin[] == 0] <- NA
      AdeqBin <- round(AdeqBin, 6)
      AdeqBin <- raster::clump(AdeqBin)
      AdeqPoints <- data.frame(raster::rasterToPoints(AdeqBin)[, 1:2])
      AdeqPoints <-
        cbind(AdeqPoints, ID = raster::extract(AdeqBin, AdeqPoints))
      # Find the patches that contain presences records
      polypoint <- as.numeric(unique(raster::extract(AdeqBin, pts1)))
      AdeqBin2 <- AdeqBin
      AdeqBin2[!AdeqBin2[] %in% polypoint] <- NA
      AdeqBin3 <- !is.na(AdeqBin2)
      if (cutoff == "PRES") {
        Mask <- AdeqBin3
        Mask[is.na(Mask)] <- 0
        Mask[is.na(Adeq[])] <- NA
        Mask2 <- Adeq * Mask
        
        
        raster::writeRaster(
          Mask2,
          paste(DirSave, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
          format = "GTiff",
          NAflag = -9999,
          overwrite = TRUE
        )
        raster::writeRaster(
          Mask,
          paste(foldCat, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
          format = "GTiff",
          NAflag = -9999,
          overwrite = TRUE
        )
        
      } else{
        # Create a vector wich contain the number (e.i. ID) of the patches
        # with presences
        filter1 <- unique(stats::na.omit(raster::values(AdeqBin2)))
        # In this step are created two data.frame one with the patches coordinates
        # that contain presences and another with patches coordinates without presences
        CoordPathP <-
          as.data.frame(AdeqPoints[AdeqPoints[, 3] %in% filter1, ])
        CoordPathNP <-
          as.data.frame(AdeqPoints[!AdeqPoints[, 3] %in% filter1, ])
        # Here are created a matrix wiht all combination between ID of patches
        # with and without presences
        
        if (ncol(CoordPathP) == 1) {
          CoordPathP <- data.frame(t(CoordPathNP))
          rownames(CoordPathP) <- NULL
        }
        
        if (ncol(CoordPathNP) == 1) {
          CoordPathNP <- data.frame(t(CoordPathNP))
          rownames(CoordPathNP) <- NULL
        }
        npatch1 <- unique(CoordPathP[, 3])
        npatch2 <- unique(CoordPathNP[, 3])
        
        DistBetweenPoly0 <- expand.grid(npatch1, npatch2)
        DistBetweenPoly0$Distance <- NA
        DistBetweenPoly0 <- as.matrix(DistBetweenPoly0)
        # Euclidean Distance between patches wiht and without presences
        for (i in 1:nrow(DistBetweenPoly0)) {
          comb <- (DistBetweenPoly0[i, 1:2])
          A <- CoordPathP[CoordPathP[, 3] == comb[1], 1:2]
          B <- CoordPathNP[CoordPathNP[, 3] == comb[2], 1:2]
          
          if (nrow(A) >= 40) {
            SEQ <- round(seq(0, nrow(A), by = (nrow(A)) / 20))
            dist <- rep(NA, length(SEQ))
            for (j in 2:length(SEQ)) {
              SEQ2 <- (SEQ[(j - 1)] + 1):SEQ[j]
              dist[j] <-
                min(flexclust::dist2(A[SEQ2,], B, method = 'euclidean', p = 2), na.rm = T)
            }
            eucdist <- min(dist[2:length(SEQ)], na.rm = T)
          } else{
            eucdist <- min(flexclust::dist2(A, B, method = 'euclidean', p = 2))
          }
          DistBetweenPoly0[i, 3] <- eucdist
        }
        
        DistBetweenPoly0 <-
          DistBetweenPoly0[order(DistBetweenPoly0[, 2]), ]
        if (is.null(nrow(DistBetweenPoly0))) {
          DistBetweenPoly0 <- t(as.matrix(DistBetweenPoly0))
        }
        # Minimum Euclidean Distance between patches wiht and without presences
        DistBetweenPoly <-
          tapply(X = DistBetweenPoly0[, 3], DistBetweenPoly0[, 2], min)
        
        # CUTOFF------
        if (cutoff == 'OBR') {
          # Cutoff based on the maximum value of the minimum distance
          spraster <- raster::rasterize(pts1, Adeq, field = 1)
          sps <- methods::as(spraster, 'SpatialPixels')@coords
          dist <- flexclust::dist2(sps, sps, method = 'euclidean', p = 2)
          dist[dist == 0] <- NA
          distmin <- apply(dist, 1, function(x)
            min(x, na.rm = T))#
          CUT <- max(distmin)
        }
        if (cutoff == "LR") {
          # Cutoff based the lower quartile distance
          CUT <- c(summary(DistBetweenPoly0[, 3]))[2]
        }
        
        #Chosen patches
        Mask <- DistBetweenPoly0[DistBetweenPoly0[, 3] <= CUT, 2]
        Mask <-
          raster::match(AdeqBin,
                        table = c(Mask, npatch1),
                        nomatch = 0)
        Mask <- Mask != 0
        # Mask <- AdeqBin%in%c(Mask,npatch1)
        Mask[is.na(Adeq)] <- NA
        Mask2 <- Adeq * Mask
        
        # Save results as raster object
        raster::writeRaster(
          Mask2,
          paste(DirSave, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
          format = "GTiff",
          NAflag = -9999,
          overwrite = TRUE
        )
        raster::writeRaster(
          Mask,
          paste(foldCat, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
          format = "GTiff",
          NAflag = -9999,
          overwrite = TRUE
        )
      }
    }
  }
}
