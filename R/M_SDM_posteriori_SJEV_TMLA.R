#written by Santiago Velazco

MSDM_Posterior <- function(RecordsData,
                           Threshold=Thr,
                           cutoff=c('MaxMin','Pres','LowerQ','MCP','MCPBuffer'),
                           PredictType=ENS,
                           DirSave=NULL,
                           DirRaster=NULL){
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
  foldCat <- paste(DirSave,"BIN",sep="/")
  for(i in 1:length(foldCat)){
    dir.create(foldCat[i])
  }
  
  # Vector with species names
  SpNames <- as.character(unique(RecordsData[, "sp"]))
  
  # Data.frame wiht two columns 1-names of the species 
  # 2-the directory of raster of each species
  if (is.null(DirRaster) == TRUE) {
    stop("Give a directory in the DirRaster argument")
  }else{
    RasterList <- list.files(DirRaster, pattern = '.tif$')
    sp <- gsub('.tif$','',RasterList)
    RasterList <- list.files(DirRaster, pattern = '.tif$', full.names = T)
    RasterList <- data.frame(sp, RasterList, stringsAsFactors = F)
    colnames(RasterList) <- c("sp",'RasterList')
  }
  
  if(cutoff=="MCPBuffer"){
    cat("Select buffer distance (km):")
    CUT_Buf <- (as.numeric(readLines(n = 1)))*1000
  }
  
  # loop to process each species
  for(s in 1:length(SpNames)){
    print(paste(s, "from", length(SpNames),":", SpNames[s]))
    # Read the raster of the species
    Adeq <- raster(RasterList[RasterList[,"sp"]==SpNames[s],'RasterList'])

    # Extract values for one species and calculate the threshold
    SpData <- data.frame(RecordsData[RecordsData[, "sp"] == SpNames[s], ])
    PredPoint <- extract(Adeq, SpData[, c("x", "y")])
    PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
    Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                     PredPoint[PredPoint$PresAbse == 0, 2])
    Thr <- unlist(c(threshold(Eval))[Threshold])
    print("First Extract")
    #### Cutoff MCP----
    if(cutoff=="MCP"){
      hull <- convHull(SpData[SpData[,"PresAbse"]==1, c("x", "y")], lonlat=TRUE)
      hull <- predict(Adeq, hull, mask=TRUE)
      Adeq[(hull[]==0)] <- 0
      writeRaster(Adeq, paste(DirSave, paste(SpNames[s],'.tif', sep = ""), sep = "/"),
                  format = "GTiff",
                  NAflag = -9999,
                  overwrite = TRUE)
        Mask <- Adeq>Thr
        writeRaster(Mask, paste(foldCat, paste(SpNames[s],'.tif', sep = ""), sep = "/"),
                    format = "GTiff",
                    NAflag = -9999,
                    overwrite = TRUE)
    }
    if(cutoff=="MCPBuffer"){
      hull <- convHull(SpData[SpData[,"PresAbse"]==1, c("x", "y")], lonlat=TRUE)
      hull <- predict(Adeq, hull, mask=TRUE)
      
      pts1 <- SpData[SpData[, "PresAbse"] == 1, c("x", "y")]
      spraster<-rasterize(pts1,Adeq,field=1)
      sps<-as(spraster,'SpatialPixels')@coords
      dist<-dist2(sps,sps,method='euclidean',p=2)
      dist[dist==0] <- NA
      distmin<-apply(dist,1,function(x) min(x, na.rm=T))
      
      hull2 <- hull
      hull2[hull2[]==0] <- NA
      hull2 <- boundaries(hull2)
      hull2[hull2[]==0] <- NA
      df <- rasterToPoints(hull2)
      df <- df[df[,3]==1,-3]
      buf <- circles(df, lonlat=T,d=CUT_Buf)
      buf <- predict(Adeq, buf,  mask=TRUE)
      buf[(hull[]==1)] <- 1
      buf[(!is.na(Adeq[])&is.na(buf[]))] <- 0
      Adeq[which(buf[]!=1)] <- 0

      writeRaster(Adeq, paste(DirSave, paste(SpNames[s],'.tif', sep = ""), sep = "/"),
                  format = "GTiff",
                  NAflag = -9999,
                  overwrite = TRUE)
        Mask <- (Adeq>Thr)
        writeRaster(Mask, paste(foldCat, paste(SpNames[s],'.tif', sep = ""), sep = "/"),
                    format = "GTiff",
                    NAflag = -9999,
                    overwrite = TRUE)
    }
    
    if(cutoff%in%c("MaxMin","LowerQ","Pres")){
      
      # Transform coordinate in a SpatialPoints object
      pts1 <- SpData[SpData[, "PresAbse"] == 1, c("x", "y")]
      coordinates(pts1) <- ~ x + y
      crs(pts1) <- crs(Adeq)
      
      # Raster with areas equal or grater than the threshold
        AdeqBin <- Adeq >= as.numeric(Thr)
        AdeqBin[AdeqBin[] == 0] <- NA
        # A "SpatialPolygonsDataFrame" which each adequability patch is a feature
        AdeqBin2 <- rasterToPolygons(AdeqBin, fun=NULL, n=8, na.rm=TRUE, digits=12, dissolve=TRUE)
        AdeqBin2 <- disaggregate(AdeqBin2)
        AdeqBin2$layer <- NULL
        print("RasterToPolygons...OK")
        # Individualize each patch with a number
        AdeqBin2$ID <- 1:length(AdeqBin2)
        # create a data.frame wiht coordinate and patch number
        AdeqPoints <- rasterToPoints(AdeqBin)[,1:2]
        AdeqPoints <- cbind(AdeqPoints, ID=extract(AdeqBin2, AdeqPoints)[,'ID'])
        # Find the patches that contain presences records
        polypoint<-intersect(AdeqBin2,pts1)
        print("Second Extract and Intersect")
        if(cutoff=="Pres"){
          Mask2 <- Adeq
          Msk <- rasterize(polypoint,Adeq,background=0)
          Msk[is.na(Adeq[])] <- NA
          Mask2[Msk!=1] <- 0
          Mask <- Mask2>=Thr
          
          writeRaster(Mask2,
                      paste(DirSave, paste(SpNames[s],'.tif', sep = ""), sep = "/"),
                      format = "GTiff",
                      NAflag = -9999,
                      overwrite = TRUE)
          writeRaster(Mask,
                      paste(foldCat, paste(SpNames[s],'.tif', sep = ""), sep = "/"),
                      format = "GTiff",
                      NAflag = -9999,
                      overwrite = TRUE)
          
        }else{
        # Create a vector wich contain the number (e.i. ID) of the patches 
        # with presences
        filter1 <- unique(polypoint$ID)
        # In this step are created two data.frame one with the patches coordinates 
        # that contain presences and another with patches coordinates without presences
        CoordPathP <- as.data.frame(AdeqPoints[AdeqPoints[,3]%in%filter1,])
        CoordPathNP <- as.data.frame(AdeqPoints[!AdeqPoints[,3]%in%filter1,])
        # Here is created a matrix wiht all combination between ID of patches 
        # with and without presences
        
        if(ncol(CoordPathP)==1){
          CoordPathP <- data.frame(t(CoordPathNP))
        rownames(CoordPathP) <- NULL
          }
        
        if(ncol(CoordPathNP)==1){
          CoordPathNP <- data.frame(t(CoordPathNP))
          rownames(CoordPathNP) <- NULL
          }
          npatch1 <- unique(CoordPathP[,3])
          npatch2 <- unique(CoordPathNP[,3])
                             
        DistBetweenPoly0 <- expand.grid(npatch1, npatch2)
        DistBetweenPoly0$Distance <- NA
        DistBetweenPoly0 <- as.matrix(DistBetweenPoly0)
        print("Distance...OK")
        # Euclidean Distance between patches wiht and without presences 
        for(i in 1:nrow(DistBetweenPoly0)){
          comb <- (DistBetweenPoly0[i,1:2])
          A <- CoordPathP[CoordPathP[,3]==comb[1],1:2]
          B <- CoordPathNP[CoordPathNP[,3]==comb[2],1:2]
          
          if(nrow(A)>=40) {
            SEQ <- round(seq(0, nrow(A), by = (nrow(A)) / 20))
            dist <- rep(NA, length(SEQ))
            for (j in 2:length(SEQ)) {
              SEQ2 <- (SEQ[(j - 1)] + 1):SEQ[j]
              dist[j] <- min(dist2(A[SEQ2, ], B, method = 'euclidean', p = 2), na.rm = T)
            }
            eucdist <- min(dist[2:length(SEQ)], na.rm = T)
          }else{
            eucdist <- min(dist2(A, B, method = 'euclidean', p = 2))
          }
          DistBetweenPoly0[i,3] <- eucdist
        }
        
        DistBetweenPoly0 <- DistBetweenPoly0[order(DistBetweenPoly0[,2]),]
        # Minimum Euclidean Distance between patches wiht and without presences
        DistBetweenPoly <- tapply(X = DistBetweenPoly0[,3], DistBetweenPoly0[,2], min)
        # Adding value of distance patches to cells
        AdeqBin2$Eucldist <- 0
        AdeqBin2$Eucldist[!AdeqBin2$ID%in%filter1] <- round(DistBetweenPoly, 4)
        print("Second Distance....Ok")
        # CUTOFF------
        if(cutoff=='MaxMin'){
        # Cutoff based on the maximum value of the minimum distance
          spraster<-rasterize(pts1,Adeq,field=1)
          sps<-as(spraster,'SpatialPixels')@coords
          dist<-dist2(sps,sps,method='euclidean',p=2)
          dist[dist==0] <- NA
          distmin<-apply(dist,1,function(x) min(x, na.rm=T))# 
          CUT<-max(distmin)
        }
        if(cutoff=="LowerQ"){
        # Cutoff based the lower quartile distance
        CUT <- c(summary(DistBetweenPoly0[,3]))[2]
        }
  
        AdeqPoints <- rasterToPoints(AdeqBin)[,1:2]
        AdeqPoints <- extract(AdeqBin2, AdeqPoints)[,'Eucldist']
        fist <- AdeqBin
        fist[fist[]==1] <- AdeqPoints
        # Threshold based on Maximum value of minimum ditance between ocurrences
        final<-fist<=CUT
        final[final==0] <- NA
        Mask2 <- Adeq
        Mask <- Adeq>=as.numeric(Thr)
        Mask[Mask==1] <- 0
        Mask[!is.na(final[])] <- 1
        Mask2[Mask!=1] <- 0
        # Save results as raster object
        writeRaster(Mask2,
                    paste(DirSave, paste(SpNames[s],'.tif', sep = ""), sep = "/"),
                    format = "GTiff",
                    NAflag = -9999,
                    overwrite = TRUE)
        writeRaster(Mask,
                    paste(foldCat, paste(SpNames[s],'.tif', sep = ""), sep = "/"),
                    format = "GTiff",
                    NAflag = -9999,
                    overwrite = TRUE)
      }
    }
  }
}
