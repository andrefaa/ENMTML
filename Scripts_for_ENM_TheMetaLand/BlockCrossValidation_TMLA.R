## Written by Santiago Velazco

BlockPartition_TMLA <- function(evnVariables = NA,
           RecordsData,
           N,
           FromTo = FromTo,
           pseudoabsencesMethod = PabM,
           PrAbRatio = PabR,
           DirSave = dirO) {
    
  # RecordsData: matrix or data frame with presences records 
  # N: 2 (dafault). interger  Number of group for data  paritioning 
  # pseudoabsences: logical, TRUE (dafault). 
  # mask: Raster object. Preferible on wiht the same resolution and extent than 
  #       variable used in the future 
  # pseudoabsencesMethod: a character string indicating which pseudo-absences
  #                       method is to be computed
  # PrAbRatio: numeric. value of PrAbRatio to be computed.
  # plot: TRUE (default). Make graphs of data paritioning process
  # evnVariables: Raster object. Variable set to be used in pseusoabsences
  # cellSize: numeric vector. a vector of values with different cell grid sizes

  #Cellsize
  cellSize = seq(0.5, 10, by = .5)
  
  # Mask
  mask <- evnVariables[[1]]
  if(class(mask)!="brick"){
    mask <- brick(mask)
  }
  names(mask) <- "Group"
  mask[!is.na(mask[,])] <- 1
  
  #Extent
  e<-extent(mask)
  
  # Loop for each species-----
  ResultList <- list()
  SpNames <- names(RecordsData)
  BestGridList <- list()
  
  # LOOP----
  for(s in FromTo:length(RecordsData)){
  print(paste(s, SpNames[s]))
  # Extract coordinates----
  presences <- RecordsData[[s]]
  
  mask2 <- mask
  mask2[] <- 0
  
  # presences2 <- data.frame(pa=c(rep(1, nrow(presences)), rep(0, nrow(absences))),rbind(presences, absences))
  presences2 <- data.frame(pa=c(rep(1, nrow(presences))),presences)
  
  # Transform the presences points in a DataFrameSpatialPoints
  coordinates(presences2)=presences2[,c("x","y")]
  crs(presences2) <-projection(mask)
  
  #### Data partitioning using a grid approach ####
  
  # Create a list of grids based on different raster resolution
  grid <- list() #List of grids
  
  # raster resolution
  DIM<-matrix(0,length(cellSize),2) # the number of rows and columns of each grid
  colnames(DIM)<-c("R","C")
  
  for(i in 1:length(cellSize)) {
    mask3 <- mask2
    res(mask3) <-cellSize[i]
    DIM[i, ] <- dim(mask3)[1:2]
    values(mask3) <- 1 # Add values to cells /
    NAS <- c(extract(mask3,presences2)) # Extract values to test if exist NAs
    if(any(is.na(NAS))) {
      while (any(is.na(NAS))) {
        extent(mask3) <- extent(mask3)+cellSize[i]
        res(mask3) <- cellSize[i] # Give to cells a size
        DIM[i, ] <- dim(mask3)[1:2]
        values(mask3) <- 1
        NAS <- extract(mask3, presences2)
      }
    }
    grid[[i]] <- rasterToPolygons(mask3)
  }
  rm(mask3)
  
  # In this section is assigned the group of each cell
  for(i in 1:length(grid)) {
    if(any(N==c(2,4,6,8,10))){
      # odds number of partition
      group<-rep(c(rep(1:N, DIM[i,2])[1:DIM[i,2]],
                   rep(c((N/2+1):N, 1:(N/2)),DIM[i,2])[1:DIM[i,2]])
                 ,nrow(grid[[i]]@data)/DIM[i,2])[1:nrow(grid[[i]]@data)]    
      grid[[i]]@data <-
        data.frame(group, expand.grid(C = 1:DIM[i, 2], R = 1:DIM[i, 1]))
    }
  }
  
  # Matrix within each columns represent the partitions of points
  # for each grid resolution
  part <- data.frame(matrix(0, nrow(presences2@data), length(grid)))
  for (i in 1:length(grid)) {
    part[, i] <- over(presences2, grid[[i]][,1])
  }
  
  part2 <- list()
  for (i in 1:length(grid)) {
    part2[[i]] <- data.frame(over(presences2, grid[[i]]), presences)
  }
  
  # Here will be deleted grids that assigned partitions less than the number
  # of groups
  pp <- sapply(part[1:nrow(presences),], function(x)
    length(unique(range(x))))
  pp <- ifelse(pp == N, TRUE, FALSE)
  # Elimination of those partition that have one record in some group
  pf <- sapply(part[1:nrow(presences),], table)
  if(is.list(pf)==TRUE){
    pf <- which(sapply(pf, min)==1)
  }else{
  pf <- which(apply(pf, 2, min)==1)
    }
  pp[pf] <- FALSE 
  grid <- grid[pp]
  part <- part[,pp]
  part2 <- part2[pp]
  
  # Performace of cells ----
  # SD of number of records per cell size-----
  pa <- presences2@data[, 1] # Vector wiht presences and absences
  Sd.Grid.P <- rep(NA, length(grid))
  for (i in 1:ncol(part)) {
    Sd.Grid.P[i] <- sd(table(part[pa == 1, i])) /
      mean(table(part[pa == 1, i]))
  }
  
  # MESS -----
  Mess.Grid.P <- rep(NA, length(grid))
  Env.P <- extract(evnVariables, presences)
  for (i in 1:ncol(part)) {
    Env.P1 <- cbind(part[i], Env.P)
    Env.P2 <- split(Env.P1[, -1], Env.P1[, 1])
    mess1 <- MESS(Env.P2[[1]], Env.P2[[2]])
    Mess.Grid.P[i] <- mean(mess1$TOTAL, na.rm = TRUE)
    rm(Env.P1)
  }
  
  # Imoran-----
  Imoran.Grid.P <- rep(NA, length(grid)) 
  for(p in 1:length(grid)) {
    part3 <- part2[[p]]
    # rownames(part3) <-
    #   paste(part3$group, part3$C, part3$R, part3$lon, part3$lat)
    mineucli <- list()
    for (r in 1:DIM[p, 'R']) {
      part4 <- part3[part3$R == r, ]
      if (nrow(part4) >= 2) {
        for (c in 1:DIM[p, 'C']) {
          A <- part4[part4$C == c,]
          B <- part4[part4$C == (c + 1),]
          euclilist <- matrix(0, (nrow(A) * nrow(B)), 3)
          if ((nrow(A) >= 1 & nrow(B) >= 1)) {
            coord <- data.frame(expand.grid(rownames(A), rownames(B)),
                                cbind(expand.grid(A[, 5], B[, 5]),
                                      expand.grid(A[, 4], B[, 4])))
            colnames(coord) <- c("A", "B", "xa", "xb", "ya", "yb")
            coord$eucli <-
              sqrt((coord$xa - coord$xb) ^ 2 + (coord$ya - coord$yb) ^ 2)
            mineucli[[length(mineucli) + 1]] <-
              coord[which.min(coord$eucli), ]
          }
        }
      }
    }
    for (r in 1:DIM[p, 'C']) {
      part4 <- part3[part3$C == r, ]
      if (nrow(part4) >= 2) {
        for (c in 1:DIM[p, 'R']) {
          A <- part4[part4$R == c,]
          B <- part4[part4$R == (c + 1),]
          if ((nrow(A) >= 1 & nrow(B) >= 1)) {
            coord <- data.frame(expand.grid(rownames(A), rownames(B)),
                                cbind(expand.grid(A[, 5], B[, 5]),
                                      expand.grid(A[, 4], B[, 4])))
            colnames(coord) <- c("A", "B", "xa", "xb", "ya", "yb")
            coord$eucli <-
              sqrt((coord$xa - coord$xb) ^ 2 + (coord$ya - coord$yb) ^ 2)
            mineucli[[length(mineucli) + 1]] <-
              coord[which.min(coord$eucli), ]
          }
        }
      }
    }
    
    euclida <- ldply(mineucli, data.frame)
    euclida$A <- as.character(euclida$A)
    euclida$B <- as.character(euclida$B)
    species2 <- presences[c(euclida$A,euclida$B),]
    dist <- as.matrix(dist(species2))
    dist <- 1/dist
    diag(dist) <- 0
    dist[which(dist==Inf)] <- 0
    species2$pc1 <- extract(evnVariables[[1]], species2)
    if(nrow(species2)==0){
      Imoran.Grid.P[p] <- NA
    }else{
      Imoran.Grid.P[p] <- 
        Moran.I(species2$pc1, dist, na.rm = T)$observed  
    }
    
    # plot(presences)
    # plot(grid[[p]], add=T)
    # points(species2, col="red", pch=19)
  }
  
  Imoran.Grid.P <- abs(Imoran.Grid.P) # OJO estamos dejando todos los valores positivos
  N.grid <- 1:length(cellSize[pp])
  
  Opt <- data.frame(N.grid, cellSize[pp], round(data.frame(Imoran.Grid.P, Mess.Grid.P, Sd.Grid.P),3))
  # Cleaning those variances based in data divided in a number of partition less than 
  # the number of groups
  
  # SELLECTION OF THE BEST CELL SIZE----
  Opt2 <- Opt
  while (nrow(Opt2) > 1) {
    # I MORAN
    if (nrow(Opt2) == 1) break
    Opt2 <- Opt2[which(Opt2$Imoran.Grid.P <= summary(Opt2$Imoran.Grid.P)[2]),]
    if (nrow(Opt2) == 1) break
    # MESS
    Opt2 <- Opt2[which(Opt2$Mess.Grid.P >= summary(Opt2$Mess.Grid.P)[5]),]
    if (nrow(Opt2) == 1) break
    # SD
    Opt2 <- Opt2[which(Opt2$Sd.Grid.P <= summary(Opt2$Sd.Grid.P)[2]),]
    if(nrow(Opt2)==2) break
    }
  
  if(nrow(Opt2)>1){
    Opt2 <- Opt2[nrow(Opt2),]
  }
  
  # Optimum size for presences
  print(Opt2)
  Optimum.Grid <- grid[[Opt2$N.grid]]
  Optimum.Grid@data[,c("C", "R")] <- NULL
  presences <- data.frame(Partition=part[,Opt2$N.grid], presences)

  # Pseudoabsences allocation-----

    # Random-----
    if(pseudoabsencesMethod=="rnd"){
      # Clip the mask raster to generate rando pseudoabsences
      pseudo.mask <- mask
      pseudo.mask2 <- list()
      
      RtoP <- data.frame(rasterToPoints(mask)[,-3])
      coordinates(RtoP)=c("x","y")
      crs(RtoP) <-projection(mask)
      
      FILTER <- over(RtoP, Optimum.Grid)
      pseudo.mask[which(pseudo.mask[]==1)] <- as.matrix(FILTER)
      writeRaster(pseudo.mask, paste(DirSave, paste(SpNames[s],'.tif',sep=""),sep='/'),
                    format = 'GTiff', NAflag = -9999, overwrite = TRUE)

      for(i in 1:N){
        mask3 <- pseudo.mask
        mask3[!mask3[]==i] <- NA 
        pseudo.mask2[[i]] <- mask3
      }
      
      pseudo.mask <- brick(pseudo.mask2)
      rm(pseudo.mask2)
      
      # plot(pseudo.mask, legend = F, col = "red")
      
      # Random allocation of pseudoabsences 
      absences <- list()
      for (i in 1:N) {
        absences.0 <- randomPoints(pseudo.mask[[i]], (1 / PrAbRatio)*sum(presences[,1]==i),
                                   ext = e,
                                   prob = FALSE)
        colnames(absences.0) <- c(x, y)
        absences[[i]] <- as.data.frame(absences.0)
      }
      names(absences) <- 1:N
    }
    # constrain
    if(pseudoabsencesMethod=="constrain"){
      
      
        pseudo.mask <- mask
        RtoP <- data.frame(rasterToPoints(pseudo.mask)[,-3])
        coordinates(RtoP)=c("x","y")
        crs(RtoP) <-projection(mask)
        FILTER <- over(RtoP, Optimum.Grid)
        pseudo.mask[which(pseudo.mask[]==1)] <- as.matrix(FILTER)
        writeRaster(pseudo.mask, paste(DirSave, paste(SpNames[s],'.tif',sep=""),sep='/'),
                    format = 'GTiff', NAflag = -9999, overwrite = TRUE)

      Model <- bioclim(evnVariables, presences[,-1])
      pseudo.mask <- dismo::predict(Model, evnVariables, ext=e)
      names(pseudo.mask) <- "Group"
      pseudo.mask <- round(pseudo.mask, 5)
      pseudo.mask <-(pseudo.mask-minValue(pseudo.mask))/
        (maxValue(pseudo.mask)-minValue(pseudo.mask))
      pseudo.mask <-(1-pseudo.mask)>=0.99 #environmental constrain
      pseudo.mask[which(pseudo.mask[,]==FALSE)] <- NA
      
      # Split the raster of environmental layer with grids
      pseudo.mask2 <- list()
      RtoP <- data.frame(rasterToPoints(pseudo.mask)[,-3])
      coordinates(RtoP)=c("x","y")
      crs(RtoP) <-projection(mask)
      FILTER <- over(RtoP, Optimum.Grid)
      pseudo.mask[which(pseudo.mask[]==1)] <- as.matrix(FILTER)
      
      for(i in 1:N){
        mask3 <- pseudo.mask
        mask3[!mask3[]==i] <- NA 
        pseudo.mask2[[i]] <- mask3
      }
      pseudo.mask <- brick(pseudo.mask2)
      rm(pseudo.mask2)
      
      absences <- list()
      for (i in 1:N) {
        absences.0 <- randomPoints(pseudo.mask[[i]], (1 / PrAbRatio)*sum(presences[,1]==i),
                                   ext = e,
                                   prob = FALSE)
        colnames(absences.0) <- c("lon", "lat")
        absences[[i]] <- as.data.frame(absences.0)
      }
      
      names(absences) <- 1:N
    }
    absences <- plyr::ldply(absences, data.frame)
    names(absences) <- c("Partition", x, y)
    absences[,c(x,y)] <- round(absences[,c(x,y)],4)
    colnames(absences) <- colnames(presences)
  # Final data.frame result----
    PresAbse <- rep(c(1, 0), sapply(list(presences, absences), nrow))
    result <- data.frame(Sp=SpNames[s], PresAbse, rbind(presences, absences), stringsAsFactors = F)
    result <- result[,c("Sp","x","y","Partition","PresAbse")]

  # Final data.frame result2----
  ResultList[[s]]<- result
  
  BestGridList[[s]]<- Opt2
  }
  
  FinalResult <- ldply(ResultList)
  write.table(FinalResult,paste(DirSave,"OccBlocks.txt",sep="\\"),sep="\t",row.names=F)
  FinalInfoGrid <- ldply(BestGridList)
  write.table(FinalInfoGrid, paste(DirSave, "BestPartitions,txt", sep = '/'), sep="\t",
              col.names = T, row.names = F)
  
  return(FinalResult)
}