## Written by Santiago Velazco

BlockPartition_TMLA <- function(evnVariables = NA,
           RecordsData,
           N,
           pseudoabsencesMethod = PabM,
           PrAbRatio = PabR,
           DirSave = DirB,
           DirM = DirM,
           MRst = sp_accessible_area,
           type = TipoMoran) {
    
  # RecordsData: matrix or data frame with presences records 
  # N: 2 (dafault). interger  Number of group for data  paritioning 
  # pseudoabsences: logical, TRUE (dafault). 
  # mask: Raster object. Preferible on wiht the same resolution and extent than 
  #       variable used in the future 
  # pseudoabsencesMethod: a character string indicating which pseudo-absences
  #                       method is to be computed
  # PrAbRatio: numeric. value of PrAbRatio to be computed.
  # evnVariables: Raster object. Variable set to be used in pseusoabsences
  # cellSize: numeric vector. a vector of values with different cell grid sizes
  
  
  
  KM <- function(cell_coord, variable, NumAbsence) {
    # cell_env = cell environmental data
    # variable = a stack raster with variables
    # NumAbsence = number of centroids sellected as absences
    # This function will be return a list whith the centroids sellected
    # for the clusters
    var <- extract(variable, cell_coord)
    Km <- kmeans(cbind(cell_coord, var), centers = NumAbsence)
    return(list(
      Centroids = Km$centers[, 1:2],
      Clusters = Km$cluster
    ))
    
    
  }
  
  # Inverse bioclim
  inv_bio <- function(e, p){
    Model <- dismo::bioclim(e, p)
    r <- dismo::predict(Model, e)
    names(r) <- "Group"
    r <- round(r, 5)
    r <- (r - minValue(r)) /
      (maxValue(r) - minValue(r))
    r <-(1-r)>=0.99 #environmental constrain
    r[which(r[,]==FALSE)] <- NA
    return(r)
  }
  
  # Inverse geo
  inv_geo <- function(e, p, d){
    Model <- dismo::circles(p, lonlat=T, d=d)
    r <- mask(e[[1]], Model@polygons, inverse=T)
    names(r) <- "Group"
    r[is.na(r)==F] <- 1 
    r[which(r[,]==FALSE)] <- NA
    return(r)
  }
  
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
  # ResultList <- rep(list(NULL),length(RecordsData))
  SpNames <- names(RecordsData)
  # BestGridList <- rep(list(NULL),length(RecordsData))
  
  #Start Cluster
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  
  # LOOP----
  results <- foreach(s = 1:length(RecordsData), .packages = c("raster","modEvA","ape","dismo")) %dopar% {
  # for(s in 1:length(RecordsData)){
  print(paste(s, SpNames[s]))
  # Extract coordinates----
  presences <- RecordsData[[s]]
  mask2 <- mask
  mask2[] <- 0
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
    if(type=="nearest"){
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
    }
    
    if(type=="all"){
      odd<-which((part3$group == 1))
      even<-which((part3$group == 2))
      dist <- as.matrix(dist(presences))
      dist <- 1/dist
      diag(dist) <- 0
      dist[which(dist==Inf)] <- 0
      dist[odd,odd] <- 0
      dist[even,even] <- 0
      mins <- apply(dist,2,max)
      for(i in 1:length(mins)){
        dist[,i] <- ifelse(dist[,i]==mins[i],mins[i],0)
      }
      species2 <- cbind(presences,pc1=extract(evnVariables[[1]], presences))
    }
    
    if(nrow(species2)<3){
      Imoran.Grid.P[p] <- NA
    }else{
      Imoran.Grid.P[p] <- 
        Moran.I(species2$pc1, dist, na.rm = T,scaled=T)$observed  
    }
  }
  
  Imoran.Grid.P <- abs(Imoran.Grid.P) # OJO estamos dejando todos los valores positivos
  N.grid <- 1:length(cellSize[pp])
  
  Opt <- data.frame(N.grid, cellSize[pp], round(data.frame(Imoran.Grid.P, Mess.Grid.P, Sd.Grid.P),3))
  # Cleaning those variances based in data divided in a number of partition less than 
  # the number of groups
  
  # SELLECTION OF THE BEST CELL SIZE----
  Opt2 <- Opt
  Dup <- !duplicated(Opt2[c("Imoran.Grid.P","Mess.Grid.P","Sd.Grid.P")])
  Opt2 <- Opt2[Dup,]
  
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
    
    if(unique(Opt2$Imoran.Grid.P) && unique(Opt2$Mess.Grid.P) && unique(Opt2$Sd.Grid.P)){
      Opt2 <- Opt2[nrow(Opt2),]
    }
  }
  
   if(nrow(Opt2)>1){
     Opt2 <- Opt2[nrow(Opt2),]
   }
  
  # Optimum size for presences
  print(Opt2)
  Optimum.Grid <- grid[[Opt2$N.grid]]
  Optimum.Grid@data[,c("C", "R")] <- NULL
  presences <- data.frame(Partition=part[,Opt2$N.grid], presences)

  #Save blocks raster
  pseudo.mask <- mask
  pseudo.mask2 <- list()
  RtoP <- data.frame(rasterToPoints(mask)[,-3])
  coordinates(RtoP)=c("x","y")
  crs(RtoP) <-projection(mask)
  FILTER <- over(RtoP, Optimum.Grid)
  pseudo.mask[which(pseudo.mask[]==1)] <- as.matrix(FILTER)
  # writeRaster(pseudo.mask, paste(DirSave, paste(SpNames[s],'.tif',sep=""),sep='/'),
  #             format = 'GTiff', NAflag = -9999, overwrite = TRUE)
  
  for(i in 1:N){
    mask3 <- pseudo.mask
    mask3[!mask3[]==i] <- NA 
    pseudo.mask2[[i]] <- mask3
  }
  
  pseudo.mask <- brick(pseudo.mask2)
  rm(pseudo.mask2)
  
#                                                          #
#              Pseudoabsences allocation-                  ####
#                                                          #

    # Pseudo-Absences with Random allocation-----
    if(pseudoabsencesMethod=="RND"){
      pseudo.mask_p <- pseudo.mask
      
      writeRaster(pseudo.mask_p, paste(DirSave, paste(SpNames[s],'.tif',sep=""),sep='/'),
                  format = 'GTiff', NAflag = -9999, overwrite = TRUE)
      
      # Random allocation of Pseudo-Absences 
      absences <- list()
      for (i in 1:N) {
        set.seed(s)
        if(MRst=="Y"){
          SpMask <- raster(file.path(DirM,paste0(SpNames[s],".tif")))
          pseudo.mask_p[[i]] <- pseudo.mask[[i]]*SpMask
          # if(sum(is.na(SpMask[])==F)<(PrAbRatio*nrow(RecordsData[[s]]))){
            # warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences") 
            # stop("Please try again with another restriction type or without restricting the extent")
          # }
        }
        absences.0 <- randomPoints(pseudo.mask[[i]], (1 / PrAbRatio)*sum(presences[,1]==i),
                                   ext = e,
                                   prob = FALSE)
        colnames(absences.0) <- c("x", "y")
        absences[[i]] <- as.data.frame(absences.0)
      }
      names(absences) <- 1:N
    }
  
    # Pseudo-Absences allocation with Environmental constrain ----
    if(pseudoabsencesMethod=="ENV_CONST"){
      
      pseudo.mask_p <- inv_bio(evnVariables, presences[,-1])
      
      # Split the raster of environmental layer with grids
      pseudo.mask_p <- mask(pseudo.mask, pseudo.mask_p)
      
      writeRaster(pseudo.mask_p, paste(DirSave, paste(SpNames[s],'.tif',sep=""),sep='/'),
                  format = 'GTiff', NAflag = -9999, overwrite = TRUE)
      
      absences <- list()
      for (i in 1:N) {
        set.seed(s)
        if(MRst=="Y"){
          SpMask <- raster(file.path(DirM, paste0(SpNames[s], ".tif")))
          pseudo.mask_p[[i]] <- pseudo.mask_p[[i]] * SpMask
          if(sum(is.na(SpMask[]) == F)<(PrAbRatio * nrow(RecordsData[[s]]))) {
            warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences") 
            stop("Please try again with another restriction type or without restricting the extent")
          }
        }
        absences.0 <- randomPoints(pseudo.mask_p[[i]], (1 / PrAbRatio)*sum(presences[,1]==i),
                                   ext = e,
                                   prob = FALSE)
        colnames(absences.0) <- c("lon", "lat")
        absences[[i]] <- as.data.frame(absences.0)
      }
      
      names(absences) <- 1:N
    }
  
    # Pseudo-Absences allocation with Geographical constrain-----
    if(pseudoabsencesMethod=="GEO_CONST"){
      
      pseudo.mask_p <- inv_geo(e=evnVariables, p=presences[,-1], d=Geo_Buf)
      
      # Split the raster of environmental layer with grids
      pseudo.mask_p <- mask(pseudo.mask, pseudo.mask_p)
      
      writeRaster(pseudo.mask_p, paste(DirSave, paste(SpNames[s],'.tif',sep=""),sep='/'),
                  format = 'GTiff', NAflag = -9999, overwrite = TRUE)
      
      absences <- list()
      for (i in 1:N) {
        set.seed(s)
        if (MRst == "Y") {
          SpMask <- raster(file.path(DirM, paste0(SpNames[s], ".tif")))
          pseudo.mask_p[[i]] <- pseudo.mask_p[[i]] * SpMask
          if(sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(RecordsData[[s]]))){
            warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences") 
            stop("Please try again with a smaller geographical buffer or without restricting the accessible area")
          }
        }
        absences.0 <- randomPoints(pseudo.mask_p[[i]], (1 / PrAbRatio)*sum(presences[,1]==i),
                                   ext = e,
                                   prob = FALSE)
        colnames(absences.0) <- c("lon", "lat")
        absences[[i]] <- as.data.frame(absences.0)
      }
      
      names(absences) <- 1:N
    }
    
    # Pseudo-Absences allocation with Environmentla and Geographical  constrain-----
    if(pseudoabsencesMethod=="GEO_ENV_CONST"){
      
    pseudo.mask_p <- inv_bio(evnVariables, presences[,-1])
    pseudo.mask_pg <- inv_geo(e=evnVariables, p=presences[,-1], d=Geo_Buf)
    pseudo.mask_p <- pseudo.mask_p*pseudo.mask_pg
    
    # Split the raster of environmental layer with grids
    pseudo.mask_p <- mask(pseudo.mask, pseudo.mask_p)
    
    writeRaster(pseudo.mask_p, paste(DirSave, paste(SpNames[s],'.tif',sep=""),sep='/'),
                format = 'GTiff', NAflag = -9999, overwrite = TRUE)
    
    absences <- list()
    for (i in 1:N) {
      set.seed(s)
      if(MRst=="Y"){
        SpMask <- raster(file.path(DirM, paste0(SpNames[s], ".tif")))
        pseudo.mask[[i]] <- pseudo.mask[[i]] * SpMask
        if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(RecordsData[[s]]))) {
          warning(
            "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
          )
          stop("Please try again with another restriction type or without restricting the extent")
        }
      }
      absences.0 <- randomPoints(pseudo.mask[[i]], (1 / PrAbRatio)*sum(presences[,1]==i),
                                 prob = FALSE)
      colnames(absences.0) <- c("lon", "lat")
      absences[[i]] <- as.data.frame(absences.0)
    }
    
    names(absences) <- 1:N
    }
  
    # Pseudo-Absences allocation with Environmentla and Geographical and k-mean constrain-----
    if(pseudoabsencesMethod=="GEO_ENV_KM_CONST"){
      
      pseudo.mask_p <- inv_bio(evnVariables, presences[,-1])
      pseudo.mask_pg <- inv_geo(e=evnVariables, p=presences[,-1], d=Geo_Buf)
      pseudo.mask_p <- pseudo.mask_p*pseudo.mask_pg
      
      # Split the raster of environmental layer with grids
      pseudo.mask_p <- mask(pseudo.mask, pseudo.mask_p)
      
      writeRaster(pseudo.mask_p, paste(DirSave, paste(SpNames[s],'.tif',sep=""),sep='/'),
                  format = 'GTiff', NAflag = -9999, overwrite = TRUE)
      
      absences <- list()
      for (i in 1:N) {
        set.seed(s)
        if(MRst=="Y"){
          SpMask <- raster(file.path(DirM,paste0(SpNames[s],".tif")))
          pseudo.mask[[i]] <- pseudo.mask[[i]]*SpMask
          if(sum(is.na(SpMask[])==F)<(PrAbRatio*nrow(RecordsData[[s]]))){
            warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences") 
            stop("Please try again with another restriction type or without restricting the extent")
          }
        }
        
        absences.0 <-
          KM(rasterToPoints(pseudo.mask_p[[i]])[,-3],
             mask(evnVariables, pseudo.mask_p[[i]]),
             (1 / PrAbRatio)*sum(presences[,1]==i))
        colnames(absences.0) <- c("lon", "lat")
        
        absences[[i]] <- as.data.frame(absences.0)
      }
      names(absences) <- 1:N
    }
  
    absences <- plyr::ldply(absences, data.frame)
    names(absences) <- c("Partition", "x", "y")
    absences[,c("x","y")] <- round(absences[,c("x","y")],4)
    colnames(absences) <- colnames(presences)
  # Final data.frame result----
    PresAbse <- rep(c(1, 0), sapply(list(presences, absences), nrow))
    result <- data.frame(Sp=SpNames[s], PresAbse, rbind(presences, absences), stringsAsFactors = F)
    result <- result[,c("Sp","x","y","Partition","PresAbse")]
    
    Opt2 <- data.frame(Sp=SpNames[s],Opt2)

  # Final data.frame result2----
    out <- list(ResultList= result,
                BestGridList = Opt2)
    return(out)
  }
  
  stopCluster(cl)
  FinalResult <- data.frame(data.table::rbindlist(do.call(c,lapply(results, "[", "ResultList"))))
  FinalInfoGrid <- data.frame(data.table::rbindlist(do.call(c,lapply(results, "[", "BestGridList"))))
  
  colnames(FinalResult) <- c("sp","x","y","Partition","PresAbse")
  write.table(FinalResult,paste(DirSave,"OccBlocks.txt",sep="\\"),sep="\t",row.names=F)
  write.table(FinalInfoGrid, paste(DirSave, "BestPartitions.txt", sep = '/'), sep="\t",
              col.names = T, row.names = F)
  
  return(FinalResult)
}
