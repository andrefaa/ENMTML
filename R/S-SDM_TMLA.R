S_SDM <- function(DirENM,
                  DirMSDM,
                  DirProj,
                  spN){
  #Function to create S-SDMs
  #Written by Andre Andrade
  
  #Parameters:
    #DirR: Result Folder
    #Fut: Projections Folder
    #MSDM: MSDM exists
  
  #Initialization
  
  #Start Cluster
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  
  #Create S-SDM Folders
  sapply(DirT,function(x) dir.create(file.path(x,"S_SDM")))
  DirSSDM <- file.path(DirT,"S_SDM")
  DirENM <- file.path(DirENM,"BIN")

  #Calculate S-SDM
  foreach(s = 1:length(DirSSDM), .packages = "raster") %dopar% {
    ENM <- sum(stack(list.files(DirENM[s],full.names=T)))
    writeRaster(ENM,file.path(DirSSDM[s],"S-SDM.tif"),format="GTiff",NAflag=-9999,overwrite=T)
  }
  
  #M-SDM folders
  if(!is.null(DirMSDM)){
    sapply(DirMSDM,function(x) dir.create(file.path(x,"S_SDM")))
    DirSSDM2 <- file.path(DirMSDM,"S_SDM")
    DirMSDM <- file.path(DirMSDM,"BIN")
    
    #Calculate S-SDM
    foreach(s = 1:length(DirSSDM2), .packages = "raster") %dopar% {
      ENM <- sum(stack(list.files(DirMSDM[s],full.names=T)))
      writeRaster(ENM,file.path(DirSSDM2[s],"S-SDM.tif"),format="GTiff",NAflag=-9999,overwrite=T)
    }
  }
  
  #Projection folders
  if(!is.null(DirProj)){
    sapply(DirProj,function(x) dir.create(file.path(x,"S_SDM")))
    DirSSDM3 <- file.path(DirProj,"S_SDM")
    DirProj <- file.path(DirProj,"BIN")
    
    #Calculate S-SDM
    foreach(s = 1:length(DirSSDM3), .packages = "raster") %dopar% {
      ENM <- sum(stack(list.files(DirProj[s],full.names=T)))
      writeRaster(ENM,file.path(DirSSDM3[s],"S-SDM.tif"),format="GTiff",NAflag=-9999,overwrite=T)
    }
  }
}
