# Function to save and project models in other places or times
PredictENM <- function(DirRDA,
                       DirSave = NA,
                       Validation = NA,
                       ThrPresent=NA,
                       PredictType = c('NoEnsemble', 'Mean', 'PCA'),
                       Variables,
                       RecordsData,
                       PresAbse,
                       Species,
                       x,
                       y,
                       Threshold = 'equal_sens_spec') {
  

  pkg<-c("dismo", "rgdal", "maptools", "raster", "XML","rJava", "kernlab", "SDMTools", "MASS",
         "randomForest", "maxlike","mgcv", "plyr", "GRaF", 'maxnet')
  sapply(pkg, require, character.only = TRUE)
  rm(pkg)
  
  rasterOptions(maxmemory = 1e+09)
  if(class(Variables)=='RasterStack'){
    Variables <- brick(Variables)
  }
  
  options(warn=-1)
  
  if(is.na(Validation)==TRUE){
    cat('Select the 1-Validation.txt archive')
    Validation <- read.table(file.choose(), header = T, stringsAsFactors = F)
  }
  
  if (any(PredictType == c('Mean', 'PCA'))) {
    if (is.na(ThrPresent) == TRUE) {
      cat('Select the 1-ThrPresent.txt archive')
      ThrPresent <- read.table(file.choose(), header = T, stringsAsFactors = F)
    }
  }
  
  if(is.na(DirSave) == TRUE){
    cat('Select a directory to saver output ')
    DirSave <- choose.dir()
  }
  
  if((length(PredictType)>1)==TRUE){
    stop(paste('PredictType argument must have one of these values:', 'NoEnsemble', 'Mean', 'PCA'))
  }
  
  SpNames <-unique(Validation[,Species])
  
  # Prediction of models------
  for (s in 1:length(SpNames)) {
    print(paste(SpNames[s], Sys.time()))
    # Model_List data per species 
    load(paste(DirRDA, '/', SpNames[s], '.rda', sep=""))
    # Algorithms names
    
    # Without Ensemble-----
    if(PredictType=='NoEnsemble'){
    
      ListRaster <- PREDICT(Variables = Variables, Models_List = Models_List)
      
      SpValidation <- Validation[Validation$Sp==SpNames[s],]
      
      for(i in 1:length(ListRaster)){
        writeRaster(ListRaster[[i]], 
                    paste(DirSave, '/',SpNames[s],"_",names(ListRaster[i]),'_Con',".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
        Thr <- SpValidation[SpValidation$Algorithm==names(ListRaster[i]), 'Thr']
        writeRaster(ListRaster[[i]]>=Thr, 
                    paste(DirSave, '/',SpNames[s],"_",names(ListRaster[i]),'_Cat',".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
      }
    }
    
    # With Mean Ensemble----
    if(PredictType=='Mean'){
      
      SpData <- RecordsData[RecordsData[,Species]==SpNames[s],]
      
      # Selection of best algorithms based on TSS
      SpValidation <- Validation[Validation$Sp==SpNames[s],]
      SpValidation$Algorithm <- as.character(SpValidation$Algorithm)
      
      if(any(SpValidation$Algorithm==c('MAXENTD','MAXENTS','MAXENTD_NEW','MAXENTS_NEW'))){
        SpValidation$Family <- gsub("D", "", gsub("S", "", 
                                                  gsub("_NEW", "", SpValidation$Algorithm)))
        Max <- max(SpValidation[SpValidation$Family=="MAXENT", 'TSS'])
        Del <- which(SpValidation[SpValidation$Family=="MAXENT", 'TSS']<Max)
        SpValidation[SpValidation$Family=="MAXENT", ][Del,] <- NA
        SpValidation <- SpValidation[complete.cases(SpValidation),] 
      }
      
      Best <- which(SpValidation$TSS>=mean(SpValidation$TSS))
      Best <- SpValidation$Algorithm[Best]
      print(Best)
      ListRaster <- PREDICT(Variables = Variables, Models_List = Models_List[names(Models_List)%in%Best])
      Final <- brick(ListRaster)
      Final <- STANDAR(Final)
      Final <- round(mean(Final), 4)
      # Threshold
      Thr <- ThrPresent[ThrPresent[,1]==SpNames[s],2]
      # Save raster
      writeRaster(Final, 
                  paste(DirSave, '/',SpNames[s],"_","Mean",'_Con',".tif", sep=""),
                  format='GTiff',
                  overwrite=TRUE)
      writeRaster(Final>=Thr, 
                  paste(DirSave, '/',SpNames[s],"_","Mean",'_Cat',".tif", sep=""),
                  format='GTiff',
                  overwrite=TRUE)
    }
    
    # With PCA Ensemble------
    if (PredictType == 'PCA') {
      SpData <- RecordsData[RecordsData[,Species]==SpNames[s],]
      
      # Selection of best algorithms based on TSS
      SpValidation <- Validation[Validation$Sp==SpNames[s],]
      SpValidation$Algorithm <- as.character(SpValidation$Algorithm)
      
      if(any(SpValidation$Algorithm==c('MAXENTD','MAXENTS','MAXENTD_NEW','MAXENTS_NEW'))){
        SpValidation$Family <- gsub("D", "", gsub("S", "", 
                                                  gsub("_NEW", "", SpValidation$Algorithm)))
        Max <- max(SpValidation[SpValidation$Family=="MAXENT", 'TSS'])
        Del <- which(SpValidation[SpValidation$Family=="MAXENT", 'TSS']<Max)
        SpValidation[SpValidation$Family=="MAXENT", ][Del,] <- NA
        SpValidation <- SpValidation[complete.cases(SpValidation),] 
      }
      
      Best <- which(SpValidation$TSS>=mean(SpValidation$TSS))
      Best <- SpValidation$Algorithm[Best]
      print(Best)
      ListRaster <- PREDICT(Variables = Variables, Models_List = Models_List[names(Models_List)%in%Best])
      Final <- brick(ListRaster)
      
      #Extracting values
      df<-rasterToPoints(Final)
      df<-na.omit(df)
      df<-df[,-c(1:2)]
      #Scale transform 
      df <- data.frame(apply(df,2,scale))
      #PCA
      data.pca <- prcomp(df, retx=TRUE)
      Values <-round(data.pca$x[,1], 4)
      Final <- Final[[1]]
      Final[complete.cases(Final[])] <- Values
      Final <- STANDAR(Final)
      
      # Threshold
      Thr <- ThrPresent[ThrPresent[,1]==SpNames[s],2]
      # Save raster
      writeRaster(Final, 
                  paste(DirSave, '/',SpNames[s],"_","PCA",'_Con',".tif", sep=""),
                  format='GTiff',
                  overwrite=TRUE)
      writeRaster(Final>=Thr, 
                  paste(DirSave, '/',SpNames[s],"_","PCA",'_Cat',".tif", sep=""),
                  format='GTiff',
                  overwrite=TRUE)
    }
  }
}
