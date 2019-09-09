Ensemble_TMLA <- function(DirR,
                          ValTable,
                          ThrTable,
                          PredictType,
                          RecordsData,
                          Threshold,
                          sensV,
                          Proj,
                          cores){
  
  
  #Start Cluster----
  cl <- makeCluster(cores,outfile="")
  registerDoParallel(cl)
  
  #Create Folders----
  DirENS <- paste(DirR,"Ensemble",sep="/")
  dir.create(DirENS)
    
  ensF <- paste(DirENS,PredictType[PredictType!="N"],sep="/")
  for(i in 1:length(ensF)){
    dir.create(ensF[i])
    assign(paste("Dir",PredictType[PredictType!="N"][i],sep=""),ensF[i])
  }
    
    #Binary ensemble directories
    ensFCat <- file.path(sort(rep(ensF,length(Threshold))),Threshold)
    for(i in 1:length(ensFCat)){
      dir.create(ensFCat[i])
      assign(paste("Dir",PredictType[PredictType!="N"][i],"Cat",sep=""),ensFCat[i])
    }
    
    #Projection directories
    if(!is.null(Proj)){
      ProjN <- list.files(file.path(DirR,"Projection"),full.names = T)
      sapply(ProjN, function(x) dir.create(file.path(x,"Ensemble"),recursive = T))
      ModFut <- file.path(ProjN,"Ensemble")
      for(i in 1:length(ModFut)){
        DirENS <- file.path(ModFut[i],PredictType)
        sapply(DirENS,function(x) dir.create(x,recursive = T))
        #Binary ensemble projection
        DirENSFCat <- file.path(sort(rep(DirENS,length(Threshold))),Threshold)
        sapply(DirENSFCat,function(x) dir.create(x,recursive = T))
      }
      DirENSFCat <- file.path(sort(rep(ModFut,length(PredictType))),PredictType,Threshold)
    }
  
  #Load Rasters----
  Folders <- list.files(file.path(DirR,"Algorithm"),full.names = T)
  spN <- substr(list.files(Folders[1],pattern="\\.tif$"),1,nchar(list.files(Folders[1],pattern="\\.tif$"))-4)
  if(!is.null(Proj)){
    ProjN <- as.list(ProjN)
    names(ProjN) <- list.files(file.path(DirR,"Projection"))
    ProjN <- lapply(ProjN,function(x) grep(list.files(x,full.names = T), pattern='Ensemble', inv=T, value=T))
  }
  
  #Thresholds List----
  ListSummary <- as.list(ensemble)
  names(ListSummary) <- ensemble
  
  
  #Evaluation Metric----
  if (any(PredictType%in%c("SUP","W_MEAN","PCA_SUP"))){
    cat("Calculate Superior Models based on which evaluation metric?\n(AUC/Kappa/TSS/Jaccard/Sorensen/Fpb)")
    Q1 <- readLines(n=1)
    while(!Q1%in%c("AUC","Kappa","TSS","Jaccard","Sorensen","Fpb")){
      cat("Please select a valid response,\n\n")
      cat("Calculate Superior Models based on which evaluation metric?\n(AUC/Kappa/TSS/Jaccard/Sorensen/Fpb)")
      Q1 <- readLines(n=1)
    }
  }
  
  
  #Species Loop----
  result <- foreach(s = 1:length(spN), 
                     .packages = c("raster", "dismo",
                                   "plyr", "RStoolbox" ),
                     .export = c( "PCA_ENS_TMLA", "STANDAR", "STANDAR_FUT", "Eval_Jac_Sor_TMLA", 
                                  "Thresholds_TMLA" )) %dopar% {
                                    
    #Species Occurrence Data & Evaluation Table----
    SpData <- RecordsData[RecordsData["sp"] == spN[s], ]
    SpVal <- ValTable[ValTable["Sp"] == spN[s], ]
    SpThr <- ThrTable[ThrTable["Sp"] == spN[s],] 
    
    #Raster Stack----
    ListRaster <- stack(file.path(Folders,paste0(spN[s],".tif")))
    names(ListRaster) <- list.files(file.path(DirR,"Algorithm"))
    
    ListFut <- lapply(ProjN, function(x) stack(file.path(x,paste0(spN[s],".tif"))))
    ListFut <- lapply(seq(ListFut), function(i) {
        y <- ListFut[[i]]
        names(y) <- list.files(file.path(DirR,"Algorithm"))
        return(y)
    })
    names(ListFut) <- list.files(file.path(DirR,"Projection"))
    
  ##%######################################################%##
  #                                                          #
  ####                  Ensemble Methods                  ####
  #                                                          #
  ##%######################################################%##
  
    # Mean Ensemble----
    if(any(PredictType=="MEAN")){
      Final <- brick(ListRaster)
      FinalT <- calc(Final,mean)
      Final <- STANDAR(FinalT)
      PredPoint <- extract(Final, SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Eval_JS <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                   a=PredPoint[PredPoint$PresAbse == 0, 2])
        
      #Final Thresholds
      Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
      ListSummary[["MEAN"]] <- data.frame(Sp=spN[s], Algorithm="MEA", Thr)
      writeRaster(Final, 
                  paste(DirMEAN, '/',spN[s],".tif", sep=""),
                  format='GTiff',
                  overwrite=TRUE)
      Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
      DirMEANCat <- grep(pattern="/MEAN/",x=ensFCat,value=T)
      for(t in 1:length(Thr_Alg)){
        writeRaster(Final>=Thr_Alg[t], 
                    paste(DirMEANCat[t], '/',spN[s],".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
          
        #Future Projection
        if(!is.null(Proj)){
          for(p in 1:length(ListFut)){
            Final <- brick(ListFut[[p]])
            Final <- calc(Final,mean)
            Final <- STANDAR_FUT(Final,FinalT)
            
            writeRaster(Final, 
                        file.path(ModFut[p],"MEAN",spN[s]),
                        format='GTiff',
                        overwrite=TRUE)
            DirMEANFCat <- grep(pattern=paste0(names(ListFut)[p],"/Ensemble/MEAN/",Threshold),x=DirENSFCat,value=T)
            for(t in 1:length(DirMEANFCat)){
              writeRaster(Final>=Thr_Alg[t], 
                          file.path(DirMEANFCat,paste(spN[s],sep="_")),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
    }
    
    # Weighted Mean Ensemble----
    if(any(PredictType=="W_MEAN")){
      
      ThResW <- unlist(SpVal[Q1])
      Final <- brick(ListRaster)
      Final <- calc(Final, function(x) x*ThResW)
      FinalT <- calc(Final,mean)
      Final <- STANDAR(FinalT)
      
      PredPoint <- extract(Final, SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Eval_JS <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                   a=PredPoint[PredPoint$PresAbse == 0, 2])
      
      #Final Thresholds
      Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
      ListSummary[["W_MEAN"]] <- data.frame(Sp=spN[s], Algorithm="WMEA", Thr)
      
      #Save Maps
      writeRaster(Final, 
                  paste(DirW_MEAN, '/',spN[s],".tif", sep=""),
                  format='GTiff',
                  overwrite=TRUE)
      Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
      DirMEANCat <- grep(pattern="/W_MEAN/",x=ensFCat,value=T)
      for(t in 1:length(Thr_Alg)){
        writeRaster(Final>=Thr_Alg[t], 
                    paste(DirMEANCat[t], '/',spN[s],".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
      }
     
      #Future Projection
      if(!is.null(Proj)){
        for(p in 1:length(ListFut)){
          Final <- brick(stack(ListFut[[p]]))
          Final <- calc(Final, function(x) x*ThResW)
          Final <- calc(Final,mean)
          Final <- STANDAR_FUT(Final,FinalT)
          
          writeRaster(Final, 
                      file.path(ModFut[p],"W_MEAN",spN[s]),
                      format='GTiff',
                      overwrite=TRUE)
          DirMEANFCat <- grep(pattern=paste0(names(ListFut)[p],"/Ensemble/W_MEAN/",Threshold),x=DirENSFCat,value=T)
          for(t in 1:length(Thr_Alg)){
            writeRaster(Final>=Thr_Alg[t], 
                        file.path(DirMEANFCat,paste0(spN[s],".tif")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }
    
    # Superior Ensemble----
    if(any(PredictType=='SUP')){
  
      Best <- SpVal[which(unlist(SpVal[Q1])>=mean(unlist(SpVal[Q1]))),"Algorithm"]
      W <- names(ListRaster)[names(ListRaster)%in%Best]
      
      Final <- brick(subset(ListRaster,subset=W))
      FinalT <- calc(Final,mean)
      Final <- STANDAR(FinalT)
      PredPoint <- raster::extract(Final, SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Eval_JS <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                   a=PredPoint[PredPoint$PresAbse == 0, 2])
      
      #Final Thresholds
      Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
      ListSummary[["SUP"]] <- data.frame(Sp=spN[s], Algorithm="SUP",Thr)
      
      #Save Final Maps  
      writeRaster(Final, 
                  paste(DirSUP, '/',spN[s],".tif", sep=""),
                  format='GTiff',
                  overwrite=TRUE)
      Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
      DirMEANCat <- grep(pattern="/SUP/",x=ensFCat,value=T)
      for(t in 1:length(Thr_Alg)){
        writeRaster(Final>=Thr_Alg[t], 
                    paste(DirMEANCat[t], '/',spN[s],".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
      }
        
      #Projection
      if(!is.null(Proj)){
        for(p in 1:length(ListFut)){
          Final <- brick(subset(ListFut[[p]],subset=W))
          Final <- calc(Final,mean)
          Final <- STANDAR_FUT(Final,FinalT)
          
          writeRaster(Final, 
                      file.path(ModFut[p],"SUP",paste(spN[s],sep="_")),
                      format='GTiff',
                      overwrite=TRUE)
          DirMEANFCat <- grep(pattern=paste0(names(ListFut)[p],"/Ensemble/SUP/",Threshold),x=DirENSFCat,value=T)
          for(t in 1:length(Thr_Alg)){
            writeRaster(Final>=Thr_Alg[t], 
                        file.path(DirMEANFCat,paste(spN[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }
    
    # With PCA ------
    if (any(PredictType == 'PCA')) {
  
      Final <- brick(ListRaster)
      Final <- PCA_ENS_TMLA(Final)
      PredPoint <- extract(Final, SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                   PredPoint[PredPoint$PresAbse == 0, 2])
      #Final Thresholds
      Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
      ListSummary[["PCA"]] <- data.frame(Sp=spN[s], Algorithm="PCA", Thr)
      
      #Save Final Maps
      writeRaster(Final, 
                  paste(DirPCA, '/',spN[s],".tif", sep=""),
                  format='GTiff',
                  overwrite=TRUE)
      Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
      DirMEANCat <- grep(pattern="/PCA/",x=ensFCat,value=T)
      for(t in 1:length(Thr_Alg)){
        writeRaster(Final>=Thr_Alg[t], 
                    paste(DirMEANCat[t], '/',spN[s],".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
      }
        
      #Future Projection
      if(!is.null(Proj)){
        for(p in 1:length(ListFut)){
          Final <- brick(ListFut[[p]])
          Final <- PCA_ENS_TMLA(Final)
          
          writeRaster(Final, 
                      file.path(ModFut[p],"PCA",spN[s]),
                      format='GTiff',
                      overwrite=TRUE)
          DirMEANFCat <- grep(pattern=paste0(names(ListFut)[p],"/Ensemble/PCA/",Threshold),x=DirENSFCat,value=T)
          for(t in 1:length(Thr_Alg)){
            writeRaster(Final>=Thr_Alg[t], 
                        file.path(DirMEANFCat,paste(spN[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }
    
    # With PCA over the Mean(Superior) Ensemble----
    if (any(PredictType == 'PCA_SUP')) {
      
      Best <- SpVal[which(unlist(SpVal[Q1])>=mean(unlist(SpVal[Q1]))),"Algorithm"]
      W <- names(ListRaster)[names(ListRaster)%in%Best]
      
      Final <- brick(subset(ListRaster,subset=W))
      Final <- PCA_ENS_TMLA(Final)
      PredPoint <- extract(Final, SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                   PredPoint[PredPoint$PresAbse == 0, 2])
      #Final Thresholds
      Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
      ListSummary[["PCS"]] <- data.frame(Sp=spN[s], Algorithm="PCS", Thr)
      
      #Save Final Maps
      writeRaster(Final, 
                  paste(DirPCA_SUP, '/',spN[s],".tif", sep=""),
                  format='GTiff',
                  overwrite=TRUE)
      Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
      DirMEANCat <- grep(pattern="/PCA_SUP/",x=ensFCat,value=T)
      for(t in 1:length(Thr_Alg)){
        writeRaster(Final>=Thr_Alg[t], 
                    paste(DirMEANCat[t], '/',spN[s],".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
      }
        
      #Future Projection
      if(!is.null(Proj)){
        for(p in 1:length(ListFut)){
          Final <- brick(subset(ListFut[[p]],subset=W))
          Final <- PCA_ENS_TMLA(Final)
          
          writeRaster(Final, 
                      file.path(ModFut[p],"PCA_SUP",paste(spN[s],sep="_")),
                      format='GTiff',
                      overwrite=TRUE)
          DirMEANFCat <- grep(pattern=paste0(names(ListFut)[p],"/Ensemble/PCA_SUP/",Threshold),x=DirENSFCat,value=T)
          for(t in 1:length(Thr_Alg)){
            writeRaster(Final>=Thr_Alg[t], 
                        file.path(DirMEANFCat,paste(spN[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }
    
    #With PCA over the threshold Ensemble----
    if (any(PredictType == 'PCA_THR')) {
      Final <- brick(ListRaster)
      ValidTHR <- SpThr[grepl(Threshold,SpThr[,"THR"],ignore.case = T),"THR_VALUE"]
      for(k in 1:nlayers(Final)){
        FinalSp <- Final[[k]]
        FinalSp[FinalSp<ValidTHR[k]] <- 0
        Final[[k]] <- FinalSp
      }
      Final <- PCA_ENS_TMLA(Final)
      PredPoint <- extract(Final, SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                   PredPoint[PredPoint$PresAbse == 0, 2])
      
      #Final Thresholds
      Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
      ListSummary[["PCT"]] <- data.frame(Sp=spN[s], Algorithm="PCT", Threshold=Thr)
        
      writeRaster(Final, 
                  paste(DirPCA_THR, '/',paste(spN[s],sep="_"),".tif", sep=""),
                  format='GTiff',
                  overwrite=TRUE)
      Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
      DirMEANCat <- grep(pattern="/PCA_THR/",x=ensFCat,value=T)
      for(t in 1:length(Thr_Alg)){
        writeRaster(Final>=Thr_Alg[t], 
                    paste(DirMEANCat[t], '/',spN[s],".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
      }
        
      #Future Projection
      if(!is.null(Proj)){
        for(p in 1:length(ListFut)){
          Final <- brick(ListFut[[p]])
          
          #Select only values above the Threshold
          for(k in 1:nlayers(ListRaster)){
            FinalSp <- Final[[k]]
            FinalSp[FinalSp<ValidTHR[k]] <- 0
            if(all(na.omit(FinalSp[])==0)){
              Final[Final[[k]]] <- NULL
            }else{
              Final[[k]] <- FinalSp
            }
          }
          Final <- PCA_ENS_TMLA(Final)
          
          writeRaster(Final, 
                      file.path(ModFut[p],"PCA_THR",paste(spN[s],sep="_")),
                      format='GTiff',
                      overwrite=TRUE)
          DirMEANFCat <- grep(pattern=paste0(names(ListFut)[p],"/Ensemble/PCA_THR/",Threshold),x=DirENSFCat,value=T)
          for(t in 1:length(Thr_Alg)){
            writeRaster(Final>=Thr_Alg[t], 
                        file.path(DirMEANFCat,paste(spN[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }
    
    #Final Data Frame Results
    result <- ldply(ListSummary,data.frame,.id=NULL)
    return(result)
  }
  # Save .txt with the models performance---- 
  FinalSummary <- result
  write.table(FinalSummary,paste(DirR,"Thresholds_Ensemble.txt", sep = '/'),sep="\t",
              col.names = T,row.names=F)
  stopCluster(cl)
}