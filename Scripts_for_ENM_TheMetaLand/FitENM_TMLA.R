## Written by Santiago Velazco & Andre Andrade

FitENM_TMLA <- function(RecordsData,
                   Variables,
                   Fut=NULL,
                   Part,
                   Algorithm,
                   PredictType,
                   spN,
                   Tst,
                   Threshold = Thresh,
                   DirSave=DirR, 
                   DirMask=NULL,
                   DirMSDM=DirPRI,
                   DirProj=NULL,
                   per=NULL,
                   repl=NULL,
                   Save="N") {
  
  Ti <- Sys.time()
  options(warn = -1)

  # Directory to save----
  folders <- paste(DirSave,Algorithm,sep="/")
  for(i in 1:length(folders)){
    dir.create(folders[i])
  }
  
  #Binary directories
  foldCat <- paste(folders,"BIN",sep="/")
  for(i in 1:length(foldCat)){
    dir.create(foldCat[i])
  }
  
  #Partition directories
  if(Save=="Y"){
    foldPart <- paste(folders,"PART",sep="/")
    for(i in 1:length(foldPart)){
      dir.create(foldPart[i])
      
    }
    #Binary directories
    PartCat <- paste(foldPart,"BIN",sep="/")
    for(i in 1:length(PartCat)){
      dir.create(PartCat[i])
    }
  }
  
  #Ensemble directory
  if(any(PredictType!="N")){
    DirENS <- paste(DirSave,"ENS",sep="/")
    dir.create(DirENS)
    
    ensF <- paste(DirENS,PredictType[PredictType!="N"],sep="/")
    for(i in 1:length(ensF)){
      dir.create(ensF[i])
      assign(paste("Dir",PredictType[PredictType!="N"][i],sep=""),ensF[i])
    }
  
    #Binary ensemble directories
    ensFCat <- paste(ensF,"BIN",sep="/")
    for(i in 1:length(ensFCat)){
      dir.create(ensFCat[i])
      assign(paste("Dir",PredictType[PredictType!="N"][i],"Cat",sep=""),ensFCat[i])
    }
    
    #Partitions directories
    if(Save=="Y"){
      ensPart <- paste(ensF,"PART",sep="/")
      for(i in 1:length(ensPart)){
        dir.create(ensPart[i])
      }
      #Binary directories
      ensCat <- paste(ensPart,"BIN",sep="/")
      for(i in 1:length(ensCat)){
        dir.create(ensCat[i])
      }
    }
  }
  
  #Projection directories
  if(is.null(Fut)==F){
    ProjN <- names(Fut)
    dir.create(file.path(DirSave,"FUT"))
    ModFut <- file.path(DirSave,"FUT",ProjN)
    for(i in 1:length(ModFut)){
      dir.create(ModFut[i])
      for(h in Algorithm){
        dir.create(file.path(ModFut[i],h))
        dir.create(file.path(ModFut[i],h,"BIN"))
      }
    }
    if(any(PredictType!="N")){
      for(j in 1:length(ModFut)){
        DirENS <- file.path(ModFut[j],"ENS")
        dir.create(DirENS)
        ensP <-file.path(DirENS,PredictType[PredictType!="N"])
        for(k in 1:length(ensP)){
          dir.create(ensP[k])
          dir.create(file.path(ensP[k],"BIN"))
        }
      }
    }
  }
  
  # Extracting enviromental variables----
  Ncol <- ncol(RecordsData) + 1
  RecordsData <- na.omit(data.frame(RecordsData, raster::extract(Variables, RecordsData[, c("x", "y")])))
  Ncol2 <- ncol(RecordsData)
  VarCol <- colnames(RecordsData[,Ncol:Ncol2])
  
  #Project to different geographical area or time period----
  VariablesP <- list()
  if(is.null(Fut)){
    VariablesP[[1]] <- Variables
  }else{
    VariablesP <- Fut
  }
  
  # Species names
  SpNames <- as.character(unique(RecordsData[, "sp"]))
  if(is.null(repl)==F){
    SpNames2 <- substr(SpNames,1,nchar(SpNames)-2)
  }else{
    SpNames2 <- SpNames
  }
  print(paste("Total species to be modeled", length(SpNames)))
  
  # Extent to predict the models
    e <- extent(Variables)
  
  # Number of partition
    N <- as.numeric(max(RecordsData[, "Partition"]))
      
  # Lists for validation-----
  for(i in 1:length(Algorithm)){
    assign(paste("Validation_", Algorithm[i], sep=""),list())
    assign(paste("Summary_", Algorithm[i], sep=""),list())
  }
  if(any(PredictType!="N")){
    for(i in 1:length(PredictType)){
      assign(paste("Validation_", PredictType[i], sep=""),list())
    assign(paste("Summary_", PredictType[i], sep=""),list())
    }
  }  
  
  #Txt of Final tables    
  if(is.null(repl)==F){
    VALNAME <- paste('Validation_Partition','_',sep="",repl,'.txt' )
  }else{
    VALNAME <- paste('Validation_Partition.txt' )
  }
  VALNAMEII <- paste('Thresholds_Complete.txt' )
  
  
  # Backqround point for Maxent_new algorithm----
  if (is.null(DirMask) == FALSE) {
     if ((any(Algorithm == "MXD") || any(Algorithm == "MXS"))) {
       RecordsDataM <- split(RecordsData,f=RecordsData$sp)
       RecordsDataMt <- list()
       
       for(i in 1:length(RecordsDataM)){
         RecordsDataMt[[i]] <- RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1,]
       }
       
       RecordsDataM <- RecordsDataMt
       names(RecordsDataM) <- SpNames
       rm(RecordsDataMt)
       
       for(i in 1:length(names(RecordsDataM))){
        msk <- raster(paste(DirMask,paste(spN[i],".tif",sep=""),sep="/"))
         if(Part=="boot"||Part=="cross"){
           NM <- 1
           RecordsDataM[[i]] <- rbind(RecordsDataM[[i]],RecordsData[RecordsData$sp==i & RecordsData$Partition==2 & RecordsData$PresAbse==0 ,])
         }else{
           NM <- max(RecordsDataM[[i]][,"Partition"])
         }
         
         for(x in 1:NM){
           msk2 <- msk
           msk2[!msk[]==x] <- NA 
           NCell <- sum(!is.na(msk2[]))
           if (NCell > 10000) {
             ab.0 <- data.frame(randomPoints(msk2,p=as.numeric(as.character(RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1, c("x", "y")])),10000))
             var.0 <- data.frame(extract(Variables,ab.0))
           }else{
             ab.0 <-
               randomPoints(msk2, p=as.numeric(as.character(RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1, c("x", "y")])),(NCell - nrow(RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1,])))
             var.0 <- data.frame(extract(Variables, ab.0))
           }
           ab.0 <- cbind(rep(spN[i],nrow(ab.0)),ab.0,rep(x,nrow(ab.0)),rep(0,nrow(ab.0)),var.0)
           ab.0 <- na.omit(ab.0)
           rm(var.0)
           colnames(ab.0) <- colnames(RecordsData)
           RecordsDataM[[i]] <- rbind(RecordsDataM[[i]],ab.0)
           rm(ab.0)
         }
       }
       rm(list=c("msk","msk2"))
       RecordsDataM <- ldply(RecordsDataM,data.frame,.id=NULL)
       cols <-  c("x","y","Partition","PresAbse",names(Variables))
       RecordsDataM[,cols] <-  apply(RecordsDataM[,cols], 2, function(x) as.numeric(as.character(x)))
     }
   }else{
    if ((any(Algorithm == "MXD") || any(Algorithm == "MXS"))) {
      RecordsDataM <- split(RecordsData,f=RecordsData$sp)
      RecordsDataMt <- list()

      for(i in 1:length(RecordsDataM)){
        RecordsDataMt[[i]] <- RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1,]
      }

      RecordsDataM <- RecordsDataMt
      names(RecordsDataM) <- SpNames
      rm(RecordsDataMt)

      for(i in names(RecordsDataM)){
        msk <- Variables[[1]]
        msk[is.na(msk[])==FALSE] <- 1

        if(Part=="boot"||Part=="cross"){
          NM <- 1
          RecordsDataM[[i]] <- rbind(RecordsDataM[[i]],RecordsData[RecordsData$sp==i & RecordsData$Partition==2 & RecordsData$PresAbse==0 ,])
        }else{
          NM <- max(RecordsDataM[[i]][,"Partition"])
        }

        for(x in 1:NM){
          msk2 <- msk
          msk2[!msk[]==x] <- NA 
          NCell <- sum(!is.na(msk2[]))
          if (NCell > 10000) {
            ab.0 <- data.frame(randomPoints(msk2,p=as.numeric(as.character(RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1, c("x", "y")])),10000))
            var.0 <- extract(Variables,ab.0)
          }else{
            ab.0 <-
              randomPoints(msk2, p=as.numeric(as.character(RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1, c("x", "y")])),(NCell - nrow(RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1,])))
            var.0 <- extract(Variables, ab.0)
          }
          ab.0 <- cbind(rep(i,nrow(ab.0)),ab.0,rep(x,nrow(ab.0)),rep(0,nrow(ab.0)),var.0)
          ab.0 <- na.omit(ab.0)
          rm(var.0)
          colnames(ab.0) <- colnames(RecordsData)
          RecordsDataM[[i]] <- rbind(RecordsDataM[[i]],ab.0)
          rm(ab.0)
        }
      }
      rm(list=c("msk","msk2"))
      RecordsDataM <- ldply(RecordsDataM,data.frame,.id=NULL)
      cols <-  c("x","y","Partition","PresAbse",names(Variables))  
      RecordsDataM[,cols] = apply(RecordsDataM[,cols], 2, function(x) as.numeric(as.character(x)))
    }
  }
  
  # List of threshold
  if (any(PredictType%in%c('Mean','Sup','PCA','PCASup','PCAThr'))) {
    ThresholdPresent <- list()
    THRNAME <- paste('Thresholds.txt', sep="")
  }

  # Construction of models LOOP-----
  for(s in 1:length(SpNames)){
    print(paste(s, SpNames[s], Sys.time()))
    
    ListRaster <- as.list(Algorithm)
    names(ListRaster) <- Algorithm
    
    #Create lists for the future
    if(is.null(Fut)==F){
       ListFut <- as.list(ProjN)
       names(ListFut) <- ProjN
       
       for(p in 1:length(ListFut)){
         ListFut[[p]] <- ListRaster
       }
    }

    # Ocurrences filtter by species and splited by partition----
    SpData <- RecordsData[RecordsData[, "sp"] == SpNames[s], ]
    if ((any(Algorithm == "MXD") | any(Algorithm == "MXS"))) {
      SpDataM <- RecordsDataM[RecordsDataM[, "sp"] == SpNames[s], ]
    }
    
    #Include MSDM----
    if(is.null(DirMSDM)==F){
      if(grepl("LatLong",DirMSDM)){
        MSDM <- stack(file.path(DirMSDM,list.files(DirMSDM,pattern=".tif")))
        names(MSDM) <- c("Lat","Long")
      }else{
        MSDM <- raster(file.path(DirMSDM,paste(SpNames[s],".tif",sep="")))
        names(MSDM) <- "MSDM"
      }
      SpDataT <- cbind(SpData,extract(MSDM,SpData[c("x","y")]))
      colnames(SpDataT) <- c(colnames(SpData),names(MSDM))
      VariablesT <- stack(Variables,MSDM)
      if ((any(Algorithm == "MXD") | any(Algorithm == "MXS"))) {
        SpDataTM <- cbind(SpDataM,extract(MSDM,SpDataM[c("x","y")]))
        colnames(SpDataTM) <- c(colnames(SpDataM),names(MSDM))
      }
      VarColT <- c(VarCol,names(MSDM))
    }else{
      VariablesT <- Variables
      VarColT <- VarCol
      SpDataT <- SpData
      SpDataTM <- SpDataM
    }
    
    #Define N due to Partition Method
    if(Part=="boot" || Part=="cross"){
      N <- 1
    }else{
      N <- N
    }
    
    # List of models for partial models-----
    RastPart <- as.list(Algorithm)
    for(i in 1:length(RastPart)){
      RastPart[[i]] <- as.list(rep(RastPart[[i]],N))
    }
    names(RastPart) <- Algorithm
      
    
    #Partition of presence data
    PAtrain <- list()
    PAtest <- list()
    for (i in 1:N) {
      PAtrain[[i]] <- SpDataT[SpDataT[, "Partition"] == i, ]
      PAtest[[i]] <- SpDataT[SpDataT[, "Partition"] != i, ]
    }
    
    #Maxent Input
    if ((any(Algorithm == "MXD") | any(Algorithm == "MXS"))) {
      PAtrainM <- list()
      PAtestM <- list()
      for (i in 1:N) {
        PAtrainM[[i]] <- SpDataTM[SpDataTM[, "Partition"] == i, ]
        PAtestM[[i]] <- SpDataTM[SpDataTM[, "Partition"] != i, ]
      }
      if(Part%in%c("boot","cross")){
        PAtestM <- PAtest
      }
    }

    #BIOCLIM (BIO)----- 
    if (any(Algorithm == "BIO")) {
      Model <- list()
      #BIO model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][PAtrain[[i]][, "PresAbse"] == 1,]
        Model[[i]] <- bioclim(dataPr[, VarColT])
      }
      
      #BIO evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        for (i in 1:N) {
          RastPart[["BIO"]][[i]] <- STANDAR(predict(Model[[i]], VariablesT))
          PredPoint <- extract(RastPart[["BIO"]][[i]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["BIO"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
        }
        
        #BIO threshold
        Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- rep(Threshold,N)
        
        #Evaluation 2.0
        res20 <- Validation2_0(Eval,Thr,PredPoint,RastPart[["BIO"]],N)
        
        #BIO result 
        Boyce <- mean(unlist(Boyce))
        Jac <- mean(unlist(res20["JAC"]))
        OPR <- mean(unlist(res20["OPR"]))
        Fcpb <- mean(unlist(res20["FCPB"]))
        Validation<-SUMMRES(Eval, N, Thr)
        if(is.null(repl)){
          Validation_BIO[[s]] <- data.frame(Sp=spN[s], Algorithm="BIO", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
        }else{
          Validation_BIO[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="BIO", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
        }
      
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(RastPart[["BIO"]][[i]],paste(grep("BIO",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[["BIO"]][[i]]>=Thr[p], 
                          paste(grep("BIO",PartCat,value=T), '/',SpNames[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }else{
            writeRaster(RastPart[["BIO"]][[i]],paste(grep("BIO",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[["BIO"]][[i]]>=Thr[p], 
                          paste(grep("BIO",PartCat,value=T), '/',SpNames[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
      
      # Save final model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Model <- bioclim(SpDataT[SpDataT[,"PresAbse"]==1, VarColT]) # only presences
        ListRaster[["BIO"]] <- STANDAR(predict(Model, VariablesT))
        PredPoint <- extract(ListRaster[["BIO"]], SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                         PredPoint[PredPoint$PresAbse == 0, 2]))
        Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- Threshold
        # Validation<-SUMMRES(Eval, 1, Thr)
        Summary_BIO[[s]] <- data.frame(Sp=spN[s], Algorithm="BIO", Threshold=Thr)
        if(is.null(Fut)==F){
          for(k in 1:length(VariablesP)){
            ListFut[[ProjN[k]]][["BIO"]] <- predict(VariablesP[[k]], Model)
            if(maxValue(ListFut[[ProjN[k]]][["BIO"]])==0){
              ListFut[[ProjN[k]]][["BIO"]] <- ListFut[[ProjN[k]]][["BIO"]]
            }else{
              ListFut[[ProjN[k]]][["BIO"]] <- STANDAR(ListFut[[ProjN[k]]][["BIO"]])
            }
          }
        }
      }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["BIO"]] <- predict(VariablesP[[k]],Model[[i]])
          if(maxValue(ListFut[[ProjN[k]]][["BIO"]])==0){
            ListFut[[ProjN[k]]][["BIO"]] <- ListFut[[ProjN[k]]][["BIO"]]
          }else{
            ListFut[[ProjN[k]]][["BIO"]] <- STANDAR(ListFut[[ProjN[k]]][["BIO"]])
          }
          
          PredPoint <- extract(ListFut[[ProjN[k]]][["BIO"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["BIO"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
        
          #BIO threshold
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- rep(Threshold,N)
        
          #Evaluation 2.0
          res20 <- Validation2_0(Eval,Thr,PredPoint,ListFut[[ProjN]][["BIO"]],N)
          
          #BIO result 
          Boyce <- mean(unlist(Boyce))
          Jac <- mean(unlist(res20["JAC"]))
          OPR <- mean(unlist(res20["OPR"]))
          Fcpb <- mean(unlist(res20["FCPB"]))
          Validation<-SUMMRES(Eval, N, Thr)
          if(is.null(repl)){
            Validation_BIO[[s]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="BIO", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }else{
            Validation_BIO[[s]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="BIO", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }
        }
      }
    }

    #MAXENT DEFAULT (MXD) ----
    if(any(Algorithm == 'MXD')){
      Model <- list()
      
      #MXD model
      for (i in 1:N) {
        dataPr <- PAtrainM[[i]]
        Model[[i]] <- maxnet(dataPr[,"PresAbse"], dataPr[,VarColT], f = 
                               maxnet.formula(dataPr[,"PresAbse"], 
                                              dataPr[,VarColT], classes="default"))
      }
      #MXD evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        for (i in 1:N) {
          RastPart[["MXD"]][[i]] <- STANDAR(predict(VariablesT,Model[[i]], clamp=F, type="cloglog"))
          PredPoint <- extract(RastPart[["MXD"]][[i]], PAtestM[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtestM[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["MXD"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        
        #MXD threshold
        Thr <- unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- rep(Threshold,N)
        
        #Evaluation 2.0
        res20 <- Validation2_0(Eval,Thr,PredPoint,RastPart[["MXD"]],N)
        
        #MXD result 
        Boyce <- mean(unlist(Boyce))
        Jac <- mean(unlist(res20["JAC"]))
        OPR <- mean(unlist(res20["OPR"]))
        Fcpb <- mean(unlist(res20["FCPB"]))
  
        Validation <- SUMMRES(Eval, N, Thr)
        if(is.null(repl)){
          Validation_MXD[[s]] <- data.frame(Sp=spN[s], Algorithm="MXD", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
        }else{
          Validation_MXD[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="MXD", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
        }
        
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            if(N!=1){
              writeRaster(RastPart[["MXD"]][[i]],paste(grep("MXD",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["MXD"]][[i]]>=Thr[p], 
                            paste(grep("MXD",PartCat,value=T), '/',SpNames[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(RastPart[["MXD"]][[i]],paste(grep("MXD",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["MXD"]][[i]]>=Thr[p], 
                            paste(grep("MXD",PartCat,value=T), '/',SpNames[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
          }
        }
        
        #Save final model
        if(per!=1 && repl==1 || per==1 || N!=1){
          Model <- maxnet(SpDataTM[,"PresAbse"], SpDataTM[,VarColT], f = 
                            maxnet.formula(SpDataTM[,"PresAbse"], SpDataTM[,VarColT], classes="default"))
          
          ListRaster[["MXD"]] <- STANDAR(predict(VariablesT,Model, clamp=F, type="cloglog"))
          PredPoint <- extract(ListRaster[["MXD"]], SpDataTM[, c("x", "y")])
          PredPoint <- data.frame(PresAbse = SpDataTM[, "PresAbse"], PredPoint)
          Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                           PredPoint[PredPoint$PresAbse == 0, 2]))
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- Threshold
          # Validation<-SUMMRES(Eval, 1, Thr)
          Summary_MXD[[s]] <- data.frame(Sp = spN[s], Algorithm = "MXD", Threshold=Thr)
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["MXD"]] <- predict(VariablesP[[k]], Model,clamp=F, type="cloglog")
              if(maxValue(ListFut[[ProjN[k]]][["MXD"]])==0){
                ListFut[[ProjN[k]]][["MXD"]] <- ListFut[[ProjN[k]]][["MXD"]]
              }else{
                ListFut[[ProjN[k]]][["MXD"]] <- STANDAR(ListFut[[ProjN[k]]][["MXD"]])
              }
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["MXD"]] <- predict(VariablesP[[k]],Model[[i]])
          if(maxValue(ListFut[[ProjN[k]]][["MXD"]])==0){
            ListFut[[ProjN[k]]][["MXD"]] <- ListFut[[ProjN[k]]][["MXD"]]
          }else{
            ListFut[[ProjN[k]]][["MXD"]] <- STANDAR(ListFut[[ProjN[k]]][["MXD"]])
          }
          
          PredPoint <- extract(ListFut[[ProjN[k]]][["MXD"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["MXD"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #BIO threshold
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- rep(Threshold,N)
          
          #Evaluation 2.0
          res20 <- Validation2_0(Eval,Thr,PredPoint,ListFut[[ProjN]][["MXD"]],N)
          
          #MXD result 
          Boyce <- mean(unlist(Boyce))
          Jac <- mean(unlist(res20["JAC"]))
          OPR <- mean(unlist(res20["OPR"]))
          Fcpb <- mean(unlist(res20["FCPB"]))
          Validation<-SUMMRES(Eval, N, Thr)
          if(is.null(repl)){
            Validation_MXD[[s]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="MXD", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }else{
            Validation_MXD[[s]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="MXD", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }
        }
      }
    }
    
    #MAXENT SIMPLES (MXS) ----
    if(any(Algorithm == 'MXS')){
      Model <- list()
      #MXS model
      for (i in 1:N) {
        dataPr <- PAtrainM[[i]]
        Model[[i]] <- maxnet2(dataPr[,"PresAbse"], dataPr[,VarColT], f = 
                                maxnet.formula(dataPr[,"PresAbse"],dataPr[,VarColT], classes="lq"))
      }
      #MXS evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        for (i in 1:N) {
          RastPart[["MXS"]][[i]] <- STANDAR(predict(VariablesT,Model[[i]], clamp=F, type="cloglog"))
          PredPoint <- extract(RastPart[["MXS"]][[i]], PAtestM[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtestM[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["MXS"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        #MXS threshold
        Thr <- unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- rep(Threshold,N)
        
        #Evaluation 2.0
        res20 <- Validation2_0(Eval,Thr,PredPoint,RastPart[["MXS"]],N)
        
        #MXS result 
        Boyce <- mean(unlist(Boyce))
        Jac <- mean(unlist(res20["JAC"]))
        OPR <- mean(unlist(res20["OPR"]))
        Fcpb <- mean(unlist(res20["FCPB"]))
        
        Validation <- SUMMRES(Eval, N, Thr)
        if(is.null(repl)){
          Validation_MXS[[s]] <- data.frame(Sp=spN[s], Algorithm="MXS", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
        }else{
          Validation_MXS[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="MXS", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
        }
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            if(N!=1){
              writeRaster(RastPart[["MXS"]][[i]],paste(grep("MXS",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["MXS"]][[i]]>=Thr[p], 
                            paste(grep("MXS",PartCat,value=T), '/',SpNames[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(RastPart[["MXS"]][[i]],paste(grep("MXS",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["MXS"]][[i]]>=Thr[p], 
                            paste(grep("MXS",PartCat,value=T), '/',SpNames[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
          }
        }
        
        #Save final model
        if(per!=1 && repl==1 || per==1 || N!=1){
          Model <- maxnet2(SpDataTM[,"PresAbse"], SpDataTM[,VarColT], f = 
                            maxnet.formula(SpDataTM[,"PresAbse"], SpDataTM[,VarColT], classes="lq"))
          ListRaster[["MXS"]] <- STANDAR(predict(VariablesT,Model, clamp=F, type="cloglog"))
          PredPoint <- extract(ListRaster[["MXS"]], SpDataTM[, c("x", "y")])
          PredPoint <- data.frame(PresAbse = SpDataTM[, "PresAbse"], PredPoint)
          Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                           PredPoint[PredPoint$PresAbse == 0, 2]))
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- Threshold
          # Validation<-SUMMRES(Eval, 1, Thr)
          Summary_MXS[[s]] <- data.frame(Sp = spN[s], Algorithm = "MXS", Threshold=Thr)
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["MXS"]] <- predict(VariablesP[[k]], Model,clamp=F, type="cloglog")
              if(maxValue(ListFut[[ProjN[k]]][["MXS"]])==0){
                ListFut[[ProjN[k]]][["MXS"]] <- ListFut[[ProjN[k]]][["MXS"]]
              }else{
                ListFut[[ProjN[k]]][["MXS"]] <- STANDAR(ListFut[[ProjN[k]]][["MXS"]])
              }
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["MXS"]] <- predict(VariablesP[[k]],Model[[i]])
          if(maxValue(ListFut[[ProjN[k]]][["MXS"]])==0){
            ListFut[[ProjN[k]]][["MXS"]] <- ListFut[[ProjN[k]]][["MXS"]]
          }else{
            ListFut[[ProjN[k]]][["MXS"]] <- STANDAR(ListFut[[ProjN[k]]][["MXS"]])
          }
          
          PredPoint <- extract(ListFut[[ProjN[k]]][["MXS"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["MXS"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #BIO threshold
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- rep(Threshold,N)
          
          #Evaluation 2.0
          res20 <- Validation2_0(Eval,Thr,PredPoint,ListFut[[ProjN]][["MXS"]],N)
          
          #MXS result 
          Boyce <- mean(unlist(Boyce))
          Jac <- mean(unlist(res20["JAC"]))
          OPR <- mean(unlist(res20["OPR"]))
          Fcpb <- mean(unlist(res20["FCPB"]))
          Validation<-SUMMRES(Eval, N, Thr)
          if(is.null(repl)){
            Validation_MXS[[s]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="MXS", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }else{
            Validation_MXS[[s]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="MXS", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }
        }
      }
    }
    
    #MAXIMUM LIKELIHOOD (MLK)----
    if(any(Algorithm == 'MLK')) {
      Model <- list()
      Fmula <- paste( " ~ ", paste(c(VarColT, paste("(",VarColT, ")^2", sep = "")),
                                   collapse = " + "), sep = "")
      Fmula <- as.formula(Fmula)
      #MLK model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][, c("x", "y")]
        Model[[i]] <- maxlike(Fmula, stack(VariablesT), dataPr,
                              link=c("logit"),
                              hessian = TRUE, removeDuplicates=FALSE)
      }
      
      #Evaluate Model
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        for (i in 1:N) {
          RastPart[["MLK"]][[i]] <- STANDAR(predict(VariablesT,Model[[i]]))
          PredPoint <- extract(RastPart[["MLK"]][[i]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["MLK"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        #MLK threshold
        Thr <- unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- rep(Threshold,N)
        
        #Evaluation 2.0
        res20 <- Validation2_0(Eval,Thr,PredPoint,RastPart[["MLK"]],N)
        
        #MLK result 
        Boyce <- mean(unlist(Boyce))
        Jac <- mean(unlist(res20["JAC"]))
        OPR <- mean(unlist(res20["OPR"]))
        Fcpb <- mean(unlist(res20["FCPB"]))
        
        Validation <- SUMMRES(Eval, N, Thr)
        if(is.null(repl)){
          Validation_MLK[[s]] <- data.frame(Sp=spN[s], Algorithm="MLK", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
        }else{
          Validation_MLK[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="MLK", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
        }
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            if(N!=1){
              writeRaster(RastPart[["MLK"]][[i]],paste(grep("MLK",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["MLK"]][[i]]>=Thr[p], 
                            paste(grep("MLK",PartCat,value=T), '/',SpNames[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(RastPart[["MLK"]][[i]],paste(grep("MLK",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["MLK"]][[i]]>=Thr[p], 
                            paste(grep("MLK",PartCat,value=T), '/',SpNames[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
          }
        }
        
        #Save final model
        if(per!=1 && repl==1 || per==1 || N!=1){
          Model <- maxlike(Fmula, stack(VariablesT), SpDataT[, c("x","y")],
                           link=c("logit"),
                           hessian = TRUE, removeDuplicates=FALSE)
          
          ListRaster[["MLK"]]<- STANDAR(predict(VariablesT, Model))
          PredPoint <- extract(ListRaster[["MLK"]], SpDataT[, c("x", "y")])
          PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2]))
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- Threshold
          # Validation<-SUMMRES(Eval, 1, Thr)
          Summary_MLK[[s]] <- data.frame(Sp = spN[s], Algorithm = "MLK", Threshold=Thr)
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["MLK"]] <- predict(VariablesP[[k]], Model)
              if(maxValue(ListFut[[ProjN[k]]][["MLK"]])==0){
                ListFut[[ProjN[k]]][["MLK"]] <- ListFut[[ProjN[k]]][["MLK"]]
              }else{
                ListFut[[ProjN[k]]][["MLK"]] <- STANDAR(ListFut[[ProjN[k]]][["MLK"]])
              }
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["MLK"]] <- predict(VariablesP[[k]],Model[[i]])
          if(maxValue(ListFut[[ProjN[k]]][["MLK"]])==0){
            ListFut[[ProjN[k]]][["MLK"]] <- ListFut[[ProjN[k]]][["MLK"]]
          }else{
            ListFut[[ProjN[k]]][["MLK"]] <- STANDAR(ListFut[[ProjN[k]]][["MLK"]])
          }
          
          PredPoint <- extract(ListFut[[ProjN[k]]][["MLK"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["MLK"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #BIO threshold
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- rep(Threshold,N)
          
          #Evaluation 2.0
          res20 <- Validation2_0(Eval,Thr,PredPoint,ListFut[[ProjN]][["MLK"]],N)
          
          #MLK result 
          Boyce <- mean(unlist(Boyce))
          Jac <- mean(unlist(res20["JAC"]))
          OPR <- mean(unlist(res20["OPR"]))
          Fcpb <- mean(unlist(res20["FCPB"]))
          Validation<-SUMMRES(Eval, N, Thr)
          if(is.null(repl)){
            Validation_MLK[[s]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="MLK", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }else{
            Validation_MLK[[s]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="MLK", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }
        }
      }
    }
    
    #SUPPORT VECTOR MACHINE (SVM)-----
    if (any(Algorithm == "SVM")) {
      Model <- list()
      Fmula <- formula(paste("PresAbse", '~ .'))
      
      #SVM model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][, c("PresAbse", VarColT)]
        set.seed(0)
        Model[[i]] <- ksvm(Fmula,data = dataPr,kernel = "rbfdot",C = 1, prob.model=T)
      }
      
      #SVM evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        for (i in 1:N) {
          RastPart[["SVM"]][[i]] <- STANDAR(predict(VariablesT,Model[[i]]))
          PredPoint <- extract(RastPart[["SVM"]][[i]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["SVM"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        
        #SVM threshold
        Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- rep(Threshold,N)
        
        #Evaluation 2.0
        res20 <- Validation2_0(Eval,Thr,PredPoint,RastPart[["SVM"]],N)
        
        #SVM result 
        Boyce <- mean(unlist(Boyce))
        Jac <- mean(unlist(res20["JAC"]))
        OPR <- mean(unlist(res20["OPR"]))
        Fcpb <- mean(unlist(res20["FCPB"]))
  
        Validation<-SUMMRES(Eval, N, Thr)
        if(is.null(repl)){
          Validation_SVM[[s]] <- data.frame(Sp=spN[s], Algorithm="SVM", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
        }else{
          Validation_SVM[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="SVM", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
        }
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            if(N!=1){
              writeRaster(RastPart[["SVM"]][[i]],paste(grep("SVM",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["SVM"]][[i]]>=Thr[p], 
                            paste(grep("SVM",PartCat,value=T), '/',SpNames[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(RastPart[["SVM"]][[i]],paste(grep("SVM",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["SVM"]][[i]]>=Thr[p], 
                            paste(grep("SVM",PartCat,value=T), '/',SpNames[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
          }
        }
        
        # Save final model
        if(per!=1 && repl==1 || per==1 || N!=1){
          # if(is.null(repl)){
          #   SpDataT <- SpDataT[SpDataT$Partition==1,]
          # }
          Model <- ksvm(Fmula,data = SpDataT[,c("PresAbse", VarColT)],
                        kernel = "rbfdot",C = 1, prob.model=T)
          ListRaster[["SVM"]] <- STANDAR(predict(VariablesT,Model))
          PredPoint <- extract(ListRaster[["SVM"]], SpDataT[, c("x", "y")])
          PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2]))
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- Threshold
          # Validation<-SUMMRES(Eval, 1, Thr)
          Summary_SVM[[s]] <- data.frame(Sp = spN[s], Algorithm = "SVM", Threshold=Thr)
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["SVM"]] <- predict(VariablesP[[k]], Model)
              if(maxValue(ListFut[[ProjN[k]]][["SVM"]])==0){
                ListFut[[ProjN[k]]][["SVM"]] <- ListFut[[ProjN[k]]][["SVM"]]
              }else{
                ListFut[[ProjN[k]]][["SVM"]] <- STANDAR(ListFut[[ProjN[k]]][["SVM"]])
              }
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["SVM"]] <- predict(VariablesP[[k]],Model[[i]])
          if(maxValue(ListFut[[ProjN[k]]][["SVM"]])==0){
            ListFut[[ProjN[k]]][["SVM"]] <- ListFut[[ProjN[k]]][["SVM"]]
          }else{
            ListFut[[ProjN[k]]][["SVM"]] <- STANDAR(ListFut[[ProjN[k]]][["SVM"]])
          }
          
          PredPoint <- extract(ListFut[[ProjN[k]]][["SVM"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["SVM"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #BIO threshold
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- rep(Threshold,N)
          
          #Evaluation 2.0
          res20 <- Validation2_0(Eval,Thr,PredPoint,ListFut[[ProjN]][["SVM"]],N)
          
          #SVM result 
          Boyce <- mean(unlist(Boyce))
          Jac <- mean(unlist(res20["JAC"]))
          OPR <- mean(unlist(res20["OPR"]))
          Fcpb <- mean(unlist(res20["FCPB"]))
          Validation<-SUMMRES(Eval, N, Thr)
          if(is.null(repl)){
            Validation_SVM[[s]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="SVM", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }else{
            Validation_SVM[[s]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="SVM", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }
        }
      }
    }
    
    #RANDOM FOREST (RDF)----
    if (any(Algorithm == "RDF")) {
      Model <- list()
      #RDF model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][, c("PresAbse", VarColT)]
        set.seed(1)
        Model[[i]] <- tuneRF(dataPr[,-1], (dataPr[,1]), trace=F,
                             stepFactor=2, ntreeTry=1000, doBest=T, plot=F)
      }
      #RDF evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        for (i in 1:N) {
          RastPart[["RDF"]][[i]] <- STANDAR(predict(VariablesT,Model[[i]]))
          PredPoint <- extract(RastPart[["RDF"]][[i]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["RDF"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        
        #RDF threshold
        Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- rep(Threshold,N)
        
        #Evaluation 2.0
        res20 <- Validation2_0(Eval,Thr,PredPoint,RastPart[["RDF"]],N)
        
        #RDF result 
        Boyce <- mean(unlist(Boyce))
        Jac <- mean(unlist(res20["JAC"]))
        OPR <- mean(unlist(res20["OPR"]))
        Fcpb <- mean(unlist(res20["FCPB"]))
        
        Validation<-SUMMRES(Eval, N, Thr)
        if(is.null(repl)){
          Validation_RDF[[s]] <- data.frame(Sp=spN[s], Algorithm="RDF", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
        }else{
          Validation_RDF[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="RDF", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
        }
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            if(N!=1){
              writeRaster(RastPart[["RDF"]][[i]],paste(grep("RDF",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["RDF"]][[i]]>=Thr[p], 
                            paste(grep("RDF",PartCat,value=T), '/',SpNames[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(RastPart[["RDF"]][[i]],paste(grep("RDF",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["RDF"]][[i]]>=Thr[p], 
                            paste(grep("RDF",PartCat,value=T), '/',SpNames[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
          }
        }
        
        # Save final model
        if(per!=1 && repl==1 || per==1 || N!=1){
          set.seed(0)
          Model <- tuneRF(SpDataT[,VarColT], (SpDataT[,"PresAbse"]), trace=F,
                          stepFactor=2, ntreeTry=500, doBest=T, plot = F)    
          ListRaster[["RDF"]] <- STANDAR(predict(VariablesT,Model))
          PredPoint <- extract(ListRaster[["RDF"]], SpDataT[, c("x", "y")])
          PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2]))
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- Threshold
          # Validation<-SUMMRES(Eval, 1, Thr)
          Summary_RDF[[s]] <- data.frame(Sp = spN[s], Algorithm = "RDF", Threshold=Thr)
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["RDF"]] <- predict(VariablesP[[k]], Model)
              if(maxValue(ListFut[[ProjN[k]]][["RDF"]])==0){
                ListFut[[ProjN[k]]][["RDF"]] <- ListFut[[ProjN[k]]][["RDF"]]
              }else{
                ListFut[[ProjN[k]]][["RDF"]] <- STANDAR(ListFut[[ProjN[k]]][["RDF"]])
              }
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["RDF"]] <- predict(VariablesP[[k]],Model[[i]])
          if(maxValue(ListFut[[ProjN[k]]][["RDF"]])==0){
            ListFut[[ProjN[k]]][["RDF"]] <- ListFut[[ProjN[k]]][["RDF"]]
          }else{
            ListFut[[ProjN[k]]][["RDF"]] <- STANDAR(ListFut[[ProjN[k]]][["RDF"]])
          }
          
          PredPoint <- extract(ListFut[[ProjN[k]]][["RDF"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["RDF"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #BIO threshold
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- rep(Threshold,N)
          
          #Evaluation 2.0
          res20 <- Validation2_0(Eval,Thr,PredPoint,ListFut[[ProjN]][["RDF"]],N)
          
          #RDF result 
          Boyce <- mean(unlist(Boyce))
          Jac <- mean(unlist(res20["JAC"]))
          OPR <- mean(unlist(res20["OPR"]))
          Fcpb <- mean(unlist(res20["FCPB"]))
          Validation<-SUMMRES(Eval, N, Thr)
          if(is.null(repl)){
            Validation_RDF[[s]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="RDF", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }else{
            Validation_RDF[[s]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="RDF", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }
        }
      }
    }
    
    #GENERALISED ADDITIVE MODEL (GAM)------
    if(any(Algorithm == 'GAM')) {
      Model <- list()
      Fmula <- paste("s(", VarColT,", k=3)", sep="")
      Fmula <- paste("PresAbse", paste(Fmula, collapse = " + "), sep = " ~ ")
      Fmula <- as.formula(Fmula)
      #GAM model
        
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][, c("PresAbse", VarColT)]
        Model[[i]] <- gam(Fmula, data = dataPr, optimizer = c("outer", "newton"), 
                          select = T, family = binomial)
        }
      #GAM evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        for (i in 1:N) {
          RastPart[["GAM"]][[i]] <- STANDAR(predict(VariablesT,Model[[i]]))
          PredPoint <- extract(RastPart[["GAM"]][[i]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["GAM"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        #GAM threshold
        Thr <- unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- rep(Threshold,N)
        
        #Evaluation 2.0
        res20 <- Validation2_0(Eval,Thr,PredPoint,RastPart[["GAM"]],N)
        
        #GAM result 
        Boyce <- mean(unlist(Boyce))
        Jac <- mean(unlist(res20["JAC"]))
        OPR <- mean(unlist(res20["OPR"]))
        Fcpb <- mean(unlist(res20["FCPB"]))
        
        Validation<-SUMMRES(Eval, N, Thr)
        if(is.null(repl)){
          Validation_GAM[[s]] <- data.frame(Sp=spN[s], Algorithm="GAM", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
        }else{
          Validation_GAM[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="GAM", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
        }
        
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            if(N!=1){
              writeRaster(RastPart[["GAM"]][[i]],paste(grep("GAM",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["GAM"]][[i]]>=Thr[p], 
                            paste(grep("GAM",PartCat,value=T), '/',SpNames[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(RastPart[["GAM"]][[i]],paste(grep("GAM",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["GAM"]][[i]]>=Thr[p], 
                            paste(grep("GAM",PartCat,value=T), '/',SpNames[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
          }
        }
        
        # Save final model
        if(per!=1 && repl==1 || per==1 || N!=1){
          Model <- gam(Fmula, data = SpDataT[, c("PresAbse",VarColT)], optimizer = c("outer", "newton"), 
                       select = T, family = binomial)
          ListRaster[["GAM"]] <- STANDAR(predict(VariablesT,Model))
          PredPoint <- extract(ListRaster[["GAM"]], SpDataT[, c("x", "y")])
          PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2]))
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- Threshold
          # Validation<-SUMMRES(Eval, 1, Thr)
          Summary_GAM[[s]] <- data.frame(Sp = spN[s], Algorithm = "GAM", Threshold=Thr)
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["GAM"]] <- predict(VariablesP[[k]], Model)
              if(maxValue(ListFut[[ProjN[k]]][["GAM"]])==0){
                ListFut[[ProjN[k]]][["GAM"]] <- ListFut[[ProjN[k]]][["GAM"]]
              }else{
                ListFut[[ProjN[k]]][["GAM"]] <- STANDAR(ListFut[[ProjN[k]]][["GAM"]])
              }
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["GAM"]] <- predict(VariablesP[[k]],Model[[i]])
          if(maxValue(ListFut[[ProjN[k]]][["GAM"]])==0){
            ListFut[[ProjN[k]]][["GAM"]] <- ListFut[[ProjN[k]]][["GAM"]]
          }else{
            ListFut[[ProjN[k]]][["GAM"]] <- STANDAR(ListFut[[ProjN[k]]][["GAM"]])
          }
          
          PredPoint <- extract(ListFut[[ProjN[k]]][["GAM"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["GAM"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #BIO threshold
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- rep(Threshold,N)
          
          #Evaluation 2.0
          res20 <- Validation2_0(Eval,Thr,PredPoint,ListFut[[ProjN]][["GAM"]],N)
          
          #GAM result 
          Boyce <- mean(unlist(Boyce))
          Jac <- mean(unlist(res20["JAC"]))
          OPR <- mean(unlist(res20["OPR"]))
          Fcpb <- mean(unlist(res20["FCPB"]))
          Validation<-SUMMRES(Eval, N, Thr)
          if(is.null(repl)){
            Validation_GAM[[s]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="GAM", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }else{
            Validation_GAM[[s]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="GAM", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }
        }
      }
    }
    
    #GENERALISED LINEAR MODEL (GLM) ------
    if(any(Algorithm == 'GLM')) {
      Model <- list()
      Fmula <- paste("PresAbse", paste(c(VarColT, paste("(",VarColT, ")^2", sep = "")),
                                     collapse = " + "), sep = " ~ ")
      Fmula <- as.formula(Fmula)
      #GLM model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][, c("PresAbse", VarColT)]
        Model[[i]] <- glm(Fmula, data = dataPr, family = binomial)
      }
      
      #GLM evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        for (i in 1:N) {
          RastPart[["GLM"]][[i]] <- STANDAR(predict(VariablesT,Model[[i]]))
          PredPoint <- extract(RastPart[["GLM"]][[i]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["GLM"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        
        #GLM threshold
        Thr <- unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- rep(Threshold,N)
        
        
        #Evaluation 2.0
        res20 <- Validation2_0(Eval,Thr,PredPoint,RastPart[["GLM"]],N)
        
        #GLM result 
        Boyce <- mean(unlist(Boyce))
        Jac <- mean(unlist(res20["JAC"]))
        OPR <- mean(unlist(res20["OPR"]))
        Fcpb <- mean(unlist(res20["FCPB"]))
        
        Validation<-SUMMRES(Eval, N, Thr)
        if(is.null(repl)){
          Validation_GLM[[s]] <- data.frame(Sp=spN[s], Algorithm="GLM", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
        }else{
          Validation_GLM[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="GLM", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
        }
        
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            if(N!=1){
              writeRaster(RastPart[["GLM"]][[i]],paste(grep("GLM",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["GLM"]][[i]]>=Thr[p], 
                            paste(grep("GLM",PartCat,value=T), '/',SpNames[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(RastPart[["GLM"]][[i]],paste(grep("GLM",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              for(p in 1:length(Thr)){
                writeRaster(RastPart[["GLM"]][[i]]>=Thr[p], 
                            paste(grep("GLM",PartCat,value=T), '/',SpNames[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
          }
        }
        
        # Save final model
        if(per!=1 && repl==1 || per==1 || N!=1){
          Model <- glm(Fmula, data = SpDataT[, c("PresAbse",VarColT)], family = binomial)
          ListRaster[["GLM"]] <- STANDAR(predict(VariablesT,Model))
          PredPoint <- extract(ListRaster[["GLM"]], SpDataT[, c("x", "y")])
          PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2]))
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- Threshold
          # Validation<-SUMMRES(Eval, 1, Thr)
          Summary_GLM[[s]] <- data.frame(Sp = spN[s], Algorithm = "GLM", Threshold=Thr)
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["GLM"]] <- predict(VariablesP[[k]], Model)
              if(maxValue(ListFut[[ProjN[k]]][["GLM"]])==0){
                ListFut[[ProjN[k]]][["GLM"]] <- ListFut[[ProjN[k]]][["GLM"]]
              }else{
                ListFut[[ProjN[k]]][["GLM"]] <- STANDAR(ListFut[[ProjN[k]]][["GLM"]])
              }
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["GLM"]] <- predict(VariablesP[[k]],Model[[i]])
          if(maxValue(ListFut[[ProjN[k]]][["GLM"]])==0){
            ListFut[[ProjN[k]]][["GLM"]] <- ListFut[[ProjN[k]]][["GLM"]]
          }else{
            ListFut[[ProjN[k]]][["GLM"]] <- STANDAR(ListFut[[ProjN[k]]][["GLM"]])
          }
          
          PredPoint <- extract(ListFut[[ProjN[k]]][["GLM"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["GLM"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #BIO threshold
          Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
          names(Thr) <- rep(Threshold,N)
          
          #Evaluation 2.0
          res20 <- Validation2_0(Eval,Thr,PredPoint,ListFut[[ProjN]][["GLM"]],N)
          
          #GLM result 
          Boyce <- mean(unlist(Boyce))
          Jac <- mean(unlist(res20["JAC"]))
          OPR <- mean(unlist(res20["OPR"]))
          Fcpb <- mean(unlist(res20["FCPB"]))
          Validation<-SUMMRES(Eval, N, Thr)
          if(is.null(repl)){
            Validation_GLM[[s]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="GLM", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }else{
            Validation_GLM[[s]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="GLM", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
          }
        }
      }
    }
    
    #GAUSSIAN (GAU) ----
    if(any(Algorithm == 'GAU')){
      Model <- list()
    #GAU model
    for (i in 1:N) {
      dataPr <- PAtrain[[i]]
      Model[[i]] <- graf(dataPr[,"PresAbse"], dataPr[,VarColT])
    }
      
    #GAU evaluation
    if((is.null(Fut)==F && Tst=="Y")==F){
      Eval <- list()
      Boyce <- list()
      for (i in 1:N) {
        RastPart[["GAU"]][[i]] <- STANDAR(predict(VariablesT,Model[[i]]))
        PredPoint <- extract(RastPart[["GAU"]][[i]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1,2],
                              PredPoint[PredPoint$PresAbse == 0,2])
        #Boyce Index
        Boyce[[i]] <- ecospat.boyce(RastPart[["GAU"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      #GAU threshold
      Thr <- unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      
      #Evaluation 2.0
      res20 <- Validation2_0(Eval,Thr,PredPoint,RastPart[["GAU"]],N)
      
      #GAU result 
      Boyce <- mean(unlist(Boyce))
      Jac <- mean(unlist(res20["JAC"]))
      OPR <- mean(unlist(res20["OPR"]))
      Fcpb <- mean(unlist(res20["FCPB"]))
      
      Validation<-SUMMRES(Eval, N, Thr)
      if(is.null(repl)){
        Validation_GAU[[s]] <- data.frame(Sp=spN[s], Algorithm="GAU", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
      }else{
        Validation_GAU[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="GAU", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
      }
      
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(RastPart[["GAU"]][[i]],paste(grep("GAU",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[["GAU"]][[i]]>=Thr[p], 
                          paste(grep("GAU",PartCat,value=T), '/',SpNames[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }else{
            writeRaster(RastPart[["GAU"]][[i]],paste(grep("GAU",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[["GAU"]][[i]]>=Thr[p], 
                          paste(grep("GAU",PartCat,value=T), '/',SpNames[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
      
      #Save final model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Model <- graf(SpDataT[,"PresAbse"], SpDataT[,VarColT])
        ListRaster[["GAU"]] <- STANDAR(predict.graf.raster(Model, VariablesT, type = "response", 
                                                           CI = 0.95, maxn = NULL)$posterior.mode)
        PredPoint <- extract(ListRaster[["GAU"]], SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2]))
        Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- Threshold
        # Validation<-SUMMRES(Eval, 1, Thr)
        Summary_GAU[[s]] <- data.frame(Sp = spN[s], Algorithm = "GAU", Threshold=Thr)
        if(is.null(Fut)==F){
          for(k in 1:length(VariablesP)){
            ListFut[[ProjN[k]]][["GAU"]] <- predict.graf.raster(Model, VariablesP[[k]], type = "response", 
                                                                CI = 0.95, maxn = NULL)$posterior.mode
            if(maxValue(ListFut[[ProjN[k]]][["GAU"]])==0){
              ListFut[[ProjN[k]]][["GAU"]] <- ListFut[[ProjN[k]]][["GAU"]]
            }else{
              ListFut[[ProjN[k]]][["GAU"]] <- STANDAR(ListFut[[ProjN[k]]][["GAU"]])
            }
          }
        }
      }
    }else{
      Eval <- list()
      Boyce <- list()
      for(k in 1:length(VariablesP)){
        ListFut[[ProjN[k]]][["GAU"]] <- predict(VariablesP[[k]],Model[[i]])
        if(maxValue(ListFut[[ProjN[k]]][["GAU"]])==0){
          ListFut[[ProjN[k]]][["GAU"]] <- ListFut[[ProjN[k]]][["GAU"]]
        }else{
          ListFut[[ProjN[k]]][["GAU"]] <- STANDAR(ListFut[[ProjN[k]]][["GAU"]])
        }
        
        PredPoint <- extract(ListFut[[ProjN[k]]][["GAU"]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        #Boyce Index
        Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["GAU"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        
        
        #GAU threshold
        Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- rep(Threshold,N)
        
        #Evaluation 2.0
        res20 <- Validation2_0(Eval,Thr,PredPoint,ListFut[[ProjN]][["GAU"]],N)
        
        #GAU result 
        Boyce <- mean(unlist(Boyce))
        Jac <- mean(unlist(res20["JAC"]))
        OPR <- mean(unlist(res20["OPR"]))
        Fcpb <- mean(unlist(res20["FCPB"]))
        Validation<-SUMMRES(Eval, N, Thr)
        if(is.null(repl)){
          Validation_GAU[[s]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="GAU", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
        }else{
          Validation_GAU[[s]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="GAU", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
        }
      }
    }
    }
    
    # Models performance----
    Obj <- ls(pattern = 'Validation_')
    ObjII <- ls(pattern = 'Summary_')
    result <- list()
    resultII <- list()
    for(i in 1:length(Obj)){
      result[[i]] <- ldply(get(Obj[i]))
      resultII[[i]] <- ldply(get(ObjII[i]))}
    result <- ldply(result,data.frame,.id=NULL)
    resultII <- ldply(resultII,data.frame,.id=NULL)
    
    # Ensemble-----
    # Without Ensemble
    SpValidation <- result[result$Sp==spN[s],]
    SpThr <- resultII[resultII$Sp==spN[s],]
    
    #Save final models
    if(per!=1 && repl==1 || per==1 || N!=1){
      if((is.null(Fut)==F && Tst=="Y")==F){
        for(i in 1:length(ListRaster)){
          writeRaster(round(ListRaster[[i]], 4), 
                      paste(folders[i], '/',spN[s],".tif", sep=""),
                      format='GTiff',
                      overwrite=TRUE)
          Thr <- SpThr[SpThr$Algorithm==names(ListRaster[i]), 'Threshold']
          writeRaster(ListRaster[[i]]>=Thr, 
                      paste(foldCat[i], '/',spN[s],".tif", sep=""),
                      format='GTiff',
                      overwrite=TRUE)
        }
      }
      
      #Save Projections
      if(is.null(Fut)==F){
        for(p in 1:length(ListFut)){
          for(o in 1:length(ListFut[[p]])){
            Thr <- SpThr[SpThr$Algorithm==names(ListFut[[p]])[o], 'Threshold']
              writeRaster(ListFut[[p]][[o]]>=Thr, 
                          file.path(ModFut[p],Algorithm[o],"BIN",paste(spN[s],sep="_")),
                          format='GTiff',
                          overwrite=TRUE)
            writeRaster(ListFut[[p]][[o]],file.path(ModFut[p],Algorithm[o],spN[s]),
                        format='GTiff',overwrite=TRUE)
          }
        }
      }
    }
    
    #Adjust invasion cenario for ensemble
    if((is.null(Fut)==F && Tst=="Y")){
      RastPart <- list(ListFut)
    }
    
    # Mean Ensemble----
    if(any(PredictType=="Mean")){
      
      #Partial Models Ensemble
      Final <- do.call(Map, c(c,RastPart))
      if(length(Final)>1){
        Final <- lapply(Final,function(x) brick(stack(x)))
        Final <- lapply(Final, function(x) STANDAR(round(mean(x),4)))
      }else{
        Final <- STANDAR(round(mean(brick(unlist(Final))),4))
      }
      
      # Threshold
      Eval <- list()
      Boyce <- list()
      for(i in 1:N){
        PredPoint <- extract(Final[[i]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      
      #Evaluation 2.0
      res20 <- Validation2_0(as.list(Eval),Thr,PredPoint,as.list(Final),N)
      
      #MEA result
      Boyce <- mean(unlist(Boyce))
      Jac <- mean(unlist(res20[names(res20)=="JAC"]))
      OPR <- mean(unlist(res20[names(res20)=="OPR"]))
      Fcpb <- mean(unlist(res20[names(res20)=="FCPB"]))
      
      Validation <- SUMMRES(Eval,N=1,Thr)
      if(is.null(repl)){
        Validation_Mean[[s]] <- data.frame(Sp=spN[s], Algorithm="MEA", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
      }else{
        Validation_Mean[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="MEA", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
      }
      
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(Final[[i]],paste(grep("\\bMean\\b",ensPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
              writeRaster(Final[[i]]>=Thr, 
                          paste(grep("\\bMean\\b",ensCat,value=T), '/',SpNames[s],"_",i,"_",".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
          }else{
            writeRaster(Final[[i]],paste(grep("\\bMean\\b",ensPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
              writeRaster(Final[[i]]>=Thr, 
                          paste(grep("\\bMean\\b",ensCat,value=T), '/',SpNames[s],"_",".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
          }
        }
      }
      
      #Final Model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Final <- brick(ListRaster)
        Final <- STANDAR(round(mean(Final),4))
        PredPoint <- extract(Final, SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2]))
        Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- Threshold
        
        Summary_Mean[[s]] <- data.frame(Sp=spN[s], Algorithm="MEA", Threshold=Thr)
        
        writeRaster(Final, 
                    paste(DirMean, '/',spN[s],".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
          writeRaster(Final>=as.numeric(Thr), 
                      paste(DirMeanCat, '/',spN[s],"_",".tif", sep=""),
                      format='GTiff',
                      overwrite=TRUE)
        
          #Future Projection
        if(is.null(Fut)==F && Tst!="Y"){
          for(p in 1:length(ListFut)){
            Final <- brick(ListFut[[p]])
            Final <- STANDAR(round(mean(Final),4))
            
            writeRaster(Final, 
                        file.path(ModFut[p],"ENS","Mean",SpNames2[s]),
                        format='GTiff',
                        overwrite=TRUE)
              writeRaster(Final>=as.numeric(Thr), 
                          file.path(ModFut[p],"ENS","Mean","BIN",paste(SpNames2[s],sep="_")),
                          format='GTiff',
                          overwrite=TRUE)
          }
        }
      }
    }

    # With Over the Mean(Superior) Ensemble----
    if(any(PredictType=='Sup')){
      SpValidation <- result[result$Sp==spN[s],]
      SpValidation$Algorithm <- as.character(SpValidation$Algorithm)
      
      Best <- which(SpValidation$TSS>=mean(SpValidation$TSS))
      Best <- SpValidation$Algorithm[Best]
        
      W <- names(ListRaster)%in%Best
      
      #Partial Models
      Final <- do.call(Map, c(c,RastPart))
      if(length(Final)>1){
        Final <- lapply(Final,function(x) brick(stack(x[W])))
        Final <- lapply(Final, function(x) STANDAR(round(mean(x),4)))
      }else{
        Final <- STANDAR(round(mean(brick(unlist(Final)[W])),4))
      }
      
      # Threshold
      Eval <- list()
      Boyce <- list()
      for(i in 1:N){
        PredPoint <- extract(Final[[i]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      
      #Evaluation 2.0
      res20 <- Validation2_0(as.list(Eval),Thr,PredPoint,as.list(Final),N)
      
      #SUP result
      Boyce <- mean(unlist(Boyce))
      Jac <- mean(unlist(res20[names(res20)=="JAC"]))
      OPR <- mean(unlist(res20[names(res20)=="OPR"]))
      Fcpb <- mean(unlist(res20[names(res20)=="FCPB"]))
      
      Validation <- SUMMRES(Eval,N=1,Thr)
      if(is.null(repl)){
        Validation_Sup[[s]] <- data.frame(Sp=spN[s], Algorithm="SUP", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
      }else{
        Validation_Sup[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="SUP", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
      }
      
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(Final[[i]],paste(grep("\\bSup\\b",ensPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
              writeRaster(Final[[i]]>=Thr, 
                          paste(grep("\\bSup\\b",ensCat,value=T), '/',SpNames[s],"_",i,"_",".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
          }else{
            writeRaster(Final[[i]],paste(grep("\\bSup\\b",ensPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
              writeRaster(Final[[i]]>=Thr, 
                          paste(grep("\\bSup\\b",ensCat,value=T), '/',SpNames[s],"_",".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
          }
        }
      }
      
      #Final Model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Final <- brick(ListRaster[W])
        Final <- STANDAR(round(mean(Final),4))
        PredPoint <- extract(Final, SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2]))
        Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- Threshold
        
        Summary_Sup[[s]] <- data.frame(Sp=spN[s], Algorithm="SUP", Threshold=Thr)
        
        writeRaster(Final, 
                    paste(DirSup, '/',paste(spN[s],sep="_"),".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
        writeRaster(Final>=as.numeric(Thr), 
                    paste(DirSupCat, '/',spN[s],"_",".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
        
        #Future Projection
        if(is.null(Fut)==F && Tst!="Y"){
          for(p in 1:length(ListFut)){
            Final <- brick(ListFut[[p]][W])
            Final <- STANDAR(round(mean(Final),4))
            
            writeRaster(Final, 
                        file.path(ModFut[p],"ENS","Sup",paste(SpNames2[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
            writeRaster(Final>=as.numeric(Thr), 
                        file.path(ModFut[p],"ENS","Sup","BIN",paste(SpNames2[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }

    # With PCA ------
    if (any(PredictType == 'PCA')) {

      
      #Partial Models Ensemble
      Final <- do.call(Map, c(c,RastPart))
      if(length(Final)>1){
        Final <- lapply(Final,function(x) brick(stack(x)))
        Final <- lapply(Final,function(x) PCA_ENS_TMLA(x))
      }else{
        Final <- brick(unlist(Final))
        Final <- PCA_ENS_TMLA(Final)
      }
      
      # Threshold
      Eval <- list()
      Boyce <- list()
      for(i in 1:N){
        PredPoint <- extract(Final[[i]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      
      #Evaluation 2.0
      res20 <- Validation2_0(as.list(Eval),Thr,PredPoint,as.list(Final),N)
      
      #PCA result
      Boyce <- mean(unlist(Boyce))
      Jac <- mean(unlist(res20[names(res20)=="JAC"]))
      OPR <- mean(unlist(res20[names(res20)=="OPR"]))
      Fcpb <- mean(unlist(res20[names(res20)=="FCPB"]))
      
      Validation <- SUMMRES(Eval,N=1,Thr)
      if(is.null(repl)){
        Validation_PCA[[s]] <- data.frame(Sp=spN[s], Algorithm="PCA", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
      }else{
        Validation_PCA[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="PCA", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
      }
      
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(Final[[i]],paste(grep("\\bPCA\\b",ensPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
              writeRaster(Final[[i]]>=Thr, 
                          paste(grep("\\bPCA\\b",ensCat,value=T), '/',SpNames[s],"_",i,"_",".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
          }else{
            writeRaster(Final[[i]],paste(grep("\\bPCA\\b",ensPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
              writeRaster(Final[[i]]>=Thr, 
                          paste(grep("\\bPCA\\b",ensCat,value=T), '/',SpNames[s],"_",".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
          }
        }
      }
      
      #Final Model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Final <- brick(ListRaster)
        Final <- PCA_ENS_TMLA(Final)
        PredPoint <- extract(Final, SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2]))
        Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- Threshold
        
        Summary_PCA[[s]] <- data.frame(Sp=spN[s], Algorithm="PCA", Threshold=Thr)
        
        writeRaster(Final, 
                    paste(DirPCA, '/',spN[s],".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
          writeRaster(Final>=as.numeric(Thr), 
                      paste(DirPCACat, '/',spN[s],"_",".tif", sep=""),
                      format='GTiff',
                      overwrite=TRUE)
        
        #Future Projection
        if(is.null(Fut)==F && Tst!="Y"){
          for(p in 1:length(ListFut)){
            Final <- brick(ListFut[[p]])
            Final <- PCA_ENS_TMLA(Final)
            
            writeRaster(Final, 
                        file.path(ModFut[p],"ENS","PCA",SpNames2[s]),
                        format='GTiff',
                        overwrite=TRUE)
              writeRaster(Final>=as.numeric(Thr), 
                          file.path(ModFut[p],"ENS","PCA","BIN",paste(SpNames2[s],sep="_")),
                          format='GTiff',
                          overwrite=TRUE)
          }
        }
      }
    }
        
    # With PCA over the Mean(Superior) Ensemble----
    if (any(PredictType == 'PCA_Sup')) {

      # Selection of best algorithms based on TSS
      SpValidation <- result[result$Sp==spN[s],]
      SpValidation$Algorithm <- as.character(SpValidation$Algorithm)
        
      Best <- which(SpValidation$TSS>=mean(SpValidation$TSS))
      Best <- SpValidation$Algorithm[Best]
        
      W <- names(ListRaster)%in%Best
      
      #Partial Models
      Final <- do.call(Map, c(c,RastPart))
      if(length(Final)>1){
        Final <- lapply(Final,function(x) brick(stack(x[W])))
        Final <- lapply(Final, function(x) PCA_ENS_TMLA(x))
      }else{
        Final <- brick(unlist(Final)[W])
        Final <- PCA_ENS_TMLA(Final)
      }
      
      # Threshold
      Eval <- list()
      Boyce <- list()
      for(i in 1:N){
        PredPoint <- extract(Final[[i]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      
      #Evaluation 2.0
      res20 <- Validation2_0(as.list(Eval),Thr,PredPoint,as.list(Final),N)
      
      #PCA_SUP result
      Boyce <- mean(unlist(Boyce))
      Jac <- mean(unlist(res20[names(res20)=="JAC"]))
      OPR <- mean(unlist(res20[names(res20)=="OPR"]))
      Fcpb <- mean(unlist(res20[names(res20)=="FCPB"]))
      
      Validation <- SUMMRES(Eval,N=1,Thr)
      if(is.null(repl)){
        Validation_PCA_Sup[[s]] <- data.frame(Sp=spN[s], Algorithm="PCS", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
      }else{
        Validation_PCA_Sup[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="PCS", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
      }
      
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(Final[[i]],paste(grep("\\bPCA_Sup\\b",ensPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            writeRaster(Final[[i]]>=Thr, 
                        paste(grep("\\bPCA_Sup\\b",ensCat,value=T), '/',SpNames[s],"_",i,"_",".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
          }else{
            writeRaster(Final[[i]],paste(grep("\\bPCA_Sup\\b",ensPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            writeRaster(Final[[i]]>=Thr, 
                        paste(grep("\\bPCA_Sup\\b",ensCat,value=T), '/',SpNames[s],"_",".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
      
      #Final Model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Final <- brick(ListRaster[W])
        Final <- PCA_ENS_TMLA(Final)
        PredPoint <- extract(Final, SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2]))
        Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- Threshold
        
        Summary_PCA_Sup[[s]] <- data.frame(Sp=spN[s], Algorithm="PCS", Threshold=Thr)
  
        writeRaster(Final, 
                    paste(DirPCA_Sup, '/',spN[s],"_",".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
        writeRaster(Final>=as.numeric(Thr), 
                    paste(DirPCA_SupCat, '/',spN[s],"_",".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
        
        #Future Projection
        if(is.null(Fut)==F && Tst!="Y"){
          for(p in 1:length(ListFut)){
            Final <- brick(ListFut[[p]][W])
            Final <- PCA_ENS_TMLA(Final)
            
            writeRaster(Final, 
                        file.path(ModFut[p],"ENS","PCA_Sup",paste(SpNames2[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
            writeRaster(Final>=as.numeric(Thr), 
                        file.path(ModFut[p],"ENS","PCA_Sup","BIN",paste(SpNames2[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }
    
    #With PCA over the threshold Ensemble----
    if (any(PredictType == 'PCA_Thr')) {
      
      SpValidation <- result[result$Sp==spN[s],]
      
      #Partial Models
      Final <- do.call(Map, c(c,RastPart))
      if(length(Final)>1){
        Final <- lapply(Final,function(x) brick(stack(x)))
        for(i in 1:length(Final)){
          for(k in Algorithm){
            FinalSp <- Final[[i]][[k]]
            FinalSp[FinalSp<SpValidation[SpValidation$Algorithm==k,"THR"]] <- 0
            Final[[i]][[k]] <- FinalSp
          }
        }
        Final <- lapply(Final, function(x) PCA_ENS_TMLA(x))
      }else{
        Final <- brick(unlist(Final))
        if(Tst=="Y"){
          names(Final) <- Algorithm
        }
        for(k in Algorithm){
          FinalSp <- Final[[k]]
          FinalSp[FinalSp<SpValidation[SpValidation$Algorithm==k,"THR"]] <- 0
          Final[[k]] <- FinalSp
        }
        Final <- PCA_ENS_TMLA(Final)
      }

      # Threshold
      Eval <- list()
      Boyce <- list()
      for(i in 1:N){
        PredPoint <- extract(Final[[i]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      
      #Evaluation 2.0
      res20 <- Validation2_0(as.list(Eval),Thr,PredPoint,as.list(Final),N)
      
      #PCA_THR result
      Boyce <- mean(unlist(Boyce))
      Jac <- mean(unlist(res20[names(res20)=="JAC"]))
      OPR <- mean(unlist(res20[names(res20)=="OPR"]))
      Fcpb <- mean(unlist(res20[names(res20)=="FCPB"]))
      
      Validation <- SUMMRES(Eval,N=1,Thr)
      if(is.null(repl)){
        Validation_PCA_Thr[[s]] <- data.frame(Sp=spN[s], Algorithm="PCT", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)
      }else{
        Validation_PCA_Thr[[s]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="PCT", Validation,Boyce=Boyce,Jack_ppac=Jac,OPR_ppac=OPR,Fcpb=Fcpb)          
      }
      
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(Final[[i]],paste(grep("\\bPCA_Thr\\b",ensPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            writeRaster(Final[[i]]>=Thr, 
                        paste(grep("\\bPCA_Thr\\b",ensCat,value=T), '/',SpNames[s],"_",i,"_",".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
          }else{
            writeRaster(Final[[i]],paste(grep("\\bPCA_Thr\\b",ensPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            writeRaster(Final[[i]]>=Thr, 
                        paste(grep("\\bPCA_Thr\\b",ensCat,value=T), '/',SpNames[s],"_",".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
      
      #Final Model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Final <- brick(ListRaster)
        for(k in Algorithm){
          FinalSp <- Final[[k]]
          FinalSp[FinalSp<SpValidation[SpValidation$Algorithm==k,"THR"]] <- 0
          Final[[k]] <- FinalSp
        }
        Final <- PCA_ENS_TMLA(Final)
        PredPoint <- extract(Final, SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2]))
        Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
        names(Thr) <- Threshold
        
        Summary_PCA_Thr[[s]] <- data.frame(Sp=spN[s], Algorithm="PCT", Threshold=Thr)
        
        writeRaster(Final, 
                    paste(DirPCA_Thr, '/',paste(spN[s],sep="_"),".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
        writeRaster(Final>=unlist(Thr), 
                      paste(DirPCA_ThrCat, '/',spN[s],"_",".tif", sep=""),
                      format='GTiff',
                      overwrite=TRUE)
          
        #Future Projection
        if(is.null(Fut)==F && Tst!="Y"){
          for(p in 1:length(ListFut)){
              Final <- brick(ListFut[[p]])
              
              #Select only values above the Threshold
              for(k in Algorithm){
                FinalSp <- Final[[k]]
                FinalSp[FinalSp<SpValidation[SpValidation$Algorithm==k,"THR"]] <- 0
                if(all(na.omit(FinalSp[])==0)){
                  Final[Final[[k]]] <- NULL
                }else{
                Final[[k]] <- FinalSp
                }
              }
              
              Final <- PCA_ENS_TMLA(Final)
  
              writeRaster(Final, 
                          file.path(ModFut[p],"ENS","PCA_Thr",paste(SpNames2[s],sep="_")),
                          format='GTiff',
                          overwrite=TRUE)
                writeRaster(Final>=unlist(Thr), 
                            file.path(ModFut[p],"ENS","PCA_Thr","BIN",paste(SpNames2[s],sep="_")),
                            format='GTiff',
                            overwrite=TRUE)
            }
        }
      }
    }
    
    # Save .txt with the models performance---- 
    
    #Partial Models
    Obj <- ls(pattern = 'Validation_')
    result <- list()
    for(i in 1:length(Obj)){
      result[[i]] <- ldply(get(Obj[i]))}
    result <- ldply(result)
    write.table(result, paste(DirSave, VALNAME, sep = '/'), sep="\t",
                col.names = T, row.names = F)
    
    #Full Models
    if(per!=1 && repl==1 || per==1 || N!=1){
    ObjII <- ls(pattern = 'Summary_')
    resultII <- list()
    for(i in 1:length(ObjII)){
      resultII[[i]] <- ldply(get(ObjII[i]))}
    resultII <- ldply(resultII)
    write.table(resultII, paste(DirSave, VALNAMEII, sep = '/'), sep="\t",
                col.names = T, row.names = F)
    }
  }#Fecha loop Especie
  
  # Save additional information and retuls----
  InfoModeling <- list(c("###########################################################"),
       paste('Start date :',Ti),
       paste('End date :',Sys.time()),
       c("Algorithm:", Algorithm),
       c("Ensemble:" , PredictType),
       c("Partition Method:" , Part),
       c("Train percentage (random partition only):", per),
       paste("PA Mask:" , DirMask),
       paste("MSDM:" , DirMSDM),
       paste("Resultados em:" , DirSave),
       paste('No_species:',length(SpNames)),
       paste("Threshold:",Threshold),
       matrix(spN))
  lapply(InfoModeling, write, 
         paste(DirSave, "/InfoModeling.txt", sep=""), append=TRUE, 
         ncolumns=20, sep='\t')
}
