## Written by Santiago Velazco & Andre Andrade

FitENM_TMLA <- function(RecordsData,
                   Variables,
                   Fut=NULL,
                   Part,
                   Algorithm,
                   PredictType,
                   Threshold = Thresh,
                   DirSave=DirR, 
                   DirMask=NULL,
                   DirMSDM=DirPRI,
                   DirProj=NULL,
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
    assign(paste("Summary_", PredictType[i], sep=""),list())
    }
  }  
  
  # List of models for trasferability-----
    for(i in 1:length(Algorithm)){
      Models_List <- as.list(c(Algorithm,PredictType))
      names(Models_List) <- c(Algorithm,PredictType)
    }
    
  if(is.null(repl)==F){
    VALNAME <- paste('Validation_Partition','_',sep="",repl,'.txt' )
    VALNAMEII <- paste('Validation_Complete','_',sep="",repl,'.txt' )
  }else{
    VALNAME <- paste('Validation_Partition.txt' )
    VALNAMEII <- paste('Validation_Complete.txt' )
  }
  
  
  # Backqround point for Maxent_new algorithm
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
       
       for(i in names(RecordsDataM)){
         if(is.null(repl)==F){
           iT <- substr(i,1,nchar(i)-2)
           msk <- raster(paste(DirMask,paste(iT,".tif",sep=""),sep="/"))
         }else{
          msk <- raster(paste(DirMask,paste(i,".tif",sep=""),sep="/"))
         }
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
    
    #Include MSDM
    if(is.null(DirMSDM)==F){
      msdm <- raster(paste(DirMSDM,paste(SpNames[s],".tif",sep=""),sep="/"))
      SpData <- cbind(SpData,extract(msdm,SpData[c("x","y")]))
      colnames(SpData)[ncol(SpData)] <- "MSDM"
      if ((any(Algorithm == "MXD") | any(Algorithm == "MXS"))) {
        SpDataM <- cbind(SpDataM,extract(msdm,SpDataM[c("x","y")]))
        colnames(SpDataM)[ncol(SpDataM)] <- "MSDM"
      }
    }
    
    #Define N due to Partition Method
    if(Part=="boot" || Part=="cross"){
      N <- 1
    }else{
      N <- N
    }
    
    #Partition of presence data
    PAtrain <- list()
    PAtest <- list()
    for (i in 1:N) {
      PAtrain[[i]] <- SpData[SpData[, "Partition"] == i, ]
      PAtest[[i]] <- SpData[SpData[, "Partition"] != i, ]
    }
    
    #Maxent Input
    if ((any(Algorithm == "MXD") | any(Algorithm == "MXS"))) {
      PAtrainM <- list()
      PAtestM <- list()
      for (i in 1:N) {
        PAtrainM[[i]] <- SpDataM[SpDataM[, "Partition"] == i, ]
        PAtestM[[i]] <- SpDataM[SpDataM[, "Partition"] != i, ]
      } 
    }

    #BIOCLIM (BIO)----- 
    if (any(Algorithm == "BIO")) {
      Model <- list()
      #BIO model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][PAtrain[[i]][, "PresAbse"] == 1,]
        Model[[i]] <- bioclim(dataPr[, VarCol])
      }
      
      #BIO evaluation
      Eval <- list()
      RastPart <- list()
      for (i in 1:N) {
        RastPart[[i]] <- STANDAR(predict(Model[[i]], Variables))
        PredPoint <- extract(RastPart[[i]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      }
      
      #BIO threshold
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      #BIO result 
      Validation<-SUMMRES(Eval, N, Thr)
      Validation_BIO[[s]] <- data.frame(Sp=SpNames[s], Algorithm="BIO", Validation)
      
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(RastPart[[i]],paste(grep("BIO",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("BIO",PartCat,value=T), '/',SpNames[s],"_",i,"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }else{
            writeRaster(RastPart[[i]],paste(grep("BIO",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("BIO",PartCat,value=T), '/',SpNames[s],"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
      
      # Save final model
      Model <- bioclim(SpData[SpData[,"PresAbse"]==1, VarCol]) # only presences
      ListRaster[["BIO"]] <- STANDAR(predict(Model, Variables))
      PredPoint <- extract(ListRaster[["BIO"]], SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                       PredPoint[PredPoint$PresAbse == 0, 2]))
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- Threshold
      Validation<-SUMMRES(Eval, 1, Thr)
      Summary_BIO[[s]] <- data.frame(Sp=SpNames[s], Algorithm="BIO", Validation)
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

    #MAXENT DEFAULT (MXD) ----
    if(any(Algorithm == 'MXD')){
      Model <- list()
      
      #MXD model
      for (i in 1:N) {
        dataPr <- PAtrainM[[i]]
        Model[[i]] <- maxnet(dataPr[,"PresAbse"], dataPr[,VarCol], f = 
                               maxnet.formula(dataPr[,"PresAbse"], 
                                              dataPr[,VarCol], classes="default"))
      }
      #MXD evaluation
      Eval <- list()
      RastPart <- list()
      for (i in 1:N) {
        RastPart[[i]] <- STANDAR(predict(Variables,Model[[i]], clamp=F, type="cloglog"))
        PredPoint <- extract(RastPart[[i]], PAtestM[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtestM[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      }
      #MXD threshold
      Thr <- unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      #MXD result
      Validation <- SUMMRES(Eval, N, Thr)
      Validation_MXD[[s]] <- data.frame(Sp=SpNames[s], Algorithm="MXD", Validation)
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(RastPart[[i]],paste(grep("MXD",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("MXD",PartCat,value=T), '/',SpNames[s],"_",i,"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }else{
            writeRaster(RastPart[[i]],paste(grep("MXD",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("MXD",PartCat,value=T), '/',SpNames[s],"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
      
      #Save final model
      Model <- maxnet(SpDataM[,"PresAbse"], SpDataM[,VarCol], f = 
                        maxnet.formula(SpDataM[,"PresAbse"], SpDataM[,VarCol], classes="default"))
      
      ListRaster[["MXD"]] <- STANDAR(predict(Variables,Model, clamp=F, type="cloglog"))
      PredPoint <- extract(ListRaster[["MXD"]], SpDataM[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpDataM[, "PresAbse"], PredPoint)
      Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                       PredPoint[PredPoint$PresAbse == 0, 2]))
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- Threshold
      Validation<-SUMMRES(Eval, 1, Thr)
      Summary_MXD[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "MXD", Validation)
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
    
    #MAXENT SIMPLES (MXS) ----
    if(any(Algorithm == 'MXS')){
      Model <- list()
      #MXS model
      for (i in 1:N) {
        dataPr <- PAtrainM[[i]]
        Model[[i]] <- maxnet2(dataPr[,"PresAbse"], dataPr[,VarCol], f = 
                                maxnet.formula(dataPr[,"PresAbse"],dataPr[,VarCol], classes="lq"))
      }
      #MXS evaluation
      Eval <- list()
      RastPart <- list()
      for (i in 1:N) {
        RastPart[[i]] <- STANDAR(predict(Variables,Model[[i]], clamp=F, type="cloglog"))
        PredPoint <- extract(RastPart[[i]], PAtestM[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtestM[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
      }
      #MXS threshold
      Thr <- unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      #MXS result
      Validation <- SUMMRES(Eval, N, Thr)
      Validation_MXS[[s]] <- data.frame(Sp=SpNames[s], Algorithm="MXS", Validation)
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(RastPart[[i]],paste(grep("MXS",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("MXS",PartCat,value=T), '/',SpNames[s],"_",i,"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }else{
            writeRaster(RastPart[[i]],paste(grep("MXS",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("MXS",PartCat,value=T), '/',SpNames[s],"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
      
      #Save final model
      Model <- maxnet2(SpDataM[,"PresAbse"], SpDataM[,VarCol], f = 
                        maxnet.formula(SpDataM[,"PresAbse"], SpDataM[,VarCol], classes="lq"))
      ListRaster[["MXS"]] <- STANDAR(predict(Variables,Model, clamp=F, type="cloglog"))
      PredPoint <- extract(ListRaster[["MXS"]], SpDataM[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpDataM[, "PresAbse"], PredPoint)
      Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                       PredPoint[PredPoint$PresAbse == 0, 2]))
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- Threshold
      Validation<-SUMMRES(Eval, 1, Thr)
      Summary_MXS[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "MXS", Validation)
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
    
    #MAXIMUM LIKELIHOOD (MLK)----
    if(any(Algorithm == 'MLK')) {
      Model <- list()
      Fmula <- paste( " ~ ", paste(c(VarCol, paste("(",VarCol, ")^2", sep = "")),
                                   collapse = " + "), sep = "")
      Fmula <- as.formula(Fmula)
      #MLK model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][, c("x", "y")]
        Model[[i]] <- maxlike(Fmula, stack(Variables), dataPr,
                              link=c("logit"),
                              hessian = TRUE, removeDuplicates=FALSE)
      }
      
      #Evaluate Model
      Eval <- list()
      RastPart <- list()
      for (i in 1:N) {
        RastPart[[i]] <- STANDAR(predict(Variables,Model[[i]]))
        PredPoint <- extract(RastPart[[i]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      }
      #MLK threshold
      Thr <- unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      #MLK result
      Validation <- SUMMRES(Eval, N, Thr)
      Validation_MLK[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "MLK", Validation)
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(RastPart[[i]],paste(grep("MLK",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("MLK",PartCat,value=T), '/',SpNames[s],"_",i,"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }else{
            writeRaster(RastPart[[i]],paste(grep("MLK",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("MLK",PartCat,value=T), '/',SpNames[s],"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
      
      #Save final model
      Model <- maxlike(Fmula, stack(Variables), SpData[, c("x","y")],
                       link=c("logit"),
                       hessian = TRUE, removeDuplicates=FALSE)
      
      ListRaster[["MLK"]]<- STANDAR(predict(Variables, Model))
      PredPoint <- extract(ListRaster[["MLK"]], SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                   PredPoint[PredPoint$PresAbse == 0, 2]))
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- Threshold
      Validation<-SUMMRES(Eval, 1, Thr)
      Summary_MLK[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "MLK", Validation)
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
    
    #SUPPORT VECTOR MACHINE (SVM)-----
    if (any(Algorithm == "SVM")) {
      Model <- list()
      Fmula <- formula(paste("PresAbse", '~ .'))
      
      #SVM model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][, c("PresAbse", VarCol)]
        set.seed(0)
        Model[[i]] <- ksvm(Fmula,data = dataPr,kernel = "rbfdot",C = 1, prob.model=T)
      }
      
      #SVM evaluation
      Eval <- list()
      RastPart <- list()
      for (i in 1:N) {
        RastPart[[i]] <- STANDAR(predict(Variables,Model[[i]]))
        PredPoint <- extract(RastPart[[i]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      }
      
      #SVM threshold
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      #SVM result 
      Validation<-SUMMRES(Eval, N, Thr)
      Validation_SVM[[s]] <- data.frame(Sp=SpNames[s], Algorithm="SVM", Validation)
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(RastPart[[i]],paste(grep("SVM",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("SVM",PartCat,value=T), '/',SpNames[s],"_",i,"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }else{
            writeRaster(RastPart[[i]],paste(grep("SVM",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("SVM",PartCat,value=T), '/',SpNames[s],"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
      
      # Save final model
      Model <- ksvm(Fmula,data = SpData[,c("PresAbse", VarCol)],
                    kernel = "rbfdot",C = 1, prob.model=T)
      ListRaster[["SVM"]] <- STANDAR(predict(Variables,Model))
      PredPoint <- extract(ListRaster[["SVM"]], SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                   PredPoint[PredPoint$PresAbse == 0, 2]))
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- Threshold
      Validation<-SUMMRES(Eval, 1, Thr)
      Summary_SVM[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "SVM", Validation)
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
    
    #RANDOM FOREST (RDF)----
    if (any(Algorithm == "RDF")) {
      Model <- list()
      #RDF model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][, c("PresAbse", VarCol)]
        set.seed(1)
        Model[[i]] <- tuneRF(dataPr[,-1], (dataPr[,1]), trace=F,
                             stepFactor=2, ntreeTry=1000, doBest=T, plot=F)
      }
      #RDF evaluation
      Eval <- list()
      RastPart <- list()
      for (i in 1:N) {
        RastPart[[i]] <- STANDAR(predict(Variables,Model[[i]]))
        PredPoint <- extract(RastPart[[i]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      }
      
      #RDF threshold
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      #RDF result 
      Validation<-SUMMRES(Eval, N, Thr)
      Validation_RDF[[s]] <- data.frame(Sp=SpNames[s], Algorithm="RDF", Validation)
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(RastPart[[i]],paste(grep("RDF",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("RDF",PartCat,value=T), '/',SpNames[s],"_",i,"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }else{
            writeRaster(RastPart[[i]],paste(grep("RDF",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("RDF",PartCat,value=T), '/',SpNames[s],"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
      
      # Save final model
      set.seed(0)
      Model <- tuneRF(SpData[,VarCol], (SpData[,"PresAbse"]), trace=F,
                      stepFactor=2, ntreeTry=500, doBest=T, plot = F)    
      ListRaster[["RDF"]] <- STANDAR(predict(Variables,Model))
      PredPoint <- extract(ListRaster[["RDF"]], SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                   PredPoint[PredPoint$PresAbse == 0, 2]))
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- Threshold
      Validation<-SUMMRES(Eval, 1, Thr)
      Summary_RDF[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "RDF", Validation)
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
    
    #GENERALISED ADDITIVE MODEL (GAM)------
    if(any(Algorithm == 'GAM')) {
      if((sum(sapply(PAtrain, function(x) nrow(x[x[,"PresAbse"]==1,])>length(VarCol)))==N)==T){
      Model <- list()
      Fmula <- paste("s(", VarCol,", k=3)", sep="")
      Fmula <- paste("PresAbse", paste(Fmula, collapse = " + "), sep = " ~ ")
      Fmula <- as.formula(Fmula)
      #GAM model
        
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][, c("PresAbse", VarCol)]
        Model[[i]] <- gam(Fmula, data = dataPr, optimizer = c("outer", "newton"), 
                          select = T, family = binomial)
        }
      #GAM evaluation
      Eval <- list()
      RastPart <- list()
      for (i in 1:N) {
        RastPart[[i]] <- STANDAR(predict(Variables,Model[[i]]))
        PredPoint <- extract(RastPart[[i]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      }
      #GAM threshold
      Thr <- unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      #GAM result
      Validation <- SUMMRES(Eval, N, Thr)
      Validation_GAM[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "GAM", Validation)
      
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(RastPart[[i]],paste(grep("GAM",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("GAM",PartCat,value=T), '/',SpNames[s],"_",i,"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }else{
            writeRaster(RastPart[[i]],paste(grep("GAM",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("GAM",PartCat,value=T), '/',SpNames[s],"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
      
      # Save final model
      Model <- gam(Fmula, data = SpData[, c("PresAbse",VarCol)], optimizer = c("outer", "newton"), 
                   select = T, family = binomial)
      ListRaster[["GAM"]] <- STANDAR(predict(Variables,Model))
      PredPoint <- extract(ListRaster[["GAM"]], SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                   PredPoint[PredPoint$PresAbse == 0, 2]))
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- Threshold
      Validation<-SUMMRES(Eval, 1, Thr)
      Summary_GAM[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "GAM", Validation)
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
      }else{
        Validation2 <- data.frame(matrix(rep(NA,(dim(Validation)[1]*dim(Validation)[2])),nrow=dim(Validation)[1],ncol=dim(Validation)[2]))
        colnames(Validation2) <- colnames(Validation)
        Validation_GAM[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "GAM", Validation2)
        Summary_GAM[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "GAM", Validation2)
        ListRaster[["GAM"]] <- Variables[[1]]*0
        if(is.null(Fut)==F){
          for(k in 1:length(VariablesP)){
            ListFut[[ProjN[k]]][["GAM"]] <- Variables[[1]]*0
          }
        }
      }
    }
    
    #GENERALISED LINEAR MODEL (GLM) ------
    if(any(Algorithm == 'GLM')) {
      if((sum(sapply(PAtrain, function(x) nrow(x[x[,"PresAbse"]==1,])>length(VarCol)))==N)==T){
      Model <- list()
      Fmula <- paste("PresAbse", paste(c(VarCol, paste("(",VarCol, ")^2", sep = "")),
                                     collapse = " + "), sep = " ~ ")
      Fmula <- as.formula(Fmula)
      #GLM model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][, c("PresAbse", VarCol)]
        Model[[i]] <- glm(Fmula, data = dataPr, family = binomial)
      }
      
      #GLM evaluation
      Eval <- list()
      RastPart <- list()
      for (i in 1:N) {
        RastPart[[i]] <- STANDAR(predict(Variables,Model[[i]]))
        PredPoint <- extract(RastPart[[i]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      }
      
      #GLM threshold
      Thr <- unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- rep(Threshold,N)
      #GLM result
      Validation <- SUMMRES(Eval, N, Thr)
      Validation_GLM[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "GLM", Validation)
      
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          if(N!=1){
            writeRaster(RastPart[[i]],paste(grep("GLM",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("GLM",PartCat,value=T), '/',SpNames[s],"_",i,"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }else{
            writeRaster(RastPart[[i]],paste(grep("GLM",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            for(p in 1:length(Thr)){
              writeRaster(RastPart[[i]]>=Thr[p], 
                          paste(grep("GLM",PartCat,value=T), '/',SpNames[s],"_",Validation$TYPE[p],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
      
      # Save final model
      Model <- glm(Fmula, data = SpData[, c("PresAbse",VarCol)], family = binomial)
      ListRaster[["GLM"]] <- STANDAR(predict(Variables,Model))
      PredPoint <- extract(ListRaster[["GLM"]], SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                   PredPoint[PredPoint$PresAbse == 0, 2]))
      Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
      names(Thr) <- Threshold
      Validation<-SUMMRES(Eval, 1, Thr)
      Summary_GLM[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "GLM", Validation)
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
      }else{
        Validation2 <- data.frame(matrix(rep(NA,(dim(Validation)[1]*dim(Validation)[2])),nrow=dim(Validation)[1],ncol=dim(Validation)[2]))
        colnames(Validation2) <- colnames(Validation)
        Validation_GLM[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "GLM", Validation2)
        Summary_GLM[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "GLM", Validation2)
        ListRaster[["GLM"]] <- Variables[[1]]*0
        if(is.null(Fut)==F){
          for(k in 1:length(VariablesP)){
            ListFut[[ProjN[k]]][["GLM"]] <- Variables[[1]]*0
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
      Model[[i]] <- graf(dataPr[,"PresAbse"], dataPr[,VarCol])
    }
      
    #GAU evaluation
    Eval <- list()
    RastPart <- list()
    for (i in 1:N) {
      RastPart[[i]] <- STANDAR(predict(Variables,Model[[i]]))
      PredPoint <- extract(RastPart[[i]], PAtest[[i]][, c("x", "y")])
      PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
      Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1,2],
                            PredPoint[PredPoint$PresAbse == 0,2])
    }
    #GAU threshold
    Thr <- unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
    names(Thr) <- rep(Threshold,N)
    #GAU result
    Validation <- SUMMRES(Eval, N, Thr)
    Validation_GAU[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "GAU", Validation)
    
    #Save Partition Predictions
    if(Save=="Y"){
      for(i in 1:N){
        if(N!=1){
          writeRaster(RastPart[[i]],paste(grep("GAU",foldPart,value=T),"/",SpNames[s],"_",i,".tif", sep=""),
                      format='GTiff',
                      overwrite=TRUE)
          for(p in 1:length(Thr)){
            writeRaster(RastPart[[i]]>=Thr[p], 
                        paste(grep("GAU",PartCat,value=T), '/',SpNames[s],"_",i,"_",Validation$TYPE[p],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }else{
          writeRaster(RastPart[[i]],paste(grep("GAU",foldPart,value=T),"/",SpNames[s],".tif", sep=""),
                      format='GTiff',
                      overwrite=TRUE)
          for(p in 1:length(Thr)){
            writeRaster(RastPart[[i]]>=Thr[p], 
                        paste(grep("GAU",PartCat,value=T), '/',SpNames[s],"_",Validation$TYPE[p],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }
    
    #Save final model
    Model <- graf(SpData[,"PresAbse"], SpData[,VarCol])
    ListRaster[["GAU"]] <- STANDAR(predict.graf.raster(Model, Variables, type = "response", 
                                                       CI = 0.95, maxn = NULL)$posterior.mode)
    PredPoint <- extract(ListRaster[["GAU"]], SpData[, c("x", "y")])
    PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
    Eval <- list(dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                 PredPoint[PredPoint$PresAbse == 0, 2]))
    Thr<-unlist(sapply(Eval, function(x) threshold(x)[Threshold]))
    names(Thr) <- Threshold
    Validation<-SUMMRES(Eval, 1, Thr)
    Summary_GAU[[s]] <- data.frame(Sp = SpNames[s], Algorithm = "GAU", Validation)
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
    
    # Models performance----
    Obj <- ls(pattern = 'Summary_')
    result <- list()
    for(i in 1:length(Obj)){
      result[[i]] <- ldply(get(Obj[i]))}
    result <- ldply(result)

    # Ensemble-----
    # Without Ensemble
      SpValidation <- result[result$Sp==SpNames[s],]
      for(i in 1:length(ListRaster)){
        writeRaster(round(ListRaster[[i]], 4), 
                    paste(folders[i], '/',SpNames2[s],".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
        Thr <- SpValidation[SpValidation$Algorithm==names(ListRaster[i]), 'THR']
        Type <- SpValidation[SpValidation$Algorithm==names(ListRaster[i]), 'TYPE']
        for (h in 1:length(Thr)){
          writeRaster(ListRaster[[i]]>=Thr[h], 
                      paste(foldCat[i], '/',SpNames2[s],"_",Type[h],".tif", sep=""),
                      format='GTiff',
                      overwrite=TRUE)
        }
      }
    
    #Save Projections
    if(is.null(Fut)==F){
      for(p in 1:length(ListFut)){
        for(o in 1:length(ListFut[[p]])){
          Thr <- SpValidation[SpValidation$Algorithm==names(ListFut[[p]])[o], 'THR']
          Type <- SpValidation[SpValidation$Algorithm==names(ListFut[[p]])[o], 'TYPE']
          for(m in 1:length(Thr)){
            writeRaster(ListFut[[p]][[o]]>=Thr[m], 
                        file.path(ModFut[p],Algorithm[o],"BIN",paste(SpNames2[s],Type[m],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
          writeRaster(ListFut[[p]][[o]],file.path(ModFut[p],Algorithm[o],SpNames2[s]),
                      format='GTiff',overwrite=TRUE)
        }
      }
    }
    
    # Mean Ensemble----
    if(any(PredictType=="Mean")){
      Final <- brick(ListRaster)
      Final <- STANDAR(round(mean(Final),4))
      
      # Threshold
      PredPoint <- extract(Final, SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Thr <- threshold(Eval)[Threshold]
      names(Thr) <- Threshold
      Validation <- SUMMRES(list(Eval),N=1,Thr)
      Summary_Mean[[s]] <- data.frame(Sp=SpNames[s], Algorithm="MEA", Validation)

      writeRaster(Final, 
                  paste(DirMean, '/',SpNames2[s],".tif", sep=""),
                  format='GTiff',
                  overwrite=TRUE)
      TYPE <- Summary_Mean[[s]]["TYPE"]
      for(h in 1:length(Thr)){
        writeRaster(Final>=as.numeric(Thr)[h], 
                    paste(DirMeanCat, '/',SpNames2[s],"_",as.character(TYPE[h,]),".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
      }
      if(is.null(Fut)==F){
        for(p in 1:length(ListFut)){
          Final <- brick(ListFut[[p]])
          Final <- STANDAR(round(mean(Final),4))
          
          writeRaster(Final, 
                      file.path(ModFut[p],"ENS","Mean",SpNames2[s]),
                      format='GTiff',
                      overwrite=TRUE)
          TYPE <- Summary_Mean[[s]]["TYPE"]
          for(h in 1:length(Thr)){
            writeRaster(Final>=as.numeric(Thr)[h], 
                        file.path(ModFut[p],"ENS","Mean","BIN",paste(SpNames2[s],as.character(TYPE[h,]),sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }

    # With Over the Mean(Superior) Ensemble----
    if(any(PredictType=='Sup')){
      SpValidation <- result[result$Sp==SpNames[s],]
      Nom <- NULL
      if("no_omission"%in%Threshold){
        Nom <- c(Nom,"LPT")
      }
      if("spec_sens"%in%Threshold){
        Nom <- c(Nom,"MAX")
      }
      Validation <- NULL
      for(h in 1:length(Nom)){
        SpValidationT <- SpValidation[SpValidation$TYPE==Nom[h],]
        SpValidationT$Algorithm <- as.character(SpValidationT$Algorithm)
      
        Best <- which(SpValidationT$TSS>=mean(SpValidationT$TSS))
        Best <- SpValidationT$Algorithm[Best]
        
        W <- names(ListRaster)%in%Best
        Final <- brick(ListRaster[W])
        Final <- STANDAR(round(mean(Final),4))
        
        # Threshold
        PredPoint <- extract(Final, SpData[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                         PredPoint[PredPoint$PresAbse == 0, 2])
        Thr <- threshold(Eval)[Threshold[h]]
        names(Thr) <- Threshold[h]
        Validation <- rbind(Validation,SUMMRES(list(Eval),N=1,Thr))

        writeRaster(Final, 
                    paste(DirSup, '/',paste(SpNames2[s],Nom[h],sep="_"),".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
        writeRaster(Final>=as.numeric(Thr), 
                    paste(DirSupCat, '/',SpNames2[s],"_",Nom[h],".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
        if(is.null(Fut)==F){
          for(p in 1:length(ListFut)){
            Final <- brick(ListFut[[p]][W])
            Final <- STANDAR(round(mean(Final),4))
            
            writeRaster(Final, 
                        file.path(ModFut[p],"ENS","Sup",paste(SpNames2[s],Nom[h],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
            writeRaster(Final>=as.numeric(Thr), 
                        file.path(ModFut[p],"ENS","Sup","BIN",paste(SpNames2[s],Nom[h],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
      Summary_Sup[[s]] <- data.frame(Sp=SpNames[s], Algorithm="SUP", Validation)
    }

    # With PCA ------
    if (any(PredictType == 'PCA')) {

      # Selection of best algorithms based on TSS
      Final <- brick(ListRaster)
      
      #PCA
      Final <- PCA_ENS_TMLA(Final)
      
      # Threshold
      PredPoint <- extract(Final, SpData[, c("x", "y")])
      PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Thr <- threshold(Eval)[Threshold]
      names(Thr) <- Threshold
      Validation <- SUMMRES(list(Eval),N=1,Thr)
      Summary_PCA[[s]] <- data.frame(Sp=SpNames[s], Algorithm="PCA", Validation)
      
      writeRaster(Final, 
                  paste(DirPCA, '/',SpNames2[s],".tif", sep=""),
                  format='GTiff',
                  overwrite=TRUE)
      TYPE <- Summary_PCA[[s]]["TYPE"]
      for(h in 1:length(Thr)){
        writeRaster(Final>=as.numeric(Thr)[h], 
                    paste(DirPCACat, '/',SpNames2[s],"_",as.character(TYPE[h,]),".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
      }
      if(is.null(Fut)==F){
        for(p in 1:length(ListFut)){
          Final <- brick(ListFut[[p]])
          Final <- PCA_ENS_TMLA(Final)
          
          writeRaster(Final, 
                      file.path(ModFut[p],"ENS","PCA",SpNames2[s]),
                      format='GTiff',
                      overwrite=TRUE)
          TYPE <- Summary_PCA[[s]]["TYPE"]
          for(h in 1:length(Thr)){
            writeRaster(Final>=as.numeric(Thr)[h], 
                        file.path(ModFut[p],"ENS","PCA","BIN",paste(SpNames2[s],as.character(TYPE[h,]),sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }
        
    # With PCA over the Mean(Superior) Ensemble----
    if (any(PredictType == 'PCA_Sup')) {

      # Selection of best algorithms based on TSS
      SpValidation <- result[result$Sp==SpNames[s],]
      
      Nom <- NULL
      if("no_omission"%in%Threshold){
        Nom <- c(Nom,"LPT")
      }
      if("spec_sens"%in%Threshold){
        Nom <- c(Nom,"MAX")
      }
      
      Validation <- NULL
      for(h in 1:length(Nom)){
        SpValidationT <- SpValidation[SpValidation$TYPE==Nom[h],]
        SpValidationT$Algorithm <- as.character(SpValidationT$Algorithm)
        
        Best <- which(SpValidationT$TSS>=mean(SpValidationT$TSS))
        Best <- SpValidationT$Algorithm[Best]
        
        W <- names(ListRaster)%in%Best
        Final <- brick(ListRaster[W])
      
        #PCA
        Final <- PCA_ENS_TMLA(Final)
        
        # Threshold
        PredPoint <- extract(Final, SpData[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
        Thr <- threshold(Eval)[Threshold[h]]
        names(Thr) <- Threshold[h]
        Validation <- rbind(Validation,SUMMRES(list(Eval),N=1,Thr))

        writeRaster(Final, 
                    paste(DirPCA_Sup, '/',SpNames2[s],"_",Nom[h],".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
        writeRaster(Final>=as.numeric(Thr), 
                    paste(DirPCA_SupCat, '/',SpNames2[s],"_",Nom[h],".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
        if(is.null(Fut)==F){
          for(p in 1:length(ListFut)){
            Final <- brick(ListFut[[p]][W])
            Final <- PCA_ENS_TMLA(Final)
            
            writeRaster(Final, 
                        file.path(ModFut[p],"ENS","PCA_Sup",paste(SpNames2[s],Nom[h],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
            writeRaster(Final>=as.numeric(Thr), 
                        file.path(ModFut[p],"ENS","PCA_Sup","BIN",paste(SpNames2[s],Nom[h],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
      Summary_PCA_Sup[[s]] <- data.frame(Sp=SpNames[s], Algorithm="PCS", Validation)
    }
    
    #With PCA over the threshold Ensemble----
    if (any(PredictType == 'PCA_Thr')) {
      
      SpValidation <- result[result$Sp==SpNames[s],]
      TYPE <- NULL
      if("no_omission"%in%Threshold){
        TYPE <- c(TYPE,"LPT")
      }
      if("spec_sens"%in%Threshold){
        TYPE <- c(TYPE,"MAX")
      }
      Validation <- NULL
      for(n in 1:length(TYPE)){
       Final <- brick(ListRaster)
       SpValType <-SpValidation[SpValidation$TYPE==TYPE[n],]
       
        #Select only values above the Threshold
        for(k in Algorithm){
          FinalSp <- Final[[k]]
          FinalSp[FinalSp<SpValType[SpValType$Algorithm==k,"THR"]] <- 0
          Final[[k]] <- FinalSp
        }
        
        #PCA
        Final <- PCA_ENS_TMLA(Final)
        
        # Threshold
        PredPoint <- extract(Final, SpData[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
        Thr <- c(threshold(Eval))[Threshold[n]]
        names(Thr) <- Threshold[n]
        Validation <- rbind(Validation,SUMMRES(list(Eval),N=1,Thr))
        
        writeRaster(Final, 
                    paste(DirPCA_Thr, '/',paste(SpNames2[s],TYPE[n],sep="_"),".tif", sep=""),
                    format='GTiff',
                    overwrite=TRUE)
        writeRaster(Final>=unlist(Thr), 
                      paste(DirPCA_ThrCat, '/',SpNames2[s],"_",TYPE[n],".tif", sep=""),
                      format='GTiff',
                      overwrite=TRUE)
      }
      Summary_PCA_Thr[[s]] <- data.frame(Sp=SpNames[s], Algorithm="PCT", Validation)
      
      if(is.null(Fut)==F){
        for(p in 1:length(ListFut)){
          for(n in 1:length(TYPE)){
            Final <- brick(ListFut[[p]])
            SpValType <-SpValidation[SpValidation$TYPE==TYPE[n],]  
            
            #Select only values above the Threshold
            for(k in Algorithm){
              FinalSp <- Final[[k]]
              FinalSp[FinalSp<SpValType[SpValType$Algorithm==k,"THR"]] <- 0
              if(all(na.omit(FinalSp[])==0)){
                Final[Final[[k]]] <- NULL
              }else{
              Final[[k]] <- FinalSp
              }
            }
            
            Final <- PCA_ENS_TMLA(Final)

            writeRaster(Final, 
                        file.path(ModFut[p],"ENS","PCA_Thr",paste(SpNames2[s],TYPE[n],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
            for(h in 1:length(Thr)){
              writeRaster(Final>=unlist(Thr), 
                          file.path(ModFut[p],"ENS","PCA_Thr","BIN",paste(SpNames2[s],TYPE[n],sep="_")),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
    }
    
    # Save .txt with the models performance---- 
    Obj <- ls(pattern = 'Validation_')
    ObjII <- ls(pattern = 'Summary_')
    result <- list()
    resultII <- list()
    for(i in 1:length(Obj)){
      result[[i]] <- ldply(get(Obj[i]))
      resultII[[i]] <- ldply(get(ObjII[i]))}
    
    result <- ldply(result)
    resultII <- ldply(resultII)
    
    write.table(result, paste(DirSave, VALNAME, sep = '/'), sep="\t",
                col.names = T, row.names = F)
    write.table(resultII, paste(DirSave, VALNAMEII, sep = '/'), sep="\t",
                col.names = T, row.names = F)
    
  }#Fecha loop Especie
  
  # Save additional information and retuls----
  InfoModeling <- list(c("###########################################################"),
       paste('Start date :',Ti),
       paste('End date :',Sys.time()),
       c("Algorithm:", Algorithm),
       c("Ensemble:" , PredictType),
       c("Partition Method:" , Part),
       paste("PA Mask:" , DirMask),
       paste("MSDM:" , DirMSDM),
       paste("Resultados em:" , DirSave),
       paste('No_species:',length(SpNames)),
       matrix(SpNames))
  lapply(InfoModeling, write, 
         paste(DirSave, "/InfoModeling.txt", sep=""), append=TRUE, 
         ncolumns=20, sep='\t')
}
