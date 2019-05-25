## Written by Santiago Velazco & Andre Andrade

FitENM_TMLA_Parallel <- function(RecordsData,
                   Variables,
                   Fut=NULL,
                   Part,
                   Algorithm,
                   PredictType,
                   VarImP,
                   spN,
                   Tst,
                   Threshold = Thresh,
                   DirSave=DirR, 
                   DirMask=NULL,
                   DirMSDM=DirPRI,
                   DirProj=NULL,
                   per=NULL,
                   repl=NULL,
                   Save="N",
                   SaveFinal=SaveFinal,
                   sensV) {
  
  Ti <- Sys.time()
  options(warn = -1)
  
  #Start Cluster
  cl <- makeCluster(detectCores()-1,outfile="")
  registerDoParallel(cl)

  # Directory to save----
  folders <- paste(DirSave,"Algorithm",Algorithm,sep="/")
  for(i in 1:length(folders)){
    dir.create(folders[i],recursive=T)
  }
  
  #Binary directories
  foldCat <- file.path(sort(rep(folders,length(Threshold))),Threshold)
  for(i in 1:length(foldCat)){
    dir.create(foldCat[i])
  }
  
  #Partition directories
  if(Save=="Y"){
    foldPart <- paste(folders,"PartialModels",sep="/")
    for(i in 1:length(foldPart)){
      dir.create(foldPart[i])
      
    }
    #Binary directories
    PartCat <- file.path(sort(rep(foldPart,length(Threshold))),Threshold)
    for(i in 1:length(PartCat)){
      dir.create(PartCat[i])
    }
  }
  
  #Ensemble directory
  if(any(PredictType!="N")){
    DirENS <- paste(DirSave,"Ensemble",sep="/")
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
  }
  
  #Projection directories
  if(is.null(Fut)==F){
    ProjN <- names(Fut)
    dir.create(file.path(DirSave,"Projection"))
    ModFut <- file.path(DirSave,"Projection",ProjN)
    for(i in 1:length(ModFut)){
      dir.create(ModFut[i])
      for(h in Algorithm){
        dir.create(file.path(ModFut[i],h))
        FutCat <- file.path(sort(rep(file.path(ModFut[i],h),length(Threshold))),Threshold)
        sapply(FutCat,function(x) dir.create(x))
        if(any(PredictType!="N")){
          DirENS <- file.path(ModFut[i],"Ensemble",PredictType)
          sapply(DirENS,function(x) dir.create(x,recursive = T))
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
  print(paste("Total species to be modeled", length(spN)))

  # Number of partition
  N <- as.numeric(max(RecordsData[, "Partition"]))
      
  #Txt of Final tables    
  if(is.null(repl)==F){
    VALNAME <- paste('Validation_Partition','_',sep="",repl,'.txt' )
  }else{
    VALNAME <- paste('Validation_Partition.txt' )
  }
  VALNAMEII <- paste('Thresholds_Complete.txt' )
  
  
  # Backqround points----
  if (!is.null(DirMask)) {
     if ((any(Algorithm == "MXD") || any(Algorithm == "MXS") || any(Algorithm == "MLK"))) {
       RecordsDataM <- split(RecordsData,f=RecordsData$sp)
       RecordsDataMt <- list()
       
       RecordsDataMt <- foreach (i=1:length(RecordsDataM))%dopar%{
         RecordsDataMt <- RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1,]
         return(RecordsDataMt)
       }
       
       RecordsDataM <- RecordsDataMt
       names(RecordsDataM) <- spN
       rm(RecordsDataMt)
       
       ab.0 <- foreach (i= 1:length(names(RecordsDataM)),.packages=c("dismo","raster","plyr"))%dopar%{
        msk <- raster(paste(DirMask,paste(spN[i],".tif",sep=""),sep="/"))
         if(Part=="BOOT"||Part=="KFOLD"){
           NM <- 1
           RecordsDataM[[i]] <- rbind(RecordsDataM[[i]],RecordsData[RecordsData$sp==i & RecordsData$Partition==2 & RecordsData$PresAbse==0 ,])
         }else{
           NM <- max(RecordsDataM[[i]]$Partition)
         }
         ab0L <- list()
         for(x in 1:NM){
           msk2 <- msk
           msk2[!msk[]==x] <- NA 
           NCell <- sum(!is.na(msk2[]))
           if (NCell > 10000) {
             ab.0 <- data.frame(randomPoints(msk2,p=RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1, c("x", "y")],10000))
             var.0 <- data.frame(extract(Variables,ab.0))
           }else{
             ab.0 <-
               randomPoints(msk2, p=RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1, c("x", "y")],abs(NCell - nrow(RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1,])))
             var.0 <- data.frame(extract(Variables, ab.0))
           }
           ab.0 <- cbind(rep(spN[i],nrow(ab.0)),ab.0,rep(x,nrow(ab.0)),rep(0,nrow(ab.0)),var.0)
           colnames(ab.0) <- colnames(RecordsData)
           ab0L[[x]] <- na.omit(ab.0)
           rm(var.0)
           # RecordsDataM <- rbind(RecordsDataM[[i]],ab.0)
           # rm(ab.0)
         }
        ab.0 <- ldply(ab0L,data.frame,.id=NULL)
        return(ab.0)
       }
       
       RecordsDataM <- lapply(seq_along(RecordsDataM), function(x) rbind(RecordsDataM[[x]], ab.0[[x]]))
       RecordsDataM <- ldply(RecordsDataM,data.frame,.id=NULL)
       cols <-  c("x","y","Partition","PresAbse",names(Variables))
       RecordsDataM[,cols] <-  apply(RecordsDataM[,cols], 2, function(x) as.numeric(as.character(x)))
     }
   }else{
    if ((any(Algorithm == "MXD") || any(Algorithm == "MXS") || any(Algorithm == "MLK"))) {
      RecordsDataM <- split(RecordsData,f=RecordsData$sp)
      RecordsDataMt <- list()

      RecordsDataMt <- foreach (i=1:length(RecordsDataM))%dopar%{
        RecordsDataMt <- RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1,]
        return(RecordsDataMt)
      }

      RecordsDataM <- RecordsDataMt
      names(RecordsDataM) <- spN
      rm(RecordsDataMt)

      ab.0 <- foreach (i=1:length(names(RecordsDataM)),.packages=c("dismo","raster","plyr"))%dopar%{
        msk <- Variables[[1]]
        msk[is.na(msk[])==FALSE] <- 1

        if(Part=="BOOT"||Part=="KFOLD"){
          NM <- 1
          RecordsDataM[[i]] <- rbind(RecordsDataM[[i]],RecordsData[RecordsData$sp==i & RecordsData$Partition==2 & RecordsData$PresAbse==0 ,])
        }else{
          NM <- max(RecordsDataM[[i]][,"Partition"])
        }
        
        ab0L <- list()
        for(x in 1:NM){
          msk2 <- msk
          msk2[!msk[]==x] <- NA 
          NCell <- sum(!is.na(msk2[]))
          if (NCell > 10000) {
            ab.0 <- data.frame(randomPoints(msk2,p=RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1, c("x", "y")],10000))
            var.0 <- extract(Variables,ab.0)
          }else{
            ab.0 <-
              randomPoints(msk2, p=RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1, c("x", "y")],(NCell - nrow(RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1,])))
            var.0 <- extract(Variables, ab.0)
          }
          ab.0 <- cbind(rep(i,nrow(ab.0)),ab.0,rep(x,nrow(ab.0)),rep(0,nrow(ab.0)),var.0)
          colnames(ab.0) <- colnames(RecordsData)
          ab0L[[x]] <- na.omit(ab.0)
          rm(var.0)
          # RecordsDataM <- rbind(RecordsDataM,ab.0)
          # rm(ab.0)
        }
        ab.0 <- ldply(ab0L,data.frame,.id=NULL)
        return(ab.0)
      }
      
      RecordsDataM <- lapply(seq_along(RecordsDataM), function(x) rbind(RecordsDataM[[x]], ab.0[[x]]))
      RecordsDataM <- ldply(RecordsDataM,data.frame,.id=NULL)
      cols <-  c("x","y","Partition","PresAbse",names(Variables))  
      RecordsDataM[,cols] = apply(RecordsDataM[,cols], 2, function(x) as.numeric(as.character(x)))
    }
   }
  if(!exists("RecordsDataM")){
    RecordsDataM <- NULL
  }

  #MESS & MOPA Calculation----
  #Within the extent (for M-Restriction)
  cat("Calculate extrapolation for the current extent?(Y/N)")
  ansE <- readLines(n=1)
  while(!ansE%in%c("Y","N")){
    cat("Calculate extrapolation for the current extent?(Y/N)")
    ansE <- readLines(n=1)
  }
  if(ansE=="Y"){
    dir.create(file.path(DirSave,"Extrapolation"))
    DirProj <- file.path(DirSave,"Extrapolation")
    MESS_and_MOP(Variables=list(Variables),RecordsData=RecordsData,RecordsDataM=RecordsDataM,algorithm=Algorithm,
                  VarCol=VarCol,DirProj=DirProj,Methods=c("MESS","MOP"))
  }
  #For projections
  if(!is.null(Fut)){
    cat("Calculate extrapolation for the projected extent?(Y/N)")
    ansE <- readLines(n=1)
    while(!ansE%in%c("Y","N")){
      cat("Calculate extrapolation for the current extent?(Y/N)")
      ansE <- readLines(n=1)
    }
    if(ansE=="Y"){
      for(i in 1:length(ModFut)){
        dir.create(file.path(ModFut[i],"Extrapolation"))
      }
      DirProj <- file.path(ModFut,"Extrapolation")
      MESS_and_MOP(Variables=Fut,RecordsData=RecordsData,RecordsDataM=RecordsDataM,algorithm=Algorithm,
                   VarCol=VarCol,DirProj=DirProj,Methods=c("MESS","MOP"))
    }
  }
  
  #Define N due to Partition Method
  if(Part=="BOOT" || Part=="KFOLD"){
    N <- 1
  }else{
    N <- N
  }
  
  # Construction of models LOOP-----
  results <- foreach(s = 1:length(spN), .packages = c("raster","dismo","kernlab","randomForest",
                                                      "maxnet","maxlike","GRaF","plyr","gam","RStoolbox",
                                                      "adehabitatHS","caret","visreg","glmnet","gbm"),
                     .export=c("Validation2_0","STANDAR","maxnet2","predict.graf.raster","PCA_ENS_TMLA","predict.maxnet","boycei",
                               "Eval_Jac_Sor_TMLA","Validation_Table_TMLA","Thresholds_TMLA","VarImp_RspCurv","hingeval","ecospat.boyce")) %dopar% {

    #Results Lists
    ListRaster <- as.list(Algorithm)
    names(ListRaster) <- Algorithm
    
    #Validation List
    ListValidation <- as.list(Algorithm)
    names(ListValidation) <- Algorithm
    
    #Summary List
    ListSummary <- as.list(Algorithm)
    names(ListSummary) <- Algorithm
    
    #Create lists for the future
    if(is.null(Fut)==F){
       ListFut <- as.list(ProjN)
       names(ListFut) <- ProjN
       
       for(p in 1:length(ListFut)){
         ListFut[[p]] <- ListRaster
       }
    }

    # Ocurrences filtter by species and splited by partition----
    SpData <- RecordsData[RecordsData[, "sp"] == spN[s], ]
    if ((any(Algorithm == "MXD") | any(Algorithm == "MXS"))) {
      SpDataM <- RecordsDataM[RecordsDataM[, "sp"] == spN[s], ]
    }
    
    #Include MSDM----
    if(is.null(DirMSDM)==F){
      if(grepl("LatLong",DirMSDM)){
        MSDM <- stack(file.path(DirMSDM,list.files(DirMSDM,pattern=".tif")))
        names(MSDM) <- c("Lat","Long")
      }else{
        MSDM <- raster(file.path(DirMSDM,paste(spN[s],".tif",sep="")))
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
      if ((any(Algorithm == "MXD") | any(Algorithm == "MXS"))) {
        SpDataTM <- SpDataM
      }
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
        PAtestM[[i]] <- SpDataT[SpDataT[, "Partition"] != i, ]
      }
      if(Part%in%c("BOOT","KFOLD")){
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
        Eval_JS <- list()
        Boyce <- list()
        for (i in 1:N) {
          RastPart[["BIO"]][[i]] <- predict(Model[[i]], PAtest[[i]][, VarColT])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["BIO"]][[i]])
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                           a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["BIO"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }

        #BIO Validation 
        Boyce <- mean(unlist(Boyce))
        Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N)
        if(is.null(repl)){
          ListValidation[["BIO"]] <- data.frame(Sp=spN[s], Algorithm="BIO",Partition=Part, Validation,Boyce=Boyce)
        }else{
          ListValidation[["BIO"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="BIO",Partition=Part, Validation,Boyce=Boyce)
        }
      
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          #Partial Thresholds
          Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
          PartRas <- STANDAR(predict(Model[[i]], VariablesT))
          if(N!=1){
            writeRaster(PartRas,paste(grep("BIO",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern="BIO",x=foldCat,value=T)
            for(t in 1:length(Thr_Alg)){
              writeRaster(PartRas>=Thr_Alg[t], 
                          paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
          if(is.null(repl)==F){
            writeRaster(PartRas,paste(grep("BIO",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern="BIO",x=foldCat,value=T)
            for(t in 1:length(Thr_Alg)){
              writeRaster(PartRas>=Thr_Alg[t], 
                          paste(foldCatAlg[t], '/',spN[s],"_",repl,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }else{
            writeRaster(PartRas,paste(grep("BIO",foldPart,value=T),"/",paste0(spN[s],repl),".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern="BIO",x=foldCat,value=T)
            for(t in 1:length(Thr_Alg)){
              writeRaster(PartRas>=Thr_Alg[t],
                          paste(grep("BIO",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
      
      # Save final model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Model <- bioclim(SpDataT[SpDataT[,"PresAbse"]==1, VarColT]) # only presences
        #Final Model "Evaluation"
        PredPoint <- predict(Model, SpDataT[, VarColT])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                         PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        #Final Thresholds
        Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
        
        #Variable Importance & Response Curves
        if(VarImP=="Y"){
          VarImp_RspCurv(Model=Model,Algorithm='BIO',spN=spN[s],SpDataT = SpDataT,
                         VarColT=VarColT,Outcome=PredPoint$PredPoint)
        }

        #Final Model Rasters
        ListSummary[["BIO"]] <- data.frame(Sp=spN[s], Algorithm="BIO", Thr)
        if(SaveFinal=="Y"){
          ListRaster[["BIO"]] <- STANDAR(predict(Model, VariablesT))
          names(ListRaster[["BIO"]]) <- spN[s]
        }
        #Future Projections
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
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["BIO"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
        
          #BIO Validation 
          Boyce <- mean(unlist(Boyce))
          Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
          if(is.null(repl)){
            ListValidation[["BIO"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="BIO", Validation,Boyce=Boyce)
          }else{
            ListValidation[["BIO"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="BIO", Validation,Boyce=Boyce)
          }
        }
      }
    }
    
    #DOMAIN (DOM)----- 
    if (any(Algorithm == "DOM")) {
      Model <- list()
      #DOM model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][PAtrain[[i]][, "PresAbse"] == 1,]
        Model[[i]] <- dismo::domain(dataPr[, VarColT])
      }
      
      #DOM evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        Eval_JS <- list()
        for (i in 1:N) {
          RastPart[["DOM"]][[i]] <- predict(Model[[i]], PAtest[[i]][, VarColT])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["DOM"]][[i]])
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["DOM"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
          
        #DOM Validation 
        Boyce <- mean(unlist(Boyce))
        Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
        if(is.null(repl)){
          ListValidation[["DOM"]] <- data.frame(Sp=spN[s], Algorithm="DOM",Partition=Part, Validation,Boyce=Boyce)
        }else{
          ListValidation[["DOM"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="DOM",Partition=Part, Validation,Boyce=Boyce)
        }
        
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            #Partial Thresholds
            Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
            PartRas <- STANDAR(predict(Model[[i]], VariablesT))
            if(N!=1){
              writeRaster(PartRas,paste(grep("DOM",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="DOM",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              writeRaster(PartRas,paste(grep("DOM",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="DOM",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",repl,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(PartRas,paste(grep("DOM",foldPart,value=T),"/",paste0(spN[s],repl),".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="DOM",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t],
                            paste(grep("DOM",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
          }
        }
        
        # Save final model
        if(per!=1 && repl==1 || per==1 || N!=1){
          Model <- dismo::domain(SpDataT[SpDataT[,"PresAbse"]==1, VarColT]) # only presences
          PredPoint <- predict(Model, SpDataT[, VarColT])
          PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                       a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
          
          #Variable Importance & Response Curves
          if(VarImP=="Y"){
            VarImp_RspCurv(Model=Model,Algorithm='DOM',spN=spN[s],SpDataT = SpDataT,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }
          
          #Final Model Rasters
          ListSummary[["DOM"]] <- data.frame(Sp=spN[s], Algorithm="DOM", Thr)
          if(SaveFinal=="Y"){
            ListRaster[["DOM"]] <- STANDAR(predict(Model, VariablesT))
            names(ListRaster[["DOM"]]) <- spN[s]
          }
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["DOM"]] <- predict(VariablesP[[k]], Model)
              if(maxValue(ListFut[[ProjN[k]]][["DOM"]])==0){
                ListFut[[ProjN[k]]][["DOM"]] <- ListFut[[ProjN[k]]][["DOM"]]
              }else{
                ListFut[[ProjN[k]]][["DOM"]] <- STANDAR(ListFut[[ProjN[k]]][["DOM"]])
              }
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["DOM"]] <- predict(VariablesP[[k]],Model[[i]])
          if(maxValue(ListFut[[ProjN[k]]][["DOM"]])==0){
            ListFut[[ProjN[k]]][["DOM"]] <- ListFut[[ProjN[k]]][["DOM"]]
          }else{
            ListFut[[ProjN[k]]][["DOM"]] <- STANDAR(ListFut[[ProjN[k]]][["DOM"]])
          }
          
          PredPoint <- extract(ListFut[[ProjN[k]]][["DOM"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["DOM"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #DOM Validation 
          Boyce <- mean(unlist(Boyce))
          Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
          if(is.null(repl)){
            ListValidation[["DOM"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="DOM", Validation,Boyce=Boyce)
          }else{
            ListValidation[["DOM"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="DOM", Validation,Boyce=Boyce)
          }
        }
      }
    }
    
    #MAHALANOBIS (MAH)----- 
    if (any(Algorithm == "MAH")) {
      Model <- list()
      #MAH model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][PAtrain[[i]][, "PresAbse"] == 1,]
        Model[[i]] <- mahal(dataPr[, VarColT])
      }
      
      #MAH evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        Eval_JS <- list()
        for (i in 1:N) {
          Ras <- STANDAR(predict(Model[[i]],VariablesT))
          RastPart[["MAH"]][[i]] <- extract(Ras,PAtest[[i]][,c("x","y")])
          rm(Ras)
          # RastPart[["MAH"]][[i]] <- predict(Model[[i]], PAtest[[i]][, VarColT])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["MAH"]][[i]])
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["MAH"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        
        #MAH Validation 
        Boyce <- mean(unlist(Boyce))
        Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
        if(is.null(repl)){
          ListValidation[["MAH"]] <- data.frame(Sp=spN[s], Algorithm="MAH",Partition=Part, Validation,Boyce=Boyce)
        }else{
          ListValidation[["MAH"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="MAH",Partition=Part, Validation,Boyce=Boyce)
        }
        
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            #Partial Thresholds
            Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
            PartRas <- STANDAR(predict(Model[[i]], VariablesT))
            if(N!=1){
              writeRaster(PartRas,paste(grep("MAH",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MAH",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              writeRaster(PartRas,paste(grep("MAH",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MAH",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",repl,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(PartRas,paste(grep("MAH",foldPart,value=T),"/",paste0(spN[s],repl),".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MAH",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t],
                            paste(grep("MAH",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
          }
        }
        
        # Save final model
        if(per!=1 && repl==1 || per==1 || N!=1){
          Model <- mahal(SpDataT[SpDataT[,"PresAbse"]==1, VarColT]) # only presences
          Ras <- STANDAR(predict(Model,VariablesT))
          PredPoint <- extract(Ras,SpDataT[,c("x","y")])
          PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                       a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
          
          #Variable Importance & Response Curves
          if(VarImP=="Y"){
            VarImp_RspCurv(Model=Model,Algorithm='MAH',spN=spN[s],SpDataT = SpDataT,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }
          
          #Final Model Rasters
          ListSummary[["MAH"]] <- data.frame(Sp=spN[s], Algorithm="MAH", Thr)
          
          if(SaveFinal=="Y"){
            ListRaster[["MAH"]] <- Ras
            names(ListRaster[["MAH"]]) <- spN[s]
          }
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["MAH"]] <- predict(VariablesP[[k]], Model)
              if(maxValue(ListFut[[ProjN[k]]][["MAH"]])==0){
                ListFut[[ProjN[k]]][["MAH"]] <- ListFut[[ProjN[k]]][["MAH"]]
              }else{
                ListFut[[ProjN[k]]][["MAH"]] <- STANDAR(ListFut[[ProjN[k]]][["MAH"]])
              }
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["MAH"]] <- predict(VariablesP[[k]],Model[[i]])
          if(maxValue(ListFut[[ProjN[k]]][["MAH"]])==0){
            ListFut[[ProjN[k]]][["MAH"]] <- ListFut[[ProjN[k]]][["MAH"]]
          }else{
            ListFut[[ProjN[k]]][["MAH"]] <- STANDAR(ListFut[[ProjN[k]]][["MAH"]])
          }
          
          PredPoint <- extract(ListFut[[ProjN[k]]][["MAH"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["MAH"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #MAH Validation 
          Boyce <- mean(unlist(Boyce))
          Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
          if(is.null(repl)){
            ListValidation[["MAH"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="MAH", Validation,Boyce=Boyce)
          }else{
            ListValidation[["MAH"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="MAH", Validation,Boyce=Boyce)
          }
        }
      }
    }
    
    #ENFA (ENF)----- 
    if (any(Algorithm == "ENF")) {
      Model <- list()
      #ENF model
      for (i in 1:N) {
        dataPr <- PAtrainM[[i]][, c("PresAbse", VarColT)]
        dudi <- dudi.pca(dataPr[, VarColT],scannf = FALSE)
        Model[[i]] <- adehabitatHS::madifa(dudi,dataPr$PresAbse,scannf = FALSE)
      }
      
      #ENF evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        Eval_JS <- list()
        for (i in 1:N) {
          PredRas <- values(VariablesT)
          POS <- which(!is.na(PredRas[,1]))
          Zli <- as.matrix(na.omit(values(VariablesT)) %*% as.matrix(Model[[i]]$co))
          POSPRE <- cellFromXY(VariablesT[[1]],PAtrainM[[i]][PAtrainM[[i]]$PresAbse==1,c("x","y")])
          ZER <- rep(0,length(POS))
          ZER[POSPRE] <- 1
          f1 <- function(x) rep(x, ZER)
          Sli <- apply(Zli, 2, f1)
          m <- apply(Sli, 2, mean)
          cov <- t(as.matrix(Zli)) %*% as.matrix(Zli)/nrow(Zli)
          PredRas <- (data.frame(MD = mahalanobis(Zli, center = m, cov = cov,inverted = F)))*-1
          XY <- xyFromCell(VariablesT[[1]],1:ncell(VariablesT[[1]]))
          PredRAS <- data.frame(cbind(XY,ENF=NA))
          PredRAS[POS,"ENF"] <- PredRas
          gridded(PredRAS) <- ~x+y
          PredRAS <- STANDAR(raster(PredRAS))
          RastPart[["ENF"]][[i]] <- extract(PredRAS,PAtestM[[i]][c("x","y")])
          PredPoint <- data.frame(PresAbse = PAtestM[[i]][, "PresAbse"], RastPart[["ENF"]][[i]])
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["ENF"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        
        #ENF Validation 
        Boyce <- mean(unlist(Boyce))
        Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
        if(is.null(repl)){
          ListValidation[["ENF"]] <- data.frame(Sp=spN[s], Algorithm="ENF",Partition=Part, Validation,Boyce=Boyce)
        }else{
          ListValidation[["ENF"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="ENF",Partition=Part, Validation,Boyce=Boyce)
        }
        
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            #Partial Thresholds
            Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
            PredRas <- values(VariablesT)
            POS <- which(!is.na(PredRas[,1]))
            Zli <- as.matrix(scale(na.omit(values(VariablesT))) %*% as.matrix(Model[[i]]$co))
            POSPRE <- cellFromXY(VariablesT[[1]],PAtrainM[[i]][PAtrainM[[i]]$PresAbse==1,c("x","y")])
            ZER <- rep(0,length(POS))
            ZER[POSPRE] <- 1
            f1 <- function(x) rep(x, ZER)
            Sli <- apply(Zli, 2, f1)
            m <- apply(Sli, 2, mean)
            cov <- t(as.matrix(Zli)) %*% as.matrix(Zli)/nrow(Zli)
            PredRas <- data.frame(MD = mahalanobis(Zli, center = m, cov = cov))
            XY <- xyFromCell(VariablesT[[1]],1:ncell(VariablesT[[1]]))
            PredRAS <- data.frame(cbind(XY,ENF=NA))
            PredRAS[POS,"ENF"] <- PredRas
            gridded(PredRAS) <- ~x+y
            PartRas <- STANDAR(raster(PredRAS))
            rm(list=c("PredRas","POS",'Zli',"POSPRE","ZER","f1","Sli","m","cov","PredRAS"))
            if(N!=1){
              writeRaster(PartRas,paste(grep("ENF",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="ENF",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              writeRaster(PartRas,paste(grep("ENF",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="ENF",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",repl,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(PartRas,paste(grep("ENF",foldPart,value=T),"/",paste0(spN[s],repl),".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="ENF",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t],
                            paste(grep("ENF",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
          }
        }
        
        # Save final model
        if(per!=1 && repl==1 || per==1 || N!=1){
          dudi <- dudi.pca(SpDataT[, VarColT],scannf = FALSE)
          Model <- adehabitatHS::madifa(dudi,SpDataT$PresAbse,scannf = FALSE)
          PredRas <- values(VariablesT)
          POS <- which(!is.na(PredRas[,1]))
          Zli <- as.matrix(na.omit(values(VariablesT)) %*% as.matrix(Model$co))
          # Zli <- sweep((sweep(Zli,2,apply(Zli,2,min),"-")),2,(apply(Zli,2,max)-apply(Zli,2,min)),"/")
          POSPRE <- cellFromXY(VariablesT[[1]],PAtrainM[[i]][PAtrainM[[i]]$PresAbse==1,c("x","y")])
          ZER <- rep(0,length(POS))
          ZER[POSPRE] <- 1
          f1 <- function(x) rep(x, ZER)
          Sli <- apply(Zli, 2, f1)
          m <- apply(Sli, 2, mean)
          cov <- t(as.matrix(Zli)) %*% as.matrix(Zli)/nrow(Zli)
          PredRas <- (data.frame(MD = mahalanobis(Zli, center = m, cov = cov,inverted = F)))*-1
          XY <- xyFromCell(VariablesT[[1]],1:ncell(VariablesT[[1]]))
          PredRAS <- data.frame(cbind(XY,ENF=NA))
          PredRAS[POS,"ENF"] <- PredRas
          gridded(PredRAS) <- ~x+y
          ListRaster[["ENF"]] <- STANDAR(raster(PredRAS))
          names(ListRaster[["ENF"]]) <- spN[s]
          PredPoint <- extract(ListRaster[["ENF"]],SpDataT[,c("x","y")])
          
          # Zli <- (as.matrix(SpDataT[, VarColT]) %*% as.matrix(Model$co))*-1
          # # Zli <- sweep((sweep(Zli,2,apply(Zli,2,min),"-")),2,(apply(Zli,2,max)-apply(Zli,2,min)),"/")
          # f1 <- function(x) rep(x, SpDataT$PresAbse)
          # Sli <- apply(Zli, 2, f1)
          # m <- apply(Sli, 2, mean)
          # cov <- t(as.matrix(Sli)) %*% as.matrix(Sli)/nrow(Sli)
          # PredPoint <- data.frame(MD = mahalanobis(Zli, center = m, cov = cov,inverted = T))
          PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Variable Importance & Response Curves
          if(VarImP=="Y"){
            VarImp_RspCurv(Model=Model,Algorithm='ENF',spN=spN[s],SpDataT = SpDataT,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
          
          #Final Model Rasters
          ListSummary[["ENF"]] <- data.frame(Sp=spN[s], Algorithm="ENF", Thr)
          rm(list=c("PredRas","POS",'Zli',"POSPRE","ZER","f1","Sli","m","cov","PredRAS"))
          # if(SaveFinal=="Y"){
            # PredRas <- values(VariablesT)
            # POS <- which(!is.na(PredRas[,1]))
            # Zli <- as.matrix(na.omit(values(VariablesT)) %*% as.matrix(Model$co))
            # # Zli <- sweep((sweep(Zli,2,apply(Zli,2,min),"-")),2,(apply(Zli,2,max)-apply(Zli,2,min)),"/")
            # POSPRE <- cellFromXY(VariablesT[[1]],PAtrainM[[i]][PAtrainM[[i]]$PresAbse==1,c("x","y")])
            # ZER <- rep(0,length(POS))
            # ZER[POSPRE] <- 1
            # f1 <- function(x) rep(x, ZER)
            # Sli <- apply(Zli, 2, f1)
            # m <- apply(Sli, 2, mean)
            # cov <- t(as.matrix(Zli)) %*% as.matrix(Zli)/nrow(Zli)
            # PredRas <- data.frame(MD = mahalanobis(Zli, center = m, cov = cov,inverted = T))
            # XY <- xyFromCell(VariablesT[[1]],1:ncell(VariablesT[[1]]))
            # PredRAS <- data.frame(cbind(XY,ENF=NA))
            # PredRAS[POS,"ENF"] <- PredRas
            # gridded(PredRAS) <- ~x+y
            # ListRaster[["ENF"]] <- STANDAR(raster(PredRAS))
            # # ListRaster[["ENF"]] <- (((STANDAR(raster(PredRAS)) - max(STANDAR(raster(PredRAS)))) * -1) + min(STANDAR(raster(PredRAS))))
            # rm(list=c("PredRas","POS",'Zli',"POSPRE","ZER","f1","Sli","m","cov","PredRAS"))
            # names(ListRaster[["ENF"]]) <- spN[s]
          # }
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              PredRas <- values(VariablesP[[k]])
              POS <- which(!is.na(PredRas[,1]))
              Zli <- as.matrix(na.omit(values(VariablesP[[k]])) %*% as.matrix(Model$co))
              Zli <- sweep((sweep(Zli,2,apply(Zli,2,min),"-")),2,(apply(Zli,2,max)-apply(Zli,2,min)),"/")
              POSPRE <- cellFromXY(VariablesP[[k]],PAtrainM[[i]][PAtrainM[[i]]$PresAbse==1,c("x","y")])
              ZER <- rep(0,length(POS))
              ZER[POSPRE] <- 1
              f1 <- function(x) rep(x, ZER)
              Sli <- apply(Zli, 2, f1)
              m <- apply(Sli, 2, mean)
              cov <- t(as.matrix(Zli)) %*% as.matrix(Zli)/nrow(Zli)
              PredRas <- data.frame(MD = mahalanobis(Zli, center = m, cov = cov))
              XY <- xyFromCell(VariablesT[[1]],1:ncell(VariablesT[[1]]))
              PredRAS <- data.frame(cbind(XY,ENF=NA))
              PredRAS[POS,"ENF"] <- PredRas
              gridded(PredRAS) <- ~x+y
              ListFut[[ProjN[k]]][["ENF"]] <- STANDAR(raster(PredRAS))
              rm(list=c("PredRas","POS",'Zli',"POSPRE","ZER","f1","Sli","m","cov","PredRas",
                        "PredRAS"))
              if(maxValue(ListFut[[ProjN[k]]][["ENF"]])==0){
                ListFut[[ProjN[k]]][["ENF"]] <- ListFut[[ProjN[k]]][["ENF"]]
              }else{
                ListFut[[ProjN[k]]][["ENF"]] <- STANDAR(ListFut[[ProjN[k]]][["ENF"]])
              }
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["ENF"]] <- predict(VariablesP[[k]],Model[[i]])
          if(maxValue(ListFut[[ProjN[k]]][["ENF"]])==0){
            ListFut[[ProjN[k]]][["ENF"]] <- ListFut[[ProjN[k]]][["ENF"]]
          }else{
            ListFut[[ProjN[k]]][["ENF"]] <- STANDAR(ListFut[[ProjN[k]]][["ENF"]])
          }
          
          PredPoint <- extract(ListFut[[ProjN[k]]][["ENF"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["ENF"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #ENF Validation 
          Boyce <- mean(unlist(Boyce))
          Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
          if(is.null(repl)){
            ListValidation[["ENF"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="ENF", Validation,Boyce=Boyce)
          }else{
            ListValidation[["ENF"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="ENF", Validation,Boyce=Boyce)
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
        Model[[i]] <- maxnet2(p=dataPr[,"PresAbse"], data=dataPr[,VarColT], f = 
                               maxnet.formula(dataPr[,"PresAbse"], 
                                              dataPr[,VarColT], classes="default"))
      }
      #MXD evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        Eval_JS <- list()
        for (i in 1:N) {
          RastPart[["MXD"]][[i]] <- c(predict(Model[[i]], PAtestM[[i]][, VarColT], clamp=F, type="cloglog"))
          PredPoint <- data.frame(PresAbse = PAtestM[[i]][, "PresAbse"], RastPart[["MXD"]][[i]])
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["MXD"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        
        #MXD Validation 
        Boyce <- mean(unlist(Boyce))
        Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
        if(is.null(repl)){
          ListValidation[["MXD"]] <- data.frame(Sp=spN[s], Algorithm="MXD",Partition=Part, Validation,Boyce=Boyce)
        }else{
          ListValidation[["MXD"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="MXD",Partition=Part, Validation,Boyce=Boyce)
        }
        
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            #Partial Thresholds
            Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
            PartRas <- STANDAR(predict(VariablesT,Model[[i]], clamp=F, type="cloglog"))
            if(N!=1){
              writeRaster(PartRas,paste(grep("MXD",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MXD",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              writeRaster(PartRas,paste(grep("MXD",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MXD",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(PartRas,paste(grep("MXD",foldPart,value=T),"/",spN[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MXD",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t],
                            paste(grep("MXD",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
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
          PredPoint <- c(predict(Model, SpDataTM[, VarColT], clamp=F, type="cloglog"))
          PredPoint <- data.frame(PresAbse = SpDataTM[, "PresAbse"], PredPoint)
          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                           PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
          
          #Variable Importance & Response Curves
          if(VarImP=="Y"){
            VarImp_RspCurv(Model=Model,Algorithm='MXD',spN=spN[s],SpDataT = SpDataTM,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }
          
          #Final Model Rasters
          ListSummary[["MXD"]] <- data.frame(Sp=spN[s], Algorithm="MXD", Thr)
          
          if(SaveFinal=="Y"){
            ListRaster[["MXD"]] <- STANDAR(predict(VariablesT,Model, clamp=F, type="cloglog"))
            names(ListRaster[["MXD"]]) <- spN[s]
          }
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
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["MXD"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #MXD Validation 
          Boyce <- mean(unlist(Boyce))
          Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
          if(is.null(repl)){
            ListValidation[["MXD"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="MXD", Validation,Boyce=Boyce)
          }else{
            ListValidation[["MXD"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="MXD", Validation,Boyce=Boyce)
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
        Eval_JS <- list()
        for (i in 1:N) {
          RastPart[["MXS"]][[i]] <- c(predict(Model[[i]],PAtestM[[i]][, VarColT],clamp=F, type="cloglog"))
          PredPoint <- data.frame(PresAbse = PAtestM[[i]][, "PresAbse"], RastPart[["MXS"]][[i]])
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["MXS"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        
        #MXS Validation 
        Boyce <- mean(unlist(Boyce))
        Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
        if(is.null(repl)){
          ListValidation[["MXS"]] <- data.frame(Sp=spN[s], Algorithm="MXS",Partition=Part, Validation,Boyce=Boyce)
        }else{
          ListValidation[["MXS"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="MXS",Partition=Part, Validation,Boyce=Boyce)
        }
        
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            #Partial Thresholds
            Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
            PartRas <- STANDAR(predict(VariablesT,Model[[i]], clamp=F, type="cloglog"))
            if(N!=1){
              writeRaster(PartRas,paste(grep("MXS",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MXS",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              writeRaster(PartRas,paste(grep("MXS",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MXS",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(PartRas,paste(grep("MXS",foldPart,value=T),"/",spN[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MXS",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t],
                            paste(grep("MXS",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
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
          PredPoint <- c(predict(Model, SpDataTM[, VarColT], clamp=F, type="cloglog"))
          PredPoint <- data.frame(PresAbse = SpDataTM[, "PresAbse"], PredPoint)
          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                           PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
          
          #Variable Importance & Response Curves
          if(VarImP=="Y"){
            VarImp_RspCurv(Model=Model,Algorithm='MXS',spN=spN[s],SpDataT = SpDataTM,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }
          
          #Final Model Rasters
          ListSummary[["MXS"]] <- data.frame(Sp=spN[s], Algorithm="MXS", Thr)
          
          if(SaveFinal=="Y"){
            ListRaster[["MXS"]] <- STANDAR(predict(VariablesT,Model, clamp=F, type="cloglog"))
            names(ListRaster[["MXS"]]) <- spN[s]
          }
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
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["MXS"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #MXS Validation 
          Boyce <- mean(unlist(Boyce))
          Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
          if(is.null(repl)){
            ListValidation[["MXS"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="MXS", Validation,Boyce=Boyce)
          }else{
            ListValidation[["MXS"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="MXS", Validation,Boyce=Boyce)
          }
        }
      }
    }
    
    #MAXIMUM LIKELIHOOD (MLK)----
    if(any(Algorithm == 'MLK')) {
      if(any(lapply(PAtrain,function(x) nrow(x[x$PresAbse==1,]))<length(VarColT)*2)){
        ListValidation[["MLK"]] <- NULL
        ListRaster[["MLK"]] <- NULL
        ListSummary[["MLK"]] <- NULL
      }else{ 
        Model <- list()
        Fmula <- paste( " ~ ", paste(c(VarColT, paste("I(",VarColT, "^2)", sep = "")),
                                     collapse = " + "), sep = "")
        Fmula <- as.formula(Fmula)
        #MLK model
        for (i in 1:N) {
          # dataPr <- PAtrain[[i]][PAtrain[[i]]$PresAbse==1, c("x", "y")]
          x <- PAtrain[[i]][PAtrain[[i]][,"PresAbse"]==1, VarColT]
          z <- PAtrainM[[i]][PAtrainM[[i]][,"PresAbse"]==0, VarColT]
          # Model[[i]] <- maxlike(Fmula, stack(Variables), dataPr,
          #                       link=c("cloglog"),method="BFGS",
          #                       hessian = TRUE, removeDuplicates=FALSE)
          Model[[i]] <- maxlike(Fmula,x=x,z=z,
                               link=c("cloglog"),method="BFGS",
                               hessian = FALSE, removeDuplicates=FALSE)
        }
        
        #Evaluate Model
        if((is.null(Fut)==F && Tst=="Y")==F){
          Eval <- list()
          Boyce <- list()
          Eval_JS <- list()
          for (i in 1:N) {
            RastPart[["MLK"]][[i]] <- c(predict(Model[[i]], PAtest[[i]][, VarColT]))
            PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["MLK"]][[i]])
            Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                  PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2])
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(RastPart[["MLK"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          }
          
          #MXS Validation 
          Boyce <- mean(unlist(Boyce))
          Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
          if(is.null(repl)){
            ListValidation[["MLK"]] <- data.frame(Sp=spN[s], Algorithm="MLK",Partition=Part, Validation,Boyce=Boyce)
          }else{
            ListValidation[["MLK"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="MLK",Partition=Part, Validation,Boyce=Boyce)
          }
          
          #Save Partition Predictions
          if(Save=="Y"){
            for(i in 1:N){
              #Partial Thresholds
              Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
              PartRas <- STANDAR(predict(VariablesT,Model[[i]]))
              if(N!=1){
                writeRaster(PartRas,paste(grep("MLK",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
                Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="MLK",x=foldCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  writeRaster(PartRas>=Thr_Alg[t], 
                              paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                              format='GTiff',
                              overwrite=TRUE)
                }
              }
              if(is.null(repl)==F){
                writeRaster(PartRas,paste(grep("MLK",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
                Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="MLK",x=foldCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  writeRaster(PartRas>=Thr_Alg[t], 
                              paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                              format='GTiff',
                              overwrite=TRUE)
                }
              }else{
                writeRaster(PartRas,paste(grep("MLK",foldPart,value=T),"/",spN[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
                Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="MLK",x=foldCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  writeRaster(PartRas>=Thr_Alg[t],
                              paste(grep("MLK",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                              format='GTiff',
                              overwrite=TRUE)
                }
              }
            }
          }
          
          #Save final model
          if(per!=1 && repl==1 || per==1 || N!=1){
            x <- SpDataTM[SpDataTM[,"PresAbse"]==1, VarColT]
            z <- SpDataTM[SpDataTM[,"PresAbse"]==0, VarColT]
            Model <- maxlike(Fmula,x=x,z=z,
                             link=c("cloglog"),hessian = FALSE, 
                             method="BFGS",removeDuplicates=FALSE)
            PredPoint <- predict(Model, SpDataT[, VarColT])
            PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
            Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            #Final Thresholds
            Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
            
            #Variable Importance & Response Curves
            if(VarImP=="Y"){
              VarImp_RspCurv(Model=Model,Algorithm='MLK',spN=spN[s],SpDataT = SpDataTM,
                             VarColT=VarColT,Outcome=PredPoint$PredPoint)
            }
            
            #Final Model Rasters
            ListSummary[["MLK"]] <- data.frame(Sp=spN[s], Algorithm="MLK", Thr)
            
            if(SaveFinal=="Y"){
              ListRaster[["MLK"]] <- STANDAR(predict(VariablesT, Model))
              names(ListRaster[["MLK"]]) <- spN[s]
            }
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
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2])
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["MLK"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
            
            #MLK Validation 
            Boyce <- mean(unlist(Boyce))
            Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
            if(is.null(repl)){
              ListValidation[["MLK"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="MLK", Validation,Boyce=Boyce)
            }else{
              ListValidation[["MLK"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="MLK", Validation,Boyce=Boyce)
            }
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
        Model[[i]] <- ksvm(Fmula,data = dataPr,type="C-svc",kernel = "rbfdot",C = 1, prob.model=T)
      }
      
      #SVM evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        Eval_JS <- list()
        for (i in 1:N) {
          RastPart[["SVM"]][[i]] <- as.numeric(kernlab::predict(object=Model[[i]], newdata=PAtest[[i]][, VarColT],type="probabilities")[,2])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["SVM"]][[i]])
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["SVM"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        
        #SVM Validation 
        Boyce <- mean(unlist(Boyce))
        Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
        if(is.null(repl)){
          ListValidation[["SVM"]] <- data.frame(Sp=spN[s], Algorithm="SVM",Partition=Part, Validation,Boyce=Boyce)
        }else{
          ListValidation[["SVM"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="SVM",Partition=Part, Validation,Boyce=Boyce)
        }
        
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            #Partial Thresholds
            Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
            PartRas <- STANDAR(predict(VariablesT,Model[[i]]))
            if(N!=1){
              writeRaster(PartRas,paste(grep("SVM",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="SVM",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              writeRaster(PartRas,paste(grep("SVM",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="SVM",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(PartRas,paste(grep("SVM",foldPart,value=T),"/",spN[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="SVM",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t],
                            paste(grep("SVM",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
          }
        }
        
        # Save final model
        if(per!=1 && repl==1 || per==1 || N!=1){
          Model <- ksvm(Fmula,data = SpDataT[, c("PresAbse", VarColT)],type="C-svc",
                        kernel = "rbfdot",C = 1, prob.model=T)
          PredPoint <- as.numeric(kernlab::predict(object=Model, newdata=SpDataT[, VarColT],type="probabilities")[,2])
          PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
          
          #Variable Importance & Response Curves
          if(VarImP=="Y"){
            VarImp_RspCurv(Model=Model,Algorithm='SVM',spN=spN[s],SpDataT = SpDataT,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }
          
          #Final Model Rasters
          ListSummary[["SVM"]] <- data.frame(Sp=spN[s], Algorithm="SVM", Thr)
          
          if(SaveFinal=="Y"){
            ListRaster[["SVM"]] <- STANDAR(predict(VariablesT,Model))
            names(ListRaster[["SVM"]]) <- spN[s]
          }
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
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["SVM"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #SVM Validation 
          Boyce <- mean(unlist(Boyce))
          Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
          if(is.null(repl)){
            ListValidation[["SVM"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="SVM", Validation,Boyce=Boyce)
          }else{
            ListValidation[["SVM"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="SVM", Validation,Boyce=Boyce)
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
        Model[[i]] <- randomForest(as.factor(PresAbse)~.,data=dataPr[,c("PresAbse",VarColT)],
                                    importance=T, type="regression")
        # Model[[i]] <- tuneRF(dataPr[,-1], (dataPr[,1]), trace=F,
        #                      stepFactor=2, ntreeTry=1000, doBest=T, plot=F)
      }
      #RDF evaluation
      if((is.null(Fut)==F && Tst=="Y")==F){
        Eval <- list()
        Boyce <- list()
        Eval_JS <- list()
        for (i in 1:N) {
          RastPart[["RDF"]][[i]] <- as.vector(predict(Model[[i]], PAtest[[i]][, VarColT],type="prob")[,2])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["RDF"]][[i]])
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["RDF"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        }
        
        #RDF Validation 
        Boyce <- mean(unlist(Boyce))
        Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
        if(is.null(repl)){
          ListValidation[["RDF"]] <- data.frame(Sp=spN[s], Algorithm="RDF",Partition=Part, Validation,Boyce=Boyce)
        }else{
          ListValidation[["RDF"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="RDF",Partition=Part, Validation,Boyce=Boyce)
        }
        
        #Save Partition Predictions
        if(Save=="Y"){
          for(i in 1:N){
            #Partial Thresholds
            Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
            PartRas <- STANDAR(predict(VariablesT,Model[[i]]))
            if(N!=1){
              writeRaster(PartRas,paste(grep("RDF",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="RDF",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              writeRaster(PartRas,paste(grep("RDF",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="RDF",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t], 
                            paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }else{
              writeRaster(PartRas,paste(grep("RDF",foldPart,value=T),"/",spN[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
              Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="RDF",x=foldCat,value=T)
              for(t in 1:length(Thr_Alg)){
                writeRaster(PartRas>=Thr_Alg[t],
                            paste(grep("RDF",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
              }
            }
          }
        }
        
        # Save final model
        if(per!=1 && repl==1 || per==1 || N!=1){
          set.seed(0)
          Model <- randomForest(as.factor(PresAbse)~.,data=SpDataT[,c("PresAbse",VarColT)],
                                importance=T, type="classification")
          # Model <- tuneRF(SpDataT[,VarColT], (SpDataT[,"PresAbse"]), trace=F,
          #                 stepFactor=2, ntreeTry=500, doBest=T, plot = F)

          PredPoint <- predict(Model, SpDataT[, VarColT],type="prob")[,2]
          PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
          
          #Variable Importance & Response Curves
          if(VarImP=="Y"){
            VarImp_RspCurv(Model=Model,Algorithm='RDF',spN=spN[s],SpDataT = SpDataT,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }
          
          #Final Model Rasters
          ListSummary[["RDF"]] <- data.frame(Sp=spN[s], Algorithm="RDF", Thr)
          
          if(SaveFinal=="Y"){
            ListRaster[["RDF"]] <- STANDAR(predict(VariablesT,Model))
            names(ListRaster[["RDF"]]) <- spN[s]
          }
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
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["RDF"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          
          #RDF Validation 
          Boyce <- mean(unlist(Boyce))
          Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
          if(is.null(repl)){
            ListValidation[["RDF"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="RDF", Validation,Boyce=Boyce)
          }else{
            ListValidation[["RDF"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="RDF", Validation,Boyce=Boyce)
          }
        }
      }
    }
    
    #GENERALISED ADDITIVE MODEL (GAM)----
    if(any(Algorithm == 'GAM')) {
      if(any(lapply(PAtrain,function(x) nrow(x))<length(VarColT)*2)){
        ListValidation[["GAM"]] <- NULL
        ListRaster[["GAM"]] <- NULL
        ListSummary[["GAM"]] <- NULL
      }else{
        Model <- list()
        Fmula <- paste("s(", VarColT,",k=3)", sep="")
        Fmula <- paste("PresAbse", paste(Fmula, collapse = " + "), sep = " ~ ")
        Fmula <- as.formula(Fmula)
        
        #GAM model
        for (i in 1:N) {
          dataPr <- PAtrain[[i]][, c("PresAbse", VarColT)]
          Model[[i]] <- mgcv::gam(Fmula, data = dataPr, optimizer = c("outer", "newton"), 
                            select = T, family = binomial)
          }
        #GAM evaluation
        if((is.null(Fut)==F && Tst=="Y")==F){
          Eval <- list()
          Boyce <- list()
          Eval_JS <- list()
          for (i in 1:N) {
            RastPart[["GAM"]][[i]] <- as.vector(predict(Model[[i]], PAtest[[i]][, VarColT],type="response"))
            PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["GAM"]][[i]])
            Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                  PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2])
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(RastPart[["GAM"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          }
          
          #GAM Validation 
          Boyce <- mean(unlist(Boyce))
          Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
          if(is.null(repl)){
            ListValidation[["GAM"]] <- data.frame(Sp=spN[s], Algorithm="GAM",Partition=Part, Validation,Boyce=Boyce)
          }else{
            ListValidation[["GAM"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="GAM",Partition=Part, Validation,Boyce=Boyce)
          }
          
          #Save Partition Predictions
          if(Save=="Y"){
            for(i in 1:N){
              #Partial Thresholds
              Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
              PartRas <- STANDAR(predict(VariablesT,Model[[i]],type="response"))
              if(N!=1){
                writeRaster(PartRas,paste(grep("GAM",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
                Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="GAM",x=foldCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  writeRaster(PartRas>=Thr_Alg[t], 
                              paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                              format='GTiff',
                              overwrite=TRUE)
                }
              }
              if(is.null(repl)==F){
                writeRaster(PartRas,paste(grep("GAM",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
                Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="GAM",x=foldCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  writeRaster(PartRas>=Thr_Alg[t], 
                              paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                              format='GTiff',
                              overwrite=TRUE)
                }
              }else{
                writeRaster(PartRas,paste(grep("GAM",foldPart,value=T),"/",spN[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
                Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="GAM",x=foldCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  writeRaster(PartRas>=Thr_Alg[t],
                              paste(grep("GAM",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                              format='GTiff',
                              overwrite=TRUE)
                }
              }
            }
          }
          
          # Save final model
          if(per!=1 && repl==1 || per==1 || N!=1){
            Model <- mgcv::gam(formula=Fmula, data = SpDataT[, c("PresAbse",VarColT)], optimizer = c("outer", "newton"), 
                         select = T, family = binomial)
            PredPoint <- c(predict(Model, SpDataT[, VarColT],type="response"))
            PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
            Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            #Final Thresholds
            Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
            
            #Variable Importance & Response Curves
            if(VarImP=="Y"){
              VarImp_RspCurv(Model=Model,Algorithm='GAM',spN=spN[s],SpDataT = SpDataT,
                             VarColT=VarColT,Outcome=PredPoint$PredPoint)
            }
            
            #Final Model Rasters
            ListSummary[["GAM"]] <- data.frame(Sp=spN[s], Algorithm="GAM", Thr)
            if(SaveFinal=="Y"){
              ListRaster[["GAM"]] <- STANDAR(predict(VariablesT,Model,type="response"))
              names(ListRaster[["GAM"]]) <- spN[s]
            }
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
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2])
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["GAM"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
            
            
            #GAM Validation 
            Boyce <- mean(unlist(Boyce))
            Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
            if(is.null(repl)){
              ListValidation[["GAM"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="GAM", Validation,Boyce=Boyce)
            }else{
              ListValidation[["GAM"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="GAM", Validation,Boyce=Boyce)
            }
          }
        }
      }
    }
    
    #GENERALISED LINEAR MODEL (GLM) -----
    if(any(Algorithm == 'GLM')) {
      if(any(lapply(PAtrain,function(x) nrow(x))<length(VarColT)*2)){
        ListValidation[["GLM"]] <- NULL
        ListRaster[["GLM"]] <- NULL
        ListSummary[["GLM"]] <- NULL
      }else{
        Model <- list()
        Fmula <- paste( "PresAbse ~ ", paste(c(VarColT, paste("I(",VarColT, "^2)", sep = "")),
                                     collapse = " + "), sep = "")
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
          Eval_JS <- list()
          for (i in 1:N) {
            RastPart[["GLM"]][[i]] <- as.vector(predict(Model[[i]], PAtest[[i]][, VarColT],type="response"))
            PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["GLM"]][[i]])
            Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                  PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2])
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(RastPart[["GLM"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          }
          
          #GLM Validation 
          Boyce <- mean(unlist(Boyce))
          Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
          if(is.null(repl)){
            ListValidation[["GLM"]] <- data.frame(Sp=spN[s], Algorithm="GLM",Partition=Part, Validation,Boyce=Boyce)
          }else{
            ListValidation[["GLM"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="GLM",Partition=Part, Validation,Boyce=Boyce)
          }
          
          #Save Partition Predictions
          if(Save=="Y"){
            for(i in 1:N){
              #Partial Thresholds
              Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
              PartRas <- STANDAR(predict(VariablesT,Model[[i]],type="response"))
              if(N!=1){
                writeRaster(PartRas,paste(grep("GLM",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
                Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="GLM",x=foldCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  writeRaster(PartRas>=Thr_Alg[t], 
                              paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                              format='GTiff',
                              overwrite=TRUE)
                }
              }
              if(is.null(repl)==F){
                writeRaster(PartRas,paste(grep("GLM",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
                Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="GLM",x=foldCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  writeRaster(PartRas>=Thr_Alg[t], 
                              paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                              format='GTiff',
                              overwrite=TRUE)
                }
              }else{
                writeRaster(PartRas,paste(grep("GLM",foldPart,value=T),"/",spN[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
                Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="GLM",x=foldCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  writeRaster(PartRas>=Thr_Alg[t],
                              paste(grep("GLM",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                              format='GTiff',
                              overwrite=TRUE)
                }
              }
            }
          }
          
          # Save final model
          if(per!=1 && repl==1 || per==1 || N!=1){
            Model <- glm(Fmula, data = SpDataT[, c("PresAbse",VarColT)], family = binomial)
            PredPoint <- c(predict(Model, SpDataT[, VarColT],type="response"))
            PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
            Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            #Final Thresholds
            Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
            
            #Variable Importance & Response Curves
            if(VarImP=="Y"){
              VarImp_RspCurv(Model=Model,Algorithm='GLM',spN=spN[s],SpDataT = SpDataT,
                             VarColT=VarColT,Outcome=PredPoint$PredPoint)
            }
            
            #Final Model Rasters
            ListSummary[["GLM"]] <- data.frame(Sp=spN[s], Algorithm="GLM", Thr)
            if(SaveFinal=="Y"){
              ListRaster[["GLM"]] <- STANDAR(predict(VariablesT,Model,type="response"))
              names(ListRaster[["GLM"]]) <- spN[s]
            }
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
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2])
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["GLM"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
            
            
            #GLM Validation 
            Boyce <- mean(unlist(Boyce))
            Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
            if(is.null(repl)){
              ListValidation[["GLM"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="GLM", Validation,Boyce=Boyce)
            }else{
              ListValidation[["GLM"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="GLM", Validation,Boyce=Boyce)
            }
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
      Eval_JS <- list()
      for (i in 1:N) {
        RastPart[["GAU"]][[i]] <- predict(Model[[i]], PAtest[[i]][, VarColT])
        RastPart[["GAU"]][[i]] <- as.vector(RastPart[["GAU"]][[i]][,"posterior mode"])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["GAU"]][[i]])
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1,2],
                              PredPoint[PredPoint$PresAbse == 0,2])
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2])
        #Boyce Index
        Boyce[[i]] <- ecospat.boyce(RastPart[["GAU"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      
      #GAU Validation 
      Boyce <- mean(unlist(Boyce))
      Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
      if(is.null(repl)){
        ListValidation[["GAU"]] <- data.frame(Sp=spN[s], Algorithm="GAU",Partition=Part, Validation,Boyce=Boyce)
      }else{
        ListValidation[["GAU"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="GAU",Partition=Part, Validation,Boyce=Boyce)
      }
      
      #Save Partition Predictions
      if(Save=="Y"){
        for(i in 1:N){
          #Partial Thresholds
          Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
          PartRas <- STANDAR(predict(VariablesT,Model[[i]]))
          if(N!=1){
            writeRaster(PartRas,paste(grep("GAU",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern="GAU",x=foldCat,value=T)
            for(t in 1:length(Thr_Alg)){
              writeRaster(PartRas>=Thr_Alg[t], 
                          paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
          if(is.null(repl)==F){
            writeRaster(PartRas,paste(grep("GAU",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern="GAU",x=foldCat,value=T)
            for(t in 1:length(Thr_Alg)){
              writeRaster(PartRas>=Thr_Alg[t], 
                          paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }else{
            writeRaster(PartRas,paste(grep("GAU",foldPart,value=T),"/",spN[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern="GAU",x=foldCat,value=T)
            for(t in 1:length(Thr_Alg)){
              writeRaster(PartRas>=Thr_Alg[t],
                          paste(grep("GAU",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
      
      #Save final model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Model <- graf(SpDataT[,"PresAbse"], SpDataT[,VarColT])
        PredPoint <- data.frame(predict(Model, SpDataT[, VarColT]))$posterior.mode
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        #Final Thresholds
        Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
        
        #Variable Importance & Response Curves
        if(VarImP=="Y"){
          VarImp_RspCurv(Model=Model,Algorithm='GAU',spN=spN[s],SpDataT = SpDataT,
                         VarColT=VarColT,Outcome=PredPoint$PredPoint)
        }
        
        #Final Model Rasters
        ListSummary[["GAU"]] <- data.frame(Sp=spN[s], Algorithm="GAU", Thr)
        
        if(SaveFinal=="Y"){
          ListRaster[["GAU"]] <- STANDAR(predict.graf.raster(Model, VariablesT, type = "response", 
                                                             CI = 0.95, maxn = NULL)$posterior.mode)
          names(ListRaster[["GAU"]]) <- spN[s]
        }
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
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2])
        #Boyce Index
        Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["GAU"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        
        
        #GAU Validation 
        Boyce <- mean(unlist(Boyce))
        Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
        if(is.null(repl)){
          ListValidation[["GAU"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="GAU", Validation,Boyce=Boyce)
        }else{
          ListValidation[["GAU"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="GAU", Validation,Boyce=Boyce)
        }
      }
    }
    }
    
    #BOOSTED REGRESSION TREE (BRT) ----
    if(any(Algorithm == 'BRT')){
      if(any(lapply(PAtrain, function(x) nrow(x[x$PresAbse==1,])*0.75)<=21)){
        ListValidation[["BRT"]] <- NULL
        ListRaster[["BRT"]] <- NULL
        ListSummary[["BRT"]] <- NULL
      }else{
        Model <- list()
        #BRT model
        for (i in 1:N) {
          dataPr <- PAtrain[[i]]
          
          learn.rate <- 0.01
          ModelT <- NULL
          while(is.null(ModelT)){
            print(i)
            print(learn.rate)
            try(ModelT <- gbm.step(data=dataPr, gbm.x=VarColT, gbm.y="PresAbse", family = "bernoulli", 
                                       tree.complexity= 5, learning.rate=learn.rate, bag.fraction= 0.75,silent=T,
                                       plot.main = F))
            learn.rate <- learn.rate-0.0005
          }
          Model[[i]] <- ModelT
        }
        
        #BRT evaluation
        if((is.null(Fut)==F && Tst=="Y")==F){
          Eval <- list()
          Boyce <- list()
          Eval_JS <- list()
          for (i in 1:N) {
            RastPart[["BRT"]][[i]] <- predict.gbm(Model[[i]], PAtest[[i]][, VarColT],
                                                  n.trees=Model[[i]]$gbm.call$best.trees,type="response")
            PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["BRT"]][[i]])
            Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1,2],
                                         PredPoint[PredPoint$PresAbse == 0,2])
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2])
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(RastPart[["BRT"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          }
          
          #BRT Validation 
          Boyce <- mean(unlist(Boyce))
          Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
          if(is.null(repl)){
            ListValidation[["BRT"]] <- data.frame(Sp=spN[s], Algorithm="BRT",Partition=Part, Validation,Boyce=Boyce)
          }else{
            ListValidation[["BRT"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="BRT",Partition=Part, Validation,Boyce=Boyce)
          }
          
          #Save Partition Predictions
          if(Save=="Y"){
            for(i in 1:N){
              #Partial Thresholds
              Thr <- Thresholds_TMLA(Eval[[i]],Eval_JS,sensV)
              PartRas <- STANDAR(predict(VariablesT,Model[[i]],
                                             n.trees=Model[[i]]$gbm.call$best.trees,type="response"))
              if(N!=1){
                writeRaster(PartRas,paste(grep("BRT",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
                Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="BRT",x=foldCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  writeRaster(PartRas>=Thr_Alg[t], 
                              paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                              format='GTiff',
                              overwrite=TRUE)
                }
              }
              if(is.null(repl)==F){
                writeRaster(PartRas,paste(grep("BRT",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
                Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="BRT",x=foldCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  writeRaster(PartRas>=Thr_Alg[t], 
                              paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                              format='GTiff',
                              overwrite=TRUE)
                }
              }else{
                writeRaster(PartRas,paste(grep("BRT",foldPart,value=T),"/",spN[s],".tif", sep=""),
                            format='GTiff',
                            overwrite=TRUE)
                Thr_Alg <- Thr[Thr$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="BRT",x=foldCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  writeRaster(PartRas>=Thr_Alg[t],
                              paste(grep("BRT",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                              format='GTiff',
                              overwrite=TRUE)
                }
              }
            }
          }
        
          #Save final model
          if(per!=1 && repl==1 || per==1 || N!=1){
            learn.rate <- 0.005
            Model <- NULL
            while(is.null(Model)){
              try(Model <- gbm.step(data=SpDataT, gbm.x=VarColT, gbm.y="PresAbse", family = "bernoulli", 
                                         tree.complexity= 5, learning.rate=learn.rate, bag.fraction= bagFr,silent=T,
                                         plot.main = F))
              learn.rate <- learn.rate-0.0005
            }
            PredPoint <- predict.gbm(Model, SpDataT[, VarColT],
                                     n.trees=Model$gbm.call$best.trees)
            PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
            Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            #Final Thresholds
            Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
            
            #Variable Importance & Response Curves
            if(VarImP=="Y"){
              VarImp_RspCurv(Model=Model,Algorithm='BRT',spN=spN[s],SpDataT = SpDataT,
                             VarColT=VarColT,Outcome=PredPoint$PredPoint)
            }
            
            #Final Model Rasters
            ListSummary[["BRT"]] <- data.frame(Sp=spN[s], Algorithm="BRT", Thr)
            
            if(SaveFinal=="Y"){
              ListRaster[["BRT"]] <- STANDAR(predict(VariablesT,Model,
                                                 n.trees=Model$gbm.call$best.trees,type="response"))
              names(ListRaster[["BRT"]]) <- spN[s]
            }
            if(is.null(Fut)==F){
              for(k in 1:length(VariablesP)){
                ListFut[[ProjN[k]]][["BRT"]] <- predict(VariablesP[[k]],Model,
                                                            n.trees=Model$gbm.call$best.trees,type="response")
                if(maxValue(ListFut[[ProjN[k]]][["BRT"]])==0){
                  ListFut[[ProjN[k]]][["BRT"]] <- ListFut[[ProjN[k]]][["BRT"]]
                }else{
                  ListFut[[ProjN[k]]][["BRT"]] <- STANDAR(ListFut[[ProjN[k]]][["BRT"]])
                }
              }
            }
          }
        }else{
          Eval <- list()
          Boyce <- list()
          for(k in 1:length(VariablesP)){
            ListFut[[ProjN[k]]][["BRT"]] <- predict(Model,VariablesP[[k]],
                                                        n.trees=Model$gbm.call$best.trees,type="response")
            if(maxValue(ListFut[[ProjN[k]]][["BRT"]])==0){
              ListFut[[ProjN[k]]][["BRT"]] <- ListFut[[ProjN[k]]][["BRT"]]
            }else{
              ListFut[[ProjN[k]]][["BRT"]] <- STANDAR(ListFut[[ProjN[k]]][["BRT"]])
            }
            
            PredPoint <- extract(ListFut[[ProjN[k]]][["BRT"]], PAtest[[i]][, c("x", "y")])
            PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
            Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2])
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["BRT"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
            
            
            #BRT Validation 
            Boyce <- mean(unlist(Boyce))
            Validation<-Validation_Table_TMLA(Eval,Eval_JS,N)
            if(is.null(repl)){
              ListValidation[["BRT"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="BRT", Validation,Boyce=Boyce)
            }else{
              ListValidation[["BRT"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="BRT", Validation,Boyce=Boyce)
            }
          }
        }
      }
    }

    #Final models----
    if(per!=1 && repl==1 || per==1 || N!=1){
      if(SaveFinal=="Y"){
        if((is.null(Fut)==F && Tst=="Y")==F){
          Thr <- lapply(ListSummary, '[', c('THR','THR_VALUE'))
          for(i in 1:length(ListRaster)){
            writeRaster(round(ListRaster[[i]], 4),
                        paste(folders[i], '/',spN[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            Thr_Alg <- Thr[[i]][Thr[[i]]$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern=names(Thr)[i],x=foldCat,value=T)
            for(t in 1:length(Thr_Alg)){
              writeRaster(ListRaster[[i]]>=Thr_Alg[t], 
                          paste(foldCatAlg[t], '/',spN[s],".tif", sep=""),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }    
      #Save Projections
      if(is.null(Fut)==F){
        Thr <- lapply(ListSummary, '[', c('THR','THR_VALUE'))
        for(p in 1:length(ListFut)){
          for(o in 1:length(ListFut[[p]])){
            writeRaster(ListFut[[p]][[o]],file.path(ModFut[p],Algorithm[o],spN[s]),
                        format='GTiff',overwrite=TRUE)
            Thr_Alg <- Thr[[o]][Thr[[o]]$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern=Algorithm[o],x=foldCat,value=T)
            for(t in 1:length(Thr_Alg)){
              writeRaster(ListFut[[p]][[o]]>=Thr_Alg[t], 
                          file.path(ModFut[p],Algorithm[o],Threshold[t],paste0(spN[s],".tif")),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
    }
    
    #Adjust invasion cenario for ensemble
    if((is.null(Fut)==F && Tst=="Y")){
      RastPart <- list(ListFut)
    }
    
    # Ensemble-----

    # Mean Ensemble----
    if(any(PredictType=="MEAN")){
      
      #Partial Models Ensemble
      Final <- do.call(Map, c(rbind,RastPart))
      Final <- lapply(Final, function (x) colMeans(x))

      # Threshold
      Eval <- list()
      Boyce <- list()
      Eval_JS <- list()
      for(i in 1:N){
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], Final[[i]])
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2])
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      
      #MEAN Validation 
      Boyce <- mean(unlist(Boyce))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N)
      
      if(is.null(repl)){
        ListValidation[["MEAN"]] <- data.frame(Sp=spN[s], Algorithm="MEA", Validation,Boyce=Boyce)
      }else{
        ListValidation[["MEAN"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="MEA", Validation,Boyce=Boyce)          
      }
 
      #Final Model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Final <- brick(ListRaster)
        Final <- calc(Final,mean)
        PredPoint <- extract(Final, SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2])
          
        #Final Thresholds
        Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
        ListSummary[["MEAN"]] <- data.frame(Sp=spN[s], Algorithm="MEA", Thr)
        if(SaveFinal=="Y"){  
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
          }
        }
        
        #Future Projection
        if(is.null(Fut)==F && Tst!="Y"){
          for(p in 1:length(ListFut)){
            Final <- brick(ListFut[[p]])
            Final <- calc(Final,mean)
            
            writeRaster(Final, 
                        file.path(ModFut[p],"Ensemble","MEAN",spN[s]),
                        format='GTiff',
                        overwrite=TRUE)
              writeRaster(Final>=as.numeric(Thr), 
                          file.path(ModFut[p],"Ensemble","MEAN","BIN",paste(spN[s],sep="_")),
                          format='GTiff',
                          overwrite=TRUE)
          }
        }
      }
    }
    
    # Weighted Mean Ensemble----
    if(any(PredictType=="W_MEAN")){
      ListValidationT <- ldply(ListValidation,data.frame,.id=NULL)
      ListValidationT <- ListValidationT[ListValidationT$Algorithm%in%Algorithm,]
      colnames(ListValidationT) <- c("Sp","Algorithm","Partition","AUC","MAX_KAPPA","MAX_TSS","JACCARD",
                                     "SORENSEN","FPB","BOYCE")
      
      #Partial Models Ensemble
      Final <- do.call(Map, c(rbind,RastPart))
      ThResW <- unlist(ListValidationT[Threshold])
      Final <- lapply(Final, function(x) sweep(x, 2, ThResW, '*'))
      Final <- lapply(Final, function (x) colMeans(x))
      
      # Threshold
      Eval <- list()
      Boyce <- list()
      Eval_JS <- list()
      for(i in 1:N){
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], Final[[i]])
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2])
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      
      #W_MEAN Validation 
      Boyce <- mean(unlist(Boyce))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N)
      if(is.null(repl)){
        ListValidation[["W_MEAN"]] <- data.frame(Sp=spN[s], Algorithm="WMEA", Validation,Boyce=Boyce)
      }else{
        ListValidation[["W_MEAN"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="WMEA", Validation,Boyce=Boyce)          
      }
      
      #Final Model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Final <- brick(ListRaster)
        Final <- calc(Final, function(x) x*ThResW)
        Final <- calc(Final,mean)
        PredPoint <- extract(Final, SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                     a=PredPoint[PredPoint$PresAbse == 0, 2])
          
        #Final Thresholds
        Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
        ListSummary[["W_MEAN"]] <- data.frame(Sp=spN[s], Algorithm="WMEA", Thr)
        
        if(SaveFinal=="Y"){  
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
        }
        
        #Future Projection
        if(is.null(Fut)==F && Tst!="Y"){
          for(p in 1:length(ListFut)){
            Final <- brick(ListFut[[p]])
            Final <- calc(Final,mean)
            
            writeRaster(Final, 
                        file.path(ModFut[p],"Ensemble","W_MEAN",spN[s]),
                        format='GTiff',
                        overwrite=TRUE)
            writeRaster(Final>=as.numeric(Thr), 
                        file.path(ModFut[p],"Ensemble","W_MEAN","BIN",paste(spN[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }

    # Superior Ensemble----
    if(any(PredictType=='SUP')){
      ListValidationT <- ldply(ListValidation,data.frame,.id=NULL)
      ListValidationT <- ListValidationT[ListValidationT$Algorithm%in%Algorithm,]
      if(is.null(repl)){
        colnames(ListValidationT) <- c("Sp","Algorithm","Partition","AUC","MAX_KAPPA","MAX_TSS","JACCARD",
                                     "SORENSEN","FPB","BOYCE")
      }else{
        colnames(ListValidationT) <- c("Sp","Replicate","Algorithm","Partition","AUC","MAX_KAPPA","MAX_TSS","JACCARD",
                                       "SORENSEN","FPB","BOYCE")
      }
      Best <- ListValidationT[which(unlist(ListValidationT[Threshold])>=mean(unlist(ListValidationT[Threshold]))),"Algorithm"]
      W <- names(ListRaster)%in%Best
      
      #Partial Models
      Final <- do.call(Map, c(rbind,RastPart[W]))
      Final <- lapply(Final, function (x) colMeans(x))

      # Threshold
      Eval <- list()
      Boyce <- list()
      Eval_JS <- list()
      for(i in 1:N){
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], Final[[i]])
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      
      #SUP Validation 
      Boyce <- mean(unlist(Boyce))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N)
      if(is.null(repl)){
        ListValidation[["SUP"]] <- data.frame(Sp=spN[s], Algorithm="SUP", Validation,Boyce=Boyce)
      }else{
        ListValidation[["SUP"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="SUP", Validation,Boyce=Boyce)          
      }
      
      #Final Model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Final <- brick(ListRaster[W])
        Final <- calc(Final,mean)
        PredPoint <- raster::extract(Final, SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                     a=PredPoint[PredPoint$PresAbse == 0, 2])
          
        #Final Thresholds
        Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
        ListSummary[["SUP"]] <- data.frame(Sp=spN[s], Algorithm="SUP",Thr)
        
        if(SaveFinal=="Y"){
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
        }

        #Future Projection
        if(is.null(Fut)==F && Tst!="Y"){
          for(p in 1:length(ListFut)){
            Final <- brick(ListFut[[p]][W])
            Final <- calc(Final,mean)
            
            writeRaster(Final, 
                        file.path(ModFut[p],"Ensemble","SUP",paste(spN[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
            writeRaster(Final>=as.numeric(Thr), 
                        file.path(ModFut[p],"Ensemble","SUP","BIN",paste(spN[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }

    # With PCA ------
    if (any(PredictType == 'PCA')) {

      #Partial Models Ensemble
      if(any(lapply(RastPart, function(x) length(x))>1)){
        Final <- do.call(Map, c(cbind, RastPart))
        Final <- lapply(Final, function(x) as.numeric(princomp(x)$scores[,1]))
        Final <- lapply(Final, function(x) (x-min(x))/(max(x)-min(x)))
      }else{
        Final <- do.call(cbind,lapply(RastPart, function(x) do.call(cbind,x)))
        Final <- as.numeric(princomp(Final)$scores[,1])
        Final <- list((Final-min(Final))/(max(Final)-min(Final)))
      }

      # Threshold
      Eval <- list()
      Boyce <- list()
      Eval_JS <- list()
      for(i in 1:N){
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], Final[[i]])
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      
      #PCA Validation 
      Boyce <- mean(unlist(Boyce))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N)
      if(is.null(repl)){
        ListValidation[["PCA"]] <- data.frame(Sp=spN[s], Algorithm="PCA", Validation,Boyce=Boyce)
      }else{
        ListValidation[["PCA"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="PCA", Validation,Boyce=Boyce)          
      }
      
      #Final Model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Final <- brick(ListRaster)
        Final <- PCA_ENS_TMLA(Final)
        PredPoint <- extract(Final, SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                   PredPoint[PredPoint$PresAbse == 0, 2])
        #Final Thresholds
        Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
        ListSummary[["PCA"]] <- data.frame(Sp=spN[s], Algorithm="PCA", Thr)
        if(SaveFinal=="Y"){ 
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
        }

        #Future Projection
        if(is.null(Fut)==F && Tst!="Y"){
          for(p in 1:length(ListFut)){
            Final <- brick(ListFut[[p]])
            Final <- PCA_ENS_TMLA(Final)
            
            writeRaster(Final, 
                        file.path(ModFut[p],"Ensemble","PCA",spN[s]),
                        format='GTiff',
                        overwrite=TRUE)
              writeRaster(Final>=as.numeric(Thr), 
                          file.path(ModFut[p],"Ensemble","PCA","BIN",paste(spN[s],sep="_")),
                          format='GTiff',
                          overwrite=TRUE)
          }
        }
      }
    }
        
    # With PCA over the Mean(Superior) Ensemble----
    if (any(PredictType == 'PCA_SUP')) {
      ListValidationT <- ldply(ListValidation,data.frame,.id=NULL)
      ListValidationT <- ListValidationT[ListValidationT$Algorithm%in%Algorithm,]
      if(is.null(repl)){
        colnames(ListValidationT) <- c("Sp","Algorithm","Partition","AUC","MAX_KAPPA","MAX_TSS","JACCARD",
                                       "SORENSEN","FPB","BOYCE")
      }else{
        colnames(ListValidationT) <- c("Sp","Replicate","Algorithm","Partition","AUC","MAX_KAPPA","MAX_TSS","JACCARD",
                                       "SORENSEN","FPB","BOYCE")
      }
      Best <- ListValidationT[which(unlist(ListValidationT[Threshold])>=mean(unlist(ListValidationT[Threshold]))),"Algorithm"]
      W <- names(ListRaster)%in%Best
      
      #Partial Models
      if(any(lapply(RastPart, function(x) length(x))>1)){
        Final <- do.call(Map, c(cbind, RastPart[W]))
        Final <- lapply(Final, function(x) as.numeric(princomp(x)$scores[,1]))
        Final <- lapply(Final, function(x) (x-min(x))/(max(x)-min(x)))
      }else{
        Final <- do.call(cbind,lapply(RastPart[W], function(x) do.call(cbind,x)))
        Final <- as.numeric(princomp(Final)$scores[,1])
        Final <- list((Final-min(Final))/(max(Final)-min(Final)))
      }

      # Threshold
      Eval <- list()
      Boyce <- list()
      Eval_JS <- list()
      for(i in 1:N){
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], Final[[i]])
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      
      #PCA_SUP Validation 
      Boyce <- mean(unlist(Boyce))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N)
      if(is.null(repl)){
        ListValidation[["PCS"]] <- data.frame(Sp=spN[s], Algorithm="PCS", Validation,Boyce=Boyce)
      }else{
        ListValidation[["PCS"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="PCS", Validation,Boyce=Boyce)          
      }
      
      #Final Model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Final <- brick(ListRaster[W])
        Final <- PCA_ENS_TMLA(Final)
        PredPoint <- extract(Final, SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                                                                                         PredPoint[PredPoint$PresAbse == 0, 2])
        #Final Thresholds
        Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
        ListSummary[["PCS"]] <- data.frame(Sp=spN[s], Algorithm="PCS", Thr)
        
        if(SaveFinal=="Y"){
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
        }
        
        #Future Projection
        if(is.null(Fut)==F && Tst!="Y"){
          for(p in 1:length(ListFut)){
            Final <- brick(ListFut[[p]][W])
            Final <- PCA_ENS_TMLA(Final)
            
            writeRaster(Final, 
                        file.path(ModFut[p],"Ensemble","PCA_SUP",paste(spN[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
            writeRaster(Final>=as.numeric(Thr), 
                        file.path(ModFut[p],"Ensemble","PCA_SUP","BIN",paste(spN[s],sep="_")),
                        format='GTiff',
                        overwrite=TRUE)
          }
        }
      }
    }
    
    #With PCA over the threshold Ensemble----
    if (any(PredictType == 'PCA_THR')) {
      ListValidationT <- ldply(ListSummary,data.frame,.id=NULL)
      ListValidationT <- ListValidationT[ListValidationT$Algorithm%in%Algorithm,]
      
      #Partial Models
      Final <- do.call(Map, c(cbind, RastPart))
      ValidTHR <- ListValidationT[grepl(Threshold,ListValidationT["THR"],ignore.case = T),"THR_VALUE"]
      for (p in 1:length(Final)){
        Final[[p]] <- sapply(seq(1:length(ValidTHR)),function(x){ifelse(Final[[p]][,x]>=ValidTHR[x],Final[[p]][,x],0)})
        Final[[p]] <- as.numeric(princomp(Final[[p]])$scores[,1])
        Final[[p]] <- (Final[[p]]-min(Final[[p]]))/(max(Final[[p]])-min(Final[[p]]))
      }
      
      # Evaluation
      Eval <- list()
      Boyce <- list()
      Eval_JS <- list()
      for(i in 1:N){
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], Final[[i]])
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
      }
      
      #PCA_SUP Validation 
      Boyce <- mean(unlist(Boyce))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N)
      
      if(is.null(repl)){
        ListValidation[["PCT"]] <- data.frame(Sp=spN[s], Algorithm="PCT", Validation,Boyce=Boyce)
      }else{
        ListValidation[["PCT"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="PCT", Validation,Boyce=Boyce)          
      }
      
      #Final Model
      if(per!=1 && repl==1 || per==1 || N!=1){
        Final <- brick(ListRaster)
        ValidTHR <- ListValidationT[grepl(Threshold,ListValidationT["THR"],ignore.case = T),"THR_VALUE"]
        for(k in 1:nlayers(Final)){
          FinalSp <- Final[[k]]
          FinalSp[FinalSp<ValidTHR[k]] <- 0
          Final[[k]] <- FinalSp
        }
        Final <- PCA_ENS_TMLA(Final)
        PredPoint <- extract(Final, SpDataT[, c("x", "y")])
        PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        
        #Final Thresholds
        Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)
        ListSummary[["PCT"]] <- data.frame(Sp=spN[s], Algorithm="PCT", Threshold=Thr)
        
        if(SaveFinal=="Y"){  
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
        }
        
        #Future Projection
        if(is.null(Fut)==F && Tst!="Y"){
          for(p in 1:length(ListFut)){
              Final <- brick(ListFut[[p]])
              
              #Select only values above the Threshold
              for(k in Algorithm){
                FinalSp <- Final[[k]]
                FinalSp[FinalSp<ListValidationT[ListValidationT$Algorithm==k,"THR"]] <- 0
                if(all(na.omit(FinalSp[])==0)){
                  Final[Final[[k]]] <- NULL
                }else{
                Final[[k]] <- FinalSp
                }
              }
              
              Final <- PCA_ENS_TMLA(Final)
  
              writeRaster(Final, 
                          file.path(ModFut[p],"Ensemble","PCA_THR",paste(spN[s],sep="_")),
                          format='GTiff',
                          overwrite=TRUE)
                writeRaster(Final>=unlist(Thr), 
                            file.path(ModFut[p],"Ensemble","PCA_THR","BIN",paste(spN[s],sep="_")),
                            format='GTiff',
                            overwrite=TRUE)
            }
        }
      }
    }
    
    #Final Data Frame Results
    result <- ldply(ListValidation,data.frame,.id=NULL)
    resultII <- ldply(ListSummary,data.frame,.id=NULL)

    out <- list(Validation = result,
                Summary = resultII)
    return(out)
  }#Fecha loop Especie

# Save .txt with the models performance---- 
FinalValidation <- data.frame(data.table::rbindlist(do.call(rbind,lapply(results, "[", "Validation"))))
FinalSummary <- data.frame(data.table::rbindlist(do.call(rbind,lapply(results, "[", "Summary"))))

write.table(FinalValidation,paste(DirSave, VALNAME, sep = '/'),sep="\t",
            col.names = T,row.names=F)
if(per!=1 && repl==1 || per==1 || N!=1){
  write.table(FinalSummary,paste(DirSave, VALNAMEII, sep = '/'),sep="\t",
              col.names = T,row.names=F)
}

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
       paste('No_species:',length(spN)),
       paste("Threshold:",Threshold),
       matrix(spN))
  lapply(InfoModeling, write, 
         paste(DirSave, "/InfoModeling.txt", sep=""), append=TRUE, 
         ncolumns=20, sep='\t')
  stopCluster(cl)
}
