## Written by Santiago Velazco & Andre Andrade

FitENM_TMLA_Parallel <- function(RecordsData,
                   Variables,
                   Fut=NULL,
                   Part,
                   Algorithm,
                   PredictType,
                   VarImp,
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
                   extrapolation=extrapolation,
                   ensemble_metric=ensemble_metric,
                   sensV,
                   cores) {

  Ti <- Sys.time()
  options(warn = -1)

  #Start Cluster
  cl <- parallel::makeCluster(cores,outfile="")
  doParallel::registerDoParallel(cl)

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
          #Binary ensemble projection
          DirENSFCat <- file.path(sort(rep(DirENS,length(Threshold))),Threshold)
          sapply(DirENSFCat,function(x) dir.create(x,recursive = T))
        }
      }
    }
    DirENSFCat <- file.path(sort(rep(ModFut,length(PredictType))),"Ensemble",PredictType,Threshold)
  }

  # raster::extracting enviromental variables----
  Ncol <- ncol(RecordsData) + 1
  RecordsData <- stats::na.omit(data.frame(RecordsData, raster::extract(Variables, RecordsData[, c("x", "y")])))
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
  message(paste("Total species to be modeled", length(spN)))

  # Number of partition
  N <- as.numeric(max(RecordsData[, "Partition"]))

  #Txt of Final tables
  if(is.null(repl)==F){
    VALNAME <- paste('Evaluation_Table','_',sep="",repl,'.txt' )
  }else{
    VALNAME <- paste('Evaluation_Table.txt' )
  }
  VALNAMEII <- paste('Thresholds_Algorithms.txt' )


  # Backqround points----
  if (!is.null(DirMask)) {
    if (any(c("MXD","MXS","MLK","ENF") %in% Algorithm)) {
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
             ab.0 <- data.frame(dismo::randomPoints(msk2,p=RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1, c("x", "y")],10000))
             var.0 <- data.frame(raster::extract(Variables,ab.0))
           }else{
             ab.0 <-
               dismo::randomPoints(msk2, p=RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1, c("x", "y")],abs(NCell - nrow(RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1,])))
             var.0 <- data.frame(raster::extract(Variables, ab.0))
           }
           ab.0 <- cbind(rep(spN[i],nrow(ab.0)),ab.0,rep(x,nrow(ab.0)),rep(0,nrow(ab.0)),var.0)
           colnames(ab.0) <- colnames(RecordsData)
           ab0L[[x]] <- stats::na.omit(ab.0)
           rm(var.0)
           # RecordsDataM <- rbind(RecordsDataM[[i]],ab.0)
           # rm(ab.0)
         }
        ab.0 <- plyr::ldply(ab0L,data.frame,.id=NULL)
        return(ab.0)
       }

       RecordsDataM <- lapply(seq_along(RecordsDataM), function(x) rbind(RecordsDataM[[x]], ab.0[[x]]))
       RecordsDataM <- plyr::ldply(RecordsDataM,data.frame,.id=NULL)
       cols <-  c("x","y","Partition","PresAbse",names(Variables))
       RecordsDataM[,cols] <-  apply(RecordsDataM[,cols], 2, function(x) as.numeric(as.character(x)))
     }
   }else{
     if (any(c("MXD","MXS","MLK","ENF") %in% Algorithm)) {
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
            ab.0 <- data.frame(dismo::randomPoints(msk2,p=RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1, c("x", "y")],10000))
            var.0 <- raster::extract(Variables,ab.0)
          }else{
            ab.0 <-
              dismo::randomPoints(msk2, p=RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1, c("x", "y")],(NCell - nrow(RecordsDataM[[i]][RecordsDataM[[i]]$PresAbse==1,])))
            var.0 <- raster::extract(Variables, ab.0)
          }
          ab.0 <- cbind(rep(spN[i],nrow(ab.0)),ab.0,rep(x,nrow(ab.0)),rep(0,nrow(ab.0)),var.0)
          colnames(ab.0) <- colnames(RecordsData)
          ab0L[[x]] <- stats::na.omit(ab.0)
          rm(var.0)
          # RecordsDataM <- rbind(RecordsDataM,ab.0)
          # rm(ab.0)
        }
        ab.0 <- plyr::ldply(ab0L,data.frame,.id=NULL)
        return(ab.0)
      }

      RecordsDataM <- lapply(seq_along(RecordsDataM), function(x) rbind(RecordsDataM[[x]], ab.0[[x]]))
      RecordsDataM <- plyr::ldply(RecordsDataM,data.frame,.id=NULL)
      cols <-  c("x","y","Partition","PresAbse",names(Variables))
      RecordsDataM[,cols] = apply(RecordsDataM[,cols], 2, function(x) as.numeric(as.character(x)))
    }
   }

  if(!exists("RecordsDataM")){
    RecordsDataM <- NULL
  }

  #MOP Calculation----
  #Within the extent (for M-Restriction)
  if(is.null(repl)||repl==1){
    ansE <- extrapolation
    if(ansE==TRUE){
      dir.create(file.path(DirSave,"Extrapolation"))
      DirProj <- file.path(DirSave,"Extrapolation")
      MOP(
        Variables = list(Variables),
        RecordsData = RecordsData,
        VarCol = VarCol,
        DirProj = DirProj,
        DirMask=DirMask
      )
    }
    #For projections
    if(!is.null(Fut)){
      ansE <- extrapolation
      if(ansE==TRUE){
        for(i in 1:length(ModFut)){
          dir.create(file.path(ModFut[i],"Extrapolation"))
        }
        DirProj <- file.path(ModFut,"Extrapolation")
        MOP(
          Variables = Fut,
          RecordsData = RecordsData,
          VarCol = VarCol,
          DirProj = DirProj,
          DirMask=DirMask
        )
      }
    }
  }

  #Define N due to Partition Method
  if(Part=="BOOT" || Part=="KFOLD"){
    N <- 1
  }else{
    N <- N
  }

  cat("Fitting Models....\n")

  # Construction of models LOOP-----
  results <- foreach(s = 1:length(spN),
                   .packages = c("raster", "dismo",
                                 "kernlab", "randomForest", "maxnet", "maxlike",
                                 "plyr", "mgcv", "RStoolbox", "adehabitatHS",
                                 "caret", "visreg", "glmnet", "gbm","dplyr"),
                     .export = c( "Validation2_0", "maxnet2",
                                  "PCA_ENS_TMLA", "predict.maxnet", "boycei"
                                  ,"partial_roc","trap_roc",
                                  "STANDAR", "STANDAR_FUT", "Eval_Jac_Sor_TMLA",
                                  "Validation_Table_TMLA",
                                  "Thresholds_TMLA", "VarImp_RspCurv",
                                  "hingeval", "ecospat.boyce",
                                  "PREDICT_DomainMahal","rem_out","PREDICT_ENFA",
                                  "graf","capture.all","graf.fit.laplace","graf.fit.ep",
                                  "cov.SE","d0","d1","d2","d3","psiline","psi",
                                  "cov.SE.d1","predict.graf.raster","predict.graf","pred")) %dopar% {

    #Results Lists
    ListRaster <- as.list(Algorithm)
    names(ListRaster) <- Algorithm
    
    #Partial Models List
    RasT <- as.list(Algorithm)
    names(RasT) <- Algorithm
    
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
    if (any(c("MXD","MXS","MLK","ENF") %in% Algorithm)) {
      SpDataM <- RecordsDataM[RecordsDataM[, "sp"] == spN[s], ]
    }

    #Include MSDM----
    if(is.null(DirMSDM)==F){
      if(grepl("XY",DirMSDM)){
        MSDM <- raster::stack(file.path(DirMSDM,list.files(DirMSDM,pattern=".tif")))
        names(MSDM) <- c("Lat","Long")
      }else{
        MSDM <- raster(file.path(DirMSDM,paste(spN[s],".tif",sep="")))
        names(MSDM) <- "MSDM"
      }
      SpDataT <- cbind(SpData,raster::extract(MSDM,SpData[c("x","y")]))
      colnames(SpDataT) <- c(colnames(SpData),names(MSDM))
      VariablesT <- raster::brick(raster::stack(Variables,MSDM))
      if (any(c("MXD","MXS","MLK","ENF") %in% Algorithm)) {
        SpDataTM <- cbind(SpDataM,raster::extract(MSDM,SpDataM[c("x","y")]))
        colnames(SpDataTM) <- c(colnames(SpDataM),names(MSDM))
      }
      VarColT <- c(VarCol,names(MSDM))
    }else{
      VariablesT <- Variables
      VarColT <- VarCol
      SpDataT <- SpData
      if (any(c("MXD","MXS","MLK","ENF") %in% Algorithm)) {
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
    if (any(c("MXD","MXS","MLK","ENF") %in% Algorithm)) {
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
        Model[[i]] <- dismo::bioclim(dataPr[, VarColT])
      }

      #BIO evaluation
      if((is.null(Fut)==F && !is.null(Tst))==F){
        Eval <- list()
        Eval_JS <- list()
        Boyce <- list()
        pROC <- list()
        Area <- list()
        for (i in 1:N) {
          RastPart[["BIO"]][[i]] <- raster::predict(Model[[i]], PAtest[[i]][, VarColT])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["BIO"]][[i]])
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                           a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["BIO"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          #Percentage of Predicted Area
          RasT[["BIO"]] <- dismo::predict(Model[[i]], VariablesT)

          ArT <- NULL
          for (j in Thr){
            RasL <- RasT[["BIO"]]>=j
            ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=RasT[["BIO"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,percentage=50)$pROC_summary
          
          #Save Partition Predictions
          if(Save=="Y"){
            #Partial Thresholds
            if(N!=1){
              raster::writeRaster(RasT[["BIO"]],paste(grep("BIO",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="BIO",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["BIO"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              raster::writeRaster(RasT[["BIO"]],paste(grep("BIO",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="BIO",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["BIO"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",repl,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }else{
              raster::writeRaster(RasT[["BIO"]],paste(grep("BIO",foldPart,value=T),"/",paste0(spN[s],repl),".tif", 
                                                sep=""),format='GTiff',overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="BIO",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["BIO"]]>=Thr_Alg[t],
                                    paste(grep("BIO",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
          }
        }

        #BIO Validation
        pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
        pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
        pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
        pROC <- mean(unlist(lapply(pROC, `[`, 1)))
        BoyceSD <- stats::sd(unlist(Boyce))
        Boyce <- mean(unlist(Boyce))
        AreaSD <- apply(do.call("rbind",Area),2,sd)
        Area <- colMeans(do.call("rbind",Area))
        Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
        if(is.null(repl)){
          ListValidation[["BIO"]] <- data.frame(Sp=spN[s], Algorithm="BIO",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }else{
          ListValidation[["BIO"]] <- data.frame(Sp=spN[s],Replicate=repl,Algorithm="BIO",Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }


      # Save final model
      if(repl==1 || is.null(repl)){
        if(is.null(repl) && N==1){
          Model <- dismo::bioclim(SpDataT[SpDataT[,"PresAbse"]==1 & SpDataT[,"Partition"]==1, VarColT]) # only presences
          FinalModelT <- dismo::predict(Model, VariablesT)
          FinalModel <- STANDAR(FinalModelT)
          PredPoint <- raster::extract(FinalModel,SpDataT[SpDataT[,"Partition"]==1,2:3])
          PredPoint <- data.frame(PresAbse = SpDataT[SpDataT[,"Partition"]==1, "PresAbse"], PredPoint)
        }else{
          Model <- dismo::bioclim(SpDataT[SpDataT[,"PresAbse"]==1, VarColT]) # only presences
          FinalModelT <- dismo::predict(Model, VariablesT)
          FinalModel <- STANDAR(FinalModelT)
          PredPoint <- raster::extract(FinalModel,SpDataT[,2:3])
          PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        }
        #Final Model "Evaluation"
        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                         PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        #Final Thresholds
        Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

        #Variable Importance & Response Curves
        if(VarImp==TRUE){
          VarImp_RspCurv(Model=Model,Algorithm='BIO',folders=folders,spN=spN[s],SpDataT = SpDataT,
                         VarColT=VarColT,Outcome=PredPoint$PredPoint)
        }

        #Final Model Rasters
        ListSummary[["BIO"]] <- data.frame(Sp=spN[s], Algorithm="BIO", Thr)
        if(SaveFinal=="Y"){
          ListRaster[["BIO"]] <- FinalModel
          names(ListRaster[["BIO"]]) <- spN[s]
        }
        #Future Projections
        if(is.null(Fut)==F){
          for(k in 1:length(VariablesP)){
            ListFut[[ProjN[k]]][["BIO"]] <- STANDAR_FUT(dismo::predict(VariablesP[[k]], Model),FinalModelT)
          }
        }
      }
      }else{
        Eval <- list()
        Boyce <- list()
        Eval_JS <- list()
        pROC <- list()
        Area <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["BIO"]] <- STANDAR(dismo::predict(VariablesP[[k]],Model[[i]]))
          PredPoint <- raster::extract(ListFut[[ProjN[k]]][["BIO"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["BIO"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          #Percentae of Predicted Area
          ArT <- NULL
          for (j in Thr){
            RasL <- ListFut[[ProjN[k]]][["BIO"]]>=j
            ArT <- c(ArT,sum(na.omit(values(RasL)))/length(na.omit(values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=ListFut[[ProjN[k]]][["BIO"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,
                                              percentage=50)$pROC_summary


          #BIO Validation
          pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
          pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
          pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
          pROC <- mean(unlist(lapply(pROC, `[`, 1)))
          BoyceSD <- stats::sd(unlist(Boyce))
          Boyce <- mean(unlist(Boyce))
          AreaSD <- apply(do.call("rbind",Area),2,sd)
          Area <- colMeans(do.call("rbind",Area))
          Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
          if(is.null(repl)){
            ListValidation[["BIO"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="BIO", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }else{
            ListValidation[["BIO"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="BIO", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
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
      if((is.null(Fut)==F && !is.null(Tst))==F){
        Eval <- list()
        Boyce <- list()
        Eval_JS <- list()
        pROC <- list()
        Area <- list()
        for (i in 1:N) {
          RastPart[["DOM"]][[i]] <- raster::predict(Model[[i]], PAtest[[i]][VarColT])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["DOM"]][[i]])
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Thr2Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["DOM"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          #Percentage of Predicted Area
          RasT[["DOM"]] <- PREDICT_DomainMahal(mod = Model[[i]], variables = VariablesT)
          ArT <- NULL
          for (j in Thr){
            RasL <- RasT[["DOM"]]>=j
            ArT <- c(ArT,sum(na.omit(values(RasL)))/length(na.omit(values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=RasT[["DOM"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,percentage=50)$pROC_summary
          #Save Partition Predictions
          if(Save=="Y"){
            if(N!=1){
              raster::writeRaster(
                RasT[["DOM"]],
                paste(grep("DOM", foldPart, value = T), "/", spN[s], "_", i, sep = ""),
                format = 'GTiff',
                overwrite = TRUE )
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="DOM",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["DOM"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              raster::writeRaster(RasT[["DOM"]],paste(grep("DOM",foldPart,value=T),"/",spN[s],"_",repl,sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="DOM",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["DOM"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",repl,sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }else{
              raster::writeRaster(RasT[["DOM"]],paste(grep("DOM",foldPart,value=T),"/",paste0(spN[s],repl),sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="DOM",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["DOM"]]>=Thr_Alg[t],
                                    paste(grep("DOM",foldCatAlg[t],value=T), '/',spN[s],sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
          }
        }

        #DOM Validation
        pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
        pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
        pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
        pROC <- mean(unlist(lapply(pROC, `[`, 1)))
        BoyceSD <- stats::sd(unlist(Boyce))
        Boyce <- mean(unlist(Boyce))
        AreaSD <- apply(do.call("rbind",Area),2,sd)
        Area <- colMeans(do.call("rbind",Area))
        Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
        if(is.null(repl)){
          ListValidation[["DOM"]] <- data.frame(Sp=spN[s], Algorithm="DOM",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }else{
          ListValidation[["DOM"]] <- data.frame(Sp=spN[s],Replicate=repl,Algorithm="DOM",Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }


        # Save final model
        if(repl==1 || is.null(repl)){
          if(is.null(repl)){
            Model <-
              dismo::domain(x = VariablesT, p = SpDataT[SpDataT[, "PresAbse"] == 1 &
                                                  SpDataT[, "Partition"] == 1, c("x", "y")]) # only presences
            FinalModelT <- PREDICT_DomainMahal(mod = Model, variables = VariablesT)
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel, SpDataT[SpDataT[,"Partition"]==1,c("x","y")])
            PredPoint <- data.frame(PresAbse = SpDataT[SpDataT[,"Partition"]==1, "PresAbse"], PredPoint)
          }else{
            Model <- dismo::domain(SpDataT[SpDataT[,"PresAbse"]==1, VarColT]) # only presences
            FinalModelT <- PREDICT_DomainMahal(mod = Model, variables = VariablesT)
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel,SpDataT[,c("x","y")])
            PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          }
          #Final Model Thresholds
          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                  PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                       a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

          #Variable Importance & Response Curves
          if(VarImp==TRUE){
            VarImp_RspCurv(Model=Model,Algorithm='DOM',folders=folders,spN=spN[s],SpDataT = SpDataT,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }

          #Final Model Rasters
          ListSummary[["DOM"]] <- data.frame(Sp=spN[s], Algorithm="DOM", Thr)

          if(SaveFinal=="Y"){
            ListRaster[["DOM"]] <- FinalModel
            names(ListRaster[["DOM"]]) <- spN[s]
          }
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              PreFut <- PREDICT_DomainMahal(Model,VariablesP[[k]])
              ListFut[[ProjN[k]]][["DOM"]] <- STANDAR_FUT(PreFut, FinalModelT)
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        Eval_JS <- list()
        pROC <- list()
        Area <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["DOM"]] <- raster::predict(Model[[i]], PAtest[[i]][VarColT])
          PredPoint <- raster::extract(ListFut[[ProjN[k]]][["DOM"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["DOM"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          #Percentae of Predicted Area
          RasT[["DOM"]] <- STANDAR(PREDICT_DomainMahal(mod = Model[[i]], variables = VariablesP[[k]]))
          ArT <- NULL
          for (j in Thr){
            RasL <- RasT[["DOM"]]>=j
            ArT <- c(ArT,sum(na.omit(values(RasL)))/length(na.omit(values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=RasT[["DOM"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,percentage=50)$pROC_summary

          #DOM Validation
          pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
          pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
          pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
          pROC <- mean(unlist(lapply(pROC, `[`, 1)))
          BoyceSD <- stats::sd(unlist(Boyce))
          Boyce <- mean(unlist(Boyce))
          AreaSD <- apply(do.call("rbind",Area),2,sd)
          Area <- colMeans(do.call("rbind",Area))
          Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
          if(is.null(repl)){
            ListValidation[["DOM"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="DOM", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }else{
            ListValidation[["DOM"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="DOM", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
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
      if((is.null(Fut)==F && !is.null(Tst))==F){
        Eval <- list()
        Eval_JS <- list()
        Boyce <- list()
        pROC <- list()
        Area <- list()
        for (i in 1:N) {
          RastPart[["MAH"]][[i]] <- raster::predict(Model[[i]], PAtest[[i]][VarColT])
          RastPart[["MAH"]][[i]][RastPart[["MAH"]][[i]][] < -10] <- -10
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["MAH"]][[i]])
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["MAH"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          #Percentage of Predicted Area
          RasT[["MAH"]] <- PREDICT_DomainMahal(mod = Model[[i]], variables = VariablesT)
          RasT[["MAH"]][RasT[["MAH"]][] < -10] <- -10
          ArT <- NULL
          for (j in Thr){
            RasL <- RasT[["MAH"]]>=j
            ArT <- c(ArT,sum(na.omit(values(RasL)))/length(na.omit(values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=RasT[["MAH"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,percentage=50)$pROC_summary
          
          #Save Partition Predictions
          if(Save=="Y"){
            if(N!=1){
              raster::writeRaster(
                RasT[["MAH"]],
                paste(grep("MAH", foldPart, value = T), "/", spN[s], "_", i, sep = ""),
                format = 'GTiff',
                overwrite = TRUE )
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MAH",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["MAH"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              raster::writeRaster(RasT[["MAH"]],paste(grep("MAH",foldPart,value=T),"/",spN[s],"_",repl,sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MAH",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["MAH"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",repl,sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }else{
              raster::writeRaster(RasT[["MAH"]],paste(grep("MAH",foldPart,value=T),"/",paste0(spN[s],repl),sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MAH",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["MAH"]]>=Thr_Alg[t],
                                    paste(grep("MAH",foldCatAlg[t],value=T), '/',spN[s],sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
          }
        }
        
        #MAH Validation
        pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
        pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
        pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
        pROC <- mean(unlist(lapply(pROC, `[`, 1)))
        BoyceSD <- stats::sd(unlist(Boyce))
        Boyce <- mean(unlist(Boyce))
        AreaSD <- apply(do.call("rbind",Area),2,sd)
        Area <- colMeans(do.call("rbind",Area))
        Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
        if(is.null(repl)){
          ListValidation[["MAH"]] <- data.frame(Sp=spN[s], Algorithm="MAH",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }else{
          ListValidation[["MAH"]] <- data.frame(Sp=spN[s],Replicate=repl,Algorithm="MAH",Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }

        # Save final model
        if(repl==1 || is.null(repl)){
          if(is.null(repl)){
            Model <-
              mahal(x = VariablesT, p = SpDataT[SpDataT[, "PresAbse"] == 1 &
                                                  SpDataT[, "Partition"] == 1, c("x", "y")]) # only presences
            FinalModelT <- PREDICT_DomainMahal(mod = Model, variables = VariablesT)
            FinalModelT[FinalModelT[] < -10] <- -10
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel, SpDataT[SpDataT[,"Partition"]==1,c("x","y")])
            PredPoint <- data.frame(PresAbse = SpDataT[SpDataT[,"Partition"]==1, "PresAbse"], PredPoint)
          }else{
            Model <- mahal(SpDataT[SpDataT[,"PresAbse"]==1, VarColT]) # only presences
            FinalModelT <- PREDICT_DomainMahal(mod = Model, variables = VariablesT)
            FinalModelT[FinalModelT[] < -10] <- -10
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel,SpDataT[,c("x","y")])
            PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          }
          #Final Model Thresholds
          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                       a=PredPoint[PredPoint$PresAbse == 0, 2])
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

          #Variable Importance & Response Curves
          if(VarImp==TRUE){
            VarImp_RspCurv(Model=Model,Algorithm='MAH',folders=folders,spN=spN[s],SpDataT = SpDataT,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }

          #Final Model Rasters
          ListSummary[["MAH"]] <- data.frame(Sp=spN[s], Algorithm="MAH", Thr)

          if(SaveFinal=="Y"){
            ListRaster[["MAH"]] <- FinalModel
            names(ListRaster[["MAH"]]) <- spN[s]
          }
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              PreFut <- PREDICT_DomainMahal(Model,VariablesP[[k]])
              PreFut[PreFut[] < -10] <- -10
              ListFut[[ProjN[k]]][["MAH"]] <-
                STANDAR_FUT(PreFut, FinalModelT)
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        Eval_JS <- list()
        pROC <- list()
        Area <- list()
        for(k in 1:length(VariablesP)){
          PredPoint <- raster::predict(Model[[i]], PAtest[[i]][VarColT])
          PredPoint[PredPoint[] < -10] <- -10
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["MAH"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          #Percentae of Predicted Area
          RasT[["MAH"]] <- RasT[["MAH"]] <- STANDAR(PREDICT_DomainMahal(mod = Model[[i]], variables = VariablesP[[k]]))
          RasT[["MAH"]][RasT[["MAH"]][] < -10] <- -10
          ArT <- NULL
          for (j in Thr){
            RasL <- RasT[["MAH"]]>=j
            ArT <- c(ArT,sum(na.omit(values(RasL)))/length(na.omit(values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=RasT[["MAH"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,
                                              percentage=50)$pROC_summary
          
          
          #MAH Validation
          pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
          pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
          pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
          pROC <- mean(unlist(lapply(pROC, `[`, 1)))
          BoyceSD <- stats::sd(unlist(Boyce))
          Boyce <- mean(unlist(Boyce))
          AreaSD <- apply(do.call("rbind",Area),2,sd)
          Area <- colMeans(do.call("rbind",Area))
          Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
          if(is.null(repl)){
            ListValidation[["MAH"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="MAH", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }else{
            ListValidation[["MAH"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="MAH", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
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
        dudi <- ade4::dudi.pca(dataPr[, VarColT],scannf = FALSE)
        Model[[i]] <- adehabitatHS::madifa(dudi,dataPr$PresAbse,scannf = FALSE)
      }

      #ENF evaluation
      if((is.null(Fut)==F && !is.null(Tst))==F){
        Eval <- list()
        Boyce <- list()
        Eval_JS <- list()
        pROC <- list()
        Area <- list()
        for (i in 1:N) {
          RastPart[["ENF"]][[i]] <- PREDICT_ENFA(Model[[i]],PAtestM[[i]][, VarColT])
          PredPoint <- data.frame(PresAbse = PAtestM[[i]][, "PresAbse"], RastPart[["ENF"]][[i]])
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["ENF"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          #Percentae of Predicted Area
          RasT[["ENF"]] <- PREDICT_ENFA(Model[[i]],VariablesT,PAtrainM[[i]])
          ArT <- NULL
          for (j in Thr){
            RasL <- RasT[["ENF"]]>=j
            ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=RasT[["ENF"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,percentage=50)$pROC_summary
          #Save Partition Predictions
          if(Save=="Y"){
            if(N!=1){
              raster::writeRaster(RasT[["ENF"]],paste(grep("ENF",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="ENF",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["ENF"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              raster::writeRaster(RasT[["ENF"]],paste(grep("ENF",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),format='GTiff',overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="ENF",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["ENF"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",repl,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }else{
              raster::writeRaster(RasT[["ENF"]],paste(grep("ENF",foldPart,value=T),"/",paste0(spN[s],repl),".tif", sep=""),format='GTiff',overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="ENF",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["ENF"]]>=Thr_Alg[t],
                                    paste(grep("ENF",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
          }
        }
        
        #ENF Validation
        pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
        pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
        pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
        pROC <- mean(unlist(lapply(pROC, `[`, 1)))
        BoyceSD <- stats::sd(unlist(Boyce))
        Boyce <- mean(unlist(Boyce))
        AreaSD <- apply(do.call("rbind",Area),2,sd)
        Area <- colMeans(do.call("rbind",Area))
        Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
        if(is.null(repl)){
          ListValidation[["ENF"]] <- data.frame(Sp=spN[s], Algorithm="ENF",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }else{
          ListValidation[["ENF"]] <- data.frame(Sp=spN[s],Replicate=repl,Algorithm="ENF",Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }

        # Save final model
        if(repl==1 || is.null(repl)){
          if(is.null(repl)){
            dudi <- ade4::dudi.pca(SpDataTM[SpDataTM[,"Partition"]==1, VarColT],scannf = FALSE)
            Model <- adehabitatHS::madifa(dudi,SpDataTM[SpDataTM[,"Partition"]==1, "PresAbse"],scannf = FALSE)
            FinalModelT <- PREDICT_ENFA(Model,VariablesT,PAtrainM[[1]])
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel,SpDataTM[,c("x","y")])
            PredPoint <- data.frame(PresAbse = SpDataTM[, "PresAbse"], PredPoint)
          }else{
            dudi <- ade4::dudi.pca(SpDataT[, VarColT],scannf = FALSE)
            Model <- adehabitatHS::madifa(dudi,SpDataT$PresAbse,scannf = FALSE)
            FinalModelT <- PREDICT_ENFA(Model,VariablesT,PAtrainM[[1]])
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel,SpDataTM[,c("x","y")])
            PredPoint <- data.frame(PresAbse = SpDataTM[, "PresAbse"], PredPoint)
          }
          #Final Model Thresholds
          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])

          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

          #Variable Importance & Response Curves
          if(VarImp==TRUE){
            VarImp_RspCurv(Model=Model,Algorithm='ENF',folders=folders,spN=spN[s],SpDataT = SpDataTM,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }

          #Final Model Rasters
          ListSummary[["ENF"]] <- data.frame(Sp=spN[s], Algorithm="ENF", Thr)
          if(SaveFinal=="Y"){
            ListRaster[["ENF"]] <- FinalModel
            names(ListRaster[["ENF"]]) <- spN[s]
          }

          #Future Projections
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              PredRas <- Predict_ENFA(Model,VariablesP[[k]],PAtrainM[[1]])
              PredRas <- STANDAR_FUT(PredRas,FinalModelT)
              if(minValue(PredRas)<0){
                PredRas <- PredRas-minValue(PredRas)
              }
              ListFut[[ProjN[k]]][["ENF"]] <- PredRas
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        Eval_JS <- list()
        pROC <- list()
        Area <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["ENF"]] <- STANDAR(Predict_ENFA(Model[[i]],VariablesP[[k]],PAtrainM[[1]]))
          PredPoint <- raster::extract(ListFut[[ProjN[k]]][["ENF"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["ENF"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          #Percentae of Predicted Area
          ArT <- NULL
          for (j in Thr){
            RasL <- ListFut[[ProjN[k]]][["ENF"]]>=j
            ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=ListFut[[ProjN[k]]][["ENF"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,
                                              percentage=50)$pROC_summary
          
          
          #ENF Validation
          pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
          pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
          pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
          pROC <- mean(unlist(lapply(pROC, `[`, 1)))
          BoyceSD <- stats::sd(unlist(Boyce))
          Boyce <- mean(unlist(Boyce))
          AreaSD <- apply(do.call("rbind",Area),2,sd)
          Area <- colMeans(do.call("rbind",Area))
          Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
          if(is.null(repl)){
            ListValidation[["ENF"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="ENF", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }else{
            ListValidation[["ENF"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="ENF", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
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
                               maxnet::maxnet.formula(dataPr[,"PresAbse"],
                                              dataPr[,VarColT], classes="default"))
      }
      #MXD evaluation
      if((is.null(Fut)==F && !is.null(Tst))==F){
        Eval <- list()
        Eval_JS <- list()
        Boyce <- list()
        pROC <- list()
        Area <- list()
        for (i in 1:N) {
          RastPart[["MXD"]][[i]] <- c(predict(Model[[i]], PAtestM[[i]][, VarColT], clamp=F, type="cloglog"))
          PredPoint <- data.frame(PresAbse = PAtestM[[i]][, "PresAbse"], RastPart[["MXD"]][[i]])
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["MXD"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          #Percentae of Predicted Area
          RasT[["MXD"]] <- raster::predict(VariablesT,Model[[i]], clamp=F, type="cloglog")
          ArT <- NULL
          for (j in Thr){
            RasL <- RasT[["MXD"]]>=j
            ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=RasT[["MXD"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,percentage=50)$pROC_summary
          
          #Save Partition Predictions
          if(Save=="Y"){
            if(N!=1){
              raster::writeRaster(RasT[["MXD"]],paste(grep("MXD",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MXD",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["MXD"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              raster::writeRaster(RasT[["MXD"]],paste(grep("MXD",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),format='GTiff',overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MXD",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["MXD"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }else{
              raster::writeRaster(RasT[["MXD"]],paste(grep("MXD",foldPart,value=T),"/",spN[s],".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MXD",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["MXD"]]>=Thr_Alg[t],
                                    paste(grep("MXD",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
          }
        }
        
        #MXD Validation
        pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
        pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
        pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
        pROC <- mean(unlist(lapply(pROC, `[`, 1)))
        BoyceSD <- stats::sd(unlist(Boyce))
        Boyce <- mean(unlist(Boyce))
        AreaSD <- apply(do.call("rbind",Area),2,sd)
        Area <- colMeans(do.call("rbind",Area))
        Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
        if(is.null(repl)){
          ListValidation[["MXD"]] <- data.frame(Sp=spN[s], Algorithm="MXD",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }else{
          ListValidation[["MXD"]] <- data.frame(Sp=spN[s],Replicate=repl,Algorithm="MXD",Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }

        #Save final model
        if(repl==1 || is.null(repl)){
          if(is.null(repl) && N==1){
            Model <- maxnet2(SpDataTM[SpDataTM$Partition==1,"PresAbse"], SpDataTM[SpDataTM$Partition==1,VarColT], f =
                               maxnet::maxnet.formula(SpDataTM[SpDataTM$Partition==1,"PresAbse"], SpDataTM[SpDataTM$Partition==1,VarColT],
                                              classes="default"))
            FinalModelT <- predict(VariablesT,Model, clamp=F, type="cloglog")
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel,SpDataTM[SpDataTM$Partition==1, 2:3])
            PredPoint <- data.frame(PresAbse = SpDataTM[SpDataTM$Partition==1, "PresAbse"], PredPoint)
          }else{
            Model <- maxnet2(SpDataTM[,"PresAbse"], SpDataTM[,VarColT], f =
                               maxnet::maxnet.formula(SpDataTM[,"PresAbse"], SpDataTM[,VarColT], classes="default"))
            FinalModelT <- predict(VariablesT,Model, clamp=F, type="cloglog")
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel,SpDataTM[, 2:3])
            PredPoint <- data.frame(PresAbse = SpDataTM[, "PresAbse"], PredPoint)
          }
          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                           PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

          #Variable Importance & Response Curves
          if(VarImp==TRUE){
            VarImp_RspCurv(Model=Model,Algorithm='MXD',folders=folders,spN=spN[s],SpDataT = SpDataTM,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }

          #Final Model Rasters
          ListSummary[["MXD"]] <- data.frame(Sp=spN[s], Algorithm="MXD", Thr)

          if(SaveFinal=="Y"){
            ListRaster[["MXD"]] <- FinalModel
            names(ListRaster[["MXD"]]) <- spN[s]
          }
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["MXD"]] <- STANDAR_FUT(predict(VariablesP[[k]], Model,clamp=F, type="cloglog"),FinalModelT)
            }
          }
        }
      }else{
        Eval <- list()
        Eval_JS <- list()
        Boyce <- list()
        pROC <- list()
        Area <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["MXD"]] <- STANDAR(predict(VariablesP[[k]],Model[[i]],clamp=F, type="cloglog"))
          PredPoint <- raster::extract(ListFut[[ProjN[k]]][["MXD"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["MXD"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          #Percentae of Predicted Area
          ArT <- NULL
          for (j in Thr){
            RasL <- ListFut[[ProjN[k]]][["MXD"]]>=j
            ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=ListFut[[ProjN[k]]][["MXD"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,
                                              percentage=50)$pROC_summary
          
          
          #MXD Validation
          pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
          pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
          pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
          pROC <- mean(unlist(lapply(pROC, `[`, 1)))
          BoyceSD <- stats::sd(unlist(Boyce))
          Boyce <- mean(unlist(Boyce))
          AreaSD <- apply(do.call("rbind",Area),2,sd)
          Area <- colMeans(do.call("rbind",Area))
          Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
          if(is.null(repl)){
            ListValidation[["MXD"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="MXD", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }else{
            ListValidation[["MXD"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="MXD", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
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
        Model[[i]] <- maxnet2(p=dataPr[,"PresAbse"], data=dataPr[,VarColT], f =
                                maxnet::maxnet.formula(dataPr[,"PresAbse"],dataPr[,VarColT], classes="lq"))
      }
      #MXS evaluation
      if((is.null(Fut)==F && !is.null(Tst))==F){
        Eval <- list()
        Eval_JS <- list()
        Boyce <- list()
        pROC <- list()
        Area <- list()
        for (i in 1:N) {
          RastPart[["MXS"]][[i]] <- c(raster::predict(Model[[i]],PAtestM[[i]][, VarColT],clamp=F, type="cloglog"))
          PredPoint <- data.frame(PresAbse = PAtestM[[i]][, "PresAbse"], RastPart[["MXS"]][[i]])
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["MXS"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          #Percentage of Predicted Area
          RasT[["MXS"]] <- raster::predict(VariablesT,Model[[i]], clamp=F, type="cloglog")
          
          ArT <- NULL
          for (j in Thr){
            RasL <- RasT[["MXS"]]>=j
            ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=RasT[["MXS"]],
                                   test_data=PredPoint[PredPoint$PresAbse==1,2],
                                   error=5,iterations=500,percentage=50)$pROC_summary
          
          
          #Save Partition Predictions
          if(Save=="Y"){
            if(N!=1){
              raster::writeRaster(RasT[["MXS"]],paste(grep("MXS",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MXS",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["MXS"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              raster::writeRaster(RasT[["MXS"]],paste(grep("MXS",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),format='GTiff',overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MXS",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["MXS"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }else{
              raster::writeRaster(RasT[["MXS"]],paste(grep("MXS",foldPart,value=T),"/",spN[s],".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="MXS",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["MXS"]]>=Thr_Alg[t],
                                    paste(grep("MXS",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
          }
        }
        
        #MXS Validation
        pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
        pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
        pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
        pROC <- mean(unlist(lapply(pROC, `[`, 1)))
        BoyceSD <- stats::sd(unlist(Boyce))
        Boyce <- mean(unlist(Boyce))
        AreaSD <- apply(do.call("rbind",Area),2,sd)
        Area <- colMeans(do.call("rbind",Area))
        Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
        if(is.null(repl)){
          ListValidation[["MXS"]] <- data.frame(Sp=spN[s], Algorithm="MXS",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }else{
          ListValidation[["MXS"]] <- data.frame(Sp=spN[s],Replicate=repl,Algorithm="MXS",Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }

        #Save final model
        if(repl==1 || is.null(repl)){
          if(is.null(repl) && N==1){
            Model <- maxnet2(SpDataTM[SpDataTM$Partition==1,"PresAbse"], SpDataTM[SpDataTM$Partition==1,VarColT], f =
                               maxnet::maxnet.formula(SpDataTM[SpDataTM$Partition==1,"PresAbse"], SpDataTM[SpDataTM$Partition==1,VarColT],
                                              classes="lq"))
            FinalModelT <- raster::predict(VariablesT,Model, clamp=F, type="cloglog")
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel,SpDataTM[SpDataTM$Partition==1, 2:3])
            PredPoint <- data.frame(PresAbse = SpDataTM[SpDataTM$Partition==1, "PresAbse"], PredPoint)
          }else{
            Model <- maxnet2(SpDataTM[,"PresAbse"], SpDataTM[,VarColT], f =
                               maxnet::maxnet.formula(SpDataTM[,"PresAbse"], SpDataTM[,VarColT], classes="lq"))
            FinalModelT <- raster::predict(VariablesT,Model, clamp=F, type="cloglog")
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel,SpDataTM[, 2:3])
            PredPoint <- data.frame(PresAbse = SpDataTM[, "PresAbse"], PredPoint)
          }

          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                           PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

          #Variable Importance & Response Curves
          if(VarImp==TRUE){
            VarImp_RspCurv(Model=Model,Algorithm='MXS',folders=folders,spN=spN[s],SpDataT = SpDataTM,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }

          #Final Model Rasters
          ListSummary[["MXS"]] <- data.frame(Sp=spN[s], Algorithm="MXS", Thr)

          if(SaveFinal=="Y"){
            ListRaster[["MXS"]] <- FinalModel
            names(ListRaster[["MXS"]]) <- spN[s]
          }
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["MXS"]] <- STANDAR_FUT(raster::predict(VariablesP[[k]], Model,clamp=F, type="cloglog"),FinalModelT)
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["MXS"]] <- STANDAR(raster::predict(VariablesP[[k]],Model[[i]],clamp=F, type="cloglog"))

          PredPoint <- raster::extract(ListFut[[ProjN[k]]][["MXS"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["MXS"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          #Percentae of Predicted Area
          ArT <- NULL
          for (j in Thr){
            RasL <- ListFut[[ProjN[k]]][["MXS"]]>=j
            ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=ListFut[[ProjN[k]]][["MXS"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,
                                              percentage=50)$pROC_summary
          
          
          #MXS Validation
          pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
          pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
          pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
          pROC <- mean(unlist(lapply(pROC, `[`, 1)))
          BoyceSD <- stats::sd(unlist(Boyce))
          Boyce <- mean(unlist(Boyce))
          AreaSD <- apply(do.call("rbind",Area),2,sd)
          Area <- colMeans(do.call("rbind",Area))
          Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
          if(is.null(repl)){
            ListValidation[["MXS"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="MXS", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }else{
            ListValidation[["MXS"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="MXS", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
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
        RastPart[["MLK"]] <- NULL
        RasT[["MLK"]] <- NULL
      }else{
        Model <- list()
        Fmula <- paste( " ~ ", paste(c(VarColT, paste("I(",VarColT, "^2)", sep = "")),
                                     collapse = " + "), sep = "")
        Fmula <- stats::as.formula(Fmula)
        EnvMLK <- raster::stack((VariablesT-colMeans(na.omit(raster::values(VariablesT))))/apply(na.omit(raster::values(VariablesT)),2,stats::sd))
        #MLK model
        for (i in 1:N) {
          dataPr <- PAtrain[[i]][PAtrain[[i]]$PresAbse==1, c("x", "y")]
          # x <- PAtrain[[i]][PAtrainM[[i]][,"PresAbse"]==1, VarColT]
          # z <- PAtrainM[[i]][PAtrainM[[i]][,"PresAbse"]==0, VarColT]
          # Model[[i]] <- maxlike::maxlike(formula=Fmula,rasters=EnvMLK
          #                                                ,points=dataPr,
          #                                                link=("cloglog"),
          #                                                method="BFGS",
          #                                                hessian = FALSE, 
          #                                                removeDuplicates=FALSE
          #                                                ,savedata=T)
          Model[[i]] <- tryCatch (expr={maxlike::maxlike(formula=Fmula,rasters=EnvMLK
                                                              ,points=dataPr,
                                                              link=("cloglog"),
                                                              method="BFGS",
                                                              hessian = FALSE,
                                                              removeDuplicates=FALSE
                                                              ,savedata=T)},
                         error=function(e) {
                           message("Trying Nelder-Mead Optimization")
                           maxlike::maxlike(formula=Fmula,
                                                          rasters=EnvMLK
                                                          ,points=dataPr,
                                                          link=c("cloglog"),
                                                          method="Nelder-Mead",
                                                          hessian = FALSE,
                                                          removeDuplicates=FALSE
                                                          ,savedata=T)})
        }

        #Evaluate Model
        if((is.null(Fut)==F && !is.null(Tst))==F){
          Eval <- list()
          Eval_JS <- list()
          Boyce <- list()
          pROC <- list()
          Area <- list()
          for (i in 1:N) {
            RastPart[["MLK"]][[i]] <- c(predict(Model[[i]], PAtest[[i]][, VarColT]))
            PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["MLK"]][[i]])
            Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                      PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                           a=PredPoint[PredPoint$PresAbse == 0, 2])
            
            #Thresholds and Final Evaluation
            Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
            Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
            Threshold <- Threshold[order(Thr)]
            Thr <- sort(Thr)
            Thr <- ifelse(Thr<0,0,Thr)
            Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
            Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
            Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(RastPart[["MLK"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
            #Percentae of Predicted Area
            RasT[["MLK"]] <- predict(Model[[i]])
            ArT <- NULL
            for (j in Thr){
              RasL <- RasT[["MLK"]]>=j
              ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
            }
            Area[[i]] <- round(ArT*100,3)
            
            #PartialROC
            pROC[[i]] <- partial_roc(prediction=RasT[["MLK"]],
                                                test_data=PredPoint[PredPoint$PresAbse==1,2],
                                                error=5,iterations=500,percentage=50)$pROC_summary
            
            #Save Partition Predictions
            if(Save=="Y"){
              if(N!=1){
                raster::writeRaster(RasT[["MLK"]],paste(grep("MLK",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
                Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="MLK",x=PartCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  raster::writeRaster(RasT[["MLK"]]>=Thr_Alg[t],
                                      paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                      format='GTiff',
                                      overwrite=TRUE)
                }
              }
              if(is.null(repl)==F){
                raster::writeRaster(RasT[["MLK"]],paste(grep("MLK",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),format='GTiff',overwrite=TRUE)
                Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="MLK",x=PartCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  raster::writeRaster(RasT[["MLK"]]>=Thr_Alg[t],
                                      paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                      format='GTiff',
                                      overwrite=TRUE)
                }
              }else{
                raster::writeRaster(RasT[["MLK"]],paste(grep("MLK",foldPart,value=T),"/",spN[s],".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
                Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="MLK",x=PartCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  raster::writeRaster(RasT[["MLK"]]>=Thr_Alg[t],
                                      paste(grep("MLK",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                                      format='GTiff',
                                      overwrite=TRUE)
                }
              }
            }
          }
          
          #MLK Validation
          pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
          pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
          pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
          pROC <- mean(unlist(lapply(pROC, `[`, 1)))
          BoyceSD <- stats::sd(unlist(Boyce))
          Boyce <- mean(unlist(Boyce))
          AreaSD <- apply(do.call("rbind",Area),2,sd)
          Area <- colMeans(do.call("rbind",Area))
          Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
          if(is.null(repl)){
            ListValidation[["MLK"]] <- data.frame(Sp=spN[s], Algorithm="MLK",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }else{
            ListValidation[["MLK"]] <- data.frame(Sp=spN[s],Replicate=repl,Algorithm="MLK",Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }

          #Save final model
          if(repl==1 || is.null(repl)){
            if(is.null(repl) && N==1){
              Model <- tryCatch (expr={maxlike::maxlike(Fmula,
                                                        points=SpDataTM[SpDataTM$Partition==1 & SpDataTM$PresAbse==1,2:3],
                                                        rasters=raster::stack((VariablesT-colMeans(na.omit(raster::values(VariablesT))))/apply(na.omit(raster::values(VariablesT)),2,stats::sd)),
                                                        link=c("cloglog"),
                                                        hessian = FALSE,savedata=TRUE,
                                                        method="BFGS",
                                                        removeDuplicates=FALSE)},
                                error=function(e) {
                                  message("Trying Nelder-Mead Optimization")
                                  maxlike::maxlike(Fmula,
                                                   points=SpDataTM[SpDataTM$Partition==1 & SpDataTM$PresAbse==1,2:3],
                                                   rasters=raster::stack((VariablesT-colMeans(na.omit(raster::values(VariablesT))))/apply(na.omit(raster::values(VariablesT)),2,stats::sd)),
                                                   link=c("cloglog"),
                                                   hessian = FALSE,savedata=TRUE,
                                                   method="Nelder-Mead",
                                                   removeDuplicates=FALSE)})

              FinalModelT <- predict(Model)
              FinalModel <- STANDAR(FinalModelT)
              PredPoint <- raster::extract(FinalModel,SpDataTM[SpDataTM$Partition==1, 2:3])
              PredPoint <- data.frame(PresAbse = SpDataTM[SpDataTM$Partition==1, "PresAbse"], PredPoint)
            }else{
              Model <- tryCatch (expr={maxlike::maxlike(Fmula,
                                                        points=SpDataTM[SpDataTM[,"PresAbse"]==1,2:3],
                                                        rasters=raster::stack((VariablesT-colMeans(na.omit(raster::values(VariablesT))))/apply(na.omit(raster::values(VariablesT)),2,stats::sd)),
                                                        link=c("cloglog"),
                                                        hessian = FALSE,savedata=TRUE,
                                                        method="BFGS",
                                                        removeDuplicates=FALSE)},
                                 error=function(e) {
                                   message("Trying Nelder-Mead Optimization")
                                   maxlike::maxlike(Fmula,
                                                    points=SpDataTM[SpDataTM[,"PresAbse"]==1,2:3],
                                                    rasters=raster::stack((VariablesT-colMeans(na.omit(raster::values(VariablesT))))/apply(na.omit(raster::values(VariablesT)),2,stats::sd)),
                                                    link=c("cloglog"),
                                                    hessian = FALSE,savedata=TRUE,
                                                    method="Nelder-Mead",
                                                    removeDuplicates=FALSE)})
              
              # Model <- maxlike::maxlike(Fmula,points=SpDataTM[SpDataTM[,"PresAbse"]==1,2:3],rasters=raster::stack((VariablesT-colMeans(na.omit(values(VariablesT))))/apply(na.omit(values(VariablesT)),2,stats::sd)),
              #                  link=c("cloglog"),hessian = FALSE,savedata=TRUE,
              #                  method="BFGS",removeDuplicates=FALSE)
              FinalModelT <- predict(Model)
              FinalModel <- STANDAR(FinalModelT)
              PredPoint <- raster::extract(FinalModel,SpDataTM[, 2:3])
              PredPoint <- data.frame(PresAbse = SpDataTM[, "PresAbse"], PredPoint)
            }
            # x <- SpDataTM[SpDataTM[,"PresAbse"]==1, VarColT]
            # z <- SpDataTM[SpDataTM[,"PresAbse"]==0, VarColT]

            Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            #Final Thresholds
            Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

            #Variable Importance & Response Curves
            if(VarImp==TRUE){
              VarImp_RspCurv(Model=Model,Algorithm='MLK',folders=folders,spN=spN[s],SpDataT = SpDataTM,
                             VarColT=VarColT,Outcome=PredPoint$PredPoint)
            }

            #Final Model Rasters
            ListSummary[["MLK"]] <- data.frame(Sp=spN[s], Algorithm="MLK", Thr)

            if(SaveFinal=="Y"){
              ListRaster[["MLK"]] <- FinalModel
              names(ListRaster[["MLK"]]) <- spN[s]
            }
            if(is.null(Fut)==F){
              for(k in 1:length(VariablesP)){
                ListFut[[ProjN[k]]][["MLK"]] <- STANDAR_FUT(predict(raster::stack((VariablesP[[k]]-colMeans(na.omit(raster::values(VariablesT))))/apply(na.omit(raster::values(VariablesT)),2,stats::sd)), Model),FinalModelT)
              }
            }
          }
        }else{
          Eval <- list()
          Boyce <- list()
          Eval_JS <- list()
          pROC <- list()
          Area <- list()
          for(k in 1:length(VariablesP)){
            ListFut[[ProjN[k]]][["MLK"]] <- STANDAR(predict(raster::stack((VariablesP[[k]]-colMeans(na.omit(raster::values(VariablesT))))/apply(na.omit(raster::values(VariablesT)),2,stats::sd)),Model[[i]]))

            PredPoint <- raster::extract(ListFut[[ProjN[k]]][["MLK"]], PAtest[[i]][, c("x", "y")])
            PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
            Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                      PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                           a=PredPoint[PredPoint$PresAbse == 0, 2])
            
            #Thresholds and Final Evaluation
            Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
            Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
            Threshold <- Threshold[order(Thr)]
            Thr <- sort(Thr)
            Thr <- ifelse(Thr<0,0,Thr)
            Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["MLK"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
            
            #Percentae of Predicted Area
            ArT <- NULL
            for (j in Thr){
              RasL <- ListFut[[ProjN[k]]][["MLK"]]>=j
              ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
            }
            Area[[i]] <- round(ArT*100,3)
            
            #PartialROC
            pROC[[i]] <- partial_roc(prediction=ListFut[[ProjN[k]]][["MLK"]],
                                                test_data=PredPoint[PredPoint$PresAbse==1,2],
                                                error=5,iterations=500,
                                                percentage=50)$pROC_summary
            
            
            #MLK Validation
            pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
            pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
            pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
            pROC <- mean(unlist(lapply(pROC, `[`, 1)))
            BoyceSD <- stats::sd(unlist(Boyce))
            Boyce <- mean(unlist(Boyce))
            AreaSD <- apply(do.call("rbind",Area),2,sd)
            Area <- colMeans(do.call("rbind",Area))
            Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
            if(is.null(repl)){
              ListValidation[["MLK"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="MLK", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
            }else{
              ListValidation[["MLK"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="MLK", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
            }
          }
        }
      }
    }

    #SUPPORT VECTOR MACHINE (SVM)-----
    if (any(Algorithm == "SVM")) {
      Model <- list()
      Fmula <- stats::formula(paste("PresAbse", '~ .'))

      #SVM model
      for (i in 1:N) {
        dataPr <- PAtrain[[i]][, c("PresAbse", VarColT)]
        set.seed(0)
        Model[[i]] <- kernlab::ksvm(Fmula,data = dataPr,type="C-svc",kernel = "rbfdot",C = 1, prob.model=T)
      }

      #SVM evaluation
      if((is.null(Fut)==F && !is.null(Tst))==F){
        Eval <- list()
        Eval_JS <- list()
        Boyce <- list()
        pROC <- list()
        Area <- list()
        for (i in 1:N) {
          RastPart[["SVM"]][[i]] <- as.numeric(kernlab::predict(object=Model[[i]], newdata=PAtest[[i]][, VarColT],type="probabilities")[,2])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["SVM"]][[i]])
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["SVM"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          #Percentae of Predicted Area
          FinalModel <- data.frame(kernlab::predict(object=Model[[i]],newdata=raster::rasterToPoints(VariablesT)[,-c(1,2)],type="probabilities"))[,2]
          RasT[["SVM"]] <- VariablesT[[1]]
          RasT[["SVM"]][!is.na(RasT[["SVM"]][])] <- FinalModel
          ArT <- NULL
          for (j in Thr){
            RasL <- RasT[["SVM"]]>=j
            ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=RasT[["SVM"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,percentage=50)$pROC_summary
          
          #Save Partition Predictions
          if(Save=="Y"){
            if(N!=1){
              raster::writeRaster(RasT[["SVM"]],paste(grep("SVM",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="SVM",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["SVM"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              raster::writeRaster(RasT[["SVM"]],paste(grep("SVM",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),format='GTiff',overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="SVM",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["SVM"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }else{
              raster::writeRaster(RasT[["SVM"]],paste(grep("SVM",foldPart,value=T),"/",spN[s],".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="SVM",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["SVM"]]>=Thr_Alg[t],
                                    paste(grep("SVM",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
          }
        }
        
        #BIO Validation
        pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
        pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
        pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
        pROC <- mean(unlist(lapply(pROC, `[`, 1)))
        BoyceSD <- stats::sd(unlist(Boyce))
        Boyce <- mean(unlist(Boyce))
        AreaSD <- apply(do.call("rbind",Area),2,sd)
        Area <- colMeans(do.call("rbind",Area))
        Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
        if(is.null(repl)){
          ListValidation[["SVM"]] <- data.frame(Sp=spN[s], Algorithm="SVM",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }else{
          ListValidation[["SVM"]] <- data.frame(Sp=spN[s],Replicate=repl,Algorithm="SVM",Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }

        # Save final model
        if(repl==1 || is.null(repl)){
          if(is.null(repl) && N==1){
            Model <- kernlab::ksvm(Fmula,data = SpDataT[SpDataT$Partition==1, c("PresAbse", VarColT)],type="C-svc",
                          kernel = "rbfdot",C = 1, prob.model=T)
            FinalModel <- data.frame(kernlab::predict(object=Model,newdata=raster::rasterToPoints(VariablesT)[,-c(1,2)],type="probabilities"))[,2]
            FinalGrid <- Variables[[1]]
            FinalGrid[!is.na(FinalGrid[])] <- FinalModel
            # FinalModel <- data.frame(cbind(rasterToPoints(VariablesT)[,1:2],FinalModel))
            # sp::gridded(FinalModel) <- ~ x+y
            # FinalModel <- (raster(FinalModel))
            FinalModelT <- FinalGrid
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel,SpDataT[SpDataT$Partition==1, 2:3])
            PredPoint <- data.frame(PresAbse = SpDataT[SpDataT$Partition==1, "PresAbse"], PredPoint)
          }else{
            Model <- kernlab::ksvm(Fmula,data = SpDataT[, c("PresAbse", VarColT)],type="C-svc",
                          kernel = "rbfdot",C = 1, prob.model=T)
            FinalModel <- data.frame(kernlab::predict(object=Model,newdata=raster::rasterToPoints(VariablesT)[,-c(1,2)],type="probabilities"))[,2]
            FinalGrid <- Variables[[1]]
            FinalGrid[!is.na(FinalGrid[])] <- FinalModel
            # FinalModel <- data.frame(cbind(rasterToPoints(VariablesT)[,1:2],FinalModel))
            # sp::gridded(FinalModel) <- ~ x+y
            FinalModelT <- FinalGrid
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel,SpDataT[, 2:3])
            PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          }
          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

          #Variable Importance & Response Curves
          if(VarImp==TRUE){
            VarImp_RspCurv(Model=Model,Algorithm='SVM',folders=folders,spN=spN[s],SpDataT = SpDataT,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }

          #Final Model Rasters
          ListSummary[["SVM"]] <- data.frame(Sp=spN[s], Algorithm="SVM", Thr)

          if(SaveFinal=="Y"){
            ListRaster[["SVM"]] <- FinalModel
            names(ListRaster[["SVM"]]) <- spN[s]
          }
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              FutureModel <- data.frame(kernlab::predict(object=Model,newdata=raster::rasterToPoints(VariablesP[[k]])[,-c(1,2)],type="probabilities"))[,2]
              RasF <- VariablesT[[1]]
              RasF[!is.na(RasF[])] <- FutureModel
              ListFut[[ProjN[k]]][["SVM"]] <- STANDAR_FUT(RasF,FinalModelT)
            }
          }
        }
      }else{
        Eval <- list()
        Boyce <- list()
        for(k in 1:length(VariablesP)){
          ProjectedModel <- data.frame(kernlab::predict(object=Model[[i]],newdata=raster::rasterToPoints(VariablesP[[k]])[,-c(1,2)],type="probabilities"))[,2]
          RasF <- VariablesT[[1]]
          RasF[!is.na(RasF[])] <- ProjectedModel
          ListFut[[ProjN[k]]][["SVM"]] <- RasF

          PredPoint <- raster::extract(ListFut[[ProjN[k]]][["SVM"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["SVM"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          #Percentae of Predicted Area
          ArT <- NULL
          for (j in Thr){
            RasL <- ListFut[[ProjN[k]]][["SVM"]]>=j
            ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=ListFut[[ProjN[k]]][["SVM"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,
                                              percentage=50)$pROC_summary


          #SVM Validation
          pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
          pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
          pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
          pROC <- mean(unlist(lapply(pROC, `[`, 1)))
          BoyceSD <- stats::sd(unlist(Boyce))
          Boyce <- mean(unlist(Boyce))
          AreaSD <- apply(do.call("rbind",Area),2,sd)
          Area <- colMeans(do.call("rbind",Area))
          Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
          if(is.null(repl)){
            ListValidation[["SVM"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="SVM", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }else{
            ListValidation[["SVM"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="SVM", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
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
        # Model[[i]] <- randomForest(as.factor(PresAbse)~.,data=dataPr[,c("PresAbse",VarColT)],
        #                             importance=T, type="regression")
        Model[[i]] <- randomForest::tuneRF(dataPr[,-1], (dataPr[,1]), trace=F,
                             stepFactor=2, ntreeTry=1000, doBest=T, plot=F)
      }
      #RDF evaluation
      if((is.null(Fut)==F && !is.null(Tst))==F){
        Eval <- list()
        Eval_JS <- list()
        Boyce <- list()
        pROC <- list()
        Area <- list()
        for (i in 1:N) {
          RastPart[["RDF"]][[i]] <- as.vector(predict(Model[[i]], PAtest[[i]][, VarColT]))
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["RDF"]][[i]])
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(RastPart[["RDF"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          #Percentage of Predicted Area
          RasT[["RDF"]] <- raster::predict(VariablesT,Model[[i]])
          ArT <- NULL
          for (j in Thr){
            RasL <- RasT[["RDF"]]>=j
            ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=RasT[["RDF"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,percentage=50)$pROC_summary
          
          #Save Partition Predictions
          if(Save=="Y"){
            if(N!=1){
              raster::writeRaster(RasT[["RDF"]],paste(grep("RDF",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="RDF",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["RDF"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
            if(is.null(repl)==F){
              raster::writeRaster(RasT[["RDF"]],paste(grep("RDF",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),format='GTiff',overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="RDF",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["RDF"]]>=Thr_Alg[t],
                                    paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }else{
              raster::writeRaster(RasT[["RDF"]],paste(grep("RDF",foldPart,value=T),"/",spN[s],".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
              Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
              foldCatAlg <- grep(pattern="RDF",x=PartCat,value=T)
              for(t in 1:length(Thr_Alg)){
                raster::writeRaster(RasT[["RDF"]]>=Thr_Alg[t],
                                    paste(grep("RDF",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
              }
            }
          }
        }

        #RDF Validation
        pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
        pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
        pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
        pROC <- mean(unlist(lapply(pROC, `[`, 1)))
        BoyceSD <- stats::sd(unlist(Boyce))
        Boyce <- mean(unlist(Boyce))
        AreaSD <- apply(do.call("rbind",Area),2,sd)
        Area <- colMeans(do.call("rbind",Area))
        Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
        if(is.null(repl)){
          ListValidation[["RDF"]] <- data.frame(Sp=spN[s], Algorithm="RDF",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }else{
          ListValidation[["RDF"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="RDF",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }

        # Save final model
        if(repl==1 || is.null(repl)){
          set.seed(0)
          if(is.null(repl) && N==1){
            # Model <- randomForest(as.factor(PresAbse)~.,data=SpDataT[SpDataT$Partition==1,c("PresAbse",VarColT)],
            #                       importance=T, type="classification")
            # FinalModelT <- 1-predict(VariablesT,Model,type="prob")
            Model <- randomForest::tuneRF(SpDataT[SpDataT$Partition==1,VarColT], (SpDataT[SpDataT$Partition==1,"PresAbse"]), trace=F,
                            stepFactor=2, ntreeTry=500, doBest=T, plot = F)
            FinalModelT <- raster::predict(VariablesT,Model)
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel,SpDataT[SpDataT$Partition==1, 2:3])
            PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          }else{
            # Model <- randomForest(as.factor(PresAbse)~.,data=SpDataT[,c("PresAbse",VarColT)],
            #                       importance=T, type="classification")
            # FinalModelT <- 1-predict(VariablesT,Model,type="prob")
            Model <- randomForest::tuneRF(SpDataT[,VarColT], (SpDataT[,"PresAbse"]), trace=F,
                            stepFactor=2, ntreeTry=500, doBest=T, plot = F)
            FinalModelT <- raster::predict(VariablesT,Model)
            FinalModel <- STANDAR(FinalModelT)
            PredPoint <- raster::extract(FinalModel,SpDataT[, 2:3])
            PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
          }

          Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2])
          #Final Thresholds
          Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

          #Variable Importance & Response Curves
          if(VarImp==TRUE){
            VarImp_RspCurv(Model=Model,Algorithm='RDF',folders=folders,spN=spN[s],SpDataT = SpDataT,
                           VarColT=VarColT,Outcome=PredPoint$PredPoint)
          }

          #Final Model Rasters
          ListSummary[["RDF"]] <- data.frame(Sp=spN[s], Algorithm="RDF", Thr)

          if(SaveFinal=="Y"){
            ListRaster[["RDF"]] <- FinalModel
            names(ListRaster[["RDF"]]) <- spN[s]
          }
          if(is.null(Fut)==F){
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["RDF"]] <- STANDAR_FUT(raster::predict(VariablesP[[k]], Model),FinalModelT)
              if(maxValue(ListFut[[ProjN[k]]][["RDF"]])==0){
                ListFut[[ProjN[k]]][["RDF"]] <- ListFut[[ProjN[k]]][["RDF"]]
              }else{
                ListFut[[ProjN[k]]][["RDF"]] <- (ListFut[[ProjN[k]]][["RDF"]])
              }
            }
          }
        }
      }else{
        Eval <- list()
        Eval_JS <- list()
        Boyce <- list()
        pROC <- list()
        Area <- list()
        for(k in 1:length(VariablesP)){
          ListFut[[ProjN[k]]][["RDF"]] <- STANDAR(raster::predict(VariablesP[[k]],Model[[i]]))
          if(maxValue(ListFut[[ProjN[k]]][["RDF"]])==0){
            ListFut[[ProjN[k]]][["RDF"]] <- ListFut[[ProjN[k]]][["RDF"]]
          }else{
            ListFut[[ProjN[k]]][["RDF"]] <- (ListFut[[ProjN[k]]][["RDF"]])
          }

          PredPoint <- raster::extract(ListFut[[ProjN[k]]][["RDF"]], PAtest[[i]][, c("x", "y")])
          PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
          Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                    PredPoint[PredPoint$PresAbse == 0, 2])
          Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                         a=PredPoint[PredPoint$PresAbse == 0, 2])
          
          #Thresholds and Final Evaluation
          Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
          Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
          Threshold <- Threshold[order(Thr)]
          Thr <- sort(Thr)
          Thr <- ifelse(Thr<0,0,Thr)
          Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                       PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
          Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                            a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
          #Boyce Index
          Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["RDF"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
          
          #Percentae of Predicted Area
          ArT <- NULL
          for (j in Thr){
            RasL <- ListFut[[ProjN[k]]][["RDF"]]>=j
            ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
          }
          Area[[i]] <- round(ArT*100,3)
          
          #PartialROC
          pROC[[i]] <- partial_roc(prediction=ListFut[[ProjN[k]]][["RDF"]],
                                              test_data=PredPoint[PredPoint$PresAbse==1,2],
                                              error=5,iterations=500,
                                              percentage=50)$pROC_summary


          #RDF Validation
          pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
          pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
          pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
          pROC <- mean(unlist(lapply(pROC, `[`, 1)))
          BoyceSD <- stats::sd(unlist(Boyce))
          Boyce <- mean(unlist(Boyce))
          AreaSD <- apply(do.call("rbind",Area),2,sd)
          Area <- colMeans(do.call("rbind",Area))
          Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
          if(is.null(repl)){
            ListValidation[["RDF"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="RDF", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }else{
            ListValidation[["RDF"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="RDF", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
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
        RastPart[["GAM"]] <- NULL
        RasT[["GAM"]] <- NULL
      }else{
        Model <- list()
        Fmula <- paste("s(", VarColT,",k=3)", sep="")
        Fmula <- paste("PresAbse", paste(Fmula, collapse = " + "), sep = " ~ ")
        Fmula <- stats::as.formula(Fmula)

        #GAM model
        for (i in 1:N) {
          dataPr <- PAtrain[[i]][, c("PresAbse", VarColT)]
          Model[[i]] <- mgcv::gam(Fmula, data = dataPr, optimizer = c("outer", "newton"),
                            select = T, family = binomial)
          }
        #GAM evaluation
        if((is.null(Fut)==F && !is.null(Tst))==F){
          Eval <- list()
          Eval_JS <- list()
          Boyce <- list()
          pROC <- list()
          Area <- list()
          for (i in 1:N) {
            RastPart[["GAM"]][[i]] <- as.vector(predict(Model[[i]], PAtest[[i]][, VarColT],type="response"))
            PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["GAM"]][[i]])
            Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                      PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                           a=PredPoint[PredPoint$PresAbse == 0, 2])
            
            #Thresholds and Final Evaluation
            Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
            Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
            Threshold <- Threshold[order(Thr)]
            Thr <- sort(Thr)
            Thr <- ifelse(Thr<0,0,Thr)
            Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
            Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
            Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(RastPart[["GAM"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
            
            #Percentage of Predicted Area
            RasT[["GAM"]] <- raster::predict(VariablesT,Model[[i]],type="response")
            ArT <- NULL
            for (j in Thr){
              RasL <- RasT[["GAM"]]>=j
              ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
            }
            Area[[i]] <- round(ArT*100,3)
            
            #PartialROC
            pROC[[i]] <- partial_roc(prediction=RasT[["GAM"]],
                                                test_data=PredPoint[PredPoint$PresAbse==1,2],
                                                error=5,iterations=500,percentage=50)$pROC_summary
            
            #Save Partition Predictions
            if(Save=="Y"){
              if(N!=1){
                raster::writeRaster(RasT[["GAM"]],paste(grep("GAM",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
                Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="GAM",x=PartCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  raster::writeRaster(RasT[["GAM"]]>=Thr_Alg[t],
                                      paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                      format='GTiff',
                                      overwrite=TRUE)
                }
              }
              if(is.null(repl)==F){
                raster::writeRaster(RasT[["GAM"]],paste(grep("GAM",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),format='GTiff',overwrite=TRUE)
                Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="GAM",x=PartCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  raster::writeRaster(RasT[["GAM"]]>=Thr_Alg[t],
                                      paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                      format='GTiff',
                                      overwrite=TRUE)
                }
              }else{
                raster::writeRaster(RasT[["GAM"]],paste(grep("GAM",foldPart,value=T),"/",spN[s],".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
                Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="GAM",x=PartCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  raster::writeRaster(RasT[["GAM"]]>=Thr_Alg[t],
                                      paste(grep("GAM",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                                      format='GTiff',
                                      overwrite=TRUE)
                }
              }
            }
          }

          #GAM Validation
          pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
          pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
          pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
          pROC <- mean(unlist(lapply(pROC, `[`, 1)))
          BoyceSD <- stats::sd(unlist(Boyce))
          Boyce <- mean(unlist(Boyce))
          AreaSD <- apply(do.call("rbind",Area),2,sd)
          Area <- colMeans(do.call("rbind",Area))
          Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
          if(is.null(repl)){
            ListValidation[["GAM"]] <- data.frame(Sp=spN[s], Algorithm="GAM",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }else{
            ListValidation[["GAM"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="GAM",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }


          # Save final model
          if(repl==1 || is.null(repl)){
            if(is.null(repl) && N==1){
              Model <- mgcv::gam(formula=Fmula, data = SpDataT[SpDataT$Partition==1, c("PresAbse",VarColT)], optimizer = c("outer", "newton"),
                                 select = T, family = binomial)
              FinalModelT <- raster::predict(VariablesT,Model,type="response")
              FinalModel <- STANDAR(FinalModelT)
              PredPoint <- raster::extract(FinalModel,SpDataT[SpDataT$Partition==1, 2:3])
              PredPoint <- data.frame(PresAbse = SpDataT[SpDataT$Partition==1, "PresAbse"], PredPoint)
            }else{
              Model <- mgcv::gam(formula=Fmula, data = SpDataT[, c("PresAbse",VarColT)], optimizer = c("outer", "newton"),
                                 select = T, family = binomial)
              FinalModelT <- raster::predict(VariablesT,Model,type="response")
              FinalModel <- STANDAR(FinalModelT)
              PredPoint <- raster::extract(FinalModel,SpDataT[, 2:3])
              PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
            }
            #Final Model Threshold
            Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            #Final Thresholds
            Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

            #Variable Importance & Response Curves
            if(VarImp==TRUE){
              VarImp_RspCurv(Model=Model,Algorithm='GAM',folders=folders,spN=spN[s],SpDataT = SpDataT,
                             VarColT=VarColT,Outcome=PredPoint$PredPoint)
            }

            #Final Model Rasters
            ListSummary[["GAM"]] <- data.frame(Sp=spN[s], Algorithm="GAM", Thr)
            if(SaveFinal=="Y"){
              ListRaster[["GAM"]] <- FinalModel
              names(ListRaster[["GAM"]]) <- spN[s]
            }
            if(is.null(Fut)==F){
              for(k in 1:length(VariablesP)){
                ListFut[[ProjN[k]]][["GAM"]] <- STANDAR_FUT(raster::predict(VariablesP[[k]], Model,type="response"),FinalModelT)
              }
            }
          }
        }else{
          Eval <- list()
          Eval_JS <- list()
          Boyce <- list()
          pROC <- list()
          Area <- list()
          for(k in 1:length(VariablesP)){
            ListFut[[ProjN[k]]][["GAM"]] <- STANDAR(raster::predict(VariablesP[[k]],Model[[i]],type="response"))
            if(maxValue(ListFut[[ProjN[k]]][["GAM"]])==0){
              ListFut[[ProjN[k]]][["GAM"]] <- ListFut[[ProjN[k]]][["GAM"]]
            }else{
              ListFut[[ProjN[k]]][["GAM"]] <- (ListFut[[ProjN[k]]][["GAM"]])
            }

            PredPoint <- raster::extract(ListFut[[ProjN[k]]][["GAM"]], PAtest[[i]][, c("x", "y")])
            PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
            Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                      PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                           a=PredPoint[PredPoint$PresAbse == 0, 2])
            
            #Thresholds and Final Evaluation
            Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
            Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
            Threshold <- Threshold[order(Thr)]
            Thr <- sort(Thr)
            Thr <- ifelse(Thr<0,0,Thr)
            Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["GAM"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
            
            #Percentae of Predicted Area
            ArT <- NULL
            for (j in Thr){
              RasL <- ListFut[[ProjN[k]]][["GAM"]]>=j
              ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
            }
            Area[[i]] <- round(ArT*100,3)
            
            #PartialROC
            pROC[[i]] <- partial_roc(prediction=ListFut[[ProjN[k]]][["GAM"]],
                                                test_data=PredPoint[PredPoint$PresAbse==1,2],
                                                error=5,iterations=500,
                                                percentage=50)$pROC_summary


            #GAM Validation
            pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
            pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
            pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
            pROC <- mean(unlist(lapply(pROC, `[`, 1)))
            BoyceSD <- stats::sd(unlist(Boyce))
            Boyce <- mean(unlist(Boyce))
            AreaSD <- apply(do.call("rbind",Area),2,sd)
            Area <- colMeans(do.call("rbind",Area))
            Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
            if(is.null(repl)){
              ListValidation[["GAM"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="GAM", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
            }else{
              ListValidation[["GAM"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="GAM", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
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
        RastPart[["GLM"]] <- NULL
        RasT[["GLM"]] <- NULL
      }else{
        Model <- list()
        Fmula <- paste( "PresAbse ~ ", paste(c(VarColT, paste("I(",VarColT, "^2)", sep = "")),
                                     collapse = " + "), sep = "")
        Fmula <- stats::as.formula(Fmula)
        #GLM model
        for (i in 1:N) {
          dataPr <- PAtrain[[i]][, c("PresAbse", VarColT)]
          Model[[i]] <- stats::glm(Fmula, data = dataPr, family = binomial)
        }

        #GLM evaluation
        if((is.null(Fut)==F && !is.null(Tst))==F){
          Eval <- list()
          Eval_JS <- list()
          Boyce <- list()
          pROC <- list()
          Area <- list()
          for (i in 1:N) {
            RastPart[["GLM"]][[i]] <- as.vector(predict(Model[[i]], PAtest[[i]][, VarColT],type="response"))
            PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["GLM"]][[i]])
            Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                      PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                           a=PredPoint[PredPoint$PresAbse == 0, 2])
            
            #Thresholds and Final Evaluation
            Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
            Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
            Threshold <- Threshold[order(Thr)]
            Thr <- sort(Thr)
            Thr <- ifelse(Thr<0,0,Thr)
            Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
            Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
            Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(RastPart[["GLM"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
            
            #Percentage of Predicted Area
            RasT[["GLM"]] <- raster::predict(VariablesT,Model[[i]],type="response")
            ArT <- NULL
            for (j in Thr){
              RasL <- RasT[["GLM"]]>=j
              ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
            }
            Area[[i]] <- round(ArT*100,3)
            
            #PartialROC
            pROC[[i]] <- partial_roc(prediction=RasT[["GLM"]],
                                                test_data=PredPoint[PredPoint$PresAbse==1,2],
                                                error=5,iterations=500,percentage=50)$pROC_summary
            
            #Save Partition Predictions
            if(Save=="Y"){
              if(N!=1){
                raster::writeRaster(RasT[["GLM"]],paste(grep("GLM",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
                Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="GLM",x=PartCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  raster::writeRaster(RasT[["GLM"]]>=Thr_Alg[t],
                                      paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                      format='GTiff',
                                      overwrite=TRUE)
                }
              }
              if(is.null(repl)==F){
                raster::writeRaster(RasT[["GLM"]],paste(grep("GLM",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),format='GTiff',overwrite=TRUE)
                Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="GLM",x=PartCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  raster::writeRaster(RasT[["GLM"]]>=Thr_Alg[t],
                                      paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                      format='GTiff',
                                      overwrite=TRUE)
                }
              }else{
                raster::writeRaster(RasT[["GLM"]],paste(grep("GLM",foldPart,value=T),"/",spN[s],".tif", sep=""),
                                    format='GTiff',
                                    overwrite=TRUE)
                Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
                foldCatAlg <- grep(pattern="GLM",x=PartCat,value=T)
                for(t in 1:length(Thr_Alg)){
                  raster::writeRaster(RasT[["GLM"]]>=Thr_Alg[t],
                                      paste(grep("GLM",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                                      format='GTiff',
                                      overwrite=TRUE)
                }
              }
            }
          }

          #GLM Validation
          pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
          pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
          pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
          pROC <- mean(unlist(lapply(pROC, `[`, 1)))
          BoyceSD <- stats::sd(unlist(Boyce))
          Boyce <- mean(unlist(Boyce))
          AreaSD <- apply(do.call("rbind",Area),2,sd)
          Area <- colMeans(do.call("rbind",Area))
          Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
          if(is.null(repl)){
            ListValidation[["GLM"]] <- data.frame(Sp=spN[s], Algorithm="GLM",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }else{
            ListValidation[["GLM"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="GLM",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
          }

          # Save final model
          if(repl==1 || is.null(repl)){
            if(is.null(repl) && N==1){
              # Model <- stats::glm(Fmula, data = SpDataT[SpDataT$Partition==1, c("PresAbse",VarColT)], family = binomial(link = "logit"))
              Model <- stats::glm(Fmula, data = SpDataT[SpDataT$Partition==1, c("PresAbse",VarColT)], family =  binomial)
              FinalModelT <- raster::predict(VariablesT,Model,type="response")
              FinalModel <- STANDAR(FinalModelT)
              PredPoint <- raster::extract(FinalModel,SpDataT[SpDataT$Partition==1, 2:3])
              PredPoint <- data.frame(PresAbse = SpDataT[SpDataT$Partition==1, "PresAbse"], PredPoint)
            }else{
              Model <- stats::glm(Fmula, data = SpDataT[, c("PresAbse",VarColT)], family =  binomial)
              FinalModelT <- raster::predict(VariablesT,Model,type="response")
              FinalModel <- STANDAR(FinalModelT)
              PredPoint <- raster::extract(FinalModel,SpDataT[, 2:3])
              PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
            }
            Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2])
            #Final Thresholds
            Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

            #Variable Importance & Response Curves
            if(VarImp==TRUE){
              VarImp_RspCurv(Model=Model,Algorithm='GLM',folders=folders,spN=spN[s],SpDataT = SpDataT,
                             VarColT=VarColT,Outcome=PredPoint$PredPoint)
            }

            #Final Model Rasters
            ListSummary[["GLM"]] <- data.frame(Sp=spN[s], Algorithm="GLM", Thr)
            if(SaveFinal=="Y"){
              ListRaster[["GLM"]] <- FinalModel
              names(ListRaster[["GLM"]]) <- spN[s]
            }
            if(is.null(Fut)==F){
              for(k in 1:length(VariablesP)){
                ListFut[[ProjN[k]]][["GLM"]] <- STANDAR_FUT(raster::predict(VariablesP[[k]], Model,type="response"),FinalModelT)
              }
            }
          }
        }else{
          Eval <- list()
          Eval_JS <- list()
          Boyce <- list()
          pROC <- list()
          Area <- list()
          for(k in 1:length(VariablesP)){
            ListFut[[ProjN[k]]][["GLM"]] <- STANDAR(raster::predict(VariablesP[[k]],Model[[i]],type="response"))
            if(maxValue(ListFut[[ProjN[k]]][["GLM"]])==0){
              ListFut[[ProjN[k]]][["GLM"]] <- ListFut[[ProjN[k]]][["GLM"]]
            }else{
              ListFut[[ProjN[k]]][["GLM"]] <- (ListFut[[ProjN[k]]][["GLM"]])
            }

            PredPoint <- raster::extract(ListFut[[ProjN[k]]][["GLM"]], PAtest[[i]][, c("x", "y")])
            PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
            Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                      PredPoint[PredPoint$PresAbse == 0, 2])
            Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                           a=PredPoint[PredPoint$PresAbse == 0, 2])
            
            #Thresholds and Final Evaluation
            Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
            Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
            Threshold <- Threshold[order(Thr)]
            Thr <- sort(Thr)
            Thr <- ifelse(Thr<0,0,Thr)
            Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                         PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
            Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                              a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
            #Boyce Index
            Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["GLM"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
            
            #Percentae of Predicted Area
            ArT <- NULL
            for (j in Thr){
              RasL <- ListFut[[ProjN[k]]][["GLM"]]>=j
              ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
            }
            Area[[i]] <- round(ArT*100,3)
            
            #PartialROC
            pROC[[i]] <- partial_roc(prediction=ListFut[[ProjN[k]]][["GLM"]],
                                                test_data=PredPoint[PredPoint$PresAbse==1,2],
                                                error=5,iterations=500,
                                                percentage=50)$pROC_summary


            #GLM Validation
            pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
            pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
            pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
            pROC <- mean(unlist(lapply(pROC, `[`, 1)))
            BoyceSD <- stats::sd(unlist(Boyce))
            Boyce <- mean(unlist(Boyce))
            AreaSD <- apply(do.call("rbind",Area),2,sd)
            Area <- colMeans(do.call("rbind",Area))
            Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
            if(is.null(repl)){
              ListValidation[["GLM"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="GLM", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
            }else{
              ListValidation[["GLM"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="GLM", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
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
      Model[[i]] <- graf(y=dataPr[,"PresAbse"], x=dataPr[,VarColT],opt.l=F,method="Laplace")
    }

    #GAU evaluation
    if((is.null(Fut)==F && !is.null(Tst))==F){
      Eval <- list()
      Eval_JS <- list()
      Boyce <- list()
      pROC <- list()
      Area <- list()
      for (i in 1:N) {
        RastPart[["GAU"]][[i]] <- predict.graf(Model[[i]], PAtest[[i]][, VarColT],type="response",CI = NULL, maxn = NULL)
        RastPart[["GAU"]][[i]] <- as.vector(RastPart[["GAU"]][[i]][,"posterior mode"])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["GAU"]][[i]])
        Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                  PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                       a=PredPoint[PredPoint$PresAbse == 0, 2])
        
        #Thresholds and Final Evaluation
        Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
        Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
        Threshold <- Threshold[order(Thr)]
        Thr <- sort(Thr)
        Thr <- ifelse(Thr<0,0,Thr)
        Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
        Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
        #Boyce Index
        Boyce[[i]] <- ecospat.boyce(RastPart[["GAU"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        
        #Percentage of Predicted Area
        RasT[["GAU"]] <- predict.graf.raster(Model[[i]], VariablesT, type = "response",
                                    CI = NULL, maxn = NULL)$posterior.mode

        ArT <- NULL
        for (j in Thr){
          RasL <- RasT[["GAU"]]>=j
          ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
        }
        Area[[i]] <- round(ArT*100,3)
        
        #PartialROC
        pROC[[i]] <- partial_roc(prediction=RasT[["GAU"]],
                                            test_data=PredPoint[PredPoint$PresAbse==1,2],
                                            error=5,iterations=500,percentage=50)$pROC_summary
        
        #Save Partition Predictions
        if(Save=="Y"){
          if(N!=1){
            raster::writeRaster(RasT[["GAU"]],paste(grep("GAU",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                                format='GTiff',
                                overwrite=TRUE)
            Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern="GAU",x=PartCat,value=T)
            for(t in 1:length(Thr_Alg)){
              raster::writeRaster(RasT[["GAU"]]>=Thr_Alg[t],
                                  paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
            }
          }
          if(is.null(repl)==F){
            raster::writeRaster(RasT[["GAU"]],paste(grep("GAU",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),
                                format='GTiff',
                                overwrite=TRUE)
            Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern="GAU",x=PartCat,value=T)
            for(t in 1:length(Thr_Alg)){
              raster::writeRaster(RasT[["GAU"]]>=Thr_Alg[t],
                                  paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
            }
          }else{
            raster::writeRaster(RasT[["GAU"]],paste(grep("GAU",foldPart,value=T),"/",spN[s],".tif", sep=""),
                                format='GTiff',
                                overwrite=TRUE)
            Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern="GAU",x=PartCat,value=T)
            for(t in 1:length(Thr_Alg)){
              raster::writeRaster(RasT[["GAU"]]>=Thr_Alg[t],
                                  paste(grep("GAU",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                                  format='GTiff',
                                  overwrite=TRUE)
            }
          }
        }
      }

      #GAU Validation
      pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
      pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
      pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
      pROC <- mean(unlist(lapply(pROC, `[`, 1)))
      BoyceSD <- stats::sd(unlist(Boyce))
      Boyce <- mean(unlist(Boyce))
      AreaSD <- apply(do.call("rbind",Area),2,sd)
      Area <- colMeans(do.call("rbind",Area))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
      if(is.null(repl)){
        ListValidation[["GAU"]] <- data.frame(Sp=spN[s], Algorithm="GAU",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }else{
        ListValidation[["GAU"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="GAU",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }

      #Save final model
      if(repl==1 || is.null(repl)){
        if(is.null(repl) && N==1){
          Model <- graf(SpDataT[SpDataT$Partition==1,"PresAbse"], SpDataT[SpDataT$Partition==1,VarColT],opt.l=F,method="Laplace")
          FinalModelT <- predict.graf.raster(Model, VariablesT, type = "response",
                                             CI = NULL, maxn = NULL)$posterior.mode
          FinalModel <- STANDAR(FinalModelT)
          PredPoint <- raster::extract(FinalModel,SpDataT[SpDataT$Partition==1, 2:3])
          PredPoint <- data.frame(PresAbse = SpDataT[SpDataT$Partition==1, "PresAbse"], PredPoint)
        }else{
          Model <- graf(SpDataT[,"PresAbse"], SpDataT[,VarColT],opt.l=F,method="Laplace")
          FinalModelT <- predict.graf.raster(Model, VariablesT, type = "response",
                                             CI = NULL, maxn = NULL)$posterior.mode
          FinalModel <- STANDAR(FinalModelT)
          PredPoint <- raster::extract(FinalModel,SpDataT[, 2:3])
          PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
        }

        Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2])
        #Final Thresholds
        Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

        #Variable Importance & Response Curves
        if(VarImp==TRUE){
          VarImp_RspCurv(Model=Model,Algorithm='GAU',folders=folders,spN=spN[s],SpDataT = SpDataT,
                         VarColT=VarColT,Outcome=PredPoint$PredPoint)
        }

        #Final Model Rasters
        ListSummary[["GAU"]] <- data.frame(Sp=spN[s], Algorithm="GAU", Thr)

        if(SaveFinal=="Y"){
          ListRaster[["GAU"]] <- FinalModel
          names(ListRaster[["GAU"]]) <- spN[s]
        }
        if(is.null(Fut)==F){
          for(k in 1:length(VariablesP)){
            ListFut[[ProjN[k]]][["GAU"]] <- STANDAR_FUT(predict.graf.raster(Model, VariablesP[[k]], type = "response",
                                                                CI = NULL, maxn = NULL)$posterior.mode,FinalModelT)
          }
        }
      }
    }else{
      Eval <- list()
      Eval_JS <- list()
      Boyce <- list()
      pROC <- list()
      Area <- list()
      for(k in 1:length(VariablesP)){
        ListFut[[ProjN[k]]][["GAU"]] <- STANDAR(predict.graf.raster(Model[[i]], VariablesP[[k]], 
                                                                    type = "response",
                                                                    CI = NULL, maxn = NULL)$posterior.mode)
        if(maxValue(ListFut[[ProjN[k]]][["GAU"]])==0){
          ListFut[[ProjN[k]]][["GAU"]] <- ListFut[[ProjN[k]]][["GAU"]]
        }else{
          ListFut[[ProjN[k]]][["GAU"]] <- (ListFut[[ProjN[k]]][["GAU"]])
        }

        PredPoint <- raster::extract(ListFut[[ProjN[k]]][["GAU"]], PAtest[[i]][, c("x", "y")])
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
        Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                  PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                       a=PredPoint[PredPoint$PresAbse == 0, 2])
        
        #Thresholds and Final Evaluation
        Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
        Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
        Threshold <- Threshold[order(Thr)]
        Thr <- sort(Thr)
        Thr <- ifelse(Thr<0,0,Thr)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
        #Boyce Index
        Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["GAU"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        
        #Percentae of Predicted Area
        ArT <- NULL
        for (j in Thr){
          RasL <- ListFut[[ProjN[k]]][["GAU"]]>=j
          ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
        }
        Area[[i]] <- round(ArT*100,3)
        
        #PartialROC
        pROC[[i]] <- partial_roc(prediction=ListFut[[ProjN[k]]][["GAU"]],
                                            test_data=PredPoint[PredPoint$PresAbse==1,2],
                                            error=5,iterations=500,
                                            percentage=50)$pROC_summary


        #GAU Validation
        pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
        pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
        pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
        pROC <- mean(unlist(lapply(pROC, `[`, 1)))
        BoyceSD <- stats::sd(unlist(Boyce))
        Boyce <- mean(unlist(Boyce))
        AreaSD <- apply(do.call("rbind",Area),2,sd)
        Area <- colMeans(do.call("rbind",Area))
        Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
        if(is.null(repl)){
          ListValidation[["GAU"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="GAU", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }else{
          ListValidation[["GAU"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="GAU", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
        }
      }
    }
    }

    #BOOSTED REGRESSION TREE (BRT) ----
    if(any(Algorithm == 'BRT')){
      if(any(lapply(PAtrain, function(x) nrow(x[x$PresAbse==1,])*0.75)<=10)){
        ListValidation[["BRT"]] <- NULL
        ListRaster[["BRT"]] <- NULL
        ListSummary[["BRT"]] <- NULL
        RastPart[["BRT"]] <- NULL
        RasT[["BRT"]] <- NULL
      }else{
        Model <- list()
        #BRT model
        for (i in 1:N) {
          dataPr <- PAtrain[[i]]
          learn.rate <- 0.005
          ModelT <- NULL
          while(is.null(ModelT)){
            print(learn.rate)
            try(ModelT <- dismo::gbm.step(data=dataPr, gbm.x=VarColT, gbm.y="PresAbse", family = "bernoulli",
                                       tree.complexity= 5, learning.rate=learn.rate, bag.fraction= 0.75,silent=T,
                                       plot.main = F))
            learn.rate <- learn.rate-0.0005
            if(learn.rate<=0){
              ListValidation[["BRT"]] <- NULL
              ListRaster[["BRT"]] <- NULL
              ListSummary[["BRT"]] <- NULL
              RastPart[["BRT"]] <- NULL
              RasT[["BRT"]] <- NULL
              break
            }
          }
          if (is.null(ModelT)){
            Model[[i]] <- NULL
          }else{
            Model[[i]] <- ModelT
          }
        }
        
        Model <- Model[sapply(Model, function(x) !is.null(x))]

        #Check BRT Models
        if (length(Model)==N){

          #BRT evaluation
          if((is.null(Fut)==F && !is.null(Tst))==F){
            Eval <- list()
            Eval_JS <- list()
            Boyce <- list()
            pROC <- list()
            Area <- list()
            for (i in 1:N) {
              RastPart[["BRT"]][[i]] <- gbm::predict.gbm(Model[[i]], PAtest[[i]][, VarColT],
                                                    n.trees=Model[[i]]$gbm.call$best.trees,type="response")
              PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], RastPart[["BRT"]][[i]])
              Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                        PredPoint[PredPoint$PresAbse == 0, 2])
              Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                             a=PredPoint[PredPoint$PresAbse == 0, 2])
              
              #Thresholds and Final Evaluation
              Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
              Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
              Threshold <- Threshold[order(Thr)]
              Thr <- sort(Thr)
              Thr <- ifelse(Thr<0,0,Thr)
              Thr2 <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
              Thr2$THR_VALUE <- ifelse(Thr2$THR_VALUE<0,0,Thr2$THR_VALUE)
              Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                           PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
              Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                                a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
              #Boyce Index
              Boyce[[i]] <- ecospat.boyce(RastPart[["BRT"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
              
              #Percentage of Predicted Area
              RasT[["BRT"]] <- raster::predict(VariablesT,Model[[i]],
                              n.trees=Model[[i]]$gbm.call$best.trees,type="response")
              ArT <- NULL
              for (j in Thr){
                RasL <- RasT[["BRT"]]>=j
                ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
              }
              Area[[i]] <- round(ArT*100,3)
              
              #PartialROC
              pROC[[i]] <- partial_roc(prediction=RasT[["BRT"]],
                                                  test_data=PredPoint[PredPoint$PresAbse==1,2],
                                                  error=5,iterations=500,percentage=50)$pROC_summary
              
              #Save Partition Predictions
              if(Save=="Y"){
                if(N!=1){
                  raster::writeRaster(RasT[["BRT"]],paste(grep("BRT",foldPart,value=T),"/",spN[s],"_",i,".tif", sep=""),
                                      format='GTiff',
                                      overwrite=TRUE)
                  Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
                  foldCatAlg <- grep(pattern="BRT",x=PartCat,value=T)
                  for(t in 1:length(Thr_Alg)){
                    raster::writeRaster(RasT[["BRT"]]>=Thr_Alg[t],
                                        paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                        format='GTiff',
                                        overwrite=TRUE)
                  }
                }
                if(is.null(repl)==F){
                  raster::writeRaster(RasT[["BRT"]],paste(grep("BRT",foldPart,value=T),"/",spN[s],"_",repl,".tif", sep=""),format='GTiff',overwrite=TRUE)
                  Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
                  foldCatAlg <- grep(pattern="BRT",x=PartCat,value=T)
                  for(t in 1:length(Thr_Alg)){
                    raster::writeRaster(RasT[["BRT"]]>=Thr_Alg[t],
                                        paste(foldCatAlg[t], '/',spN[s],"_",i,".tif", sep=""),
                                        format='GTiff',
                                        overwrite=TRUE)
                  }
                }else{
                  raster::writeRaster(RasT[["BRT"]],paste(grep("BRT",foldPart,value=T),"/",spN[s],".tif", sep=""),
                                      format='GTiff',
                                      overwrite=TRUE)
                  Thr_Alg <- Thr2[Thr2$THR%in%Threshold,2]
                  foldCatAlg <- grep(pattern="BRT",x=PartCat,value=T)
                  for(t in 1:length(Thr_Alg)){
                    raster::writeRaster(RasT[["BRT"]]>=Thr_Alg[t],
                                        paste(grep("BRT",foldCatAlg[t],value=T), '/',spN[s],".tif", sep=""),
                                        format='GTiff',
                                        overwrite=TRUE)
                  }
                }
              }
            }

            #BRT Validation
            pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
            pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
            pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
            pROC <- mean(unlist(lapply(pROC, `[`, 1)))
            BoyceSD <- stats::sd(unlist(Boyce))
            Boyce <- mean(unlist(Boyce))
            AreaSD <- apply(do.call("rbind",Area),2,sd)
            Area <- colMeans(do.call("rbind",Area))
            Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
            if(is.null(repl)){
              ListValidation[["BRT"]] <- data.frame(Sp=spN[s], Algorithm="BRT",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
            }else{
              ListValidation[["BRT"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="BRT",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
            }

            #Save final model
            if(repl==1 || is.null(repl)){
              learn.rate <- 0.005
              Model <- NULL
              if(is.null(repl) && N==1){
                while(is.null(Model)){
                  try(Model <- dismo::gbm.step(data=SpDataT[SpDataT$Partition==1,], gbm.x=VarColT, gbm.y="PresAbse", family = "bernoulli",
                                        tree.complexity= 5, learning.rate=learn.rate, bag.fraction= 0.75,silent=T,
                                        plot.main = F))
                  learn.rate <- learn.rate-0.0005
                  if(learn.rate<=0){
                    ListValidation[["BRT"]] <- NULL
                    ListRaster[["BRT"]] <- NULL
                    ListSummary[["BRT"]] <- NULL
                    RastPart[["BRT"]] <- NULL
                    break
                  }
                }
                FinalModelT <- raster::predict(VariablesT,Model,
                                       n.trees=Model$gbm.call$best.trees,type="response")
                FinalModel <- STANDAR(FinalModelT)
                PredPoint <- raster::extract(FinalModel,SpDataT[SpDataT$Partition==1, 2:3])
                PredPoint <- data.frame(PresAbse = SpDataT[SpDataT$Partition==1, "PresAbse"], PredPoint)
              }else{
                while(is.null(Model)){
                  try(Model <- dismo::gbm.step(data=SpDataT, gbm.x=VarColT, gbm.y="PresAbse", family = "bernoulli",
                                        tree.complexity= 5, learning.rate=learn.rate, bag.fraction= 0.75,silent=T,
                                        plot.main = F))
                  learn.rate <- learn.rate-0.0005
                  if(learn.rate<=0){
                    ListValidation[["BRT"]] <- NULL
                    ListRaster[["BRT"]] <- NULL
                    ListSummary[["BRT"]] <- NULL
                    RastPart[["BRT"]] <- NULL
                    RasT[["BRT"]] <- NULL
                    break
                  }
                }
                FinalModelT <- raster::predict(VariablesT,Model,
                                       n.trees=Model$gbm.call$best.trees,type="response")
                FinalModel <- STANDAR(FinalModelT)
                PredPoint <- raster::extract(FinalModel,SpDataT[, 2:3])
                PredPoint <- data.frame(PresAbse = SpDataT[, "PresAbse"], PredPoint)
              }
              Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                           PredPoint[PredPoint$PresAbse == 0, 2])
              Eval_JS <- Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                                           PredPoint[PredPoint$PresAbse == 0, 2])
              #Final Thresholds
              Thr <- Thresholds_TMLA(Eval,Eval_JS,sensV)

              #Variable Importance & Response Curves
              if(VarImp==TRUE){
                VarImp_RspCurv(Model=Model,Algorithm='BRT',folders=folders,spN=spN[s],SpDataT = SpDataT,
                               VarColT=VarColT,Outcome=PredPoint$PredPoint)
              }

              #Final Model Rasters
              ListSummary[["BRT"]] <- data.frame(Sp=spN[s], Algorithm="BRT", Thr)

              if(SaveFinal=="Y"){
                ListRaster[["BRT"]] <- FinalModel
                names(ListRaster[["BRT"]]) <- spN[s]
              }
              if(is.null(Fut)==F){
                for(k in 1:length(VariablesP)){
                  ListFut[[ProjN[k]]][["BRT"]] <- STANDAR_FUT(raster::predict(VariablesP[[k]],Model,
                                                              n.trees=Model$gbm.call$best.trees,type="response"),FinalModelT)
                }
              }
            }
          }else{
            Eval <- list()
            Eval_JS <- list()
            Boyce <- list()
            pROC <- list()
            Area <- list()
            for(k in 1:length(VariablesP)){
              ListFut[[ProjN[k]]][["BRT"]] <- STANDAR(raster::predict(Model,VariablesP[[k]],
                                                          n.trees=Model$gbm.call$best.trees,type="response"))
              if(maxValue(ListFut[[ProjN[k]]][["BRT"]])==0){
                ListFut[[ProjN[k]]][["BRT"]] <- ListFut[[ProjN[k]]][["BRT"]]
              }else{
                ListFut[[ProjN[k]]][["BRT"]] <- (ListFut[[ProjN[k]]][["BRT"]])
              }

              PredPoint <- raster::extract(ListFut[[ProjN[k]]][["BRT"]], PAtest[[i]][, c("x", "y")])
              PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], PredPoint)
              Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                        PredPoint[PredPoint$PresAbse == 0, 2])
              Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                             a=PredPoint[PredPoint$PresAbse == 0, 2])
              
              #Thresholds and Final Evaluation
              Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
              Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
              Threshold <- Threshold[order(Thr)]
              Thr <- sort(Thr)
              Thr <- ifelse(Thr<0,0,Thr)
              Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                           PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
              Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                                a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
              #Boyce Index
              Boyce[[i]] <- ecospat.boyce(ListFut[[ProjN[k]]][["BRT"]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
              
              #Percentage of Predicted Area
              ArT <- NULL
              for (j in Thr){
                RasL <- ListFut[[ProjN[k]]][["BRT"]]>=j
                ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
              }
              Area[[i]] <- round(ArT*100,3)
              
              #PartialROC
              pROC[[i]] <- partial_roc(prediction=ListFut[[ProjN[k]]][["BRT"]],
                                                  test_data=PredPoint[PredPoint$PresAbse==1,2],
                                                  error=5,iterations=500,
                                                  percentage=50)$pROC_summary


              #BRT Validation
              pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
              pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
              pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
              pROC <- mean(unlist(lapply(pROC, `[`, 1)))
              BoyceSD <- stats::sd(unlist(Boyce))
              Boyce <- mean(unlist(Boyce))
              AreaSD <- apply(do.call("rbind",Area),2,sd)
              Area <- colMeans(do.call("rbind",Area))
              Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
              if(is.null(repl)){
                ListValidation[["BRT"]] <- data.frame(Sp=spN[s],Projection=names(VariablesP)[k] ,Algorithm="BRT", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
              }else{
                ListValidation[["BRT"]] <- data.frame(Sp=spN[s],Replicate=repl,Projection=names(VariablesP)[k] ,Algorithm="BRT", Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
              }
            }
          }
        }
      }
    }

    #Final models----
    if(repl==1 || is.null(repl)){
      if(SaveFinal=="Y"){
        if((is.null(Fut)==F && !is.null(Tst))==F){
          Thr <- lapply(ListSummary, '[', c('THR','THR_VALUE'))
          for(i in 1:length(ListRaster)){
            foldAlg <- grep(pattern=names(Thr)[i],x=folders,value=T)
            raster::writeRaster(round(ListRaster[[i]], 4),
                        paste(foldAlg, '/',spN[s],".tif", sep=""),
                        format='GTiff',
                        overwrite=TRUE)
            Thr_Alg <- Thr[[i]][Thr[[i]]$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern=names(Thr)[i],x=foldCat,value=T)
            for(t in 1:length(Thr_Alg)){
              raster::writeRaster(ListRaster[[i]]>=Thr_Alg[t],
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
          ListFut[[p]] <- ListFut[[p]][unlist(lapply(ListFut[[p]],function(x) class(x)=="RasterLayer"))]
          for(o in 1:length(ListFut[[p]])){
            raster::writeRaster(ListFut[[p]][[o]],file.path(ModFut[p],names(Thr)[o],spN[s]),
                        format='GTiff',overwrite=TRUE)
            Thr_Alg <- Thr[[o]][Thr[[o]]$THR%in%Threshold,2]
            foldCatAlg <- grep(pattern=Algorithm[o],x=foldCat,value=T)
            for(t in 1:length(Thr_Alg)){
              raster::writeRaster(ListFut[[p]][[o]]>=Thr_Alg[t],
                          file.path(ModFut[p],names(Thr)[o],Threshold[t],paste0(spN[s],".tif")),
                          format='GTiff',
                          overwrite=TRUE)
            }
          }
        }
      }
    }

    #Ensemble----
    # Mean Ensemble----
    if(any(PredictType=="MEAN")){

      #Partial Models Ensemble
      Final <- do.call(Map, c(rbind,RastPart))
      Final <- lapply(Final, function (x) colMeans(x))

      # Threshold
      Eval <- list()
      Eval_JS <- list()
      Boyce <- list()
      pROC <- list()
      Area <- list()
      for(i in 1:N){
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], Final[[i]])
        Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                  PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                       a=PredPoint[PredPoint$PresAbse == 0, 2])
        
        #Thresholds and Final Evaluation
        Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
        Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
        Threshold <- Threshold[order(Thr)]
        Thr <- sort(Thr)
        Thr <- ifelse(Thr<0,0,Thr)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        
        #Percentage of Predicted Area
        RasT <- RasT[lapply(RasT,length)>1]
        ENST <- mean(stack(RasT))
        ArT <- NULL
        for (j in Thr){
          RasL <- ENST>=j
          ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
        }
        Area[[i]] <- round(ArT*100,3)
        
        #PartialROC
        pROC[[i]] <- partial_roc(prediction=ENST,
                                            test_data=PredPoint[PredPoint$PresAbse==1,2],
                                            error=5,iterations=500,percentage=50)$pROC_summary
      }

      #MEAN Validation
      pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
      pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
      pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
      pROC <- mean(unlist(lapply(pROC, `[`, 1)))
      BoyceSD <- stats::sd(unlist(Boyce))
      Boyce <- mean(unlist(Boyce))
      AreaSD <- apply(do.call("rbind",Area),2,sd)
      Area <- colMeans(do.call("rbind",Area))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)

      if(is.null(repl)){
        ListValidation[["MEAN"]] <- data.frame(Sp=spN[s], Algorithm="MEA", Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }else{
        ListValidation[["MEAN"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="MEA", Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }
    }

    # Weighted Mean Ensemble----
    if(any(PredictType=="W_MEAN")){
      ListValidationT <- plyr::ldply(ListValidation,data.frame,.id=NULL)
      ListValidationT <- ListValidationT[ListValidationT$Algorithm%in%Algorithm,]

      #Partial Models Ensemble
      Final <- do.call(Map, c(rbind,RastPart))
      ThResW <- unlist(ListValidationT[ensemble_metric])
      Final <- lapply(Final, function(x) sweep(x, 2, ThResW, '*'))
      Final <- lapply(Final, function (x) colMeans(x))

      # Threshold
      Eval <- list()
      Eval_JS <- list()
      Boyce <- list()
      Eval_JS <- list()
      pROC <- list()
      Area <- list()
      for(i in 1:N){
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], Final[[i]])
        Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                  PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                       a=PredPoint[PredPoint$PresAbse == 0, 2])
        
        #Thresholds and Final Evaluation
        Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
        Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
        Threshold <- Threshold[order(Thr)]
        Thr <- sort(Thr)
        Thr <- ifelse(Thr<0,0,Thr)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        
        #Percentage of Predicted Area
        ENST <- mean(stack(RasT)*ThResW)
        
        ArT <- NULL
        for (j in Thr){
          RasL <- ENST>=j
          ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
        }
        Area[[i]] <- round(ArT*100,3)
        
        #PartialROC
        pROC[[i]] <- partial_roc(prediction=ENST,
                                            test_data=PredPoint[PredPoint$PresAbse==1,2],
                                            error=5,iterations=500,percentage=50)$pROC_summary
      }

      #W_MEAN Validation
      pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
      pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
      pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
      pROC <- mean(unlist(lapply(pROC, `[`, 1)))
      BoyceSD <- stats::sd(unlist(Boyce))
      Boyce <- mean(unlist(Boyce))
      AreaSD <- apply(do.call("rbind",Area),2,sd)
      Area <- colMeans(do.call("rbind",Area))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
      if(is.null(repl)){
        ListValidation[["W_MEAN"]] <- data.frame(Sp=spN[s], Algorithm="WMEA",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }else{
        ListValidation[["W_MEAN"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="WMEA",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }
    }

    # Superior Ensemble----
    if(any(PredictType=='SUP')){
      ListValidationT <- plyr::ldply(ListValidation,data.frame,.id=NULL)
      ListValidationT <- ListValidationT[ListValidationT$Algorithm%in%Algorithm,]

      Best <- ListValidationT[which(unlist(ListValidationT[ensemble_metric])>=mean(unlist(ListValidationT[ensemble_metric]))),"Algorithm"]
      W <- names(ListRaster)%in%Best

      #Partial Models
      Final <- do.call(Map, c(rbind,RastPart[W]))
      Final <- lapply(Final, function (x) colMeans(x))
      

      # Threshold
      Eval <- list()
      Eval_JS <- list()
      Boyce <- list()
      Eval_JS <- list()
      pROC <- list()
      Area <- list()
      for(i in 1:N){
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], Final[[i]])
        Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                  PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                       a=PredPoint[PredPoint$PresAbse == 0, 2])
        
        #Thresholds and Final Evaluation
        Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
        Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
        Threshold <- Threshold[order(Thr)]
        Thr <- sort(Thr)
        Thr <- ifelse(Thr<0,0,Thr)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        
        #Percentage of Predicted Area
        RasT <- RasT[lapply(RasT,length)>1]
        ENST <- mean(stack(RasT[Best]))
        
        ArT <- NULL
        for (j in Thr){
          RasL <- ENST>=j
          ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
        }
        Area[[i]] <- round(ArT*100,3)
        
        #PartialROC
        pROC[[i]] <- partial_roc(prediction=ENST,
                                            test_data=PredPoint[PredPoint$PresAbse==1,2],
                                            error=5,iterations=500,percentage=50)$pROC_summary
      }

      #SUP Validation
      pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
      pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
      pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
      pROC <- mean(unlist(lapply(pROC, `[`, 1)))
      BoyceSD <- stats::sd(unlist(Boyce))
      Boyce <- mean(unlist(Boyce))
      AreaSD <- apply(do.call("rbind",Area),2,sd)
      Area <- colMeans(do.call("rbind",Area))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
      if(is.null(repl)){
        ListValidation[["SUP"]] <- data.frame(Sp=spN[s], Algorithm="SUP", Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }else{
        ListValidation[["SUP"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="SUP",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }
    }

    # With PCA ------
    if (any(PredictType == 'PCA')) {

      #Partial Models Ensemble
      if(any(lapply(RastPart, function(x) length(x))>1)){
        Final <- do.call(Map, c(cbind, RastPart))
        Final <- lapply(Final, function(x) as.numeric(stats::princomp(x)$scores[,1]))
        Final <- lapply(Final, function(x) (x-min(x))/(max(x)-min(x)))
      }else{
        Final <- do.call(cbind,lapply(RastPart, function(x) do.call(cbind,x)))
        Final <- as.numeric(stats::princomp(Final)$scores[,1])
        Final <- list((Final-min(Final))/(max(Final)-min(Final)))
      }

      # Threshold
      Eval <- list()
      Eval_JS <- list()
      Boyce <- list()
      Eval_JS <- list()
      pROC <- list()
      Area <- list()
      for(i in 1:N){
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], Final[[i]])
        Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                  PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                       a=PredPoint[PredPoint$PresAbse == 0, 2])
        
        #Thresholds and Final Evaluation
        Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
        Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
        Threshold <- Threshold[order(Thr)]
        Thr <- sort(Thr)
        Thr <- ifelse(Thr<0,0,Thr)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        
        #Percentage of Predicted Area
        RasT <- RasT[lapply(RasT,length)>1]
        ENST <- PCA_ENS_TMLA(brick(stack(RasT)))
        
        ArT <- NULL
        for (j in Thr){
          RasL <- ENST>=j
          ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
        }
        Area[[i]] <- round(ArT*100,3)
        
        #PartialROC
        pROC[[i]] <- partial_roc(prediction=ENST,
                                            test_data=PredPoint[PredPoint$PresAbse==1,2],
                                            error=5,iterations=500,percentage=50)$pROC_summary
      }

      #PCA Validation
      pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
      pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
      pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
      pROC <- mean(unlist(lapply(pROC, `[`, 1)))
      BoyceSD <- stats::sd(unlist(Boyce))
      Boyce <- mean(unlist(Boyce))
      AreaSD <- apply(do.call("rbind",Area),2,sd)
      Area <- colMeans(do.call("rbind",Area))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
      if(is.null(repl)){
        ListValidation[["PCA"]] <- data.frame(Sp=spN[s], Algorithm="PCA",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }else{
        ListValidation[["PCA"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="PCA", Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }
    }

    # With PCA over the Mean(Superior) Ensemble----
    if (any(PredictType == 'PCA_SUP')) {
      ListValidationT <- plyr::ldply(ListValidation,data.frame,.id=NULL)
      ListValidationT <- ListValidationT[ListValidationT$Algorithm%in%Algorithm,]
      Best <- ListValidationT[which(unlist(ListValidationT[ensemble_metric])>=mean(unlist(ListValidationT[ensemble_metric]))),"Algorithm"]
      W <- names(ListRaster)%in%Best

      #Partial Models
      if(any(lapply(RastPart, function(x) length(x))>1)){
        Final <- do.call(Map, c(cbind, RastPart[W]))
        Final <- lapply(Final, function(x) as.numeric(stats::princomp(x)$scores[,1]))
        Final <- lapply(Final, function(x) (x-min(x))/(max(x)-min(x)))
      }else{
        Final <- do.call(cbind,lapply(RastPart[W], function(x) do.call(cbind,x)))
        Final <- as.numeric(stats::princomp(Final)$scores[,1])
        Final <- list((Final-min(Final))/(max(Final)-min(Final)))
      }

      # Threshold
      Eval <- list()
      Eval_JS <- list()
      Boyce <- list()
      Eval_JS <- list()
      pROC <- list()
      Area <- list()
      for(i in 1:N){
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], Final[[i]])
        Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                  PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                       a=PredPoint[PredPoint$PresAbse == 0, 2])
        
        #Thresholds and Final Evaluation
        Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
        Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
        Threshold <- Threshold[order(Thr)]
        Thr <- sort(Thr)
        Thr <- ifelse(Thr<0,0,Thr)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        
        #Percentage of Predicted Area
        RasT <- RasT[lapply(RasT,length)>1]
        ENST <- PCA_ENS_TMLA(brick(stack(RasT[Best])))
        
        ArT <- NULL
        for (j in Thr){
          RasL <- ENST>=j
          ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
        }
        Area[[i]] <- round(ArT*100,3)
        
        #PartialROC
        pROC[[i]] <- partial_roc(prediction = ENST,
                                            test_data=PredPoint[PredPoint$PresAbse==1,2],
                                            error=5,iterations=500,percentage=50)$pROC_summary
      }

      #PCA_SUP Validation
      pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
      pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
      pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
      pROC <- mean(unlist(lapply(pROC, `[`, 1)))
      BoyceSD <- stats::sd(unlist(Boyce))
      Boyce <- mean(unlist(Boyce))
      AreaSD <- apply(do.call("rbind",Area),2,sd)
      Area <- colMeans(do.call("rbind",Area))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)
      if(is.null(repl)){
        ListValidation[["PCS"]] <- data.frame(Sp=spN[s], Algorithm="PCS",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }else{
        ListValidation[["PCS"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="PCS",Partition=Part, Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }
    }

    #With PCA over the threshold Ensemble----
    if (any(PredictType == 'PCA_THR')) {
      ListValidationT <- plyr::ldply(ListSummary,data.frame,.id=NULL)
      ListValidationT <- ListValidationT[ListValidationT$Algorithm%in%Algorithm,]

      #Partial Models
      Final <- do.call(Map, c(cbind, RastPart))
      ValidTHR <- ListValidationT[grepl(Threshold,as.character(ListValidationT$THR),ignore.case = T),"THR_VALUE"]
      for (p in 1:length(Final)){
        Final[[p]] <- sapply(seq(1:length(ValidTHR)),function(x){ifelse(Final[[p]][,x]>=ValidTHR[x],Final[[p]][,x],0)})
        Final[[p]] <- as.numeric(stats::princomp(Final[[p]])$scores[,1])
        Final[[p]] <- (Final[[p]]-min(Final[[p]]))/(max(Final[[p]])-min(Final[[p]]))
      }

      # Evaluation
      Eval <- list()
      Eval_JS <- list()
      Boyce <- list()
      Eval_JS <- list()
      pROC <- list()
      Area <- list()
      for(i in 1:N){
        PredPoint <- data.frame(PresAbse = PAtest[[i]][, "PresAbse"], Final[[i]])
        Eval_T <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                  PredPoint[PredPoint$PresAbse == 0, 2])
        Eval_JS_T <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                       a=PredPoint[PredPoint$PresAbse == 0, 2])
        
        #Thresholds and Final Evaluation
        Thr <- Thresholds_TMLA(Eval_T,Eval_JS_T,sensV)
        Thr <- Thr[match(Threshold,Thr$THR),"THR_VALUE"]
        Threshold <- Threshold[order(Thr)]
        Thr <- sort(Thr)
        Thr <- ifelse(Thr<0,0,Thr)
        Eval[[i]] <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                                     PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)
        Eval_JS[[i]] <- Eval_Jac_Sor_TMLA(p=PredPoint[PredPoint$PresAbse == 1, 2],
                                          a=PredPoint[PredPoint$PresAbse == 0, 2],thr=Thr)
        Boyce[[i]] <- ecospat.boyce(Final[[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
        
        #Percentage of Predicted Area
        RasT <- RasT[lapply(RasT,length)>1]
        RasT <- sapply(seq(1:length(ValidTHR)),function(x){RasT[[x]]*(RasT[[x]]>=ValidTHR[x])})
        ENST <- PCA_ENS_TMLA(brick(stack(RasT)))
        
        ArT <- NULL
        for (j in Thr){
          RasL <- ENST>=j
          ArT <- c(ArT,sum(na.omit(raster::values(RasL)))/length(na.omit(raster::values(RasL))))
        }
        Area[[i]] <- round(ArT*100,3)
        
        #PartialROC
        pROC[[i]] <- partial_roc(prediction=ENST,
                                            test_data=PredPoint[PredPoint$PresAbse==1,2],
                                            error=5,iterations=500,percentage=50)$pROC_summary
      }

      #PCA_SUP Validation
      pvalROCSD <- stats::sd(unlist(lapply(pROC, `[`, 2)))
      pvalROC <- mean(unlist(lapply(pROC, `[`, 2)))
      pROCSD <- stats::sd(unlist(lapply(pROC, `[`, 1)))
      pROC <- mean(unlist(lapply(pROC, `[`, 1)))
      BoyceSD <- stats::sd(unlist(Boyce))
      Boyce <- mean(unlist(Boyce))
      AreaSD <- apply(do.call("rbind",Area),2,sd)
      Area <- colMeans(do.call("rbind",Area))
      Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)

      if(is.null(repl)){
        ListValidation[["PCT"]] <- data.frame(Sp=spN[s], Algorithm="PCT", Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }else{
        ListValidation[["PCT"]] <- data.frame(Sp=spN[s],Replicate=repl, Algorithm="PCT", Partition=Part,Validation,pROC=pROC,pROC_SD=pROCSD,p_value_pROC=pvalROC,p_value_pROC_SD=pvalROCSD,Percentage_predicted_area=Area,Percentage_predicted_area_SD=Area,Boyce=Boyce,Boyce_SD=BoyceSD)
      }
    }

    #Final Data Frame Results

    result <- plyr::ldply(ListValidation,data.frame,.id=NULL)
    resultII <- plyr::ldply(ListSummary,data.frame,.id=NULL)

    out <- list(Validation = result,
                Summary = resultII)
    return(out)
  }#Close Species Loop


# Save .txt with the models performance----
FinalValidation <- data.frame(data.table::rbindlist(do.call(rbind,lapply(results, "[", "Validation"))))
FinalSummary <- data.frame(data.table::rbindlist(do.call(rbind,lapply(results, "[", "Summary"))))

utils::write.table(FinalValidation,paste(DirSave, VALNAME, sep = '/'),sep="\t",
            col.names = T,row.names=F)
if(repl==1 || is.null(repl)){
  utils::write.table(FinalSummary,paste(DirSave, VALNAMEII, sep = '/'),sep="\t",
              col.names = T,row.names=F)
}

cat("Models fitted!\n")

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
                     paste("Resuls in:" , DirSave),
                     paste('No species:',length(spN)),
                     paste("Threshold:",Threshold),
                     matrix(spN))
  lapply(InfoModeling, write,
         paste(DirSave, "/InfoModeling.txt", sep=""), append=TRUE,
         ncolumns=20, sep='\t')
  parallel::stopCluster(cl)
}
