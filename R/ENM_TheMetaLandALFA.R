#' Create and process Ecological Niche and Species Distribution Models
#'
ENMs_TheMetaLand<-function(dir,
                           sp,
                           x,
                           y,
                           min_occ=10,
                           thin_occ,
                           colin_var,
                           imp_var,
                           transfer,
                           eval_occ,
                           sp_accessible_area,
                           pres_abs_ratio=1,
                           pseudoabs_method,
                           part,
                           save_part="N",
                           save_final="Y",
                           algorithm,
                           thr,
                           msdm,
                           ensemble,
                           s_sdm){
  
  
#1.Check Function Arguments  
  
  er <- NULL
  if(missing(dir)){
    er <- c(er,paste("'dir' unspecified argument, specify the directory of environmental variables | "))
  }
  if(missing(sp)){
    er <- c(er,paste("'sp' unspecified argument, specify the column' name with the species name  | "))
  }
  if(missing(x)){
    er <- c(er,paste("'x' unspecified argument, specify the column with the longitude values | "))
  }
  if(missing(y)){
    er <- c(er,paste("'y' unspecified argument, specify the column with the latitude values | "))
  }
  if(missing(colin_var)){
    er <- c(er,paste("'colin_var' unspecified argument, specify whether you want to perform PCA on environmental variables | "))
  }
  if(missing(transfer)){
    er <- c(er,paste("'transfer' unspecified argument, specify whether you want to project the model for another region/time period | "))
  }
  if(missing(pres_abs_ratio)){
    er <- c(er,paste("'pres_abs_ratio' unspecified argument, specify a prevalence between train/test | "))
  }
  if(missing(pseudoabs_method)){
    er <- c(er,paste("'pseudoabs_method' unspecified argument, specify the allocation method of PseudoAusencias | "))
  }
  if(missing(part)){
    er <- c(er,paste("'part' unspecified argument, specify the method of participation | "))
  }
  if(missing(algorithm)){
    er <- c(er,paste("'algorithm' unspecified argument, specify which algorithms you want to use | "))
  }
  if(missing(thr)){
    er <- c(er,paste("'thr' unspecified argument, specify which Threshold wants to use | "))
  }
  if(missing(msdm)){
    er <- c(er,paste("'msdm' unspecified argument, specify if you want to use spatial restrictions in your models | "))
  }
  if(missing(ensemble)){
    er <- c(er,paste("'ensemble' unspecified argument, specify whether you want to perform the evaluation of the models | "))
  }
  if(!is.null((er))){
    print(er)
    stop("Argumentos faltantes, please, check if the argumentos listed above")
  }
  
  if(!(colin_var%in%c("Pearson","VIF","PCA","N"))){
    stop("'colin_var' Argument is not valid!(Pearson,VIF,PCA,N)")
  }
  if(!(transfer%in%c("Y","N"))){
    stop("'transfer' Argument is not valid!(Y/N)")
  }
  if(pres_abs_ratio<=0){
    stop("'pres_abs_ratio' Argument is not valid!(pres_abs_ratio>=0)")
  }
  if(!(pseudoabs_method%in%c("Rnd","EnvConst","GeoConst"))){
    stop("'pseudoabs_method' Argument is not valid!(Rnd/EnvConst/GeoConst)")
  }
  if(length(pseudoabs_method)>1){
    stop("Please choose only one Pseudo-absence allocation method")
  }
  if(!(part%in%c("boot","cross","band","check"))){
    stop("'part' Argument is not valid!(boot/cross/band/check)")
  }
  if(any(!algorithm%in%c("BIO","GLM","GAM","SVM","RDF","MXS","MXD","MLK","GAU"))){
    stop(paste("Algorithm",algorithm[!(algorithm%in%c("BIO","GLM","GAM","SVM","RDF","MXS","MXD","MLK","GAU"))],"is not valid"))
  }
  if(any(!thr%in%c("no_omission","spec_sens","kappa","equal_sens_spec","prevalence","sensitivity"))){
    stop("'thr' Argument is not valid!")
  }
  if(!(msdm%in%c("N","XY","MIN","CML","KER","POST"))){
    stop("'msdm' Argument is not valid!(N/XY/MIN/CML/KER/POST)")
  }
  if(length(msdm)>1){
    stop("Please choose only one 'msdm' method")
  }
  
#1.Load Packages ----
  
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  
  
  ipak(c("raster","sp","dismo","kernlab","randomForest","rgdal","gam",
         "maxnet","maptools","maxlike","mgcv", "plyr", "GRaF",
         "RStoolbox","flexclust","ape","tools","modEvA","SDMTools","SpatialEpi",
         "rgeos", "foreach", "doParallel","data.table","devtools","spThin","geoR",
         "usdm","pracma","gbm","caret","adehabitatHS"))

  #1.1. Choose.dir correction for Linux and MAC
  if(Sys.info()['sysname']!="Windows"){
    choose.dir <- function() {
      system("osascript -e 'tell app \"R\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
             intern = FALSE, ignore.stderr = TRUE)
      p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
      return(ifelse(length(p), p, NA))
    }
  }
  
#2.Adjust Algorithm Names----
  Ord <- c("BIO","DOM","MAH","ENF","MXD","MXS","MLK","SVM","RDF","GAM","GLM","GAU","BRT")
  algorithm <- Ord[Ord%in%algorithm]
  
#3.Predictors ----
  options(warn=1)
  setwd(dir)
  
  env <- unique(file_ext(list.files()))
  form <- c('bil','asc','txt','tif')
  env <- env[env%in%form]
  if(length(env)>1){
    stop("More than one file format in dir")
  }
  
  if(any(env == c('asc', 'bil', 'tif'))){
    envT<-brick(stack(list.files(pattern=paste0('\\.',env,'$'))))
  }
  if(env == 'txt'){
    envT<-read.table(list.files(pattern='\\.txt$'),h=T)
    gridded(envT)<- ~x+y
    envT<-brick(stack(envT))
  }
  
  #3.0.Check predictors consistency
  if(length(unique(colSums(!is.na(envT[]))))>1){
    envT[is.na((sum(envT))[])] <- NA
    print("Variables had differences, setting any empty cells to NA in all variables")
  }
  
  #3.1.Projection----
  if(transfer=="Y"){
    print("Select folder containing GCM folders:")
    DirP<-choose.dir(getwd())
    Pfol<-file.path(DirP,list.files(DirP))
    if(any(file_ext(list.files(DirP))%in%form)){
      stop("Select a folder containing GCM folders, NOT a folder with GCM variables!")
    }
    
    PfolN <- list.files(DirP)
    
    #Check Present/Future Names Consistency
    FutN <- list()
    for(i in 1:length(Pfol)){
      ProjT <- unique(file_ext(list.files(Pfol[[i]])))
      form <- c('bil','asc','txt','tif')
      ProjT <- ProjT[ProjT%in%form]
      FutN[[i]] <- file_path_sans_ext(list.files(Pfol[[i]],pattern=ProjT))
    }
    if(any(unlist(lapply(FutN, function(x) (names(envT)!=x))))){
      stop("Present/Future Variables Do Not Match! Make sure Present/Future Variables have the same names")
    }
  }
  
  #3.1. Variable Colinearity----
  #3.1.1.VIF----
  if(colin_var=="VIF"){
    VF <- vifstep(envT,th=nlayers(envT)*2)
    envT <- exclude(envT,VF)
    if(transfer=="Y"){
      RasM <- colMeans(na.omit(values(envT)))
      RasSTD <- apply(na.omit(values(envT)),2,std)
    }
    envT <- scale(envT)
    
    if(transfer=="Y"){
      EnvF <- list()
      for(i in 1:length(Pfol)){
        ProjT <- unique(file_ext(list.files(Pfol[i])))
        form <- c('bil','asc','txt','tif')
        ProjT <- ProjT[ProjT%in%form]
        if(length(ProjT)>1){
          stop("More than one file format in DirP")
        }
        
        if(any(ProjT == c('asc', 'bil', 'tif'))){
          EnvF[[i]]<-brick(stack(file.path(Pfol[i],list.files(Pfol[i],pattern=paste0('\\.',ProjT,'$')))))
        }
        if(ProjT == 'txt'){
          ProjTT<-read.table(file.path(Pfol[i],list.files(Pfol[i],pattern='\\.txt$'),h=T))
          gridded(ProjTT)<- ~x+y
          EnvF[[i]]<-brick(stack(ProjTT))
          rm(ProjTT)
        }
        
        EnvF[[i]] <- EnvF[[names(envT)]]
        EnvF[[i]] <- (EnvF[[i]]-RasM)/RasSTD
      }
    }
  }
  
  #3.1.2.PCA----
  
  if (colin_var=="PCA"){
      
    #Projection PCA
    if(transfer=="Y"){
      EnvF <- list()
      for(i in 1:length(Pfol)){
        EnvF[[i]] <- PCAFuturo(Env=envT,Dir=dir,DirP=Pfol[i],Save="Y")
      }
      names(EnvF) <- PfolN
      envT <- brick(stack(file.path(dir,"PCA",list.files(file.path(dir,"PCA"),pattern='PC'))))
    }else{
      envT<-PCA_env_TMLA(env = envT, Dir = dir)
    }
  }
  
  #3.3.3.Pearson----
  if(colin_var=="Pearson"){
    cat("Select correlation threshold:(0-1)")
    Cor_TH <- as.numeric(readLines(n=1))
    Pear <- layerStats(envT, 'pearson', na.rm=T)
    corr_matrix <- abs(Pear$'pearson correlation coefficient')
    corr_matrix[upper.tri(corr_matrix)] <- 0
    diag(corr_matrix) <- 0
    envT <- envT[[names(envT)[!apply(corr_matrix,2,function(x) any(x > 0.70))]]]
    if(transfer=="Y"){
      RasM <- colMeans(na.omit(values(envT)))
      RasSTD <- apply(na.omit(values(envT)),2,std)
    }
    envT <- scale(envT)
    
    if(transfer=="Y"){
      EnvF <- list()
      for(i in 1:length(Pfol)){
        ProjT <- unique(file_ext(list.files(Pfol[i])))
        form <- c('bil','asc','txt','tif')
        ProjT <- ProjT[ProjT%in%form]
        if(length(ProjT)>1){
          stop("More than one file format in DirP")
        }
        
        if(any(ProjT == c('asc', 'bil', 'tif'))){
          EnvF[[i]]<-brick(stack(file.path(Pfol[i],list.files(Pfol[i],pattern=paste0('\\.',ProjT,'$')))))
        }
        if(ProjT == 'txt'){
          ProjTT<-read.table(file.path(Pfol[i],list.files(Pfol[i],pattern='\\.txt$'),h=T))
          gridded(ProjTT)<- ~x+y
          EnvF[[i]]<-brick(stack(ProjTT))
          rm(ProjTT)
        }
        
        EnvF[[i]] <- EnvF[[names(envT)]]
        EnvF[[i]] <- (EnvF[[i]]-RasM)/RasSTD
      }
    }
  }
  
  #3.3.Erro Futuro e msdm
  if(transfer=="Y" && msdm!="N"){
    warning("msdm can not be used with future projections")
    warning("Setting msdm to N")
    msdm <- "N"
  }
  
  #3.4.Aviso caso min_occ<NPreditores
  if(min_occ<nlayers(envT)){
    warning("The minimum number of occurrences is smaller than the number of predictors.
            This may cause some issues while fitting certain algorithms!")
  }
  
  
#4.Occurrence Data ----
  
  DirR<-"Result"
  setwd("..")
  if (file.exists(file.path(getwd(),DirR))){
    DirR<-file.path(getwd(), DirR)
  }else{
    dir.create(file.path(getwd(), DirR))
    DirR<-file.path(getwd(), DirR)
  }
  
  print("Select occurrence file(.txt):")
  occ<-read.table(file.choose(getwd()),h=T,sep='\t',stringsAsFactors = F)
  occ<-occ[,c(sp,x,y)]
  colnames(occ) <- c("sp","x","y")
  occ_xy <- split(occ[,-1],f=occ$sp)
  spN<-names(occ_xy)
  
  
    #4.1.Unique Occurrences----
    occA<-Occ_Unicas_TMLA(env=envT[[1]],occ.xy=occ_xy,DirO=DirR)

    #4.2.Thining----
    if(thin_occ=="Y"){
      cat(("Select thinning method:\n1-Distance defined by Moran Variogram\n2-User defined distance\n3-Distance defined by 10x cellsize (Haversine Transformation)"))
      ThinMet <- as.integer(readLines(n=1))
      occA <- OccsThin(occA,envT,ThinMet,colin_var,DirR)
    }
  
    #4.3.Remove species with less than min_occ----
    occ <- occA[sapply(occA,function (x) nrow(x)>=min_occ)]
    spN<-names(occ)
  
    
    #4.4.Species with few records----
    if(length(occ)!=length(occ_xy)){
      print(paste("Species with less than ",min_occ, " Unique Occurrences were removed!"))
      print(names(occ_xy)[names(occ_xy)%in%spN==F])
      ndb <- ldply(occ)[,1:3]
      write.table(ndb,file.path(DirR,"Occ_Filtered.txt"),sep="\t",row.names=F)
      rm(ndb)
    }
    occ_xy <- lapply(occ,function(x) x[,c("x","y")])
    
    #4.5.GAM and GLM usage----
    if(any(sapply(occ,function(x) nrow(x))<nlayers(envT)) && any(algorithm%in%c("GAM","GLM"))){
      warning("A species has fewer records than the number of predictors, impossible to adjust GAM and GLM! GAM and GLM will be excluded")
      algorithm <- algorithm[!algorithm%in%c("GAM","GLM")]
    }
    
    
#5. Restrict Extent per Species----
    if(sp_accessible_area=="Y"){
      cat("Select restriction type (buffer / ecoregions):")
      method <- as.character(readLines(n = 1))
      while(!method%in%c("buffer","ecoregions")){
        warning("Please Select a valid restriction type (buffer / ecoregions)")
        cat("Select restriction type (buffer / ecoregions):")
        method <- as.character(readLines(n = 1))
      }
      DirM <- M_delimited(var=envT,
                  occ_xy=occ_xy,
                  method = method,
                  BufferDistanceKm=NULL,
                  EcoregionsFile=NULL,
                  Dir=dir,
                  spN=spN,
                  SaveM = TRUE)
    }
  
    
#6. Geographical Partition----
    if(part=="band" || part=="check"){
      
      if(any(grepl("PC",names(envT)))==T || any(grepl("pc",names(envT)))==T){
        colin_var<-"PCA"
      }
      
      if(eval_occ=="Y"){
        warning("Invalid combination! eval_occ can't be Y with Geographical partitions! Changing eval_occ to N")
        eval_occ <- "N"
      }
      
      if(part=="band"){  
        #6.1.Bands----
        
        DirB<-"Bands"
        if (file.exists(file.path(DirR,DirB))){
          DirB<-file.path(DirR,DirB)
        } else {
          dir.create(file.path(DirR,DirB))
          DirB<-file.path(DirR,DirB)
        }
        if(all(paste0(spN,".tif")%in%list.files(DirB,pattern=".tif"))){
          warning("Partition Already Exist! Using pre-created partitions! ")
          setwd(DirB)
          occINPUT <- read.table(file.path(DirB,"OccBands.txt"),sep="\t",header=T)
          occINPUT[,4] <- as.numeric(occINPUT[,4])
          occINPUT[,5] <- as.numeric(occINPUT[,5])
        }else{
          if(colin_var!="PCA"){
            envTT<-PCA_env_TMLA(env = envT, Dir = dir)
          }else{
            envTT<-envT
          }
          print("Use Longitudinal(1) or Latitudinal Bands(2)?")
          bands <- as.integer(readLines(n = 1))
          while(is.na(bands)||!(bands%in%c(1,2))){
            warning("Please choose bands by its number [1 or 2]")
            print("Use Longitudinal(1) or Latitudinal Bands(2)?")
            bands <- as.integer(readLines(n = 1))
          }
          print("Select Moran Calculation Type (all/nearest):")
          TipoMoran <- as.character(readLines(n = 1))
          while(is.na(TipoMoran)||!(TipoMoran%in%c("all","nearest"))){
            warning("Please choose a valid Moran Calculation Type [all/nearest]")
            print("Select Moran Calculation Type (all/nearest):")
            TipoMoran <- as.character(readLines(n = 1))
          }
          
          #Check for M-Restriction
          if(exists("DirM")){
            DirM <- DirM
          }else{
            DirM <- NULL
          }
          occINPUT <- BandsPartition_TMLA(evnVariables=envTT,RecordsData=occ_xy,N=bands,
                                      pseudoabsencesMethod=pseudoabs_method,PrAbRatio=pres_abs_ratio,DirSave=DirB,
                                      DirM=DirM,MRst=sp_accessible_area,type=TipoMoran)
          occINPUT[,4] <- as.numeric(occINPUT[,4])
          occINPUT[,5] <- as.numeric(occINPUT[,5])
          rm(envTT)
        }

      }
      if(part=="check"){
        #6.2.Block----
        
        DirB<-"Blocks"
        if (file.exists(file.path(DirR,DirB))){
          DirB<-file.path(DirR,DirB)
        } else {
          dir.create(file.path(DirR,DirB))
          DirB<-file.path(DirR,DirB)
        }
        
        if(all(paste0(spN,".tif")%in%list.files(DirB,pattern=".tif"))){
          print("Partition Already Exist! Using pre-created partitions! ")
          setwd(DirB)
          occINPUT <- read.table(file.path(DirB,"OccBlocks.txt"),sep="\t",header=T)
          occINPUT[,4] <- as.numeric(occINPUT[,4])
          occINPUT[,5] <- as.numeric(occINPUT[,5])
        }else{
          if(colin_var!="PCA"){
            envTT<-PCA_env_TMLA(env = envT, Dir = dir)
          }else{
            envTT<-envT
          }
          print("Select Moran Calculation Type (all/nearest):")
          TipoMoran <- as.character(readLines(n = 1))
          while(is.na(TipoMoran)||!(TipoMoran%in%c("all","nearest"))){
            warning("Please choose a valid Moran Calculation Type [all/nearest]")
            print("Select Moran Calculation Type (all/nearest):")
            TipoMoran <- as.character(readLines(n = 1))
          }
          
          #Check for M-Restriction
          if(exists("DirM")){
            DirM <- DirM
          }else{
            DirM <- NULL
          }
          occINPUT <- BlockPartition_TMLA(evnVariables=envTT,RecordsData=occ_xy,N=2,
                                      pseudoabsencesMethod=pseudoabs_method,PrAbRatio=pres_abs_ratio,DirSave=DirB,
                                      DirM=DirM,MRst=sp_accessible_area,type=TipoMoran)

          occINPUT[,4] <- as.numeric(occINPUT[,4])
          occINPUT[,5] <- as.numeric(occINPUT[,5])
          rm(envTT)
        }
      }
      
      #6.3.msdm A PRIORI----
      if(msdm=="N"||msdm=="POST"){
        DirPRI <- NULL
      }
      
      if(msdm%in%c("XY","MIN","CML", "KER")){
        print("Creating msdm Layers...")
        
        DirMSDM<-"msdm"
        if (file.exists(file.path(dir,DirMSDM))){
          DirMSDM<-file.path(dir,DirMSDM)
        } else {
          dir.create(file.path(dir,DirMSDM))
          DirMSDM<-file.path(dir,DirMSDM)
        }
        DirPRI <- MSDM_Priori_TMLA(Species=occ_xy,var=envT,MSDM=msdm,DirMSDM=DirMSDM)
      }

      #6.4. Future Projections ----
      if(transfer=="Y"){
        Fut <- EnvF
      }else{
        Fut <- NULL
      }
      
      #6.5.Adjust Checkerboard when using Geographical Restrictions (For Maxent Sampling)
      if(sp_accessible_area=="Y"){
        Ms <- stack(file.path(DirM,list.files(DirM)))
        Cs <- stack(file.path(DirB,list.files(DirB,pattern=".tif$")))
        Cs <- Ms*Cs
        writeRaster(Cs,file.path(DirB,names(Cs)),format="GTiff",
                    bylayer=T,overwrite=T,NAflag=-9999)
      }
      
      #6.5. Fit ENM for Geographical Partition
      FitENM_TMLA_Parallel(RecordsData=occINPUT,Variables=envT,VarImP=imp_var,Fut=Fut,Part=part,Algorithm=Alg,PredictType=ensemble,spN=spN,
                  Tst=eval_occ,Threshold=thr,DirSave=DirR,DirMask=DirB,DirMSDM=DirPRI,Save=save_part,
                  SaveFinal=save_final,repl=NULL,per=NULL)
    }
    
#7.Random Partition----
    
    if(part=="boot"||part=="cross"){
      
      #7.0.Dataset for evaluation
      if(eval_occ=="Y"){
        cat("Select the occurrence dataset for evaluation:")
        OccTst <- read.table(file.choose(),sep="\t",h=T)
        OccTst<-OccTst[,c(sp,x,y)]
        colnames(OccTst) <- c("sp","x","y")
        OccTst_xy <- split(OccTst[,-1],f=OccTst$sp)
      }
      
      #7.1.msdm A PRIORI----
        if(msdm=="N"||msdm=="POST"){
          DirPRI <- NULL
        }
      
        if(msdm%in%c("XY","MIN","CML", "KER")){
          print("Creating msdm Layers...")
          
          DirMSDM<-"msdm"
          if (file.exists(file.path(dir,DirMSDM))){
            DirMSDM<-file.path(dir,DirMSDM)
          } else {
            dir.create(file.path(dir,DirMSDM))
            DirMSDM<-file.path(dir,DirMSDM)
          }
          
          DirPRI <- MSDM_Priori_TMLA(Species=occ_xy,var=envT,MSDM=msdm,DirMSDM=DirMSDM)
        }

      #7.2. Data Partition----
      if(part=="boot"){
        cat("Select the number of replicates (>=1):")
        rep <- as.integer(readLines(n = 1))
        while(is.na(rep)||rep<1){
          warning("Please Select a valid number of replicates (>=1)")
          cat("Select the number of replicates (>=1):")
          rep <- as.integer(readLines(n = 1))
        }
        cat("Select the proportion of occurrences used for fitting the model(0-1):")
        per<-as.numeric(readLines(n = 1))
        while(is.na(per)||per<=0 || per>1){
          warning("Please Select a valid partition of data (0-1)")
          cat("Select the percentage of occurrences used for fitting the model(0-1):")
          per<-as.numeric(readLines(n = 1))
        }
      }
      if(part=="cross"){
        cat("Select the number of k-folds (>=1):")
        rep <- as.integer(readLines(n = 1))
        per<-1
        while(is.na(rep)||rep<1){
          warning("Please Select a valid number of k-folds (>=1)")
          cat("Select the number of k-folds (>=1):")
          rep <- as.integer(readLines(n = 1))
        }
        occFold<- lapply(occ_xy, function(x) cbind(x,kfold(x,rep)))
        colsK <-  c("x","y","Partition");    
        occFold <- lapply(occFold, setNames, colsK)
        write.table(ldply(occFold,data.frame,.id="sp"),file.path(DirR,"GruposCrossValidation.txt"),sep="\t",row.names=F)
      }
      
      #Adjusting for determined evaluation dataset
      if(eval_occ=="Y" && part=="boot" && per!=1 || eval_occ=="Y" && part=="cross" && rep!=1){
        if(part=="boot"){
          warning("Adjusting data partition to one!")
          per <- 1
        }
        if(part=="cross"){
          warning("Adjusting partition to bootstrap and data partition to one! 
          Replicates will be equal to the original number of folds")
          part <- "boot"
          per <- 1
        }
      }
      
      #7.3.Replicates & Model Input----  
      occINPUT <- list()
      occTREINO <- list()
      occTESTE <- list()
        
      for(k in 1:rep){
        set.seed(k)
        if(rep==1){
          k <- NULL
        }
        if(part=="boot"){
          if(rep!=1){
            print(paste("Replicate.....",k),sep="")
          }
          tr <- lapply(occ_xy, function(x) sample(1:nrow(x),round(nrow(x)*per)))
          occTR <- list()
          occTS <- list()
          for(i in 1:length(occ_xy)){
            occTR[[i]] <- occ_xy[[i]][tr[[i]],]
            if(per==1){
              occTS[[i]] <- occTR[[i]]
            }else{
              occTS[[i]] <- occ_xy[[i]][-tr[[i]],]
            }
            occTR[[i]] <- cbind(occTR[[i]], rep(1,nrow(occTR[[i]])),rep(1,nrow(occTR[[i]])))
            occTS[[i]] <- cbind(occTS[[i]], rep(2,nrow(occTS[[i]])),rep(1,nrow(occTS[[i]])))
          }
          names(occTR) <- names(occ_xy)
          names(occTS) <- names(occTR)
          if(eval_occ=="Y"){
            occTS <- OccTst_xy
            occTS <- lapply(occTS, function(x) cbind(x, rep(2,nrow(x)),rep(1,nrow(x))))
            names(occTS) <- names(occTR)
          }
        }
        if(part=="cross"){
          print(paste("Adjsuting fold....",k,sep=""))
          occFoldK <- lapply(occFold, function(x) ifelse(x$Partition!=k,1,2))
          occFoldK <- Map(cbind,occ_xy,occFoldK)
          occFoldK <- lapply(occFoldK, setNames, colsK)
          occTR <- lapply(occFoldK, function(x) split(x,f=x$Partition)[[1]])
          PresAbse <- lapply(occTR, function(x) rep(1,nrow(x)))
          occTR <- Map(cbind,occTR,PresAbse)
          occTS <- lapply(occFoldK, function(x) split(x,f=x$Partition)[[2]])
          PresAbse <- lapply(occTS, function(x) rep(1,nrow(x)))
          occTS <- Map(cbind,occTS,PresAbse)
          names(occTR) <- names(occ_xy)
          names(occTS) <- names(occTR)
        }
        
      #7.4. Generating Pseudo-Absences----
        #Random Pseudo-Absences
          if(pseudoabs_method=="Rnd"){
            if(transfer=="Y"&& eval_occ=="Y"){
              pseudo.mask <- envT[[1]]
              pseudo.maskP <- EnvF[[1]][[1]]
            }else{
              pseudo.mask <- envT[[1]]
              pseudo.maskP <- envT[[1]]
            }
            
            absencesTR <- list()
            absencesTS <- list()
            for(s in 1:length(occTR)){
              set.seed(s)
              if(sp_accessible_area=="Y"){
                SpMask <- raster(file.path(DirM,paste0(names(occTR)[s],".tif")))
                SpMask <- pseudo.mask*SpMask
                if(sum(is.na(SpMask[])==F)<(pres_abs_ratio*nrow(occTR[[s]]))){
                  warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences") 
                  stop("Please try again with another restriction type or without restricting the extent")
                  
                }
                if(eval_occ=="Y"){
                  SpMaskP <- pseudo.maskP
                }else{
                  SpMaskP <- SpMask
                }
                absencesTR[[s]] <- randomPoints(SpMask, (1 / pres_abs_ratio)*nrow(occTR[[s]]),ext = extent(SpMask),prob = FALSE)
                absencesTS[[s]] <- randomPoints(SpMaskP, (1 / pres_abs_ratio)*nrow(occTS[[s]]),ext = extent(SpMask),prob = FALSE)
              }else{
                absencesTR[[s]] <- randomPoints(pseudo.mask, (1 / pres_abs_ratio)*nrow(occTR[[s]]),ext = extent(pseudo.mask),prob = FALSE)
                absencesTS[[s]] <- randomPoints(pseudo.maskP, (1 / pres_abs_ratio)*nrow(occTS[[s]]),ext = extent(pseudo.mask),prob = FALSE)
              }
            }
            absencesTR <- lapply(absencesTR, function(x) cbind(x, rep(1,nrow(x))))
            absencesTR <- lapply(absencesTR, function(x) cbind(x, rep(0,nrow(x))))
            absencesTS <- lapply(absencesTS, function(x) cbind(x, rep(2,nrow(x))))
            absencesTS <- lapply(absencesTS, function(x) cbind(x, rep(0,nrow(x))))
            if(is.null(k) && per==1 && eval_occ=="N"){
              for(i in 1:length(absencesTS)){
                absencesTS[[i]][,c("x","y")] <- absencesTR[[i]][,c("x","y")]
              }
            }
            DirCons <- NULL
          }
          
        #Bioclimatic Constrained Pseudo-Absences
          if(pseudoabs_method=="EnvConst"){
            
            DirCons <- "EnvConstrain"
            if (file.exists(file.path(dir,DirCons))){
              DirCons<-file.path(dir,DirCons)
            } else {
              dir.create(file.path(dir,DirCons))
              DirCons<-file.path(dir,DirCons)
            }
            
            #Check for Environmental Constrain Existence
            EnvMsk <- "N"
            if(all(paste0(spN,".tif")%in%list.files(DirCons,pattern=".tif"))){
              print("Environmental constrain already exists! Using already-created masks!")
              EnvMsk <- "Y"
            }
            
            absencesTR <- list()
            absencesTS <- list()
            
            for (i in 1:length(occTR)) {
              set.seed(i)
              if(EnvMsk=="N"){
                Model <- bioclim(envT, occTR[[i]][,c("x","y")])
                pseudo.mask <- dismo::predict(Model, envT, ext=extent(envT[[1]]))
                pseudo.mask <- round(pseudo.mask, 5)
                pseudo.mask <- (pseudo.mask-minValue(pseudo.mask))/
                                    (maxValue(pseudo.mask)-minValue(pseudo.mask))
                pseudo.mask <-(1-pseudo.mask)>=0.99
                writeRaster(pseudo.mask,paste(DirCons,spN[i],sep="/"),format="GTiff",overwrite=T)
                pseudo.mask[which(pseudo.mask[,]==FALSE)] <- NA
              }else{
                pseudo.mask <- raster(file.path(DirCons,paste0(spN[i],".tif")))
              }
            
              if(transfer=="Y"&& eval_occ=="Y"){
                pseudo.maskP <- EnvF[[1]][[1]]
              }else{
                pseudo.maskP <- pseudo.mask
              }

              if(sp_accessible_area=="Y"){
                SpMask <- raster(file.path(DirM,paste0(names(occTR)[i],".tif")))
                SpMask <- pseudo.mask*SpMask
                if(sum(is.na(SpMask[])==F)<(pres_abs_ratio*nrow(occTR[[i]]))){
                  warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences") 
                  stop("Please try again with another restriction type or without restricting the extent")
                  
                }
                if(eval_occ=="Y"){
                  SpMaskP <- pseudo.maskP
                }else{
                  SpMaskP <- SpMask
                }
                absencesTR.0 <- randomPoints(SpMask, (1 / pres_abs_ratio)*nrow(occTR[[i]]),ext = extent(SpMask),prob = FALSE)
                absencesTS.0 <- randomPoints(SpMaskP, (1 / pres_abs_ratio)*nrow(occTS[[i]]),ext = extent(SpMask),prob = FALSE)
              }else{
              absencesTR.0 <- randomPoints(pseudo.mask, (1 / pres_abs_ratio)*nrow(occTR[[i]]),
                                         ext = extent(pseudo.mask),
                                         prob = FALSE)
              if(eval_occ=="Y"){
                absencesTS.0 <- randomPoints(pseudo.maskP, (1 / pres_abs_ratio)*nrow(occTS[[i]]),
                                             ext = extent(pseudo.mask),
                                             prob = FALSE)                
              }else{
                absencesTS.0 <- randomPoints(pseudo.maskP, (1 / pres_abs_ratio)*nrow(occTS[[i]]),
                                           ext = extent(pseudo.mask),
                                           prob = FALSE)
              }
              }
              absencesTR[[i]] <- as.data.frame(absencesTR.0)
              absencesTS[[i]] <- as.data.frame(absencesTS.0)
              rm(absencesTR.0,absencesTS.0)
            }
            absencesTR <- lapply(absencesTR, function(x) cbind(x, rep(1,nrow(x)), rep(0,nrow(x))))
            absencesTS <- lapply(absencesTS, function(x) cbind(x, rep(2,nrow(x)), rep(0,nrow(x))))
            if(is.null(k) && per==1 && eval_occ=="N"){
              for(i in 1:length(absencesTS)){
                absencesTS[[i]][,c("x","y")] <- absencesTR[[i]][,c("x","y")]
              }
            }
          }
        
        #Geographical Constrained Pseudo-Absences
        if(pseudoabs_method=="GeoConst"){
          
          #Define Buffer distance:
          cat("Select buffer distance(in km):")
          Geo_Buf <- as.integer(readLines(n = 1))*1000
          
          
          DirCons <- "GeoConstrain"
          if (file.exists(file.path(dir,DirCons))){
            DirCons<-file.path(dir,DirCons)
          } else {
            dir.create(file.path(dir,DirCons))
            DirCons<-file.path(dir,DirCons)
          }
          
          #Check for Environmental Constrain Existence
          EnvMsk <- "N"
          if(all(paste0(spN,".tif")%in%list.files(DirCons,pattern=".tif"))){
            print("Geographical constrain already exists! Using already-created masks!")
            EnvMsk <- "Y"
          }
          
          absencesTR <- list()
          absencesTS <- list()
          
          for (i in 1:length(occTR)) {
            set.seed(i)
            if(EnvMsk=="N"){
              Model <- circles(occTR[[i]][,c("x","y")], lonlat=T,d=Geo_Buf)
              pseudo.mask <- mask(envT[[1]],Model@polygons,inverse=T)
              pseudo.mask[is.na(pseudo.mask)==F] <- 1 
              writeRaster(pseudo.mask,paste(DirCons,spN[i],sep="/"),format="GTiff",overwrite=T)
              pseudo.mask[which(pseudo.mask[,]==FALSE)] <- NA
            }else{
              pseudo.mask <- raster(file.path(DirCons,paste0(spN[i],".tif")))
            }
            
            if(transfer=="Y"&& eval_occ=="Y"){
              pseudo.maskP <- EnvF[[1]][[1]]
            }else{
              pseudo.maskP <- pseudo.mask
            }
            
            if(sp_accessible_area=="Y"){
              SpMask <- raster(file.path(DirM,paste0(names(occTR)[i],".tif")))
              SpMask <- pseudo.mask*SpMask
              if(sum(is.na(SpMask[])==F)<(pres_abs_ratio*nrow(occTR[[i]]))){
                warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences") 
                stop("Please try again with a smaller geographical buffer or without restricting the accessible area")
                
              }
              if(eval_occ=="Y"){
                SpMaskP <- pseudo.maskP
              }else{
                SpMaskP <- SpMask
              }
              absencesTR.0 <- randomPoints(SpMask, (1 / pres_abs_ratio)*nrow(occTR[[i]]),ext = extent(SpMask),prob = FALSE)
              absencesTS.0 <- randomPoints(SpMaskP, (1 / pres_abs_ratio)*nrow(occTS[[i]]),ext = extent(SpMask),prob = FALSE)
            }else{
              absencesTR.0 <- randomPoints(pseudo.mask, (1 / pres_abs_ratio)*nrow(occTR[[i]]),
                                           ext = extent(pseudo.mask),
                                           prob = FALSE)
              if(eval_occ=="Y"){
                absencesTS.0 <- randomPoints(pseudo.maskP, (1 / pres_abs_ratio)*nrow(occTS[[i]]),
                                             ext = extent(pseudo.mask),
                                             prob = FALSE)                
              }else{
                absencesTS.0 <- randomPoints(pseudo.maskP, (1 / pres_abs_ratio)*nrow(occTS[[i]]),
                                             ext = extent(pseudo.mask),
                                             prob = FALSE)
              }
            }
            absencesTR[[i]] <- as.data.frame(absencesTR.0)
            absencesTS[[i]] <- as.data.frame(absencesTS.0)
            rm(absencesTR.0,absencesTS.0)
          }
          absencesTR <- lapply(absencesTR, function(x) cbind(x, rep(1,nrow(x)), rep(0,nrow(x))))
          absencesTS <- lapply(absencesTS, function(x) cbind(x, rep(2,nrow(x)), rep(0,nrow(x))))
          if(is.null(k) && per==1 && eval_occ=="N"){
            for(i in 1:length(absencesTS)){
              absencesTS[[i]][,c("x","y")] <- absencesTR[[i]][,c("x","y")]
            }
          }
        }
          
        #Model Input
          for(i in 1:length(occTR)){
            occTR[[i]] <- cbind(rep(names(occTR)[i],nrow(occTR[[i]])),occTR[[i]])
            colnames(occTR[[i]]) <- c("sp","x","y","Partition","PresAbse")
            absencesTR[[i]] <- data.frame(cbind(rep(names(occTR)[i],nrow(absencesTR[[i]])),absencesTR[[i]]))
            colnames(absencesTR[[i]]) <- colnames(occTR[[i]])
            occTS[[i]] <- cbind(rep(names(occTR)[i],nrow(occTS[[i]])),occTS[[i]])
            colnames(occTS[[i]]) <- colnames(occTR[[i]])
            absencesTS[[i]] <- data.frame(cbind(rep(names(occTR)[i],nrow(absencesTS[[i]])),absencesTS[[i]]))
            colnames(absencesTS[[i]]) <- colnames(occTR[[i]])
          }
          occTR <- ldply(occTR,data.frame,.id=NULL)
          absencesTR <- ldply(absencesTR,data.frame,.id=NULL)
          occTS <- ldply(occTS,data.frame,.id=NULL)
          absencesTS <- ldply(absencesTS,data.frame,.id=NULL)
          occINPUT <- rbind(occTR,absencesTR,occTS,absencesTS)
          cols = c("x","y","Partition","PresAbse");    
          occINPUT[,cols] = apply(occINPUT[,cols], 2, function(x) as.numeric(as.character(x)))
          
      #7.5. Define Projection----
        if(transfer=="Y"){
          Fut <- EnvF
        }else{
          Fut <- NULL
        }
      
      #7.6. Calculate Moran & MESS----
      if(per!=1 || part=="cross"){
        Bootstrap_Moran_e_MESS_TMLA(Env=envT,RecordsData=occINPUT,DirO=DirR,repl=k)
      }
          
      #7.7. save_part Fix----
      if(save_part=="Y"&& per==1){
        save_part <- "N"
        warning("There are no partitions to be saved!")
      }
      
      #7.8. Background restriction----
      if(exists("DirM")){
        DirB <- DirM
      }else{
        DirB <- NULL
      }
          
      #7.9. Run FitENM----
        FitENM_TMLA_Parallel(RecordsData=occINPUT,Variables=envT,VarImP=imp_var,Fut=Fut,Part=part,Algorithm=Alg,PredictType=ensemble,spN=spN,
                    Tst=eval_occ,Threshold=thr,DirSave=DirR,DirMask=DirB,DirMSDM=DirPRI,Save=save_part,
                    SaveFinal=save_final,per=per,repl=k)
        
      #7.10. Create Occurrence Table for Replicates----
        if(rep!=1 || part=="cross"){
          occTREINO[[k]] <- occINPUT[occINPUT$Partition==1,]
          occTESTE[[k]] <- occINPUT[occINPUT$Partition==2,]
        }
      }#Fechas as replicas ou kfolds
        
      #7.11.Save Final Occurrence Table & Validation File----
        if(rep!=1){
          #Save Final Occurrence Table
          occTREINO <- ldply(occTREINO,data.frame,.id=NULL)
          occTESTE <- ldply(occTESTE,data.frame,.id=NULL)
          write.table(occTREINO,file.path(DirR,"OcorrenciasTreino.txt"),sep="\t",row.names=F)
          write.table(occTESTE,file.path(DirR,"OcorrenciasTeste.txt"),sep="\t",row.names=F)
        
          #Save Final Validation File
          val <- list.files(DirR,pattern="Validation_Partition")
          valF <- list()
          for(i in 1:length(val)){
            valF[[i]] <- read.table(file.path(DirR,val[i]),sep="\t",header=T)
          }
          valF <- ldply(valF,data.frame,.id=NULL)
          valF <- valF[order(as.character(valF[,1])),]
          unlink(file.path(DirR,val))
          write.table(valF,file.path(DirR,"PartialModels_Validation.txt"),sep="\t",row.names=F)
          
          # valFII <- ldply(valFII,data.frame,.id=NULL)
          # valFII <- valFII[order(as.character(valFII[,1])),]
          # unlink(file.path(DirR,valII))
          # write.table(valFII,file.path(DirR,"FullModels_Thresholds.txt"),sep="\t",row.names=F)
          
          
          #Save Final Bootstrap File
          if(per!=1 || part=="cross"){
            Boot <- list.files(DirR,pattern="Bootstrap_Moran_MESS")
            BootF <- list()
            for(i in Boot){
              BootF[[i]] <- read.table(file.path(DirR,i),sep="\t",header=T)
            }
            BootF <- ldply(BootF,data.frame,.id=NULL)
            BootF <- BootF[order(as.character(BootF[,1])),]
            unlink(file.path(DirR,Boot))
            if(part=="boot"){
              write.table(BootF,file.path(DirR,"Bootstrap_Moran_MESS.txt"),sep="\t",row.names=F)
            }
            if(part=="cross"){
              write.table(BootF,file.path(DirR,"CrossValidation_Moran_MESS.txt"),sep="\t",row.names=F)
            }
          }
        }
    }#Fecha partition boot|jknife
    
#8.MSDM Posteriori----
    
    if(msdm=="POST"){
      
      cat("Choose L-MSDM type (OBR/LR/PRES/MCP/MCP-B)")
      Q0 <- as.character(readLines(n = 1))
      while(Q0 %in% c("OBR","LR","PRES","MCP","MCP-B")==F){
        warning("Choose a valid L-MSDM type!(OBR/LR/PRES/MCP/MCP-B)")
        Q0 <- as.character(readLines(n = 1))
      }
      
      if(any(list.files(DirR)=="ENS")){
        cat("Perform L-MSDM only on Ensemble?(Y/N)")
        Q1 <- as.character(readLines(n = 1))
        while(!(Q1 %in% c("Y","N"))){
          warning("Select a valid response (Y/N):")
          Q1 <-as.character(readLines(n = 1))
        }
        if(Q1=="Y"){
          DirT <- file.path(DirR,"ENS",ensemble[ensemble!="N"])
          DirPost <- "MSDMPosterior"
          DirPost <- file.path(DirT,DirPost)
        }
        if(Q1=="N"){
          DirT <- file.path(DirR,Alg)
          DirPost <- "MSDMPosterior"
          DirPost <- file.path(DirT,DirPost)
        }
        for(i in DirPost){
          dir.create(i)
        }
        for(i in 1:length(DirPost)){
          print(paste("Diretorio.....",i,"/",length(DirPost),sep=""))
          MSDM_Posterior(RecordsData=occINPUT,Threshold=thr,cutoff=Q0,PredictType=ensemble,
                         DirSave=DirPost[i],DirRaster=DirT[i])
        }
      }else{
        Q1 <- "N"
        DirT <- file.path(DirR,Alg)
        DirPost <- "MSDMPosterior"
        DirPost <- file.path(DirT,DirPost)
        for(i in DirPost){
          dir.create(i)
        }
        for(i in 1:length(DirPost)){
          print(paste("Diretorio.....",i,"/",length(DirPost),sep=""))
          MSDM_Posterior(RecordsData=occINPUT,Threshold=thr,cutoff=Q0,PredictType=ensemble,
                         DirSave=DirPost[i],DirRaster=DirT[i])
        }
      }
      if(Q1=="N" ||!("ensemble"%in%list.files(DirR))){
        
        cat("Perform Ensemble on Algorithm L-MSDM?(Y/N)")
        Q4 <- as.character(readLines(n = 1))
        while(!(Q4 %in% c("Y","N"))){
          warning("Select a valid response (Y/N):")
          Q4 <- as.character(readLines(n = 1))
        }
        
        if(Q4=="Y"){
          DirT <- file.path(DirR,Alg,"MSDMPosterior")
          DirPost <- file.path(DirR,"ENS",ensemble,"MSDMPosterior")
          ENS_Posterior(RecordsData=occINPUT,Algorithm=Alg,PredictType=ensemble,Threshold=thr,DirAlg=DirT,DirSave=DirR)
        }
      }
    }
    
#9.S-SDM----
    if(s_sdm=="Y"){
      
      #Where to S-SDM
      if(ensemble!="N"){
        cat("Create S-SDM only for Ensemble?(Y/N)")
        Q1 <- as.character(readLines(n=1))
        if(Q1=="Y"){
          DirT <- file.path(DirR,"ENS",ensemble)
        }else{
          DirT <- file.path(DirR,Alg)
        }
      }else{
        cat("Creating S-SDM for each Algorithm")
        DirT <- file.path(DirR,Alg)
      }
      
      #S-SDM on MSDM?
      if(msdm=="POST"){
        cat("Create S-SDM for M-SDM?(Y/N)\nS-SDM for M-SDM will be based on previous M-SDM choice!")
        Q2 <- as.character(readLines(n=1))
        if(Q2=="Y"){
          DirT2 <- file.path(DirPost)
        }else{
          DirT2 <- NULL
        }
      }
      
      #S-SDM on Projections?
      if(transfer=="Y"){
        cat("Create S-SDM for Projections?(Y/N)")
        Q3 <- as.character(readLines(n=1))
        if(Q3=="Y" && Q1=="N"){
          DirT3 <- file.path(DirR,names(EnvF),Alg)
        }else if (Q3=="Y" && Q1=="Y"){
          DirT3 <- file.path(DirR,names(EnvF),ensemble)
        }else{
          DirT3 <- NULL
        }
      }
      
      #Calculate S-SDM
      S_SDM(DirENM=DirT,DirMSDM=DirT2,DirProj=DirT3,spN)
    }
}
