## Written by Andre Andrade & Santiago Velazco

ENMs_TheMetaLand<-function(Dir,
                           Sp,
                           x,
                           y,
                           NMin=10,
                           VarColin,
                           VarImP,
                           Proj,
                           SetEval,
                           SpeciesM,
                           PabR=1,
                           PabM,
                           Part,
                           SavePart="N",
                           SaveFinal="Y",
                           Alg,
                           Thr,
                           MSDM,
                           ENS,
                           S_SDM){
  
  #Script to fit ENMs from TheMetaLand Lab!
  
  #Parametro iniciais:
    #Dir: Folder with environmental variables
    #Sp : Name of the column containing species' names
    #x : Name of the column with longitude data
    #y : Name of the column with latitude data
    #NMin:Minimum number of occurrence to fit the models; spcies that do not meet this number will be excluded
    #VarColin: Perform methods to deal with collinearity in variables
      #N: Do not perform any method to deal with collinearity
      #Pearson: Select variables by Pearson Correlation (Threshold specified by user)
      #VIF: Variance Inflation Factor (Chatterjee and Hadi 2006)
      #PCA: Perform a PCA on predictors and use PCs as environmental variables
    #Proj : Project the model onto another set of predictors? (Y/N)
    #SetEval : Use an pre-determined set of occurrences for validation? (Y/N)
    #SpeciesM : Restrict the acessible area M? (Species-specific) (Y/N)
    #PabR:Presence-Absence Ratio
    #PabM:Pseudo-absence Selection Method
      #rnd:Random
      #EnvConst: Constrained by a Bioclim Model
      #GeoConst: Constrained by a geographical buffer
    #Part: Data partition methods 
      #boot : Random partition based on a "train" percentage
      #cross: Random partition of the data in k-folds
      #band : Geographic partition in latitudinal or longitudinal bands
      #check : Geographic partition in a checkerboard
    #Alg: Algorithms list
      #BIO : Bioclim
      #MXS : Maxent Linear and Quaratic Features (MaxNet)
      #MXD : Maxent Default (MaxNet)
      #SVM : Support Vector Machine
      #GLM : Generalized Linear Model
      #GAM : Generalizes Additive Model
      #RDF : Random Forest
      #MLK : Maximum Likelihood
      #GAU : Gaussian
    #Thr : Threshold used for presence-absence maps
      #no_omission : The highest threshold at which there is no omission
      #spec_sens : Threshold at which the sum of the sensitivity and specificity is highest
      #kappa: the threshold at which kappa is highest ("max kappa")
      #prevalence: modeled prevalence is closest to observed prevalence
      #equal_sens_spec: equal sensitivity and speciﬁcity
      #sensitivty: ﬁxed (speciﬁed) sensitivity
      #Any number between 0-1
    #MSDM: Spatial restrictions
      #N: none
      #XY: Latitudinal and Longitudinal information of each cell
      #MIN: Distance to the nearest occurrence
      #CML: Cummulative distance to all occurrences
      #KER: Kernel-Gauss
      #POST: Posterior M-SDM Methods (If chosen, prefered method will be asked later)
        #OBR: Occurrence based restriciton, uses the distance between points to exclude far suitable patches (Mendes et al, in prep)
        #LR: Lower Quantile, select the nearest 25% patches (Mendes et al, in prep)
        #PRES: Select only the patches with confirmed occurrence data (Mendes et al, in prep)
        #MCP: Excludes suitable cells outside the Minimum Convex Polygon of the occurrence data (Kremen et al, 2008)
        #MCP-B: Creates a Buffer around the MCP (distance defined by user; Kremen et al, 2008)
    #ENS: Create Ensemble Model
      #N : none
      #MEAN : Simple Average
      #W_MEAN: Weighted Average
      #SUP : Average of the best models (TSS over the average)
      #PCA : PCA with all models
      #PCA_SUP : PCA of the best models (TSS over the average)
      #PCA_THR : PCA only with cells above the threshold

#1.Check Function Arguments  
  
  er <- NULL
  if(missing(Dir)){
    er <- c(er,paste("Argumento Dir não definido,defina o diretorio dos preditores | "))
  }
  if(missing(Sp)){
    er <- c(er,paste("Argumento Sp não definido,defina a coluna com nome das especies | "))
  }
  if(missing(x)){
    er <- c(er,paste("Argumento x não definido,defina a coluna com os valores de longitude | "))
  }
  if(missing(y)){
    er <- c(er,paste("Argumento y não definido, defina a coluna com os valores de latitude | "))
  }
  if(missing(VarColin)){
    er <- c(er,paste("Argumento VarColin não definido, defina se deseja realizar PCA nos preditores | "))
  }
  if(missing(Proj)){
    er <- c(er,paste("Argumento Proj não definido, defina se deseja projetar o modelo para outra extensao/periodo temporal | "))
  }
  if(missing(PabR)){
    er <- c(er,paste("Argumento PabR não definido, defina a prevalencia entre treino/teste | "))
  }
  if(missing(PabM)){
    er <- c(er,paste("Argumento PabM não definido, defina o metodo de escolha de PseudoAusencias | "))
  }
  if(missing(Part)){
    er <- c(er,paste("Argumento Part não definido, defina o tipo de particao desejada | "))
  }
  if(missing(Alg)){
    er <- c(er,paste("Argumento Alg não definido, defina quais algoritmos deseja utilizar | "))
  }
  if(missing(Thr)){
    er <- c(er,paste("Argumento Thr não definido, defina qual Threshold deseja utilizar | "))
  }
  if(missing(MSDM)){
    er <- c(er,paste("Argumento MSDM não definido, defina deseja utilizar restricoes espaciais nos seus modelos | "))
  }
  if(missing(ENS)){
    er <- c(er,paste("Argumento ENS não definido, defina se deseja realizar a avaliacao dos modelos | "))
  }
  if(!is.null((er))){
    print(er)
    stop("Argumentos faltantes, por favor cheque os argumentos listados acima")
  }
  
  if(!(VarColin%in%c("Pearson","VIF","PCA","N"))){
    stop("VarColin Argument is not valid!(Pearson,VIF,PCA,N)")
  }
  if(!(Proj%in%c("Y","N"))){
    stop("Proj Argument is not valid!(Y/N)")
  }
  if(PabR<=0){
    stop("PabR Argument is not valid!(PabR>=0)")
  }
  if(!(PabM%in%c("Rnd","EnvConst","GeoConst"))){
    stop("PabM Argument is not valid!(Rnd/EnvConst/GeoConst)")
  }
  if(length(PabM)>1){
    stop("Please choose only one Pseudo-absence allocation method")
  }
  if(!(Part%in%c("boot","cross","band","check"))){
    stop("Part Argument is not valid!(boot/cross/band/check)")
  }
  if(any(!Alg%in%c("BIO","GLM","GAM","SVM","RDF","MXS","MXD","MLK","GAU"))){
    stop(paste("Algorithm",Alg[!(Alg%in%c("BIO","GLM","GAM","SVM","RDF","MXS","MXD","MLK","GAU"))],"is not valid"))
  }
  if(any(!Thr%in%c("no_omission","spec_sens","kappa","equal_sens_spec","prevalence","sensitivity"))){
    stop("Thr Argument is not valid!")
  }
  if(!(MSDM%in%c("N","XY","MIN","CML","KER","POST"))){
    stop("MSDM Argument is not valid!(N/XY/MIN/CML/KER/POST)")
  }
  if(length(MSDM)>1){
    stop("Please choose only one M-SDM method")
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
  Alg <- Ord[Ord%in%Alg]
  
#3.Predictors ----
  options(warn=1)
  setwd(Dir)
  
  env <- unique(file_ext(list.files()))
  form <- c('bil','asc','txt','tif')
  env <- env[env%in%form]
  if(length(env)>1){
    stop("More than one file format in Dir")
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
  if(Proj=="Y"){
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
  if(VarColin=="VIF"){
    VF <- vifstep(envT,th=nlayers(envT)*2)
    envT <- exclude(envT,VF)
    if(Proj=="Y"){
      RasM <- colMeans(na.omit(values(envT)))
      RasSTD <- apply(na.omit(values(envT)),2,std)
    }
    envT <- scale(envT)
    
    if(Proj=="Y"){
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
  
  if (VarColin=="PCA"){
      
    #Projection PCA
    if(Proj=="Y"){
      EnvF <- list()
      for(i in 1:length(Pfol)){
        EnvF[[i]] <- PCAFuturo(Env=envT,Dir=Dir,DirP=Pfol[i],Save="Y")
      }
      names(EnvF) <- PfolN
      envT <- brick(stack(file.path(Dir,"PCA",list.files(file.path(Dir,"PCA"),pattern='PC'))))
    }else{
        envT<-PCA_env_TMLA(env=envT,Dir=Dir)
    }
  }
  
  #3.3.3.Pearson----
  if(VarColin=="Pearson"){
    cat("Select correlation threshold:(0-1)")
    Cor_TH <- as.numeric(readLines(n=1))
    Pear <- layerStats(envT, 'pearson', na.rm=T)
    corr_matrix <- abs(Pear$'pearson correlation coefficient')
    corr_matrix[upper.tri(corr_matrix)] <- 0
    diag(corr_matrix) <- 0
    envT <- envT[[names(envT)[!apply(corr_matrix,2,function(x) any(x > 0.70))]]]
    if(Proj=="Y"){
      RasM <- colMeans(na.omit(values(envT)))
      RasSTD <- apply(na.omit(values(envT)),2,std)
    }
    envT <- scale(envT)
    
    if(Proj=="Y"){
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
  
  #3.3.Erro Futuro e MSDM
  if(Proj=="Y" && MSDM!="N"){
    warning("MSDM can not be used with future projections")
    warning("Setting MSDM to N")
    MSDM <- "N"
  }
  
  #3.4.Aviso caso NMin<NPreditores
  if(NMin<nlayers(envT)){
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
  occ<-occ[,c(Sp,x,y)]
  colnames(occ) <- c("sp","x","y")
  occ_xy <- split(occ[,-1],f=occ$sp)
  spN<-names(occ_xy)
  
  
    #4.1.Unique Occurrences----
    occA<-Occ_Unicas_TMLA(env=envT[[1]],occ.xy=occ_xy,DirO=DirR)

    #4.2.Thining----
    if(Thin=="Y"){
      cat(("Select thinning method:\n1-Distance defined by Moran Variogram\n2-User defined distance\n3-Distance defined by 10x cellsize (Haversine Transformation)"))
      ThinMet <- as.integer(readLines(n=1))
      occA <- OccsThin(occA,envT,ThinMet,VarColin,DirR)
    }
  
    #4.3.Remove species with less than NMin----
    occ <- occA[sapply(occA,function (x) nrow(x)>=NMin)]
    spN<-names(occ)
  
    
    #4.4.Species with few records----
    if(length(occ)!=length(occ_xy)){
      print(paste("Species with less than ",NMin, " Unique Occurrences were removed!"))
      print(names(occ_xy)[names(occ_xy)%in%spN==F])
      ndb <- ldply(occ)[,1:3]
      write.table(ndb,file.path(DirR,"Occ_Filtered.txt"),sep="\t",row.names=F)
      rm(ndb)
    }
    occ_xy <- lapply(occ,function(x) x[,c("x","y")])
    
    #4.5.GAM and GLM usage----
    if(any(sapply(occA,function(x) nrow(x))<nlayers(envT)) && any(Alg%in%c("GAM","GLM"))){
      warning("A species has fewer records than the number of predictors, impossible to adjust GAM and GLM! GAM and GLM will be excluded")
      Alg <- Alg[!Alg%in%c("GAM","GLM")]
    }
    
    
#5. Restrict Extent per Species----
    if(SpeciesM=="Y"){
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
                  Dir=Dir,
                  spN=spN,
                  SaveM = TRUE)
    }
  
    
#6. Geographical Partition----
    if(Part=="band" || Part=="check"){
      
      if(any(grepl("PC",names(envT)))==T || any(grepl("pc",names(envT)))==T){
        VarColin<-"PCA"
      }
      
      if(SetEval=="Y"){
        warning("Invalid combination! SetEval can't be Y with Geographical partitions! Changing SetEval to N")
        SetEval <- "N"
      }
      
      if(Part=="band"){  
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
          if(VarColin!="PCA"){
            envTT<-PCA_env_TMLA(envT,Dir)
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
                                      pseudoabsencesMethod=PabM,PrAbRatio=PabR,DirSave=DirB,
                                      DirM=DirM,MRst=SpeciesM,type=TipoMoran)
          occINPUT[,4] <- as.numeric(occINPUT[,4])
          occINPUT[,5] <- as.numeric(occINPUT[,5])
          rm(envTT)
        }

      }
      if(Part=="check"){
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
          if(VarColin!="PCA"){
            envTT<-PCA_env_TMLA(envT,Dir)
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
                                      pseudoabsencesMethod=PabM,PrAbRatio=PabR,DirSave=DirB,
                                      DirM=DirM,MRst=SpeciesM,type=TipoMoran)

          occINPUT[,4] <- as.numeric(occINPUT[,4])
          occINPUT[,5] <- as.numeric(occINPUT[,5])
          rm(envTT)
        }
      }
      
      #6.3.MSDM A PRIORI----
      if(MSDM=="N"||MSDM=="POST"){
        DirPRI <- NULL
      }
      
      if(MSDM%in%c("XY","MIN","CML", "KER")){
        print("Creating MSDM Layers...")
        
        DirMSDM<-"MSDM"
        if (file.exists(file.path(Dir,DirMSDM))){
          DirMSDM<-file.path(Dir,DirMSDM)
        } else {
          dir.create(file.path(Dir,DirMSDM))
          DirMSDM<-file.path(Dir,DirMSDM)
        }
        DirPRI <- MSDM_Priori_TMLA(Species=occ_xy,var=envT,MSDM=MSDM,DirMSDM=DirMSDM)
      }

      #6.4. Future Projections ----
      if(Proj=="Y"){
        Fut <- EnvF
      }else{
        Fut <- NULL
      }
      
      #6.5.Adjust Checkerboard when using Geographical Restrictions (For Maxent Sampling)
      if(SpeciesM=="Y"){
        Ms <- stack(file.path(DirM,list.files(DirM)))
        Cs <- stack(file.path(DirB,list.files(DirB,pattern=".tif$")))
        Cs <- Ms*Cs
        writeRaster(Cs,file.path(DirB,names(Cs)),format="GTiff",
                    bylayer=T,overwrite=T,NAflag=-9999)
      }
      
      #6.5. Fit ENM for Geographical Partition
      FitENM_TMLA_Parallel(RecordsData=occINPUT,Variables=envT,VarImP=VarImP,Fut=Fut,Part=Part,Algorithm=Alg,PredictType=ENS,spN=spN,
                  Tst=SetEval,Threshold=Thr,DirSave=DirR,DirMask=DirB,DirMSDM=DirPRI,Save=SavePart,
                  SaveFinal=SaveFinal,repl=NULL,per=NULL)
    }
    
#7.Random Partition----
    
    if(Part=="boot"||Part=="cross"){
      
      #7.0.Dataset for evaluation
      if(SetEval=="Y"){
        cat("Select the occurrence dataset for evaluation:")
        OccTst <- read.table(file.choose(),sep="\t",h=T)
        OccTst<-OccTst[,c(Sp,x,y)]
        colnames(OccTst) <- c("sp","x","y")
        OccTst_xy <- split(OccTst[,-1],f=OccTst$sp)
      }
      
      #7.1.MSDM A PRIORI----
        if(MSDM=="N"||MSDM=="POST"){
          DirPRI <- NULL
        }
      
        if(MSDM%in%c("XY","MIN","CML", "KER")){
          print("Creating MSDM Layers...")
          
          DirMSDM<-"MSDM"
          if (file.exists(file.path(Dir,DirMSDM))){
            DirMSDM<-file.path(Dir,DirMSDM)
          } else {
            dir.create(file.path(Dir,DirMSDM))
            DirMSDM<-file.path(Dir,DirMSDM)
          }
          
          DirPRI <- MSDM_Priori_TMLA(Species=occ_xy,var=envT,MSDM=MSDM,DirMSDM=DirMSDM)
        }

      #7.2. Data Partition----
      if(Part=="boot"){
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
      if(Part=="cross"){
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
      if(SetEval=="Y" && Part=="boot" && per!=1 || SetEval=="Y" && Part=="cross" && rep!=1){
        if(Part=="boot"){
          warning("Adjusting data partition to one!")
          per <- 1
        }
        if(Part=="cross"){
          warning("Adjusting partition to bootstrap and data partition to one! 
          Replicates will be equal to the original number of folds")
          Part <- "boot"
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
        if(Part=="boot"){
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
          if(SetEval=="Y"){
            occTS <- OccTst_xy
            occTS <- lapply(occTS, function(x) cbind(x, rep(2,nrow(x)),rep(1,nrow(x))))
            names(occTS) <- names(occTR)
          }
        }
        if(Part=="cross"){
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
          if(PabM=="rnd"){
            if(Proj=="Y"&& SetEval=="Y"){
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
              if(SpeciesM=="Y"){
                SpMask <- raster(file.path(DirM,paste0(names(occTR)[s],".tif")))
                SpMask <- pseudo.mask*SpMask
                if(sum(is.na(SpMask[])==F)<(PabR*nrow(occTR[[s]]))){
                  warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences") 
                  stop("Please try again with another restriction type or without restricting the extent")
                  
                }
                if(SetEval=="Y"){
                  SpMaskP <- pseudo.maskP
                }else{
                  SpMaskP <- SpMask
                }
                absencesTR[[s]] <- randomPoints(SpMask, (1 / PabR)*nrow(occTR[[s]]),ext = extent(SpMask),prob = FALSE)
                absencesTS[[s]] <- randomPoints(SpMaskP, (1 / PabR)*nrow(occTS[[s]]),ext = extent(SpMask),prob = FALSE)
              }else{
                absencesTR[[s]] <- randomPoints(pseudo.mask, (1 / PabR)*nrow(occTR[[s]]),ext = extent(pseudo.mask),prob = FALSE)
                absencesTS[[s]] <- randomPoints(pseudo.maskP, (1 / PabR)*nrow(occTS[[s]]),ext = extent(pseudo.mask),prob = FALSE)
              }
            }
            absencesTR <- lapply(absencesTR, function(x) cbind(x, rep(1,nrow(x))))
            absencesTR <- lapply(absencesTR, function(x) cbind(x, rep(0,nrow(x))))
            absencesTS <- lapply(absencesTS, function(x) cbind(x, rep(2,nrow(x))))
            absencesTS <- lapply(absencesTS, function(x) cbind(x, rep(0,nrow(x))))
            if(is.null(k) && per==1 && SetEval=="N"){
              for(i in 1:length(absencesTS)){
                absencesTS[[i]][,c("x","y")] <- absencesTR[[i]][,c("x","y")]
              }
            }
            DirCons <- NULL
          }
          
        #Bioclimatic Constrained Pseudo-Absences
          if(PabM=="EnvConst"){
            
            DirCons <- "EnvConstrain"
            if (file.exists(file.path(Dir,DirCons))){
              DirCons<-file.path(Dir,DirCons)
            } else {
              dir.create(file.path(Dir,DirCons))
              DirCons<-file.path(Dir,DirCons)
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
            
              if(Proj=="Y"&& SetEval=="Y"){
                pseudo.maskP <- EnvF[[1]][[1]]
              }else{
                pseudo.maskP <- pseudo.mask
              }

              if(SpeciesM=="Y"){
                SpMask <- raster(file.path(DirM,paste0(names(occTR)[i],".tif")))
                SpMask <- pseudo.mask*SpMask
                if(sum(is.na(SpMask[])==F)<(PabR*nrow(occTR[[i]]))){
                  warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences") 
                  stop("Please try again with another restriction type or without restricting the extent")
                  
                }
                if(SetEval=="Y"){
                  SpMaskP <- pseudo.maskP
                }else{
                  SpMaskP <- SpMask
                }
                absencesTR.0 <- randomPoints(SpMask, (1 / PabR)*nrow(occTR[[i]]),ext = extent(SpMask),prob = FALSE)
                absencesTS.0 <- randomPoints(SpMaskP, (1 / PabR)*nrow(occTS[[i]]),ext = extent(SpMask),prob = FALSE)
              }else{
              absencesTR.0 <- randomPoints(pseudo.mask, (1 / PabR)*nrow(occTR[[i]]),
                                         ext = extent(pseudo.mask),
                                         prob = FALSE)
              if(SetEval=="Y"){
                absencesTS.0 <- randomPoints(pseudo.maskP, (1 / PabR)*nrow(occTS[[i]]),
                                             ext = extent(pseudo.mask),
                                             prob = FALSE)                
              }else{
                absencesTS.0 <- randomPoints(pseudo.maskP, (1 / PabR)*nrow(occTS[[i]]),
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
            if(is.null(k) && per==1 && SetEval=="N"){
              for(i in 1:length(absencesTS)){
                absencesTS[[i]][,c("x","y")] <- absencesTR[[i]][,c("x","y")]
              }
            }
          }
        
        #Geographical Constrained Pseudo-Absences
        if(PabM=="GeoConst"){
          
          #Define Buffer distance:
          cat("Select buffer distance(in km):")
          Geo_Buf <- as.integer(readLines(n = 1))*1000
          
          
          DirCons <- "GeoConstrain"
          if (file.exists(file.path(Dir,DirCons))){
            DirCons<-file.path(Dir,DirCons)
          } else {
            dir.create(file.path(Dir,DirCons))
            DirCons<-file.path(Dir,DirCons)
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
            
            if(Proj=="Y"&& SetEval=="Y"){
              pseudo.maskP <- EnvF[[1]][[1]]
            }else{
              pseudo.maskP <- pseudo.mask
            }
            
            if(SpeciesM=="Y"){
              SpMask <- raster(file.path(DirM,paste0(names(occTR)[i],".tif")))
              SpMask <- pseudo.mask*SpMask
              if(sum(is.na(SpMask[])==F)<(PabR*nrow(occTR[[i]]))){
                warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences") 
                stop("Please try again with a smaller geographical buffer or without restricting the accessible area")
                
              }
              if(SetEval=="Y"){
                SpMaskP <- pseudo.maskP
              }else{
                SpMaskP <- SpMask
              }
              absencesTR.0 <- randomPoints(SpMask, (1 / PabR)*nrow(occTR[[i]]),ext = extent(SpMask),prob = FALSE)
              absencesTS.0 <- randomPoints(SpMaskP, (1 / PabR)*nrow(occTS[[i]]),ext = extent(SpMask),prob = FALSE)
            }else{
              absencesTR.0 <- randomPoints(pseudo.mask, (1 / PabR)*nrow(occTR[[i]]),
                                           ext = extent(pseudo.mask),
                                           prob = FALSE)
              if(SetEval=="Y"){
                absencesTS.0 <- randomPoints(pseudo.maskP, (1 / PabR)*nrow(occTS[[i]]),
                                             ext = extent(pseudo.mask),
                                             prob = FALSE)                
              }else{
                absencesTS.0 <- randomPoints(pseudo.maskP, (1 / PabR)*nrow(occTS[[i]]),
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
          if(is.null(k) && per==1 && SetEval=="N"){
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
        if(Proj=="Y"){
          Fut <- EnvF
        }else{
          Fut <- NULL
        }
      
      #7.6. Calculate Moran & MESS----
      if(per!=1 || Part=="cross"){
        Bootstrap_Moran_e_MESS_TMLA(Env=envT,RecordsData=occINPUT,DirO=DirR,repl=k)
      }
          
      #7.7. SavePart Fix----
      if(SavePart=="Y"&& per==1){
        SavePart <- "N"
        warning("There are no partitions to be saved!")
      }
      
      #7.8. Background restriction----
      if(exists("DirM")){
        DirB <- DirM
      }else{
        DirB <- NULL
      }
          
      #7.9. Run FitENM----
        FitENM_TMLA_Parallel(RecordsData=occINPUT,Variables=envT,VarImP=VarImP,Fut=Fut,Part=Part,Algorithm=Alg,PredictType=ENS,spN=spN,
                    Tst=SetEval,Threshold=Thr,DirSave=DirR,DirMask=DirB,DirMSDM=DirPRI,Save=SavePart,
                    SaveFinal=SaveFinal,per=per,repl=k)
        
      #7.10. Create Occurrence Table for Replicates----
        if(rep!=1 || Part=="cross"){
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
          if(per!=1 || Part=="cross"){
            Boot <- list.files(DirR,pattern="Bootstrap_Moran_MESS")
            BootF <- list()
            for(i in Boot){
              BootF[[i]] <- read.table(file.path(DirR,i),sep="\t",header=T)
            }
            BootF <- ldply(BootF,data.frame,.id=NULL)
            BootF <- BootF[order(as.character(BootF[,1])),]
            unlink(file.path(DirR,Boot))
            if(Part=="boot"){
              write.table(BootF,file.path(DirR,"Bootstrap_Moran_MESS.txt"),sep="\t",row.names=F)
            }
            if(Part=="cross"){
              write.table(BootF,file.path(DirR,"CrossValidation_Moran_MESS.txt"),sep="\t",row.names=F)
            }
          }
        }
    }#Fecha partition boot|jknife
    
#8.MSDM Posteriori----
    
    if(MSDM=="POST"){
      
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
          DirT <- file.path(DirR,"ENS",ENS[ENS!="N"])
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
          MSDM_Posterior(RecordsData=occINPUT,Threshold=Thr,cutoff=Q0,PredictType=ENS,
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
          MSDM_Posterior(RecordsData=occINPUT,Threshold=Thr,cutoff=Q0,PredictType=ENS,
                         DirSave=DirPost[i],DirRaster=DirT[i])
        }
      }
      if(Q1=="N" ||!("ENS"%in%list.files(DirR))){
        
        cat("Perform Ensemble on Algorithm L-MSDM?(Y/N)")
        Q4 <- as.character(readLines(n = 1))
        while(!(Q4 %in% c("Y","N"))){
          warning("Select a valid response (Y/N):")
          Q4 <- as.character(readLines(n = 1))
        }
        
        if(Q4=="Y"){
          DirT <- file.path(DirR,Alg,"MSDMPosterior")
          DirPost <- file.path(DirR,"ENS",ENS,"MSDMPosterior")
          ENS_Posterior(RecordsData=occINPUT,Algorithm=Alg,PredictType=ENS,Threshold=Thr,DirAlg=DirT,DirSave=DirR)
        }
      }
    }
    
#9.S-SDM----
    if(S_SDM=="Y"){
      
      #Where to S-SDM
      if(ENS!="N"){
        cat("Create S-SDM only for Ensemble?(Y/N)")
        Q1 <- as.character(readLines(n=1))
        if(Q1=="Y"){
          DirT <- file.path(DirR,"ENS",ENS)
        }else{
          DirT <- file.path(DirR,Alg)
        }
      }else{
        cat("Creating S-SDM for each Algorithm")
        DirT <- file.path(DirR,Alg)
      }
      
      #S-SDM on MSDM?
      if(MSDM=="POST"){
        cat("Create S-SDM for M-SDM?(Y/N)\nS-SDM for M-SDM will be based on previous M-SDM choice!")
        Q2 <- as.character(readLines(n=1))
        if(Q2=="Y"){
          DirT2 <- file.path(DirPost)
        }else{
          DirT2 <- NULL
        }
      }
      
      #S-SDM on Projections?
      if(Proj=="Y"){
        cat("Create S-SDM for Projections?(Y/N)")
        Q3 <- as.character(readLines(n=1))
        if(Q3=="Y" && Q1=="N"){
          DirT3 <- file.path(DirR,names(EnvF),Alg)
        }else if (Q3=="Y" && Q1=="Y"){
          DirT3 <- file.path(DirR,names(EnvF),ENS)
        }else{
          DirT3 <- NULL
        }
      }
      
      #Calculate S-SDM
      S_SDM(DirENM=DirT,DirMSDM=DirT2,DirProj=DirT3,spN)
    }
}
