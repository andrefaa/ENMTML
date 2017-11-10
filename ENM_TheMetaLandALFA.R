## Written by Andre Andrade & Santiago Velazco

ENMs_TheMetaLand<-function(Dir,
                           Sp,
                           x,
                           y,
                           NMin=10,
                           PCA,
                           Proj,
                           PabR,
                           PabM,
                           Part,
                           SavePart="N",
                           Alg,
                           Thr,
                           MSDM,
                           ENS){
  
  #Funcao para criar modelos de distribuicao do laboratorio TheMetaLand
  
  #Parametro iniciais:
    #Dir: Diretorio com as camadas ambientais
    #Sp : Nome da coluna com os nomes das especies
    #x : Nome da coluna com os dados de longitude
    #y : Nome da coluna com os dados de latitude
    #NMin:Número Mínimo de Ocorrências por Especie, Especies com valores abaixo deste numero serao descartadas
    #PCA: Deseja realizar uma PCA nas var.ambientais e utilizar os eixos?(S/N)
    #Proj : Deseja projetar o modelo em outra area? (S/N)
    #PabR:Presence-Absence Ratio
    #PabM:Pseudo-absence Selection Method
      #rnd:Random
      #const: Constrained by a Bioclim Model
      #zoo: Constrained by ZooRegions(Holt et al 2013)
    #Part: Tipo de particao de dados 
      #boot : Particao aleatoria dos dados baseada em uma porcentagem
      #cross: Particao aleatoria dos dados em k grupos
      #band : Particao geografica dos dados em bandas latitudinais ou longitudinais
      #check : Particao geografica dos dados em quadriculas
    #Alg: Lista de Algoritmos
      #BIO : Bioclim
      #MXS : Maxent Simples (MaxNet)
      #MXD : Maxent Default (MaxNet)
      #SVM : Support Vector Machine
      #GLM : Generalized Linear Model
      #GAM : Generalizes Additive Model
      #RDF : Random Forest
      #MDA : Mixture Discriminant Analysis
      #MLK : Maximum Likelihood
      #GAU : Gaussian
      #ANN : Artificial Neural Network
    #Thr : Definir um limiar para cortar os modelos em presença-ausência
      #LPT : The highest threshold at which there is no omission
      #MAX : Threshold at which the sum of the sensitivity and specificity is highest
    #MSDM: Incluir Restricoes Espaciais
      #N: nao incluir
      #LatLong: coordenadas x e y
      #Min: Distancia ao ponto mais proximo
      #Cum: Distancia cumulativa aos pontos
      #Kern: Kernel-Gauss
      #Land: Manchas de Paisagem
    #ENS: Criar um modelo consenso dos vários algoritmos
      #N : Nao cria Consenso
      #Mean : Cria um consenso com a media simples
      #Sup : Cria um consenso medio dos modelos que possuiram TSS acima da media
      #PCA : Realiza uma PCA geral
      #PCA_Sup : Realiza uma PCA apenas com modelos com TSS acima da media
      #PCA_Thr : Realiza uma PCA ignorando celulas abaixo do Threshold

#1.Check Function Arguments  
  
  er <- NULL
  if(missing(Dir)){
    er <- c(er,paste("Argumento dir não definido,defina o diretorio dos preditores | "))
  }
  if(missing(Sp)){
    er <- c(er,paste("Argumento sp não definido,defina a coluna com nome das especies | "))
  }
  if(missing(x)){
    er <- c(er,paste("Argumento x não definido,defina a coluna com os valores de longitude | "))
  }
  if(missing(y)){
    er <- c(er,paste("Argumento y não definido, defina a coluna com os valores de latitude | "))
  }
  if(missing(PCA)){
    er <- c(er,paste("Argumento pca não definido, defina se deseja realizar PCA nos preditores | "))
  }
  if(missing(Proj)){
    er <- c(er,paste("Argumento pj não definido, defina se deseja projetar o modelo para outra extensao/periodo temporal | "))
  }
  if(missing(PabR)){
    er <- c(er,paste("Argumento prev não definido, defina a prevalencia entre treino/teste | "))
  }
  if(missing(PabM)){
    er <- c(er,paste("Argumento PabM não definido, defina o metodo de escolha de PseudoAusencias | "))
  }
  if(missing(Part)){
    er <- c(er,paste("Argumento part não definido, defina o tipo de particao desejada | "))
  }
  if(missing(Alg)){
    er <- c(er,paste("Argumento alg não definido, defina quais algoritmos deseja utilizar | "))
  }
  if(missing(Thr)){
    er <- c(er,paste("Argumento Thr não definido, defina qual Threshold deseja utilizar | "))
  }
  if(missing(MSDM)){
    er <- c(er,paste("Argumento msdm não definido, defina deseja utilizar restricoes espaciais nos seus modelos | "))
  }
  if(missing(ENS)){
    er <- c(er,paste("Argumento eval não definido, defina se deseja realizar a avaliacao dos modelos | "))
  }
  if(!is.null((er))){
    print(er)
    stop("Argumentos faltantes, por favor cheque os argumentos listados acima")
  }
  
  if(!(PCA%in%c("S","N"))){
    stop("PCA Argument is not valid!(S/N)")
  }
  if(!(Proj%in%c("S","N"))){
    stop("Proj Argument is not valid!(S/N)")
  }
  if(PabR<=0){
    stop("PabR Argument is not valid!(PabR>=0)")
  }
  if(!(PabM%in%c("rnd","const","zoo"))){
    stop("PabM Argument is not valid!(rnd/const/zoo)")
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
  if(any(!Thr%in%c("LPT","MAX"))){
    stop("Thr Argument is not valid!(LPT/MAX)")
  }
  if(!(MSDM%in%c("N","LatLong","Min","Cum","Kern","Land"))){
    stop("MSDM Argument is not valid!(N/LatLong/Min/Cum/Kern/Land)")
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
  
  
  ipak(c("raster","sp","dismo","kernlab","xlsx","randomForest","mda","rgdal","dummies",
         "MASS","ade4","gam","mvtnorm","progress","maxnet","maptools","XML","maxlike",
         "mgcv", "plyr", "GRaF","RStoolbox","flexclust","ape","tools","modEvA","XML",
         "SDMTools"))
  
#2.Load Auxiliary Functions ----
  
  if (file.exists("C:\\Scripts_for_ENM_TheMetaLand")==F){
    stop("Place the folder 'Scripts_for_ENM_TheMetaLand' in C:")
  }
  source("C:\\Scripts_for_ENM_TheMetaLand\\PCA_env_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\PCA_ENS_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\Occ_Unicas_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\BandsPartition_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\BlockPartition_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\Random_Moran_e_Schoener_TML.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\Moran_for_Quadrants_Pair_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\Evaluation_TML.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\BackZoo_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\SUMMRES.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\predict.graf.raster.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\FitENM_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\MSDM_Priori_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\AuxiliaryFuncENM_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\PredictENM_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\M_SDM_posteriori_SJEV_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\PCAFuturo_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\ENS_Posterior_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\maxnet2.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\Moran_for_Bootstrap_TMLA.R")
  source("C:\\Scripts_for_ENM_TheMetaLand\\Bootstrap_Moran_e_MESS_TMLA.R")
  
#3.Predictors ----
  
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
  
  #3.0.Consistencia entre as variaveis ambientais!
  if(length(unique(colSums(!is.na(envT[]))))>1){
    envT[is.na((sum(envT))[])] <- NA
    print("Variables had differences, setting any empty cells to NA in all variables")
  }
  
    #3.1.Projection----
  
    if (Proj=='S'){
      print("Select folder containing GCM folders:")
      DirP<-choose.dir(getwd())
      Pfol<-file.path(DirP,list.files(DirP))
      if(any(file_ext(list.files(DirP))%in%form)){
        stop("Select a folder containing GCM folders, NOT a folder with GCM variables!")
      }
      PfolN <- list.files(DirP)
      if(PCA=="S"){
        EnvF <- list()
        for(i in 1:length(Pfol)){
          EnvF[[i]] <- PCAFuturo(Env=envT,Dir=Dir,DirP=Pfol[i],Save="Y")
        }
        names(EnvF) <- PfolN
      }else{
        EnvF <- list()
        for(i in 1:length(Pfol)){
          ProjT <- unique(file_ext(list.files(Pfol[i])))
          form <- c('bil','asc','txt','tif')
          ProjT <- ProjT[ProjT%in%form]
          if(length(ProjT)>1){
            stop("More than one file format in DirP")
          }
          
          if(ProjT == 'bil'){
            EnvF[[i]]<-brick(stack(file.path(Pfol[i],list.files(Pfol[i],pattern='.bil'))))
          }
          
          if(ProjT == 'asc'){
            EnvF[[i]]<-brick(stack(file.path(Pfol[i],list.files(Pfol[i],pattern='.asc'))))
          }
          
          if(ProjT == 'txt'){
            ProjTT<-read.table(file.path(Pfol[i],list.files(Pfol[i],pattern='.txt'),h=T))
            gridded(ProjTT)<- ~x+y
            EnvF[[i]]<-brick(stack(ProjTT))
            rm(ProjTT)
          }
          
          if(ProjT == 'tif'){
            EnvF[[i]]<-brick(stack(file.path(Pfol[i],list.files(Pfol[i],pattern='.tif'))))
          }
        }
        names(EnvF) <- PfolN
      }
    }
  
    #3.2.PCA----
  
    if (PCA=="S"){
      if(Proj=="S"){
       envT <- brick(stack(list.files(pattern='PC')))
      }else{
        envT<-PCA_env_TMLA(envT,Dir)
      }
    }
  
    #3.3. Check Present/Future Consistency
    if(Proj=="S" && (!PabM%in%c("zoo","const"))){
      if(all(names(envT)!=names(EnvF[[1]]))){
        stop("Present/Future Variables Do Not Match! Make sure Present/Future Variables have the same names")
      }
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
    occ<-Occ_Unicas_TMLA(envT[[1]],occ_xy)
    occ <- occ[sapply(occ,function (x) nrow(x)>=NMin)]
    spN<-names(occ)
    
    #4.2.Species with few records----
    if(length(occ)!=length(occ_xy)){
      print(paste("Species with less than ",NMin, " Unique Occurrences were removed!"))
      print(names(occ_xy)[names(occ_xy)%in%spN==F])
      uni <- data.frame(Species=spN,UniqueOcc=sapply(occ,function(x) nrow(x)))
      write.table(uni,"N_Occ_Unicas.txt",sep="\t",row.names=F)
      ndb <- ldply(occ)[,1:3]
      write.table(ndb,paste(DirR,"Occ_Filtered.txt",sep="/"),sep="\t",row.names=F)
      rm(list=c("ndb","uni"))
    }
    occ_xy <- lapply(occ,function(x) x[,c("x","y")])
  
    
#5. Geographical Partition----
    if(Part=="band" || Part=="check"){
      
      if(any(grepl("PC",names(envT)))==T || any(grepl("pc",names(envT)))==T){
        PCA<-"S"
      }
      if(Part=="band"){  
        #5.1.Bands----
        
        DirB<-"Bands"
        if (file.exists(file.path(DirR,DirB))){
          DirB<-file.path(DirR,DirB)
        } else {
          dir.create(file.path(DirR,DirB))
          DirB<-file.path(DirR,DirB)
        }
        
        if(length(list.files(DirB,pattern=".tif")) ==(length(occ))){
          warning("Partition Already Exist! Using pre-created partitions! ")
          setwd(DirB)
          occT <- read.table(file.path(DirB,"OccBands.txt"),sep="\t",header=T)
          occT[,4] <- as.numeric(occT[,4])
          occT[,5] <- as.numeric(occT[,5])
        }else{
          if(PCA=="N"){
            envT<-PCA_env_TMLA(envT)
          }else{
            envT<-envT
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
          
          occT <- BandsPartition_TMLA(evnVariables=envT,RecordsData=occ_xy,N=bands,
                                      pseudoabsencesMethod=PabM,PrAbRatio=PabR,DirSave=DirB,
                                      type=TipoMoran)
          occT[,4] <- as.numeric(occT[,4])
          occT[,5] <- as.numeric(occT[,5])
        }

      }
      if(Part=="check"){
        #5.2.Block----
        
        DirB<-"Blocks"
        if (file.exists(file.path(DirR,DirB))){
          DirB<-file.path(DirR,DirB)
        } else {
          dir.create(file.path(DirR,DirB))
          DirB<-file.path(DirR,DirB)
        }
        
        if(length(list.files(DirB,pattern=".tif")) ==(length(occ))){
          print("Partition Already Exist! Using pre-created partitions! ")
          setwd(DirB)
          occT <- read.table(file.path(DirB,"OccBlocks.txt"),sep="\t",header=T)
          occT[,4] <- as.numeric(occT[,4])
          occT[,5] <- as.numeric(occT[,5])
        }else{
          if(PCA=="N"){
            envT<-PCA_env_TMLA(envT)
          }else{
            envT<-envT
          }
          print("Select Number of Blocks (>=2)")
          blocks <- as.integer(readLines(n = 1))
          while(is.na(blocks)||blocks<2){
            warning("Please Select a valid number of blocks (>=2)")
            print("Select Number of Blocks (>=2)")
            blocks <- as.integer(readLines(n = 1))
          }
          print("Select Moran Calculation Type (all/nearest):")
          TipoMoran <- as.character(readLines(n = 1))
          while(is.na(TipoMoran)||!(TipoMoran%in%c("all","nearest"))){
            warning("Please choose a valid Moran Calculation Type [all/nearest]")
            print("Select Moran Calculation Type (all/nearest):")
            TipoMoran <- as.character(readLines(n = 1))
          }
          
          occT <- BlockPartition_TMLA(evnVariables=envT,RecordsData=occ_xy,N=blocks,
                                      pseudoabsencesMethod=PabM,PrAbRatio=PabR,DirSave=DirB,
                                      type=TipoMoran)
          occT[,4] <- as.numeric(occT[,4])
          occT[,5] <- as.numeric(occT[,5])
        }
      }
      
      #5.3. Define Threshold----
      Thresh <- NULL
      if(any(Thr=="LPT")){
        Thresh <- c(Thresh,'no_omission')
      }
      if(any(Thr=="MAX")){
        Thresh <- c(Thresh,'spec_sens')
      }
      
      #5.4. Future Projections ----
      if(Proj=="S"){
        Fut <- EnvF
      }else{
        Fut <- NULL
      }
      
      FitENM_TMLA(RecordsData=occT,Variables=envT,Fut=Fut,Part=Part,Algorithm=Alg,PredictType=ENS,
                  Threshold=Thresh,DirSave=DirR,DirMask=DirB,DirMSDM=NULL,repl=NULL,Save=SavePart)
    }

#6.Random Partition----
    
    if(Part=="boot"||Part=="cross"){
      
      #6.1.MSDM A PRIORI----
        if(MSDM=="N"||MSDM=="Land"){
          DirPRI <- NULL
        }
      
        if(MSDM%in%c("LatLong","Min","Cum", "Kern")){
          print("Creating MSDM Layers...")
          
          DirMSDM<-"MSDM"
          if (file.exists(file.path(Dir,DirMSDM))){
            DirMSDM<-file.path(Dir,DirMSDM)
          } else {
            dir.create(file.path(Dir,DirMSDM))
            DirMSDM<-file.path(Dir,DirMSDM)
          }
          
          DirPRI <- MSDM_Priori_TMLA(occ_xy,envT,MSDM,DirMSDM)
          envM <- stack(paste(DirPRI,list.files(DirPRI),sep="/"))
        }

      #6.2. Data Partition----
      if(Part=="boot"){
        cat("Select the number of replicates (>=1):")
        rep <- as.integer(readLines(n = 1))
        while(is.na(rep)||rep<1){
          warning("Please Select a valid number of replicates (>=1)")
          cat("Select the number of replicates (>=1):")
          rep <- as.integer(readLines(n = 1))
        }
        cat("Select the percentage of occurrences used for fitting the model(0-1):")
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
      
      #6.3.Replicates & Model Input----  
      occINPUT <- list()
        
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
        
      #6.4. Generating Pseudo-Absences----
        #Random Pseudo-Absences
          if(PabM=="rnd"){
            pseudo.mask <- envT[[1]]
            absencesTR <- list()
            absencesTS <- list()
            for(s in 1:length(occTR)){
              set.seed(s)
              absencesTR[[s]] <- randomPoints(pseudo.mask, (1 / PabR)*nrow(occTR[[s]]),ext = extent(pseudo.mask),prob = FALSE)
              absencesTS[[s]] <- randomPoints(pseudo.mask, (1 / PabR)*nrow(occTS[[s]]),ext = extent(pseudo.mask),prob = FALSE)
            }
            absencesTR <- lapply(absencesTR, function(x) cbind(x, rep(1,nrow(x))))
            absencesTR <- lapply(absencesTR, function(x) cbind(x, rep(0,nrow(x))))
            absencesTS <- lapply(absencesTS, function(x) cbind(x, rep(2,nrow(x))))
            absencesTS <- lapply(absencesTS, function(x) cbind(x, rep(0,nrow(x))))
            DirCons <- NULL
          }
          
        #Bioclimatic Constrained Pseudo-Absences
          if(PabM=="const"){
            
            DirCons <- "Constrain"
            if (file.exists(file.path(Dir,DirCons))){
              DirCons<-file.path(Dir,DirCons)
            } else {
              dir.create(file.path(Dir,DirCons))
              DirCons<-file.path(Dir,DirCons)
            }
            
            Model <- lapply(occTR, function(x) bioclim(envT, x[,c("x","y")]))
            pseudo.mask <- lapply(Model, function(x) dismo::predict(x, envT, ext=extent(envT[[1]])))
            pseudo.mask <- lapply(pseudo.mask, function(x) round(x, 5))
            pseudo.mask <-lapply(pseudo.mask, function(x) (x-minValue(x))/
                                (maxValue(x)-minValue(x)))
            pseudo.mask <-lapply(pseudo.mask, function(x) (1-x)>=0.99) #environmental constrain
            writeRaster(stack(pseudo.mask),paste(DirCons,names(pseudo.mask),sep="/"),bylayer=T,format="GTiff",overwrite=T)
            for(i in 1:length(pseudo.mask)){
              pseudo.mask[[i]][which(pseudo.mask[[i]][,]==FALSE)] <- NA
            }
  
            absencesTR <- list()
            absencesTS <- list()
            for (i in 1:length(occTR)) {
              set.seed(i)
              absencesTR.0 <- randomPoints(pseudo.mask[[i]], (1 / PabR)*nrow(occTR[[i]]),
                                         ext = extent(pseudo.mask[[i]]),
                                         prob = FALSE)
              absencesTS.0 <- randomPoints(pseudo.mask[[i]], (1 / PabR)*nrow(occTS[[i]]),
                                           ext = extent(pseudo.mask[[i]]),
                                           prob = FALSE)
              absencesTR[[i]] <- as.data.frame(absencesTR.0)
              absencesTS[[i]] <- as.data.frame(absencesTS.0)
              rm(absencesTR.0,absencesTS.0)
            }
            absencesTR <- lapply(absencesTR, function(x) cbind(x, rep(1,nrow(x)), rep(0,nrow(x))))
            absencesTS <- lapply(absencesTS, function(x) cbind(x, rep(2,nrow(x)), rep(0,nrow(x))))
          }
          
        #Zoogeographical Constrained Pseudo-Absences
          if(PabM=="zoo"){
              print("Select Folder with ZooRegions Mask(.tif):")
              DirZ <- choose.dir(getwd())
              DirZO <- BackZoo_TMLA(Dir=Dir,DirZ=DirZ,occ=occ_xy)
              DirCons <- DirZO
              
              pseudo.mask <- stack(file.path(DirCons,paste(spN,".tif",sep="")))
              pseudo.mask <- unstack(pseudo.mask)
              
              for(i in 1:length(pseudo.mask)){
                pseudo.mask[[i]][which(pseudo.mask[[i]][,]==FALSE)] <- NA
              }
              
              absencesTR <- list()
              absencesTS <- list()
              for (i in 1:length(occTR)) {
                set.seed(i)
                absencesTR.0 <- randomPoints(pseudo.mask[[i]], (1 / PabR)*nrow(occTR[[i]]),
                                             ext = extent(pseudo.mask[[i]]),
                                             prob = FALSE)
                absencesTS.0 <- randomPoints(pseudo.mask[[i]], (1 / PabR)*nrow(occTS[[i]]),
                                             ext = extent(pseudo.mask[[i]]),
                                             prob = FALSE)
                absencesTR[[i]] <- as.data.frame(absencesTR.0)
                absencesTS[[i]] <- as.data.frame(absencesTS.0)
                rm(absencesTR.0,absencesTS.0)
              }
              absencesTR <- lapply(absencesTR, function(x) cbind(x, rep(1,nrow(x)), rep(0,nrow(x))))
              absencesTS <- lapply(absencesTS, function(x) cbind(x, rep(2,nrow(x)), rep(0,nrow(x))))
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
          occT <- rbind(occTR,absencesTR,occTS,absencesTS)
          cols = c("x","y","Partition","PresAbse");    
          occT[,cols] = apply(occT[,cols], 2, function(x) as.numeric(as.character(x)))
          if(rep!=1){
            occT[,"sp"] <- paste(occT[,"sp"],k,sep="_")
            occINPUT[[k]] <- occT
          }else{
            occINPUT <- occT
          }
          
          # occINPUT <- ldply(occINPUT,data.frame,.id=NULL)
          
      #6.6. Define Threshold----
        if(Thr=="LPT"){
          Thresh <- 'no_omission'
        }else{
          Thresh <- 'spec_sens'
        }
        
      #6.6. Define Projection----
        if(Proj=="S"){
          Fut <- EnvF
        }else{
          Fut <- NULL
        }
      
      #6.7. Calculate Moran & MESS
      if(per!=1 || Part=="cross"){
        Bootstrap_Moran_e_MESS_TMLA(Env=envT,RecordsData=occINPUT,DirO=DirR)
      }
          
        #6.8. Run FitENM
          FitENM_TMLA(RecordsData=occINPUT,Variables=envT,Fut=Fut,Part=Part,Algorithm=Alg,PredictType=ENS,
                      Threshold = Thresh,DirSave=DirR,DirMask=NULL,DirMSDM=DirPRI,Save=SavePart,repl=k)
        
          #6.9. Create Occurrence Table for Replicates
          if(rep!=1 || Part=="cross"){
            occTREINO[[k]] <- occTR
            occTESTE[[k]] <- occTS
          }
        }#Fechas as replicas ou kfolds
        
        #6.10.Save Final Occurrence Table & Validation File
        if(rep!=1){
          #Save Final Occurrence Table
          occTREINO <- ldply(occTREINO,data.frame,.id=NULL)
          occTESTE <- ldply(occTESTE,data.frame,.id=NULL)
          write.table(occTREINO,file.path(DirR,"OcorrenciasTreino.txt"),sep="\t",row.names=F)
          write.table(occTESTE,file.path(DirR,"OcorrenciasTeste.txt"),sep="\t",row.names=F)
        
          #Save Final Validation File
          val <- list.files(DirR,pattern="Validation")
          valF <- list()
          for(i in val){
            valF[[i]] <- read.table(file.path(DirR,i),sep="\t",header=T)
          }
          valF <- ldply(valF,data.frame,.id=NULL)
          valF <- valF[order(as.character(valF[,1])),]
          unlink(file.path(DirR,val))
          write.table(valF,file.path(DirR,"Validation.txt"),sep="\t",row.names=F)
          
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
    }
    
#7.MSDM Posteriori----
    
    if(MSDM=="Land"){
      
      cat("Choose L-MSDM type (MaxMin/LowerQ/Pres/MCP/MCPBuffer)")
      Q0 <- as.character(readLines(n = 1))
      while(Q0 %in% c("MaxMin","LowerQ","Pres","MCP","MCPBuffer")==F){
        warning("Choose a valid L-MSDM type!(MaxMin/LowerQ/Pres/MCP/MCPBuffer)")
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
          MSDM_Posterior(RecordsData=occT,Threshold=Thresh,cutoff=Q0,PredictType=ENS,
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
          ENS_Posterior(RecordsData=occT,Algorithm=Alg,PredictType=ENS,Threshold=Thresh,DirAlg=DirT,DirSave=DirR)
        }
      }
    }
}


ENMs_TheMetaLand(Dir="",
                 Sp="",x="",y="",NMin=,PCA="",Proj="",PabR=,PabM="",
                 Part="",SavePart="",Alg="",Thr="",MSDM="",ENS="")
