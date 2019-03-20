#' Create and process Ecological Niche and Species Distribution Models
#'
#' @param pred_dir character. Directory path with predictors (file formats supported are: ASC, BILL, TIFF or TXT)
#' @occ_dir character. Directory path with tab delimited TXT file with species names, latitude and longitude 
#' @param sp character. Name of the column with information about species names
#' @param x character. Name of the column with information about longitude
#' @param y character. Name of the column with information about latitude
#' @param min_occ interger. Minimum number of unique occurrences (species with less than this number will be excluded)
#' @param thin_occ chracater. Perform a spatial filtering (Thinning) on the presences? (Y/N)
#' @param colin_var character. Method to reduce variable collinearity: 
#' \itemize {
#'   \item N: Use original variables.
#'   \item PEARSON: Select variables by Pearson correlation (threshold specified by user).
#'   \item VIF: Variance Inflation Factor (Chatterjee and Hadi 2006).
#'   \item PCA: Perform a Principal Component Analysis on predictors and use Principal Componets as environmental variables
#' }
#' @param imp_var character. Perform importance of variable and curves response for selected algorithms? (Y/N)
#' @param transfer character. Project the model onto another region or time period? (Y/N)
#' @param eval_occ character. Use a pre-determined set of occurrences for validation? (Y/N)
#' @param sp_accessible_area character. Restrict for each species the accessible area,i.e. the area used to construct the model? (Y/N)
#' @param pres_abs_ratio numeric. Presence-Absence ratio (values between 0 and 1)
#' @param pseudoabs_method character. Pseudo-absence allocation method:
#' \itemize{ 
#' \item RND: Random allocation throughout area used to fit models. 
#' \item ENV_CONST: Pseudo-absences are environmentally constrained to region with lower suitability values predicted by a Bioclim model. 
#' \item GEO_CONST: Pseudo-absences are allocated far from occurrences, constrained by a geographical buffer. 
#' \item GEO_ENV_CONST: Pseudo-absences are cosntrained both environmentally (Bioclim Model) and geographically (buffer). 
#' \item GEO_ENV_KM_CONST: Pseudo-absences constrained on a three-level proccedure (Santi). 
#' }
#' @param part character. Partition method for model's validation:
#' \itemize{
#'   \item BOOT: Random bootstrap partition (e.g. 70"%" training and 30"%" test).
#'   \item KFOLD: Random partition in k-fold cross validation.
#'   \item BAND: Geographic partition structured as bands (latitudinal or longitudinal).
#'   \item BLOCK: Geographic partition structured as a checkerboard.
#' }
#' @param save_part character. Save .tif files of the partitions ? (Y/N) (Default="Y")
#' @param save_final character. Save .tif files of the final model (fitted with all the data)? (Y/N) (Default="Y")
#' @param algorithm character. Algorithm to construct ecological niche models: 
#' \itemize{
#'   \item BIO: Bioclim
#'   \item MAH: Mahalanobis
#'   \item DOM: Domain
#'   \item ENF: Ecological Niche Factor Analysis
#'   \item MXS: Maxent Simple (only linear and quadratic features, based on MaxNet package)
#'   \item MXD: Maxent Default (all features, based on MaxNet package)
#'   \item SVM: Support Vector Machine
#'   \item GLM: Generalized Linear Model
#'   \item GAM: Generalizes Additive Model
#'   \item BRT: Boosted Regression Tree
#'   \item RDF: Random Forest
#'   \item MLK: Maximum Likelihood
#'   \item GAU: Gaussian
#' }
#' 
#' @param thr character. Threshold used for presence-absence maps:
#' \itemize{
#'   \item LPT: The highest threshold at which there is no omission
#'   \item MAX_TSS: Threshold at which the sum of the sensitivity and specificity is highest
#'   \item MAX_KAPPA: the threshold at which kappa is highest ("max kappa")
#'   \item SENSITIVITY: fixed (specified) sensitivity
#'   \item JACCARD: the threshold at which Jaccard is highest
#'   \item SORENSEN: the threshold at which Sorensen is highest
#'   }
#'   
#' @param msdm character. Include spatial restrictions. These methods restrict Ecological Niche Models in order to have less potential prediction and turn ENMs closer to species distribution models (SDMs). They are classified in a Priori and a Posteriori methods:
#' 
#' a Priori methods: 
#' 
#' \itemize{
#'   \item N: Do not perform MSDM
#'   \item XY: Create two layers latitude and longitude layer (added as a predictor)
#'   \item MIN: Create a layer with information of the distance from each cell to the closest occurrence (added as a predictor)
#'   \item CML: Create a layer with information of the summed distance from each cell to ALL occurrences (added as a predictor)
#'   \item KER: Create a layer with a Gaussian-Kernel on the occurrence data (added as a predictor)
#'   }
#' a Posteriori methods
#' \itemize{
#'   \item POST: Posterior M-SDM Methods (If chosen, preferred method will be asked later) [NOT added as a predictor]
#'   \item OBR: Occurrence based restriction, uses the distance between points to exclude far suitable patches (Mendes et al, in prep)
#'   \item LR: Lower Quantile, select the nearest 25% patches (Mendes et al, in prep)
#'   \item PRES: Select only the patches with confirmed occurrence data (Mendes et al, in prep)
#'   \item MCP: Excludes suitable cells outside the Minimum Convex Polygon of the occurrence data (Kremen et al, 2008)
#'   \item MCP-B: Creates a Buffer around the MCP (distance defined by user; Kremen et al, 2008)
#'   }
#'
#' @param ensemble character. Method used to ensemble different algorithms:
#'   \itemize{
#'   \item N: No ensemble
#'   \item MEAN: Simple average of the different models
#'   \item W_MEAN: Weighted Average
#'   \item SUP: Average of the best models (TSS over the average)
#'   \item PCA: Performs a Principal Component Analysis (PCA) and returns the first axis
#'   \item PCA_SUP: PCA of the best models (TSS over the average)
#'   \item PCA_THR: PCA only with cells above the threshold
#'   }
#' 
#' @param s_sdm character. Perform a stacked of Species Distribution Model (richness map)? (Y/N)
#'
#' 
#' @references
#'\itemize{
#'\item Kremen, C., Cameron, A., Moilanen, A., Phillips, S. J., Thomas, C. D.,
#'Beentje, H., . Zjhra, M. L. (2008). Aligning Conservation Priorities Across
#'Taxa in Madagascar with High-Resolution Planning Tools. Science, 320(5873),
#'222-226. doi:10.1126/science.1155193
#'}
#'
#'@examples
#' library(ENMTheMetaLand)
#' 
#' 
#' 
#' @export
ENMs_TheMetaLand <- function(pred_dir,
                             occ_dir,
                             sp,
                             x,
                             y,
                             min_occ = 10,
                             thin_occ,
                             colin_var,
                             imp_var,
                             transfer,
                             eval_occ,
                             sp_accessible_area,
                             pres_abs_ratio = 1,
                             pseudoabs_method,
                             part,
                             save_part = "N",
                             save_final = "Y",
                             algorithm,
                             thr,
                             msdm,
                             ensemble,
                             s_sdm) {
  
#1.Check Function Arguments  
  
  er <- NULL
  if(missing(pred_dir)){
    er <- c(er,paste("'pred_dir' unspecified argument, specify the directory of environmental variables | "))
  }
  if(missing(occ_dir)){
    er <- c(er,paste("'occ_dir' unspecified argument, specify the directory of occurrence species data | "))
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
    stop("Argumentos faltantes, please, check  the argumentos listed above")
  }
  
  if(!(colin_var%in%c("PEARSON","VIF","PCA","N"))){
    stop("'colin_var' Argument is not valid!(PEARSON, VIF, PCA, N)")
  }
  if(!(transfer%in%c("Y","N"))){
    stop("'transfer' Argument is not valid!(Y/N)")
  }
  if(pres_abs_ratio<=0){
    stop("'pres_abs_ratio' Argument is not valid!(pres_abs_ratio>=0)")
  }
  if(!(pseudoabs_method%in%c("RND", "ENV_CONST", "GEO_CONST", "GEO_ENV_CONST", "GEO_ENV_KM_CONST"))){
    stop("'pseudoabs_method' Argument is not valid!(Rnd/EnvConst/GeoConst)")
  }
  if(length(pseudoabs_method)>1){
    stop("Please choose only one Pseudo-absence allocation method")
  }
  if(!(part%in%c("BOOT","KFOLD","BANDS","BLOCK"))){
    stop("'part' Argument is not valid!(BOOT/KFOLD/BANDS/BLOCK)")
  }
  if(any(!algorithm%in%c("BIO","MAH","DOM","ENF","GLM","GAM","SVM","BRT","RDF","MXS","MXD","MLK","GAU"))){
    stop(paste("Algorithm",algorithm[!(algorithm%in%c("BIO","MAH","DOM","ENF","GLM","GAM","SVM","BRT","RDF","MXS","MXD","MLK","GAU"))],"is not valid"))
  }
  if(any(!thr%in%c("LPT","MAX_TSS","MAX_KAPPA","SENSITIVITY","JACCARD","SORENSEN"))){
    stop("'thr' Argument is not valid!")
  }
  if(!(msdm%in%c("N","XY","MIN","CML","KER","POST"))){
    stop("'msdm' Argument is not valid!(N/XY/MIN/CML/KER/POST)")
  }
  if(length(msdm)>1){
    stop("Please choose only one 'msdm' method")
  }
  
#1.Load Packages ----
  
  ipak <- function(pkg) {
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
  
#2.Adjust Names----
  Ord <- c("BIO","DOM","MAH","ENF","MXD","MXS","MLK","SVM","RDF","GAM","GLM","GAU","BRT")
  algorithm <- Ord[Ord%in%algorithm]
  
  ThrNames <- c("LPT", "MAX_TSS", "MAX_KAPPA", "SENSITIVITY","JACCARD", "SORENSEN")
  thr <- ThrNames[c("LPT","MAX_TSS","MAX_KAPPA","SENSITIVITY","JACCARD", "SORENSEN")%in%thr]
  rm(ThrNames)
  
#3.Predictors ----
  options(warn = -1)
  setwd(pred_dir)
  
  env <- unique(file_ext(list.files()))
  form <- c('bil', 'asc', 'txt', 'tif')
  env <- env[env %in% form]
  if (length(env) > 1) {
    stop("More than one file format in pred_dir")
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
  if(colin_var=="VIF") {
    VF <- vifstep(envT, th = nlayers(envT) * 2)
    envT <- exclude(envT, VF)
    if (transfer == "Y") {
      RasM <- colMeans(na.omit(values(envT)))
      RasSTD <- apply(na.omit(values(envT)), 2, std)
    }
    envT <- scale(envT)
    
    if(transfer=="Y") {
      EnvF <- list()
      for (i in 1:length(Pfol)) {
        ProjT <- unique(file_ext(list.files(Pfol[i])))
        form <- c('bil', 'asc', 'txt', 'tif')
        ProjT <- ProjT[ProjT %in% form]
        if (length(ProjT) > 1) {
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
  
  if (colin_var=="PCA") {
      
    #Projection PCA
    if(transfer=="Y"){
      EnvF <- list()
      for(i in 1:length(Pfol)){
        EnvF[[i]] <- PCAFuturo(Env=envT,Dir=pred_dir,DirP=Pfol[i],Save="Y")
      }
      names(EnvF) <- PfolN
      envT <- brick(stack(file.path(pred_dir,"PCA",list.files(file.path(pred_dir,"PCA"),pattern='PC'))))
    }else{
      envT<-PCA_env_TMLA(env = envT, Dir = pred_dir)
    }
  }
  
  #3.3.3.Pearson----
  if(colin_var=="PEARSON"){
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
  
  # Read txt with occurences data
  occ <- read.table(occ_dir,
                    h = T,
                    sep = '\t',
                    stringsAsFactors = F)
  occ<-occ[,c(sp,x,y)]
  colnames(occ) <- c("sp","x","y")
  occ_xy <- split(occ[,-1],f=occ$sp)
  spN<-names(occ_xy)
  
  
    #4.1.Unique Occurrences----
    occA<-Occ_Unicas_TMLA(env=envT[[1]], occ.xy=occ_xy, DirO=DirR)

    #4.2.Thining----
    if(thin_occ=="Y"){
      cat(("Select thinning method:\n1-Distance defined by Moran Variogram\n2-User defined distance\n3-Distance defined by 2x cellsize (Haversine Transformation)"))
      ThinMet <- as.integer(readLines(n=1))
      occA <- OccsThin(occA, envT, ThinMet, colin_var, DirR)
    }
  
    #4.3.Save Thinned & Unique Occurrences
    ndb <- ldply(occA)[,1:3]
    write.table(ndb,file.path(DirR,"Occ_Cleaned.txt"),sep="\t",row.names=F)
  
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
      cat("Select restriction type (buffer / mask):")
      method <- as.character(readLines(n = 1))
      while(!method%in%c("buffer","mask")){
        warning("Please Select a valid restriction type (buffer / mask)")
        cat("Select restriction type (buffer / mask):")
        method <- as.character(readLines(n = 1))
      }
      DirM <- M_delimited(var=envT,
                  occ_xy=occ_xy,
                  method = method,
                  BufferDistanceKm=NULL,
                  EcoregionsFile=NULL,
                  Dir=pred_dir,
                  spN=spN,
                  SaveM = TRUE)
    }
  
    if(grepl("GEO", pseudoabs_method)){
      #Define Buffer distance:
      cat("Select buffer distance(in km):")
      Geo_Buf <- as.integer(readLines(n = 1))*1000
    }
    
#6. Geographical Partition----
    if(part=="BANDS" || part=="BLOCK"){
      
      if(any(grepl("PC",names(envT)))==T || any(grepl("pc",names(envT)))==T){
        colin_var<-"PCA"
      }
      
      if(eval_occ=="Y"){
        warning("Invalid combination! eval_occ can't be Y with Geographical partitions! Changing eval_occ to N")
        eval_occ <- "N"
      }
      
      if(part=="BANDS"){  
        #6.1.Bands----
        
        DirB<-"BANDS"
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
            envTT<-PCA_env_TMLA(env = envT, Dir = pred_dir)
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
          TipoMoran <- "all"
          # print("Select Moran Calculation Type (all/nearest):")
          # TipoMoran <- as.character(readLines(n = 1))
          # while(is.na(TipoMoran)||!(TipoMoran%in%c("all","nearest"))){
          #   warning("Please choose a valid Moran Calculation Type [all/nearest]")
          #   print("Select Moran Calculation Type (all/nearest):")
          #   TipoMoran <- as.character(readLines(n = 1))
          # }
          
          #Check for M-Restriction
          if(exists("DirM")){
            DirM <- DirM
          }else{
            DirM <- NULL
          }
          occINPUT <-
            BandsPartition_TMLA(
              evnVariables = envTT,
              RecordsData = occ_xy,
              N = bands,
              pseudoabsencesMethod = pseudoabs_method,
              PrAbRatio = pres_abs_ratio,
              DirSave = DirB,
              DirM = DirM,
              MRst = sp_accessible_area,
              type = TipoMoran
            )
          occINPUT[,4] <- as.numeric(occINPUT[,4])
          occINPUT[,5] <- as.numeric(occINPUT[,5])
          rm(envTT)
        }

      }
      if(part=="BLOCK"){
        #6.2.Block----
        
        DirB<-"BLOCK"
        if (file.exists(file.path(DirR, DirB))){
          DirB<-file.path(DirR,DirB)
        } else {
          dir.create(file.path(DirR,DirB))
          DirB<-file.path(DirR,DirB)
        }
        
        if(all(paste0(spN,".tif")%in%list.files(DirB,pattern=".tif"))){
          print("Partition Already Exist! Using pre-created partitions! ")
          setwd(DirB)
          occINPUT <- read.table(file.path(DirB,"OccCheck.txt"),sep="\t",header=T)
          occINPUT[,4] <- as.numeric(occINPUT[,4])
          occINPUT[,5] <- as.numeric(occINPUT[,5])
        }else{
          if(colin_var!="PCA"){
            envTT<-PCA_env_TMLA(env = envT, Dir = pred_dir)
          }else{
            envTT<-envT
          }
          TipoMoran <- "all"
          # print("Select Moran Calculation Type (all/nearest):")
          # TipoMoran <- as.character(readLines(n = 1))
          # while(is.na(TipoMoran)||!(TipoMoran%in%c("all","nearest"))){
          #   warning("Please choose a valid Moran Calculation Type [all/nearest]")
          #   print("Select Moran Calculation Type (all/nearest):")
          #   TipoMoran <- as.character(readLines(n = 1))
          # }
          
          #Check for M-Restriction
          if(exists("DirM")){
            DirM <- DirM
          }else{
            DirM <- NULL
          }
          
          occINPUT <-
            BlockPartition_TMLA(
              evnVariables = envTT,
              RecordsData = occ_xy,
              N = 2,
              pseudoabsencesMethod = pseudoabs_method,
              PrAbRatio = pres_abs_ratio,
              DirSave = DirB,
              DirM = DirM,
              MRst = sp_accessible_area,
              type = TipoMoran,
              Geo_Buf=Geo_Buf
            )
          
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
        if (file.exists(file.path(pred_dir,DirMSDM))){
          DirMSDM<-file.path(pred_dir,DirMSDM)
        } else {
          dir.create(file.path(pred_dir,DirMSDM))
          DirMSDM<-file.path(pred_dir,DirMSDM)
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
      FitENM_TMLA_Parallel(
        RecordsData = occINPUT,
        Variables = envT,
        VarImP = imp_var,
        Fut = Fut,
        Part = part,
        Algorithm = algorithm,
        PredictType = ensemble,
        spN = spN,
        Tst = eval_occ,
        Threshold = thr,
        DirSave = DirR,
        DirMask = DirB,
        DirMSDM = DirPRI,
        Save = save_part,
        SaveFinal = save_final,
        repl = NULL,
        per = NULL
      )
    }
    
#7.Random Partition----
    
    if(part=="BOOT"||part=="KFOLD"){
      
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
          if (file.exists(file.path(pred_dir,DirMSDM))){
            DirMSDM<-file.path(pred_dir,DirMSDM)
          } else {
            dir.create(file.path(pred_dir,DirMSDM))
            DirMSDM<-file.path(pred_dir,DirMSDM)
          }
          
          DirPRI <- MSDM_Priori_TMLA(Species=occ_xy,var=envT,MSDM=msdm,DirMSDM=DirMSDM)
        }

      #7.2. Data Partition----
      if(part=="BOOT"){
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
      if(part=="KFOLD"){
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
      if(eval_occ=="Y" && part=="BOOT" && per!=1 || eval_occ=="Y" && part=="KFOLD" && rep!=1){
        if(part=="BOOT"){
          warning("Adjusting data partition to one!")
          per <- 1
        }
        if(part=="KFOLD"){
          warning("Adjusting partition to bootstrap and data partition to one! 
          Replicates will be equal to the original number of folds")
          part <- "BOOT"
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
        if(part=="BOOT"){
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
        
        if(part=="KFOLD"){
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
        # Pseudo-Absences with Random allocation-----
        if(pseudoabs_method=="RND"){
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
          
        # Pseudo-Absences allocation with Environmental constrain ----
        if(pseudoabs_method=="ENV_CONST"){
            
            DirCons <- "EnvConstrain"
            if (file.exists(file.path(pred_dir,DirCons))){
              DirCons<-file.path(pred_dir,DirCons)
            } else {
              dir.create(file.path(pred_dir,DirCons))
              DirCons<-file.path(pred_dir,DirCons)
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
                pseudo.mask <- inverse_bio(envT, occTR[[i]][,c("x","y")])
                
                writeRaster(pseudo.mask,paste(DirCons,spN[i],sep="/"),format="GTiff",overwrite=T)
                
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
        
        # Pseudo-Absences allocation with Geographical constrain-----
        if(pseudoabs_method=="GEO_CONST"){
          
          DirCons <- "GeoConstrain"
          if (file.exists(file.path(pred_dir,DirCons))){
            DirCons<-file.path(pred_dir,DirCons)
          } else {
            dir.create(file.path(pred_dir,DirCons))
            DirCons<-file.path(pred_dir,DirCons)
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
              
              pseudo.mask <- inv_geo(e=envT, p=occTR[[i]][,c("x","y")], d=Geo_Buf)
              writeRaster(pseudo.mask,paste(DirCons,spN[i],sep="/"),format="GTiff",overwrite=T)
            
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
          
        # Pseudo-Absences allocation with Environmentla and Geographical  constrain-----
        if(pseudoabs_method=="GEO_ENV_CONST"){
          DirCons <- "GeoEnvConstrain"
          if (file.exists(file.path(pred_dir,DirCons))){
            DirCons<-file.path(pred_dir,DirCons)
          } else {
            dir.create(file.path(pred_dir,DirCons))
            DirCons<-file.path(pred_dir,DirCons)
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
              
              pseudo.mask <- inv_geo(e=envT, p=occTR[[i]][,c("x","y")], d=Geo_Buf)
              pseudo.mask1 <- inverse_bio(envT, occTR[[i]][,c("x","y")])
              pseudo.mask <- pseudo.mask*pseudo.mask1
              writeRaster(pseudo.mask,paste(DirCons,spN[i],sep="/"),format="GTiff",overwrite=T)
              
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
        
        # Pseudo-Absences allocation with Environmentla and Geographical  constrain-----
        if(pseudoabs_method=="GEO_ENV_KM_CONST"){
          DirCons <- "GeoEnvConstrain"
          if (file.exists(file.path(pred_dir,DirCons))){
            DirCons<-file.path(pred_dir,DirCons)
          } else {
            dir.create(file.path(pred_dir,DirCons))
            DirCons<-file.path(pred_dir,DirCons)
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
              
              pseudo.mask <- inv_geo(e=envT, p=occTR[[i]][,c("x","y")], d=Geo_Buf)
              pseudo.mask1 <- inverse_bio(envT, occTR[[i]][,c("x","y")])
              pseudo.mask <- pseudo.mask*pseudo.mask1
              writeRaster(pseudo.mask,paste(DirCons,spN[i],sep="/"),format="GTiff",overwrite=T)
              
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
              
              absencesTR.0 <- KM(rasterToPoints(SpMask)[,-3],
                 mask(evnVariables, SpMask),
                 ((1 / pres_abs_ratio)*nrow(occTR[[i]])))
              absencesTS.0 <- KM(rasterToPoints(SpMaskP)[,-3],
                 mask(evnVariables, SpMaskP),
                 ((1 / pres_abs_ratio)*nrow(occTS[[i]])))
            }else{
              absencesTR.0 <- KM(cell_coord = rasterToPoints(pseudo.mask)[,-3],
                                 variable = mask(evnVariables, pseudo.mask),
                                 NumAbsence = (1 / pres_abs_ratio)*nrow(occTR[[i]]))
              if(eval_occ=="Y"){
                absencesTS.0 <- KM(rasterToPoints(pseudo.maskP)[,-3],
                                   mask(evnVariables, pseudo.mask),
                                   (1 / pres_abs_ratio)*nrow(occTS[[i]]))
              }else{
                absencesTS.0 <- KM(rasterToPoints(pseudo.maskP)[,-3],
                                   mask(evnVariables, pseudo.mask),
                                   (1 / pres_abs_ratio)*nrow(occTS[[i]]))
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
      if(per!=1 || part=="KFOLD"){
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
        FitENM_TMLA_Parallel(RecordsData=occINPUT,Variables=envT,VarImP=imp_var,Fut=Fut,Part=part,Algorithm=algorithm,PredictType=ensemble,spN=spN,
                    Tst=eval_occ,Threshold=thr,DirSave=DirR,DirMask=DirB,DirMSDM=DirPRI,Save=save_part,
                    SaveFinal=save_final,per=per,repl=k)
        
      #7.10. Create Occurrence Table for Replicates----
        if(rep!=1 || part=="KFOLD"){
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
          if(per!=1 || part=="KFOLD"){
            Boot <- list.files(DirR,pattern="Bootstrap_Moran_MESS")
            BootF <- list()
            for(i in Boot){
              BootF[[i]] <- read.table(file.path(DirR,i),sep="\t",header=T)
            }
            BootF <- ldply(BootF,data.frame,.id=NULL)
            BootF <- BootF[order(as.character(BootF[,1])),]
            unlink(file.path(DirR,Boot))
            if(part=="BOOT"){
              write.table(BootF,file.path(DirR,"Bootstrap_Moran_MESS.txt"),sep="\t",row.names=F)
            }
            if(part=="KFOLD"){
              write.table(BootF,file.path(DirR,"CrossValidation_Moran_MESS.txt"),sep="\t",row.names=F)
            }
          }
        }
    }#Fecha partition BOOT|KFOLD
    
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
