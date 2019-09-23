#' Create and process Ecological Niche and Species Distribution Models
#'
#' @param pred_dir character. Directory path with predictors (file formats supported are: ASC, BILL, TIFF or TXT)
#'
#' @param proj_dir character. Directory path containing folders with predictors for different regions or time periods used to project models (file formats supported are: ASC, BILL, TIFF or TXT).
#'
#' @param occ_file character. Directory path with tab-delimited TXT file with species names, latitude and longitude
#'
#' @param sp character. Name of the column with information about species names
#'
#' @param x character. Name of the column with information about longitude

#' @param y character. Name of the column with information about latitude
#'
#' @param min_occ integer. Minimum number of unique occurrences (species with less than this number will be excluded)
#'
#' @param thin_occ character Default NULL. Perform spatial filtering (Thinning, based on spThin package) on the presences. For this augment it is necessary provide a vector in which its elements need to have the names 'method' or 'method' and 'distance' (more information below). Three thinning methods are available:
#' \itemize{
#' \item MORAN-Distance defined by Moran Variogram, usage thin_occ=c(method='MORAN').
#' \item CELLSIZE-Distance defined by 2x cellsize (Haversine Transformation), usage thin_occ=c(method='CELLSIZE').
#' \item USER-DEFINED-User defined distance. For this option it is neede provide a vector with two values. Usage thin_occ=c(method='USER-DEFINED', ditance='300'). The second numeric value refers to the distance in km that will be used for thinning. So distance=300 means that all records within a radius of 300 km will be deleted
#' }
#'
#' @param eval_occ character. Directory path with tab-delimited TXT file with species names, latitude and longitude, this three columns must have the same columns names than the databased used in the occ_file argument. This is an external occurrence database that will be used to external models validation. (default NULL)
#'
#' @param colin_var character. Method to reduce variable collinearity:
#' \itemize{
#'   \item PCA: Perform a Principal Component Analysis on predictors and use Principal Components as environmental variables, usage colin_var=c(method='PCA').
#'   \item VIF: Variance Inflation Factor (Chatterjee and Hadi 2006), , usage colin_var=c(method='VIF').
#'   \item PEARSON: Select variables by Pearson correlation, a threshold of maximum correlation must be specified by user, usage colin_var=c(method='PEARSON', threshold='0.7').
#' }
#'
#' @param imp_var character. Perform variable importance and curves response for selected algorithms? (Y/N)
#'
#' @param sp_accessible_area character. Restrict for each species the accessible area, i.e., the area used to construct the model. It is necessary to provide a vector for this argument. Three methods were implemented
#' \itemize{
#'   \item BUFFER and based on maximum distance among pair of occurrences for each species. Usage sp_accessible_area=c(method='BUFFER', type='1').
#'   \item BUFFER and based on A single buffer for all species expressed in km. Usage sp_accessible_area=c(method='BUFFER', type='2', width='300').
#'   \item MASK: this method consist in delimit the area used to model calibration based on the polygon where a species occurrences fall. For instance it is possible delimit the calibration area based on ecorregion shapefile. For this option it is necessary inform the path to the file that will be used as mask. Next file format can be loaded '.bil', '.asc', '.tif', '.shp', and '.txt'. Usage sp_accessible_area=c(method='MASK', filepath='C:/Users/mycomputer/ecoregion/olson.shp').
#' }
#'
#' @param pseudoabs_method character. Pseudo-absence allocation method. It is necessary to provide a vector for this argument. Only one method can be chosen. The next methods are implemented:
#' \itemize{
#' \item RND: Random allocation throughout area used to fit models. Usage pseudoabs_method=c(method='RND').
#' \item ENV_CONST: Pseudo-absences are environmentally constrained to a region with lower suitability values predicted by a Bioclim model. Usage pseudoabs_method=c(method='ENV_CONST').
#' \item GEO_CONST: Pseudo-absences are allocated far from occurrences based on a geographical buffer. For this method it is necessary provie a second value wich express the buffer width in km. Usage pseudoabs_method=c(method='GEO_CONST', width='50').
#' \item GEO_ENV_CONST: Pseudo-absences are constrained environmentally (based on Bioclim model) but distributed geographically far from occurrences based on a geographical buffer. For this method it is necessary provie a second value wich express the buffer width in km. Usage pseudoabs_method=c(method='GEO_ENV_CONST', width='50')
#' \item GEO_ENV_KM_CONST: Pseudo-absences are constrained on a three-level procedure; it is similar to the GEO_ENV_CONST with an additional step which distributes the pseudo-absences in the environmental space using k-means cluster analysis. For this method it is necessary provie a second value wich express the buffer width in km. Usage pseudoabs_method=c(method='GEO_ENV_KM_CONST', width='50')
#' }
#'
#' @param pres_abs_ratio numeric. Presence-Absence ratio (values between 0 and 1)

#' @param part character. Partition method for model's validation. Only one method can be chosen. It is necessary to provide a vector for this argument. The next methods are implemented:
#' \itemize{
#'   \item BOOT: Random bootstrap partition (e.g. 70 % training and 30 % test). Usage •	part=c(method='BOOT', replicates='2',  proportion='0.7'). 'replicate' refers to the number of replicates, it assumes a value >=1. 'proportion' refres to the proportion of occurrences used for fitting the model, it assumes a value >0 and <=1.
#'   \item KFOLD: Random partition in k-fold cross-validation. Usage part=c(method= 'KFOLD', folds='5'). 'folds' referes to the number of k-folds for patitioning, it assume value >=1.
#'   \item BANDS: Geographic partition structured as bands (latitudinal(1) or longitudinal(2)). Usage part=c(method= 'BANDS', type='1'). 'type' refers to the bands disposition
#'   \item BLOCK: Geographic partition structured as a checkerboard. Usage part=c(method= 'BLOCK').
#' }
#'
#' @param save_part logical. If TRUE, function will save .tif files of partial models, i.e. model created by each occurrence partitions. (default FALSE).
#'
#' @param save_final logical. If TRUE, function will Save .tif files of the final model, i.e. fitted with all occurrences data. (default TRUE)
#'
#' @param algorithm character. Algorithm to construct ecological niche models (it is possible to use more than one method):
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
#'   \item GAU: Gaussian Process
#' }
#'
#' @param thr character. Threshold used for presence-absence predictions. It is possible to use more than one threshol type. It is necessary to provide a vector for this argument:
#' \itemize{
#'   \item LPT: The highest threshold at which there is no omission. Usage thr=c(type='LPT').
#'   \item MAX_TSS: Threshold at which the sum of the sensitivity and specificity is highest.
#'   Usage thr=c(type='MAX_TSS').
#'   \item MAX_KAPPA: The threshold at which kappa is highest ("max kappa"). Usage thr=c(type='MAX_KAPPA').
#'   \item SENSITIVITY: Fixed (specified) sensitivity. For this type of threshold thr must be use as thr=c(type='SENSITIVITY', sens='0.6'). 'sens' refers to models will be bynarized using this suitability value. Note that this method assume 'sens' value for all algorithm.
#'   \item JACCARD: The threshold at which Jaccard is highest. Usage thr=c(type='JACCARD').
#'   \item SORENSEN: The threshold at which Sorensen is highest. Usage thr=c(type='SORENSEN').
#'   }
#' In the case of use more than one threshold type it is necessary concatenate the names of threshold types, e.g., thr=c(type=c('LPT', 'MAX_TSS', 'JACCARD')). In the case of SENSITIVITY threshold is used it is necessayr specify the desired sensitivity value, e.g. thr=c(type=c('LPT', 'MAX_TSS', 'SENSITIVITY'), sens='0.8')
#'
#' @param msdm character. Include spatial restrictions. These methods restrict ecological niche models in order to have less potential prediction and turn models closer to species distribution models. They are classified in a Priori and a Posteriori methods. The firt one encompase method taht include geographical layers as predictor of models construction, whereas a Posteriori constrain models based on occurrence and suitability pattenrs. This argument is filled only wiht a method, in the case of use MCP-B method msdm is filled in a different way se below:
#'
#' a Priori methods:
#'
#' \itemize{
#'   \item XY: Create two layers latitude and longitude layer (added as a predictor).Usage msdm=c(method='XY').
#'   \item MIN: Create a layer with information of the distance from each cell to the closest occurrence (added as a predictor).Usage msdm=c(method='MIN').
#'   \item CML: Create a layer with information of the summed distance from each cell to ALL occurrences (added as a predictor).Usage msdm=c(method='CML').
#'   \item KER: Create a layer with a Gaussian-Kernel on the occurrence data (added as a predictor).Usage msdm=c(method='KER').
#'   }
#' a Posteriori methods
#' \itemize{
#'   \item OBR: Occurrence based restriction, uses the distance between points to exclude far suitable patches (Mendes et al., in prep).Usage msdm=c(method='OBR').
#'   \item LR: Lower Quantile, select the nearest 25% patches (Mendes et al., in prep).Usage msdm=c(method='LR').
#'   \item PRES: Select only the patches with confirmed occurrence data (Mendes et al, in prep).Usage msdm=c(method='PRES').
#'   \item MCP: Excludes suitable cells outside the Minimum Convex Polygon of the occurrence data.Usage msdm=c(method='MCP').
#'   \item MCP-B: Creates a Buffer around the MCP (distance defined by user in km). Usage msdm=c(method='MCP-B', width=100).
#'   }
#'
#' @param ensemble character. Method used to ensemble different algorithms. It is possible to use more than one method. It is necessary to provide a vector for this argument.  It is possible to use more than one ensemble method. For SUP, W_MEAN or PCA_SUP method it is necesary provide an evaluation metric to ensemble arguemnts (i.e. AUC, Kappa, TSS, Jaccard, Sorensen or Fpb) see below. (default NULL):
#'   \itemize{
#'   \item MEAN: Simple average of the different models. Usage ensemble=c(method='MEAN').
#'   \item W_MEAN: Weighted average. Usage ensemble=c(method='W_MEAN').
#'   \item SUP: Average of the best models (TSS over the average). Usage ensemble=c(method='SUP').
#'   \item PCA: Performs a Principal Component Analysis (PCA) and returns the first axis. Usage ensemble=c(method='PCA').
#'   \item PCA_SUP: PCA of the best models (TSS over the average). For this metdhos it is necesary provide an evaluation metric. Usage ensemble=c(method='PCA_SUP', metric='Fpd').
#'   \item PCA_THR: PCA only with cells above the threshold. Usage ensemble=c(method='PCA_THR').
#'   }
#'
#'  In the case of use more than one ensemble method it is necessary concatenate the names of ensemble mehtods within the argument, e.g., ensemble=c(type=c('MEAN', 'PCA')), ensemble=c(method=c('MEAN, 'W_MEAN', 'PCA_SUP'), metric='Fpb')
#' @param extrapolation logical. If TRUE the function will calculate extrapolation based on Mobility-oriented parity (MOP) method, for current and future conditions.
#' @param cores numeric. Define the number number of CPU cores to run modeling procedures in parallel.
#'
#'
#'@examples
#' library(ENMTML)
#'
#'
#'
#' @export
#'
ENMTML <- function(pred_dir,
                   proj_dir=NULL,
                   occ_file,
                   sp,
                   x,
                   y,
                   min_occ = 10,
                   thin_occ=NULL,
                   eval_occ=NULL,
                   colin_var=NULL,
                   imp_var=NULL,
                   sp_accessible_area=NULL,
                   pseudoabs_method,
                   pres_abs_ratio = 1,
                   part,
                   save_part = FALSE,
                   save_final = TRUE,
                   algorithm,
                   thr,
                   msdm=NULL,
                   ensemble=NULL,
                   extrapolation=FALSE,
                   cores=1) {

  #1.Check Function Arguments
  cat("Checking for function arguments ...\n")

  er <- NULL
  if(missing(pred_dir)){
    er <- c(er,paste("'pred_dir' unspecified argument, specify the directory of environmental variables | "))
  }
  if(missing(occ_file)){
    er <- c(er,paste("'occ_file' unspecified argument, specify the directory of occurrence species data | "))
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
  if(is.null(ensemble)) {
    ensemble2 <- "N"
  } else{
    ensemble2 <- ensemble[grep('method', names(ensemble))]
    if (any(ensemble2 %in% c("SUP", "W_MEAN", "PCA_SUP")) & !any(names(ensemble) == 'metric')) {
      stop("If you used SUP, W_MEAN or PCA_SUP ensemble method it is necessay provide an evaluation metric to ensemble arguemnt (AUC, Kappa, TSS, Jaccard, Sorensen or Fpb). \n e.g., ensemble=c(method=c('W_MEAN', 'PCA_SUP'), metric='Fpb')")
    } else if (any(ensemble2 %in% c("SUP", "W_MEAN", "PCA_SUP"))) {
      ensemble_metric <- ensemble[grep('metric', names(ensemble))]
    }
    if (!any(ensemble2 %in% c("SUP", "W_MEAN", "PCA_SUP"))) {
      ensemble_metric <- NULL
    }
  }

  if(!is.null((er))){
    print(er)
    stop("Missing arguments. Please check the argumentos listed above")
  }

  if(!is.null(thin_occ)){
    if(!(thin_occ['method']%in%c('MORAN','CELLSIZE','USER-DEFINED'))){
      stop("'thin_occ' Argument is not valid! a character vector is nedeed, e.g., thin_occ=c(method='MORAN')")
    }
    if(thin_occ['method']=="USER-DEFINED"&length(thin_occ)!=2){
      stop("'thin_occ' Argument is not valid for method=USER-DEFINED! a distance in km must be provided e.g., thin_occ=c(method='USER-DEFINED', distance='300')")
    }
  }

  if(!is.null(colin_var)){
    if(!(colin_var['method']%in%c("PEARSON","VIF","PCA","N"))){
      stop("'colin_var' Argument is not valid!(PEARSON, VIF, PCA, N)")
    }
    if(colin_var['method']==c("PEARSON")&length(colin_var)!=2){
      stop("'colin_var' Argument is not valid for PEARSON method! a threshold between 0-1 is needed e.g., colin_var=c(method='PEARSON', threshold='0.7')")
    }
  }


  if(pres_abs_ratio<=0){
    stop("'pres_abs_ratio' Argument is not valid!(pres_abs_ratio>=0)")
  }
  if(!(pseudoabs_method['method']%in%c("RND", "ENV_CONST", "GEO_CONST", "GEO_ENV_CONST", "GEO_ENV_KM_CONST"))){
    stop("'pseudoabs_method' Argument is not valid!. Can be used: 'RND', 'ENV_CONST', 'GEO_CONST', 'GEO_ENV_CONST', 'GEO_ENV_KM_CONST')")
  }
  if(length(pseudoabs_method['method'])>1){
    stop("Please choose only one Pseudo-absence allocation method")
  }
  if(!(part['method']%in%c("BOOT","KFOLD","BANDS","BLOCK"))){
    stop("'part' Argument is not valid!(BOOT/KFOLD/BANDS/BLOCK)")
  }
  if(any(!algorithm%in%c("BIO","MAH","DOM","ENF","GLM","GAM","SVM","BRT","RDF","MXS","MXD","MLK","GAU"))){
    stop(paste("Algorithm", algorithm[!(algorithm%in%c("BIO","MAH","DOM","ENF","GLM","GAM","SVM","BRT","RDF","MXS","MXD","MLK","GAU"))],"is not valid"))
  }
  if(any(!thr[grep('type', names(thr))]%in%c("LPT","MAX_TSS","MAX_KAPPA","SENSITIVITY","JACCARD","SORENSEN"))){
    stop("'thr' Argument is not valid!")
  }
  if(any(thr[grep('type', names(thr))]%in%"SENSITIVITY") && !any(names(thr)%in%'sens')){
    stop("provide a sensitivity value in the vector used in 'thr' argument, see ENMTML function help")
  }
  if(!is.null(msdm)){
    if(length(msdm)==2){
      msdm <- msdm['method']
      msdm_width <- as.numeric(msdm['width'])
    }
    if(!(msdm['method']%in%c('XY','MIN','CML','KER', 'OBR', 'LR', 'PRES', 'MCP', 'MCP-B'))){
      stop("'msdm' Argument is not valid!(XY/MIN/CML/KER/OBR/LR/PRES/MCP/MCP-B)")
    }
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
         "usdm","pracma","gbm","caret","adehabitatHS", "visreg"))

  #2.Adjust Names----
  Ord <- c("BIO","DOM","MAH","ENF","MXD","MXS","MLK","SVM","RDF","GAM","GLM","GAU","BRT")
  algorithm <- Ord[Ord%in%algorithm]

  #3.Predictors ----
  cat("Loading environmental variables ...\n")

  options(warn = 1)
  setwd(pred_dir)

  env <- unique(file_ext(list.files()))
  form <- c('bil', 'asc', 'txt', 'tif')
  env <- env[env %in% form]
  if (length(env) > 1) {
    stop("More than one file format in pred_dir")
  }

  if(any(env == c('asc', 'bil', 'tif'))){
    envT<-stack(list.files(pattern=paste0('\\.',env,'$')))
    try(envT <- brick(envT))
    if(class(envT)=="RasterBrick"){
      print("RasterBrick successfully created!")
    }else{
      print("Failed to create RasterBrick, using Raster Stack instead!")
    }
  }
  if(env == 'txt'){
    envT<-read.table(list.files(pattern='\\.txt$'),h=T)
    gridded(envT)<- ~x+y
    envT<-stack(envT)
    try(envT <- brick(envT))
    if(class(envT)=="RasterBrick"){
      print("RasterBrick successfully created!")
    }else{
      print("Failed to create RasterBrick, using Raster Stack instead!")
    }
  }

  #CRS
  if(is.na(crs(envT))){
    crs(envT) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  }

  #3.0.Check predictors consistency
  if(length(unique(colSums(!is.na(envT[]))))>1){
    envT[is.na((sum(envT))[])] <- NA
    print("Variables had differences, setting any empty cells to NA in all variables")
  }

  #3.1.Projection----
  if(!is.null(proj_dir)){
    DirP<-proj_dir
    Pfol<-file.path(DirP,list.files(DirP))
    if(any(file_ext(list.files(DirP))%in%form)){
      stop("Select a folder containing folders with environment conditions for different regions or time periods, NOT a folder with this variables!")
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
  if(!is.null(colin_var)){
    #3.1.1.VIF----
    if(colin_var%in%"VIF") {
      cat("Performing a reduction of variables collinearity ...\n")
      VF <- vifstep(envT, th = 10)
      envT <- exclude(envT, VF)
      if (!is.null(proj_dir)) {
        RasM <- colMeans(na.omit(values(envT)))
        RasSTD <- apply(na.omit(values(envT)), 2, std)
      }
      envT <- raster::scale(envT)

      if(!is.null(proj_dir)) {
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
      cat("Performing a reduction of variables collinearity ...\n")
      #Projection PCA
      if(!is.null(proj_dir)){
        EnvF <- PCAFuturo(Env=envT,Dir=pred_dir,DirP=Pfol,Save="Y")
        names(EnvF) <- PfolN
        envT <- brick(stack(list.files(file.path(pred_dir,"PCA"),pattern='PC',full.names = T)))
      }else{
        envT<-PCA_env_TMLA(env = envT, Dir = pred_dir)
      }
    }

    #3.3.3.Pearson----
    if(colin_var['method']=="PEARSON"){
      cat("Performing a reduction of variables collinearity ...\n")
      cat("Select correlation threshold:(0-1)\n")
      Cor_TH <- as.numeric(colin_var["threshold"])
      Pear <- layerStats(envT, 'pearson', na.rm=T)
      corr_matrix <- abs(Pear$'pearson correlation coefficient')
      corr_matrix[upper.tri(corr_matrix)] <- 0
      diag(corr_matrix) <- 0
      envT <- envT[[names(envT)[!apply(corr_matrix,2,function(x) any(x > 0.70))]]]
      if(!is.null(proj_dir)){
        RasM <- colMeans(na.omit(values(envT)))
        RasSTD <- apply(na.omit(values(envT)),2,std)
      }
      envT <- scale(envT)

      if(!is.null(proj_dir)){
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
  }else{
  #3.3.4.colin_var='N'----
    if (!is.null(proj_dir)) {
      EnvF <- list()
      for (i in 1:length(Pfol)) {
        ProjT <- unique(file_ext(list.files(Pfol[i])))
        form <- c('bil', 'asc', 'txt', 'tif')
        ProjT <- ProjT[ProjT %in% form]
        if (length(ProjT) > 1) {
          stop("More than one file format in DirP")
        }
        if (any(ProjT == c('asc', 'bil', 'tif'))) {
          EnvF[[i]] <-
            brick(stack(file.path(
              Pfol[i], list.files(Pfol[i], pattern = paste0('\\.', ProjT, '$'))
            )))
        }
        if (ProjT == 'txt') {
          ProjTT <-
            read.table(file.path(Pfol[i], list.files(Pfol[i], pattern = '\\.txt$'), h =
                                   T))
          gridded(ProjTT) <- ~ x + y
          EnvF[[i]] <- brick(stack(ProjTT))
          rm(ProjTT)
        }
      }
      names(EnvF) <- PfolN
    }
  }

  #3.3.Erro Futuro e msdm
  if(!is.null(proj_dir) && !is.null(msdm)){
    warning("msdm can not be used with future projections")
    warning("Setting msdm to NULL")
    msdm <- NULL
  }

  #3.4.Aviso caso min_occ<NPreditores
  if(min_occ<nlayers(envT)){
    warning("The minimum number of occurrences is smaller than the number of predictors.
            This may cause some issues while fitting certain algorithms!")
  }

  #4.Occurrence Data ----
  cat("Loading and processing species occurrence data ...\n")

  DirR<-"Result"
  setwd("..")
  if (file.exists(file.path(getwd(),DirR))){
    DirR<-file.path(getwd(), DirR)
  }else{
    dir.create(file.path(getwd(), DirR))
    DirR<-file.path(getwd(), DirR)
  }

  # Read txt with occurences data
  occ <- read.table(occ_file,
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
  if(!is.null(thin_occ)){
    if(thin_occ['method']%in%c('MORAN','CELLSIZE')){
      occA <- OccsThin(occ=occA, envT, as.numeric(thin_occ['method']), colin_var['method'], DirR, pred_dir)
    }
    if(thin_occ['method']=='USER-DEFINED'){
      occA <- OccsThin(occ=occA, envT, as.numeric(thin_occ['method']), colin_var['method'], DirR, pred_dir, distance=as.numeric(thin_occ['distance']))
    }
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

  #5. Restrict Extent per Species----
  if(!is.null(sp_accessible_area)){
    cat("Generating masks for species acessible area  ...\n")

    if(sp_accessible_area['method']=="BUFFER"&
       sp_accessible_area['type']=='1') {
      DirM <- M_delimited(var=envT,
                          occ_xy=occ_xy,
                          method = sp_accessible_area['method'],
                          BufferDistanceKm=NULL,
                          EcoregionsFile=NULL,
                          Dir=DirR,
                          spN=spN,
                          SaveM = TRUE,
                          Buffer_Opt=as.numeric(sp_accessible_area['type']))
    }
    if(sp_accessible_area['method']=="BUFFER"&
       sp_accessible_area['type']=='2') {
      DirM <- M_delimited(var=envT,
                          occ_xy=occ_xy,
                          method = sp_accessible_area['method'],
                          BufferDistanceKm=as.numeric(sp_accessible_area['width']),
                          EcoregionsFile=NULL,
                          Dir=DirR,
                          spN=spN,
                          SaveM = TRUE,
                          Buffer_Opt=as.numeric(sp_accessible_area['type']))
    }
    if(sp_accessible_area['method']=="MASK") {
      DirM <- M_delimited(var=envT,
                          occ_xy=occ_xy,
                          method = sp_accessible_area['method'],
                          BufferDistanceKm=NULL,
                          EcoregionsFile=sp_accessible_area['filepath'],
                          Dir=DirR,
                          spN=spN,
                          SaveM = TRUE)
    }
    if(sp_accessible_area['method']=="USER-DEFINED") {
      DirM <- M_delimited(var=envT,
                          occ_xy=occ_xy,
                          method = sp_accessible_area['method'],
                          BufferDistanceKm=NULL,
                          EcoregionsFile=sp_accessible_area['filepath'],
                          Dir=DirR,
                          spN=spN,
                          SaveM = TRUE)
    }
  }

  if(grepl("GEO", pseudoabs_method['method'])){
    Geo_Buf <- as.numeric(pseudoabs_method['width'])*1000
  }else{
    Geo_Buf=NA
  }

  #6. Geographical Partition----
  cat('Performing partition of species occurrence data ...\n')
  if(part['method']=="BANDS" || part['method']=="BLOCK"){

    if(any(grepl("PC",names(envT)))==T || any(grepl("pc",names(envT)))==T){
      colin_var <- c(method="PCA")
    }

    if(!is.null(eval_occ)){
      warning("Invalid combination! eval_occ can't be used with Geographical partitions! Changing eval_occ to NULL")
      eval_occ <- NULL
    }

    if(part['method']=="BANDS"){
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
        occINPUT <- read.table(file.path(DirB,"OccBands.txt"),sep="\t",header=T, stringsAsFactors = F)
        occINPUT[,4] <- as.numeric(occINPUT[,4])
        occINPUT[,5] <- as.numeric(occINPUT[,5])
      }else{
        if(!is.null(colin_var)){
          if(colin_var['method']!="PCA"){
            envTT<-PCA_env_TMLA(env = envT, Dir = pred_dir)
          }else{
            envTT<-envT
          }
        }else{
          envTT<-PCA_env_TMLA(env = envT, Dir = pred_dir)
        }


        # print("Use Longitudinal(1) or Latitudinal Bands(2)?")
        bands <- as.integer(part['type'])
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
            pseudoabsencesMethod = pseudoabs_method['method'],
            PrAbRatio = pres_abs_ratio,
            DirSave = DirB,
            DirM = DirM,
            MRst = sp_accessible_area['method'],
            type = TipoMoran,
            Geo_Buf = Geo_Buf,
            cores = cores
          )
        occINPUT[,4] <- as.numeric(occINPUT[,4])
        occINPUT[,5] <- as.numeric(occINPUT[,5])
        rm(envTT)
      }

    }
    if(part['method']=="BLOCK"){
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
        occINPUT <- read.table(file.path(DirB,"OccBlocks.txt"),sep="\t",header=T)
        occINPUT[,4] <- as.numeric(occINPUT[,4])
        occINPUT[,5] <- as.numeric(occINPUT[,5])
      }else{
        if(!is.null(colin_var)){
          if(colin_var['method']!="PCA"){
            envTT<-PCA_env_TMLA(env = envT, Dir = pred_dir)
          }else{
            envTT<-envT
          }
        }else{
          envTT<-PCA_env_TMLA(env = envT, Dir = pred_dir)
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
            pseudoabsencesMethod = pseudoabs_method['method'],
            PrAbRatio = pres_abs_ratio,
            DirSave = DirB,
            DirM = DirM,
            MRst = sp_accessible_area['method'],
            type = TipoMoran,
            Geo_Buf = Geo_Buf,
            cores = cores
          )

        occINPUT[,4] <- as.numeric(occINPUT[,4])
        occINPUT[,5] <- as.numeric(occINPUT[,5])
        rm(envTT)
      }
    }

    #6.3.msdm A PRIORI----
    if(is.null(msdm)||msdm["method"]%in%c('OBR', 'LR', 'PRES', 'MCP', 'MCPB')){
      DirPRI <- NULL
    }else{
      print("Creating msdm layers...")

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
    if(!is.null(proj_dir)){
      Fut <- EnvF
    }else{
      Fut <- NULL
    }

    #6.5.Adjust Checkerboard when using Geographical Restrictions----
    if(!is.null(sp_accessible_area)){
      Ms <- stack(file.path(DirM,list.files(DirM)))
      Ms <- Ms[[spN]]
      Cs <- stack(file.path(DirB,list.files(DirB,pattern=".tif$")))
      # nomesCs <- gsub("_"," ",names(Cs))
      for(i in 1:nlayers(Ms)){
        writeRaster(Ms[[i]]*Cs[[i]], file.path(DirB,names(Cs)[i]),format="GTiff",
                    bylayer=F,overwrite=T,NAflag=-9999)
      }
    }


    #6.6.Value for Sensitivity Threshold
    if(any(thr[grep('type', names(thr))]%in%"SENSITIVITY")){
      # cat("Specify the desired sensitivity value (0-1):\n")
      sensV <- as.numeric(thr['sens'])
    }else{
      sensV <- NULL
    }

    #6.7. Fit ENM for Geographical Partition----
    FitENM_TMLA_Parallel(
      RecordsData = occINPUT,
      Variables = envT,
      VarImp = imp_var,
      Fut = Fut,
      Part = part['method'],
      Algorithm = algorithm,
      PredictType = ensemble2,
      spN = spN,
      Tst = eval_occ,
      Threshold = thr,
      DirSave = DirR,
      DirMask = DirB,
      DirMSDM = DirPRI,
      Save = ifelse(save_part, 'Y', 'N'),
      SaveFinal = ifelse(save_final, 'Y', 'N'),
      sensV=sensV,
      repl = NULL,
      per = NULL,
      extrapolation=extrapolation,
      cores=cores
    )
  }

  #7.Random Partition----

  if(part['method']=="BOOT"||part['method']=="KFOLD"){

    #7.0.Dataset for evaluation
    if(!is.null(eval_occ)){
      OccTst <- read.table(eval_occ,sep="\t",h=T)
      OccTst<-OccTst[,c(sp,x,y)]
      colnames(OccTst) <- c("sp","x","y")
      OccTst_xy <- split(OccTst[,-1],f=OccTst$sp)
    }

    #7.1.msdm A PRIORI----
    if(is.null(msdm)||msdm['method']%in%c('OBR', 'LR', 'PRES', 'MCP', 'MCPB')){
      DirPRI <- NULL
    }else{
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
    if(part['method']=="BOOT"){
      rep <- as.numeric(part['replicates'])
      per<-as.numeric(part['proportion'])
    }
    if(part['method']=="KFOLD"){
      rep <- as.integer(part['folds'])
      per<-1
      occFold<- lapply(occ_xy, function(x) cbind(x,kfold(x,rep)))
      colsK <-  c("x","y","Partition");
      occFold <- lapply(occFold, setNames, colsK)
      write.table(ldply(occFold,data.frame,.id="sp"),file.path(DirR,"GruposCrossValidation.txt"),sep="\t",row.names=F)
    }

    #Adjusting for determined evaluation dataset
    if(!is.null(eval_occ) && part['method']=="BOOT" && per!=1 || !is.null(eval_occ) && part['method']=="KFOLD" && rep!=1){
      if(part['method']=="BOOT"){
        warning("Adjusting data partition to one!")
        part=c(method='BOOT', replicates='1',  proportion='1')
        per <- 1
      }
      if(part['method']=="KFOLD"){
        warning("Adjusting partition to bootstrap and data partition to one!
          Replicates will be equal to the original number of folds")
        part=c(method='BOOT', replicates='1',  proportion='1')
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
      if(part['method']=="BOOT"){
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
        if(!is.null(eval_occ)){
          occTS <- OccTst_xy
          occTS <- lapply(occTS, function(x) cbind(x, rep(2,nrow(x)),rep(1,nrow(x))))
          names(occTS) <- names(occTR)
        }
      }

      if(part['method']=="KFOLD"){
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
      if(pseudoabs_method['method']=="RND"){
        if(!is.null(proj_dir)&& !is.null(eval_occ)){
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
          if(!is.null(sp_accessible_area)){
            SpMask <- raster(file.path(DirM,paste0(names(occTR)[s],".tif")))
            SpMask <- pseudo.mask*SpMask
            if(sum(is.na(SpMask[])==F)<(pres_abs_ratio*nrow(occTR[[s]]))){
              warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences")
              stop("Please try again with another restriction type or without restricting the extent")

            }
            if(!is.null(eval_occ)){
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
        if(is.null(k) && per==1 && is.null(eval_occ)){
          for(i in 1:length(absencesTS)){
            absencesTS[[i]][,c("x","y")] <- absencesTR[[i]][,c("x","y")]
          }
        }
        DirCons <- NULL
      }

      # Pseudo-Absences allocation with Environmental constrain ----
      if(pseudoabs_method['method']=="ENV_CONST"){

        DirCons <- "EnvConstrain"
        if (file.exists(file.path(DirR,DirCons))){
          DirCons<-file.path(DirR,DirCons)
        } else {
          dir.create(file.path(DirR,DirCons))
          DirCons<-file.path(DirR,DirCons)
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
            pseudo.mask <- inv_bio(envT, occTR[[i]][,c("x","y")])

            writeRaster(pseudo.mask,paste(DirCons,spN[i],sep="/"),format="GTiff",overwrite=T)

          }else{
            pseudo.mask <- raster(file.path(DirCons,paste0(spN[i],".tif")))
          }

          if(!is.null(proj_dir)&& !is.null(eval_occ)){
            pseudo.maskP <- EnvF[[1]][[1]]
          }else{
            pseudo.maskP <- pseudo.mask
          }

          if(!is.null(sp_accessible_area)){
            SpMask <- raster(file.path(DirM,paste0(names(occTR)[i],".tif")))
            SpMask <- pseudo.mask*SpMask
            if(sum(is.na(SpMask[])==F)<(pres_abs_ratio*nrow(occTR[[i]]))){
              warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences")
              stop("Please try again with another restriction type or without restricting the extent")

            }
            if(!is.null(eval_occ)){
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
            if(!is.null(eval_occ)){
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
        if(is.null(k) && per==1 && is.null(eval_occ)){
          for(i in 1:length(absencesTS)){
            absencesTS[[i]][,c("x","y")] <- absencesTR[[i]][,c("x","y")]
          }
        }
      }

      # Pseudo-Absences allocation with Geographical constrain-----
      if(pseudoabs_method['method']=="GEO_CONST"){

        DirCons <- "GeoConstrain"
        if (file.exists(file.path(DirR,DirCons))){
          DirCons<-file.path(DirR,DirCons)
        } else {
          dir.create(file.path(DirR,DirCons))
          DirCons<-file.path(DirR,DirCons)
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

          if(!is.null(proj_dir)&& !is.null(eval_occ)){
            pseudo.maskP <- EnvF[[1]][[1]]
          }else{
            pseudo.maskP <- pseudo.mask
          }

          if(!is.null(sp_accessible_area)){
            SpMask <- raster(file.path(DirM,paste0(names(occTR)[i],".tif")))
            SpMask <- pseudo.mask*SpMask
            if(sum(is.na(SpMask[])==F)<(pres_abs_ratio*nrow(occTR[[i]]))){
              warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences")
              stop("Please try again with a smaller geographical buffer or without restricting the accessible area")

            }
            if(!is.null(eval_occ)){
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
            if(!is.null(eval_occ)){
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
        if(is.null(k) && per==1 && is.null(eval_occ)){
          for(i in 1:length(absencesTS)){
            absencesTS[[i]][,c("x","y")] <- absencesTR[[i]][,c("x","y")]
          }
        }
      }

      # Pseudo-Absences allocation with Environmentla and Geographical  constrain-----
      if(pseudoabs_method['method']=="GEO_ENV_CONST"){
        DirCons <- "GeoEnvConstrain"
        if (file.exists(file.path(DirR,DirCons))){
          DirCons<-file.path(DirR,DirCons)
        } else {
          dir.create(file.path(DirR,DirCons))
          DirCons<-file.path(DirR,DirCons)
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
            pseudo.mask1 <- inv_bio(envT, occTR[[i]][,c("x","y")])
            pseudo.mask <- pseudo.mask*pseudo.mask1
            writeRaster(pseudo.mask,paste(DirCons,spN[i],sep="/"),format="GTiff",overwrite=T)

          }else{
            pseudo.mask <- raster(file.path(DirCons,paste0(spN[i],".tif")))
          }

          if(!is.null(proj_dir)&& !is.null(eval_occ)){
            pseudo.maskP <- EnvF[[1]][[1]]
          }else{
            pseudo.maskP <- pseudo.mask
          }

          if(!is.null(sp_accessible_area)){
            SpMask <- raster(file.path(DirM,paste0(names(occTR)[i],".tif")))
            SpMask <- pseudo.mask*SpMask
            if(sum(is.na(SpMask[])==F)<(pres_abs_ratio*nrow(occTR[[i]]))){
              warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences")
              stop("Please try again with a smaller geographical buffer or without restricting the accessible area")
            }
            if(!is.null(eval_occ)){
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
            if(!is.null(eval_occ)){
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
        if(is.null(k) && per==1 && is.null(eval_occ)){
          for(i in 1:length(absencesTS)){
            absencesTS[[i]][,c("x","y")] <- absencesTR[[i]][,c("x","y")]
          }
        }
      }

      # Pseudo-Absences allocation with Environmentla and Geographical  constrain-----
      if(pseudoabs_method['method']=="GEO_ENV_KM_CONST"){
        DirCons <- "GeoEnvConstrain_KM"
        if (file.exists(file.path(DirR,DirCons))){
          DirCons<-file.path(DirR,DirCons)
        } else {
          dir.create(file.path(DirR,DirCons))
          DirCons<-file.path(DirR,DirCons)
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
            pseudo.mask1 <- inv_bio(envT, occTR[[i]][,c("x","y")])
            pseudo.mask <- pseudo.mask*pseudo.mask1
            writeRaster(pseudo.mask,paste(DirCons,spN[i],sep="/"),format="GTiff",overwrite=T)

          }else{
            pseudo.mask <- raster(file.path(DirCons,paste0(spN[i],".tif")))
          }

          if(!is.null(proj_dir)&& !is.null(eval_occ)){
            pseudo.maskP <- EnvF[[1]][[1]]
          }else{
            pseudo.maskP <- pseudo.mask
          }

          if(!is.null(sp_accessible_area)){
            SpMask <- raster(file.path(DirM,paste0(names(occTR)[i],".tif")))
            SpMask <- pseudo.mask*SpMask
            if(sum(is.na(SpMask[])==F)<(pres_abs_ratio*nrow(occTR[[i]]))){
              warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences")
              stop("Please try again with a smaller geographical buffer or without restricting the accessible area")
            }
            if(!is.null(eval_occ)){
              SpMaskP <- pseudo.maskP
            }else{
              SpMaskP <- SpMask
            }

            absencesTR.0 <- KM(rasterToPoints(SpMask)[,-3],
                               mask(envT, SpMask),
                               ((1 / pres_abs_ratio)*nrow(occTR[[i]])))
            absencesTS.0 <- KM(rasterToPoints(SpMaskP)[,-3],
                               mask(envT, SpMaskP),
                               ((1 / pres_abs_ratio)*nrow(occTS[[i]])))
          }else{
            absencesTR.0 <- KM(cell_coord = rasterToPoints(pseudo.mask)[,-3],
                               variable = mask(envT, pseudo.mask),
                               NumAbsence = (1 / pres_abs_ratio)*nrow(occTR[[i]]))
            if(!is.null(eval_occ)){
              absencesTS.0 <- KM(rasterToPoints(pseudo.maskP)[,-3],
                                 mask(envT, pseudo.mask),
                                 (1 / pres_abs_ratio)*nrow(occTS[[i]]))
            }else{
              absencesTS.0 <- KM(rasterToPoints(pseudo.maskP)[,-3],
                                 mask(envT, pseudo.mask),
                                 (1 / pres_abs_ratio)*nrow(occTS[[i]]))
            }
          }
          absencesTR[[i]] <- data.frame(x=absencesTR.0[,1],y=absencesTR.0[,2])
          absencesTS[[i]] <- data.frame(x=absencesTS.0[,1],y=absencesTS.0[,2])
          rm(absencesTR.0,absencesTS.0)
        }
        absencesTR <- lapply(absencesTR, function(x) cbind(x, rep(1,nrow(x)), rep(0,nrow(x))))
        absencesTS <- lapply(absencesTS, function(x) cbind(x, rep(2,nrow(x)), rep(0,nrow(x))))
        if(is.null(k) && per==1 && is.null(eval_occ)){
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
      if(!is.null(proj_dir)){
        Fut <- EnvF
      }else{
        Fut <- NULL
      }

      #7.6. Calculate Moran & MESS----
      if(per!=1 || part['method']=="KFOLD"){
        Bootstrap_Moran_e_MESS_TMLA(Env=envT,RecordsData=occINPUT,DirO=DirR,repl=k)
      }

      #7.7. save_part Fix----
      if(save_part && per==1){
        save_part <- FALSE
        warning("There are no partitions to be saved!")
      }

      #7.8. Background restriction----
      if(exists("DirM")){
        DirB <- DirM
      }else{
        DirB <- NULL
      }

      #7.9.Value for Sensitivity Threshold
      if(any(thr[grep('type', names(thr))]%in%"SENSITIVITY")){
        # cat("Specify the desired sensitivity value (0-1):")
        sensV <- as.numeric(thr['sens'])
      }else{
        sensV <- NULL
      }

      #7.9. Run FitENM----
      FitENM_TMLA_Parallel(
        RecordsData = occINPUT,
        Variables = envT,
        VarImp = imp_var,
        Fut = Fut,
        Part = part['method'],
        Algorithm = algorithm,
        PredictType = ensemble2,
        spN = spN,
        Tst = eval_occ,
        Threshold = thr[grep('type', names(thr))],
        DirSave = DirR,
        DirMask = DirB,
        DirMSDM = DirPRI,
        Save = ifelse(save_part, 'Y', 'N'),
        SaveFinal = ifelse(save_final, 'Y', 'N'),
        sensV = sensV,
        per = per,
        repl = k,
        extrapolation=extrapolation,
        cores=cores
      )

      #7.10. Create Occurrence Table for Replicates----
      if(rep!=1){
        occTREINO[[k]] <- occINPUT[occINPUT$Partition==1,]
        occTESTE[[k]] <- occINPUT[occINPUT$Partition==2,]
      }
    }#Fechas as replicas ou kfolds

    #7.11.Save Final Occurrence Table & Validation File----
    if(rep!=1){
      #Save Final Occurrence Table
      occTREINO <- ldply(occTREINO,data.frame,.id=NULL)
      occTESTE <- ldply(occTESTE,data.frame,.id=NULL)
      write.table(occTREINO,file.path(DirR,"Occurrences_Fitting.txt"),sep="\t",row.names=F)
      write.table(occTESTE,file.path(DirR,"Occurrences_Evaluation.txt"),sep="\t",row.names=F)

      #Save Final Validation File
      val <- list.files(DirR,pattern="Validation_Partition")
      valF <- list()
      for(i in 1:length(val)){
        valF[[i]] <- read.table(file.path(DirR,val[i]),sep="\t",header=T)
      }
      valF <- ldply(valF,data.frame,.id=NULL)
      valF <- valF[order(as.character(valF[,1])),]
      valF <- valF[!colnames(valF) %in% "Boyce_SD"]
      valF_Mean <- aggregate(.~Sp+Algorithm, data=valF, mean)
      valF_SD <- aggregate(.~Sp+Algorithm, data=valF, sd)
      valF_SD <- valF_SD[,-c(1:4)]
      colnames(valF_SD) <- paste0(colnames(valF_SD),"_SD")
      valF <- cbind(valF_Mean,valF_SD)
      valF$Replicate <- NULL
      valF$Partition <- part['method']
      unlink(file.path(DirR,val))
      write.table(valF,file.path(DirR,"Validation_Partition.txt"),sep="\t",row.names=F)

      #Save Final Bootstrap File
      if(per!=1 || part['method']=="KFOLD"){
        Boot <- list.files(DirR,pattern="Bootstrap_Moran_MESS")
        BootF <- list()
        for(i in Boot){
          BootF[[i]] <- read.table(file.path(DirR,i),sep="\t",header=T)
        }
        BootF <- ldply(BootF,data.frame,.id=NULL)
        BootF <- BootF[order(as.character(BootF[,1])),]
        BootF_Mean <- aggregate(.~sp, data=BootF, mean)
        BootF_SD <- aggregate(.~sp, data=BootF, sd)
        BootF_SD <- BootF_SD[,-c(1,4)]
        colnames(BootF_SD) <- paste0(colnames(BootF_SD),"_SD")
        BootF <- cbind(BootF_Mean,BootF_SD)
        BootF$Replicate <- NULL
        unlink(file.path(DirR,Boot))
        if(part['method']=="BOOT"){
          write.table(BootF,file.path(DirR,"Bootstrap_Moran_MESS.txt"),sep="\t",row.names=F)
        }
        if(part['method']=="KFOLD"){
          write.table(BootF,file.path(DirR,"CrossValidation_Moran_MESS.txt"),sep="\t",row.names=F)
        }
      }
    }
  }#Fecha partition BOOT|KFOLD

  #8.Ensemble----
  if (any(ensemble2!="N")){
    cat("Performing Ensemble....\n")
    ThrTable <- read.table(file.path(DirR,"Thresholds_Complete.txt"),sep="\t",h=T)
    ValF <- read.table(file.path(DirR,"Validation_Partition.txt"),sep="\t",h=T)

    Ensemble_TMLA(DirR = DirR,
                  ValTable = ValF,
                  ThrTable = ThrTable,
                  PredictType = ensemble2,
                  RecordsData = occINPUT,
                  Threshold = thr[grep('type', names(thr))],
                  sensV = sensV,
                  Proj = proj_dir,
                  ensemble_metric=ensemble_metric,
                  cores = cores)
    cat("Ensemble created!\n")
  }

  #9.MSDM Posteriori----

  if(!is.null(msdm) && msdm['method']%in%c('OBR', 'LR', 'PRES', 'MCP', 'MCPB')){
    if(any(list.files(DirR)=="Ensemble")){
      DirT <- file.path(DirR,"Ensemble",ensemble2)
      DirPost <- "MSDMPosterior"
      DirPost <- file.path(DirT,DirPost)
      for(i in DirPost){
        dir.create(i)
      }
    }else{
      DirT <- file.path(DirR,"Algorithm",algorithm)
      DirPost <- "MSDMPosterior"
      DirPost <- file.path(DirT,DirPost)
    }
    #Create MSDM Folders
    for(i in DirPost){
      dir.create(i)
    }
    #Perform MSDM Posterior
    for(i in 1:length(DirPost)){
      print(paste("Folder.....",i,"/",length(DirPost),sep=""))
      MSDM_Posterior(
        RecordsData = occINPUT,
        Threshold = thr[grep('type', names(thr))],
        cutoff = msdm,#Aqui é o tipo de MSDM-Posterior
        CUT_Buf = msdm_width,#Aqui é a distancia do Buffer
        DirSave = DirPost[i],
        DirRaster = DirT[i]
      )
    }
  }
}
