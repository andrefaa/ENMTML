.arg_check <- function(){
  
  #ID message
  cat("Checking for function arguments ...\n")
  
  #Check missing key predictors
  er <- NULL
  if(any(sapply(missing(pred_dir,occ_file,sp,x,y,pa_ratio,pa_method,partition,algorithm,thr)){
    which(missing(pred_dir,occ_file,sp,x,y,pa_ratio,pa_method,partition,algorithm,thr))
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
      ensemble_metric <- NULL
    } else{
      ensemble2 <- ensemble[grep('method', names(ensemble))]
      if (any(ensemble2 %in% c("SUP", "W_MEAN", "PCA_SUP")) & !any(names(ensemble) == 'metric')) {
        stop("While using SUP, W_MEAN or PCA_SUP ensemble methods, provide an evaluation metric to ensemble arguemnt (AUC, Kappa, TSS, Jaccard, Sorensen or Fpb). \n e.g., ensemble=c(method=c('W_MEAN', 'SUP'), metric='Jaccard')")
      } else if (any(ensemble2 %in% c("SUP", "W_MEAN", "PCA_SUP"))) {
        ensemble_metric <- ensemble[grep('metric', names(ensemble))]
        if(!ensemble_metric%in%c('AUC', 'Kappa', 'TSS', 'Jaccard', 'Sorensen', 'Fpb', 'Boyce')){
          stop("'metric' used in the 'ensemble' argument did not match with ENMTML metrics\n",
               "Please select one of these: \n AUC\n Kappa\n TSS\n Jaccard\n Sorensen\n Fpb\n Boyce")
        }
      }
      if (!any(ensemble2 %in% c("SUP", "W_MEAN", "PCA_SUP"))) {
        ensemble_metric <- NULL
      }
    }
    
    
    
    if(!is.null((er))){
      print(er)
      stop("Missing arguments. Please check the arguments listed above")
    }
    
    if(!is.null(thin_occ)){
      if(!(thin_occ['method']%in%c('moran','cellsize','user-defined'))){
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
    if((pseudoabs_method['method']%in%c("GEO_CONST", "GEO_ENV_CONST", "GEO_ENV_KM_CONST")) & length(pseudoabs_method)!=2){
      stop("'pseudoabs_method' Argument is not valid!. Please specify a 'width' for the buffer!")
    }
    if(length(pseudoabs_method['method'])>1){
      stop("Please choose only one Pseudo-absence allocation method")
    }
    if(!(part['method']%in%c("BOOT","KFOLD","BANDS","BLOCK"))){
      stop("'part' Argument is not valid!(BOOT/KFOLD/BANDS/BLOCK)")
    }
    if(any(!algorithm%in%c("BIO","MAH","DOM","ENF","GLM","GAM","SVM","SVM-B","BRT","RDF","MXS","MXD","MLK","GAU"))){
      stop(paste("Algorithm", algorithm[!(algorithm%in%c("BIO","MAH","DOM","ENF","GLM","GAM","SVM","SVM-B","BRT","RDF","MXS","MXD","MLK","GAU"))],"is not valid"))
    }
    if(any(!thr[grep('type', names(thr))]%in%c("LPT","MAX_TSS","MAX_KAPPA","SENSITIVITY","JACCARD","SORENSEN"))){
      stop("'thr' Argument is not valid!")
    }
    if(any(thr[grep('type', names(thr))]%in%"SENSITIVITY") && !any(names(thr)%in%'sens')){
      stop("provide a sensitivity value in the vector used in 'thr' argument, see ENMTML function help")
    }
    if(!is.null(msdm)){
      if(!(msdm['method']%in%c('XY','MIN','CML','KER', 'OBR', 'LR', 'PRES', 'MCP', 'MCP-B'))){
        stop("'msdm' Argument is not valid!(XY/MIN/CML/KER/OBR/LR/PRES/MCP/MCP-B)")
      }
      
      if(length(msdm)>1){
        if(any(names(msdm)=='width')){
          msdm_width <- as.numeric(msdm['width'])
        }else{
          stop("More than one M-SDM method chosen: Please select a SINGLE M-SDM method")
        }
      }else{
        msdm_width <- NULL
      }
      
      if(any((msdm['method']%in%c('XY','MIN','CML','KER'))) & !is.null(proj_dir)){
        stop("It is not yet possible to combine priori M-SDMs with Projections!\nPlease remove the priori M-SDM or set the projection argument to NULL")
      }
    }
    
    if(!is.null(sp_accessible_area)){
      if(!(sp_accessible_area['method']%in%c('BUFFER','MASK','USER-DEFINED'))){
        stop("'sp_accessible_area' Argument is not valid!(BUFFER/MASK/USER-DEFINED)")
      }
      if(sp_accessible_area['method']=="USER-DEFINED"&length(sp_accessible_area)!=2){
        stop("'sp_accessible_area' Argument is not valid for method=USER-DEFINED! A folder containing the masks must be provided e.g., sp_accessible_area=c(method='USER-DEFINED', filepath='C:/Users/mycomputer/accessibleareafolder')")
      }
      if(sp_accessible_area['method']=="MASK"&length(sp_accessible_area)!=2){
        stop("'sp_accessible_area' Argument is not valid for method=MASK! A filepath containing the file used to generate the species-specific masks must be provided e.g., sp_accessible_area=c(method='USER-DEFINED', filepath='C:/Users/mycomputer/ecoregion/olson.shp')")
      }
    }
  
  #2.Adjust Names----
  Ord <- c("BIO","DOM","MAH","ENF","MXD","MXS","MLK","SVM","SVM-B","RDF","GAM","GLM","GAU","BRT")
  algorithm <- Ord[Ord%in%algorithm]
  
  Ord_Thr <- c("MAX_KAPPA","MAX_TSS","LPT","SENSITIVITY","JACCARD","SORENSEN")
  thr <- Ord_Thr[Ord_Thr%in%thr]
  names(thr) <- rep("type", length(thr))
}

#Arguments check from sub-functions
if (!any(c("pearson", "vif", "pca", "fa") %in% colin)) {
  stop(
    "argument 'colin' was misused, select one of the available methods: pearson, vif, pca, fa"
  )
}