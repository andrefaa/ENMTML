.occurrences <- function() {
  cat("Loading and processing species occurrence data ...\n")
  
  if(is.null(result_dir)){
    DirR <- file.path(dirname(pred_dir),"Result")
    if (file.exists(DirR)){
      warning("Result folder already exists, files may be overwritten!")
    }else{
      dir.create(DirR, recursive = T)
    }
  }else{
    DirR <- result_dir
    if (!grepl('/', DirR)){
      message("Folder with results will be created at the same level of the predictors' folder")
      if (file.exists(file.path(dirname(pred_dir),DirR))){
        message("Result folder already exists, files may be overwritten!")
      }
      DirR <- file.path(dirname(pred_dir),result_dir)
      dir.create(DirR, recursive = T)
    }else{
      if (file.exists(DirR)){
        message("Result folder already exists, files may be overwritten!")
      }
      dir.create(DirR, recursive = T)
    }
  }
  cat(paste0("Results can be found at:  ","\n",DirR,"\n"))
  
  # Read txt with occurences data
  occ <- utils::read.table(occ_file,
                           h = TRUE,
                           sep = '\t',
                           stringsAsFactors = FALSE)
  occ<-occ[,c(sp,x,y)]
  suppressWarnings(occ[[x]] <- as.numeric(occ[[x]]))
  suppressWarnings(occ[[y]] <- as.numeric(occ[[y]]))
  
  if(any(is.na(occ[,x]) | is.na(occ[,y]))){
    message("Some records were removed because data are no numeric in latitude or longitude")
    message(paste("Number of records removed:", sum(is.na(occ[,x]) | is.na(occ[,y]))))
  }
  
  colnames(occ) <- c("sp","x","y")
  #Correct issues caused by species name separated by space
  occ$sp <- gsub(" ","_",occ$sp)
  occ_xy <- split(occ[,-1],f=occ$sp)
  spN<-names(occ_xy)
  
  
  #4.1.Unique Occurrences----
  occA<-Occ_Unicas_TMLA(env=envT[[1]], occ.xy=occ_xy, DirO=DirR)
  
  #4.2.Thinning----
  if(!is.null(thin_occ)){
    cat("Thinning occurrences...\n")
    if(thin_occ['method']%in%c('MORAN','CELLSIZE')){
      occA <- OccsThin(occ=occA, envT, as.character(thin_occ['method']), colin_var['method'], DirR, pred_dir)
    }
    if(thin_occ['method']=='USER-DEFINED'){
      occA <- OccsThin(occ=occA, envT, as.character(thin_occ['method']), colin_var['method'], DirR, pred_dir, distance=as.numeric(thin_occ['distance']))
    }
  }
  
  #4.3.Save Thinned & Unique Occurrences
  ndb <- plyr::ldply(occA, .id='sp')[,1:3]
  utils::write.table(ndb,file.path(DirR,"Occurrences_Cleaned.txt"),sep="\t",row.names=F)
  
  #4.3.Remove species with less than min_occ----
  occ <- occA[sapply(occA,function (x) nrow(x)>=min_occ)]
  spN<-names(occ)
  
  
  #4.4.Species with few records----
  if(length(occ)!=length(occ_xy)){
    cat(paste("Species with less than ",min_occ, " Unique Occurrences were removed! \n"))
    print(names(occ_xy)[names(occ_xy)%in%spN==F])
    ndb <- plyr::ldply(occ)[,1:3]
    utils::write.table(ndb,file.path(DirR,"Occurrences_Filtered.txt"),sep="\t",row.names=F)
    rm(ndb)
  }
  occ_xy <- lapply(occ,function(x) x[,c("x","y")])
}