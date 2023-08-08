.predictors <- function(pred_dir,
                        proj_dir,
                        colin) {
  
  #1.Read Files in pred_dir----
  cat("Loading environmental variables...\n")
  env <- terra::rast(list.files(pred_dir,full.names = T, pattern = "\\.tif$|\\.nc$"))
  
  #1.1.CRS
  if(is.na(terra::crs(env))){
    terra::crs(envT) <- "EPSG:4326"
  }
  
  #1.2.Check predictors consistency
  env <- .syncNA(env)
  
  #2.Projection----
  if(!is.null(proj_dir)){
    
    #Check folders structure
    if(length(list.files(proj_dir,pattern = "\\.tif$|\\.nc$")>0)){
      stop("Select a parent folder that contains sub-folders with predictors for different regions or time periods, NOT a folder containing those variables!")
    }
    
    #Check names consistency

    #Current
    env_n <- list.files(pred_dir,pattern = "\\.tif$|\\.nc$")
    
    #Projections
    pfol<-list.files(proj_dir,full.names=T)
    proj_n <- lapply(pfol,function(x) list.files(x,pattern = "\\.tif$|\\.nc$"))

    if(any(unlist(lapply(proj_n, function(x) (names(env_n)!=x))))){
      stop("Current/Projections predictors names DO NOT MATCH!")
    }
    
    #Check values scale
    for(i in 1:length(pfol)){
      proj_t <- terra::rast(list.files(pfol[i],full.names=T,pattern = "\\.tif$|\\.nc$"))
      m <- terra::minmax(env-proj_t)
      if (!all(abs(m) < 1e10)){
        nom <- paste(names(env)[which(abs(m) > 10,arr.ind=T)[,2]],collapse='\n')
        warning(c('Check if current and projected variables scales match!\n',nom))
      }
    }
  }
  
  #3.Collinearity----
  if (!is.null(colin)){
    
    env <- .colin_var(env,
                      proj_dir,
                      colin)
  }else{
    
    result <- list(
      env = env
    )
  }
  
  return(env)
  
}