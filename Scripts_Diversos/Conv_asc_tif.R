Conv_asc_tif <- function(DirIn,DirOut){
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  ipak(c("raster","doParallel","foreach"))
  
  Ascs <- file.path(DirIn,list.files(DirIn,pattern="\\.asc$"))
  nomes <- substr(list.files(DirIn,pattern="\\.asc$"),1,nchar(list.files(DirIn,pattern="\\.asc$"))-4)
  
  cl <- makeCluster(detectCores()-2)    #create a cluster 
  registerDoParallel(cl)
  
  x1 <- foreach(i= 1:length(Ascs), .packages = "raster")%dopar%{
    asc <- raster(Ascs[i])
    writeRaster(asc,file.path(DirOut,nomes[i]),bylayer=T,format="GTiff",overwrite=T)
  }
  stopCluster(cl)
}

Conv_asc_tif("D:\\Doc_ENM\\Narrow\\NaoEquilibrio\\Distribuicao",
             "D:\\Doc_ENM\\Narrow\\NaoEquilibrio\\Distribuicao")
