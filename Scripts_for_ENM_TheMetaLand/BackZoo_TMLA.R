BackZoo_TMLA <- function(Dir,DirZ,occ){
  
  DirZO<-"ZooRegions"
  if (file.exists(file.path(Dir,DirZO))){
    DirZO<-file.path(Dir,DirZO)
  } else {
    dir.create(file.path(Dir,DirZO))
    DirZO<-file.path(Dir,DirZO)
  }
  
  if(length(list.files(DirZO))==length(occ)){
    print("ZooMasks Found! Using Already Created ZooMasks")
  }else{
    
    #Import Zooregions ASCs
    Zoo <- brick(raster(paste(DirZ,list.files(DirZ,pattern = "\\.tif$"),sep="/")))
    
    #Get Species Zoogeographical Locations
    Ext <- lapply(occ,function(x) extract(Zoo,x)) 
    sp.zoo <- lapply(Ext,function(x) unique(x))
    sp.zoo <- lapply(sp.zoo, function(x) x[x!=0])
    for(i in 1:length(sp.zoo)){
      ZooSp <- Zoo%in%sp.zoo[[i]]
      writeRaster(ZooSp,paste(DirZO,paste(names(occ)[i],".tif",sep=""),sep="/"),format="GTiff",overwrite=T)
    }
  }
  return(DirZO)
}

 
