Adequability.Subset<-function(dir_adeq,dir_subset,dir_output){
  library(SDMTools)
  library(raster)
  
  adeq<-list.files(dir_adeq,pattern = ".asc")
  models<-substr(adeq,1,3)
  adeq<-unique(models)
  setwd(dir_output)

  for (a in 1:length(adeq)){
    rasters<-list.files(dir_adeq,pattern= adeq[a])
    rasters<-grep(".asc",rasters,value=T)
    for (b in 1:length(rasters)){
      asc<-read.asc(paste(dir_adeq,rasters[b],sep="\\"))
      ras<-raster(asc)
      names(ras)<-substr(rasters[b],1,nchar(rasters[b])-4)
      mask.novel<-read.asc(paste(dir_subset,paste("Novel",rasters[b],sep="_"),sep="\\"))
      mask.novel<-raster(mask.novel)
      mask.over<-read.asc(paste(dir_subset,paste("Overlap",rasters[b],sep="_"),sep="\\"))
      mask.over<-raster(mask.over)
      
      ras.novel<-ras
      ras.novel[mask.novel[]!=1]<-0
      writeRaster(ras.novel,paste(paste("Novel","Adeq",rasters[b],sep="_")),format="ascii")
      
      ras.over<-ras
      ras.over[mask.over[]!=1]<-0
      writeRaster(ras.over,paste(paste("Overlap","Adeq",rasters[b],sep="_")),format="ascii")
    }
  }
}

Adequability.Subset("C:\\OneDrive\\Cap2(Vulnerabilidade)\\Resultados ENMs\\Presente\\Ensemble\\ENS",
                    "C:\\OneDrive\\Cap2(Vulnerabilidade)\\Resultados ENMs\\Overlap\\Presente-A1B",
                    "C:\\OneDrive\\Cap2(Vulnerabilidade)\\Resultados ENMs\\Overlap\\Presente-A1B(Adequabilidade)")

