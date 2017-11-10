Split.NCDF<-function(dir_in,dir_out,resample=""){
  #Function to read NCDF files and save each layer as an ASC
  #Arguments
    #dir_in:directory of the netcdf file
    #dir_out:directory to save ASC files
  library(raster)
  library(SDMTools)
  
  setwd(dir_in)
  ncd<-list.files(pattern=".nc")
  ncd<-brick(ncd,varname="v")
  if(resample=="Y"){
    msk<-read.asc(file.choose())
    msk<-raster(msk)
  }
  setwd(dir_out)
  yr<-1992
  mnt<-9
  for (x in 1:nlayers(ncd)){
    if (mnt==12){
      mnt<-0
      yr<-yr+1
    }
    mnt<-mnt+1
    sst<-raster(ncd,layer=x)
    sst<-rotate(sst)
    if (resample=="Y"){
      sst<-resample(sst,msk)
      sst<-mask(sst,msk)
    }
    writeRaster(sst,paste(paste("MerVeloc",mnt,yr,sep="_"),".asc",sep=""),format="ascii")
    print(paste(paste(mnt,yr,sep="/"),"03/2014",sep="..."))
  }
}

Split.NCDF ("C:\\OneDrive\\IBM_Tubastraea\\IBM_Layers\\Month_SST","C:\\OneDrive\\IBM_Tubastraea\\IBM_Layers\\Month_SST",resample="Y")
