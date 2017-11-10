Leave_1_Out_Aval<-function(dir_input,ocor){
  setwd(dir_input)
  
  library(SDMTools)
  library(raster)
  
  setwd(dir_input)
  
  lista<-list.files(pattern=".asc")
  asc_base<-read.asc(lista[1])
  adeq<-raster(asc_base)
  shunda<-rasterToPoints(adeq)
  nrow(shunda)
  camadas<-lista[-1]
  
  for (x in 1:length(camadas)){
    asc<-read.asc(camadas[x])
    adeq.temp<-raster(asc)
    adeq<-stack(adeq,adeq.temp)
  }
  names(adeq)<-substr(lista,1,nchar(lista)-4)
  
  ocor<-read.table(ocor,sep="\t",header=T)
  ocor.xy<-ocor[,2:3]
  ocor.adeq<-NULL
  for (y in 1:nrow(ocor.xy)){
    ocor.adeq.t<-extract(adeq[[y]],as.data.frame(ocor.xy[y,]))
    ocor.adeq<-c(ocor.adeq,ocor.adeq.t)
  }
  ocor<-cbind(ocor,ocor.adeq)
  write.table(ocor,"Aval_Th.txt",row.names=F,sep="\t")
}

Leave_1_Out_Aval("C:\\Users\\decoa\\Desktop\\Modelagem Zandim\\Andre\\Spp\\Resul_MaxS","C:\\Users\\decoa\\Desktop\\Modelagem Zandim\\Andre\\Avaliacao\\Avaliacao_Thr_Spp_MaxS.txt")
