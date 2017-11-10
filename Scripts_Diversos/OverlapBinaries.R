Overlap.Binary<-function(dir1,dir2,dir_output){
  library(SDMTools)
  library(raster)
  
  ascs1<-list.files(dir1,pattern = ".asc")
  ascs2<-list.files(dir2,pattern = ".asc")
  n.over<-NULL
  n.nov<-NULL
  n.los<-NULL
  n.ras1<-NULL
  n.ras2<-NULL
  
  setwd(dir_output)

  for (a in 1:length(ascs1)){
    asc1<-read.asc(paste(dir1,ascs1[a],sep="\\"))
    ras1<-raster(asc1)
    asc2<-read.asc(paste(dir2,ascs1[a],sep="\\"))
    ras2<-raster(asc2)
    
    names(ras1)<-substr(ascs1[a],1,nchar(ascs1[a])-4)
    names(ras2)<-substr(ascs2[a],1,nchar(ascs2[a])-4)
    
    overlap<-ras1==1 & ras2==1
    novel<-ras1==0 & ras2==1
    lost<-ras1==1 & ras2==0
    
    n.ras1<-c(n.ras1,sum(na.omit(ras1[])))
    n.ras2<-c(n.ras2,sum(na.omit(ras2[])))
    n.over<-c(n.over,sum(na.omit(overlap[])))
    n.nov<-c(n.nov,sum(na.omit(novel[])))
    n.los<-c(n.los,sum(na.omit(lost[])))
    
    writeRaster(overlap,paste("Overlap",names(ras1),sep="_"),format="ascii")
    writeRaster(novel,paste("Novel",names(ras1),sep="_"),format="ascii")
    writeRaster(lost,paste("Lost",names(ras1),sep="_"),format="ascii")
  }
  cell.n<-cbind(ascs1,n.ras1,n.ras2,n.over,n.nov,n.los)
  colnames(cell.n)<-c("sp","Presente","Futuro","Overlap","Novel","Lost")
  write.table(cell.n,"N_Celulas.txt",sep="\t",row.names=F)
}

Overlap.Binary("C:\\OneDrive\\Cap2(Vulnerabilidade)\\ENMs\\Resultados ENMs\\Presente\\Ensemble\\Binary",
               "C:\\OneDrive\\Cap2(Vulnerabilidade)\\ENMs\\Resultados ENMs\\Binarios\\Futuro(A1B)",
               "C:\\OneDrive\\Cap2(Vulnerabilidade)\\ENMs\\Resultados ENMs\\Overlap\\Novos\\Pres-A1B")
