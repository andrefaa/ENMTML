Binary.Maps.Livia <- function(dir_input,dir_output){
  
  library(SDMTools)
  library(raster)
  
  #dir_input:Diretorio com os SDM
  #dir_output: Diretorio onde serao salvos os mapas binarios e a tabela
  #Tabela de thresholds deve conter 2 colunas(Especie,Th_Roc); em formato .txt
  
    setwd(dir_input)
  
    thrs<-read.table(file.choose(),h=T,sep='\t')
    
    modelos<-list.files(pattern=".asc")
  
    asc.base<-read.asc(modelos[1])
    ras<-raster(asc.base)
    camadas<-modelos[-1]
    
    for (b in 1:length(camadas)){
      asc<-read.asc(paste(dir_input,camadas[b],sep="\\"))
      ras.t<-raster(asc)
      ras<-stack(ras,ras.t)
    }
    names(ras)<-substr(modelos,1,nchar(modelos)-4)
    setwd(dir_output)
    
    matriz.bin<-NULL
    
    for (c in 1:nrow(thrs)){
      th.sp<-thrs[c,2]
      bin.sp<-ras[[c]]>=th.sp
      matr.sp<-as.data.frame(bin.sp)
      matr.sp<-as.numeric(matr.sp[,1])
      matriz.bin<-cbind(matriz.bin,matr.sp)
      #writeRaster(bin.sp,paste(names(ras[[c]]),".asc",sep=""),format="ascii")
    }
    colnames(matriz.bin)<-substr(modelos,1,nchar(modelos)-4)
    msk<-ras[[c]]
    xy<-xyFromCell(msk,cell=1:ncell(msk))
    colnames(xy)<-c("x","y")
    matriz.bin<-as.data.frame(cbind(xy,matriz.bin))
    write.table(matriz.bin,"Matriz_Presenca_Ausencia.txt",sep="\t",row.names=F)
}

Binary.Maps.Livia("C:\\Users\\decoa\\Dropbox\\asc\\Futuro","C:\\Users\\decoa\\Dropbox\\asc\\Futuro")
