setwd("C:/Users/Andre/Google Drive/Mestrado/Capitulo 1(Transferabilidade)/Ocorrencias_20/Avaliacao Modelos Maxent/PCA/2Q/Test-train")
lista<-list.files()
length(lista)
c<-matrix(0,nrow=((length(lista)*10)+10))


Automatic.File.Import_Size.of.Predicted.Area<-function (lista){
  for (x in 1:length(lista)){
    
    setwd(paste("C:/Users/Andre/Google Drive/Mestrado/Capitulo 1(Transferabilidade)/Ocorrencias_20/Avaliacao Modelos Maxent/PCA/2Q/Test-train/",lista[x],paste("/"),sep=""))
    sublista.xls<-list.files(paste(getwd(),"/",sep=""),pattern="mod.xls")
    sublista.asc<-list.files(paste(getwd(),"/",sep=""),pattern=".asc")
    tab.avaliacao<-read.xlsx(paste(sublista.xls),sheetIndex=1)
    th.roc<-tab.avaliacao["th_roc"]
    b<-matrix(0,nrow=nrow(th.roc))
    
    for (y in 1:nrow(th.roc)){
      table <- read.asciigrid(sublista.asc[y])
      naomit<-na.exclude(table@data)
      tresh<-ifelse(naomit>=th.roc[y,],1,0)
      soma<-sum(tresh)
      b[y,]<-soma
      b
    }
    c[((x*10):(((x+1)*10)-1)),1]<-b
  }
  write.xlsx(c,"Size_of_predicted_area.xlsx")
}

Automatic.File.Import_Size.of.Predicted.Area(lista)
