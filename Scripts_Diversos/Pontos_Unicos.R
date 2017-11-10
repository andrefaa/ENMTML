Unique_Pontos <- function(diretorio_mask,output){
  #Parametros:
   #diretorio_mask:diretorio do arquivo .asc que sera usado como mascara
   #output: diretorio que serao salvos as ocorrencias unicas(.txt)
  
  library(SDMTools)
  library(raster)
  
  #Importar arquivo de ocorrencias (.txt)
  
  print("Selecione o arquivo com as ocorrencias(.txt):")
  ocorrencias.r<-read.table(file.choose(),header=T,sep="\t")
  ocorrencias<-ocorrencias.r[,which(names(ocorrencias.r) %in% c("Taxa","Long","Lat"))]
  
  #Carregar camadas ambientais
  
  setwd(diretorio_mask)
  
  layers<-list.files(pattern=".asc") #pegar todos layers e salvar como objeto  atribuindo nomes
  environmental<-read.asc(layers[1]) #ler os nomes dos arquivos
  environmental<-raster(environmental)
  
  ocorrencias.xy<-ocorrencias[,c(2,3)]
  celula<-extract(environmental,ocorrencias.xy,cellnumber=T)
  tabela<-cbind(ocorrencias,celula)
  
  setwd(output)
  
  #Arredar ocorrencias para dentro da mascara
  
  ocorrencias.var.sna<-na.omit(tabela)
  ocorrencias.na<-tabela[(is.na(tabela[,5])),]
  row.names(ocorrencias.na)<-NULL
  
  if (nrow(ocorrencias.na)!=0){
    ocorrencias.arredado<-ocorrencias.na
    mask.tabela<-rasterToPoints(environmental)
    d1<-NULL
    
    mask.tabela.x<-mask.tabela[,1]
    mask.tabela.y<-mask.tabela[,2]
    
    ocorrencias.na.x<-ocorrencias.na[,2]
    ocorrencias.na.y<-ocorrencias.na[,3]
    
    for (b in 1:nrow(ocorrencias.na)){
      print(ocorrencias.arredado[b,c(1:3)])
      for (c in 1:nrow(mask.tabela)){
        d1[c]<-sqrt((mask.tabela.x[c]-ocorrencias.na.x[b])^2+(mask.tabela.y[c]-ocorrencias.na.y[b])^2)
      }
      posicao<-match(min(d1),d1)
      menor.dist<-mask.tabela[posicao,]
      ocorrencias.arredado[b,c(2,3)]<-menor.dist[c(1,2)]
      print(ocorrencias.arredado[b,c(1:3)])
      print(paste("Linha_",b,"/",nrow(ocorrencias.na),"....OK",sep=""))
      d1<-NULL
    }
    
    ocorrencias.xy.arredado<-ocorrencias.arredado[,c(2,3)]
    ocorrencias.var.arredado<-extract(environmental,ocorrencias.xy.arredado,cellnumber=T)
    ocorrencias.var.arredado<-cbind(ocorrencias.arredado[,c(1:3)],ocorrencias.var.arredado)
    ocorrencias.finais<-rbind(ocorrencias.var.sna,ocorrencias.var.arredado)
    ocorrencias.finais.sem.var<-ocorrencias.finais[,which(names(ocorrencias.finais) %in% c("Taxa","Long","Lat"))]
  
    }else{
    ocorrencias.finais<-tabela
    ocorrencias.finais.sem.var<-ocorrencias.finais[,which(names(ocorrencias.finais) %in% c("Taxa","Long","Lat"))]
  }
  
  write.table(ocorrencias.finais.sem.var,"Ocorrencias.arredadas.txt",sep="\t",row.names=F)
  
  #Cria um vetor com as especies presentes na planilha
  especies<-unique(ocorrencias.var.arredado[,1])
  
  #Remover presencas duplicadas
  
  ocorrencias.var<-ocorrencias.finais
  
  ocorrencias.unicas<-NULL
  registros.ssp<-NULL
  registros.unicos<-NULL
  
  for (b in 1:length(especies)){
    
    ocorrencias.especie<-ocorrencias.var[ocorrencias.var[,1]==especies[b],]
    duplicadas<-which(duplicated(ocorrencias.especie[,'cells'])==T)
    registros.ssp<-nrow(ocorrencias.especie[-duplicadas])
    registros.unicos<-rbind(registros.unicos,registros.ssp)
    
    if (length(duplicadas)!=0){
      ocorrencias.unicas<-rbind(ocorrencias.unicas,ocorrencias.especie[-duplicadas,])
    }else{
      ocorrencias.unicas<-rbind(ocorrencias.unicas,ocorrencias.especie)
    }    
  }
  
  rownames(ocorrencias.unicas)<-NULL
  registros.unicos<-cbind(as.character(especies),registros.unicos)
  colnames(registros.unicos)<-c("Especie","Registros_Unicos")
  write.table(ocorrencias.unicas,"Ocorrencias_Unicas.txt",row.names=F,sep="\t")
  write.table(registros.unicos,"Registros_Unicos.txt",row.names=F,sep="\t")
}

Unique_Pontos("C:\\Users\\Andre\\Google Drive\\Mestrado\\Poliquetas\\ENM\\ASC\\25km","C:\\Users\\Andre\\Google Drive\\Mestrado\\Poliquetas\\ENM\\ASC\\25km")
