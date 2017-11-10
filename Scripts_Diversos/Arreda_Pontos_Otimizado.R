Arreda_Pontos_Otimizado<-function (diretorio_output){
 #Parametros
    #diretorio_output: diretorio onde sera salvo o arquivo com os pontos arredados

  require(SDMTools)
  require(raster)
  
#Importar a máscara
  
  print("Select the asc file that will be used as a mask to move the points into:")
  asc_base<-read.asc(file.choose())
  environmental<-raster(asc_base)

#Selecionar arquivo com os registros de ocorrencia(.txt):
  
  print("Select the occurrences file(.txt):")
  ocorrencias.r<-read.table(file.choose(),h=T,sep='\t')
  ocorrencias<-ocorrencias.r[,which(names(ocorrencias.r) %in% c("Species","Long","Lat"))]
  ocorrencias.xy<-ocorrencias[,-1]

  pontos.env<-extract(environmental,ocorrencias.xy,cellnumbers=TRUE)
  pontos.env<-cbind(ocorrencias,pontos.env)
  pontos.na<-pontos.env[(is.na(pontos.env[,5])),]
  pontos.correto<-pontos.env[(is.na(pontos.env[,5])==F),]

#Identificar a coluna e linha de cada ponto

  rows<-rowFromY(environmental,pontos.na$Lat) #posicao da celula no eixo Y
  cols<-colFromX(environmental,pontos.na$Long) #posicao da celula no eixo X
  coord<-cbind(rows,cols) #Colunas e linha da célula
  pontos.na<-cbind(pontos.na,coord)

#Arredar os pontos com base na vizinhança


for (L in 1:nrow(pontos.na)){
  
  cell.value<-NULL
  print(paste("Linha Original:"))
  print(paste(pontos.na[L,]))
  
  for (x in 1:10){
    for (L1 in (pontos.na[L,6]-x):(pontos.na[L,6]+x)){
      for (C1 in (pontos.na[L,7]-x):(pontos.na[L,7]+x)){
        cell.value<-rbind(cell.value,c(environmental[L1,C1],L1,C1))
      }
    }
    if (all(is.na(cell.value[,1]))==F){
      cell.coord<-which(is.na(cell.value[,1])==F)
      cells.validas<-NULL
      cells.validas<-rbind(cells.validas,cell.value[cell.coord,])
      #cells.validas<-as.matrix(cell.value[cell.coord,])
      cell.number<-NULL
      cell.xy<-NULL
      d<-NULL
      for (b in 1:nrow(cells.validas)){
        print(x)
        print(cells.validas)
        cell.number[b]<-cellFromRowCol(environmental, cells.validas[b,2],cells.validas[b,3])
        cell.xy<-rbind(cell.xy,xyFromCell(environmental,cell.number[b],spatial=FALSE))
        d<-rbind(d,sqrt((cell.xy[b,1]-pontos.na[L,2])^2+(cell.xy[b,2]-pontos.na[L,3])^2))
      }
      cell.xy<-cbind(cell.xy,d)
      cell.prox<-which(cell.xy[,3]==min(cell.xy[,3]))
      pontos.na[L,2:3]<-cell.xy[cell.prox,1:2]
    break}
  }
  print(paste("Linha Arredada:"))
  print(paste(pontos.na[L,]))
  print(paste("Porcentagem completa:",L/nrow(pontos.na)))
}

pontos.na.novos<-pontos.na[,1:3]
pontos.novos<-extract(environmental,pontos.na.novos[,2:3],cellnumber=T)
pontos.novos<-cbind(pontos.na.novos,pontos.novos)

pontos.arredados<-rbind(pontos.correto,pontos.novos)

setwd(diretorio_output)

write.table(pontos.arredados[,1:3],"Ocorrencias_Arredadas.txt",row.names=F,sep="\t")
}

Arreda_Pontos_Otimizado("C:\\OneDrive\\SDM_Aedes")
