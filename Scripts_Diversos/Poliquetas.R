#1.Camadas Ambientais_Poliquetas

setwd("C:/Users/Andre/Google Drive/Mestrado/Layers/BioOracle/BioOracle_05_Poliq")
rm(list=ls())

env_prep<-list.files(pattern=".asc")
library(SDMTools)
asc_base<-read.asc(env_prep[2])
environmental<-raster(asc_base)
camadas<-env_prep[-1]

for (a in 1:length(camadas)){
  asc<-read.asc(camadas[a])
  raster<-raster(asc)
  environmental<-stack(environmental,raster)
}

names(environmental)<-env_prep

extent<- c(-85,-25,-56,10)
environmental<- crop(environmental,extent)

setwd("C:/Users/Andre/Google Drive/Mestrado/Poliquetas")

env_prep2<-list.files(pattern=".asc")
asc_base2<-read.asc(env_prep2[2])
environmental2<-raster(asc_base2)
environmental2<-crop(environmental2,extent)
names(environmental2)<-env_prep2[2]

environmental<-stack(environmental,environmental2)
names(environmental)
tabela.env<-rasterToPoints(environmental)

#2.ARQUIVOS DE OCORRENCIA

ocorrencias<-read.table(file.choose(),h=T,sep=",")
ocorrencias[,3]<-paste(ocorrencias[,3],ocorrencias[,4],sep="_")
ocorrencias<-ocorrencias[,-4]
ocorrencias.rear<-ocorrencias[,1]
ocorrencias.xy<-ocorrencias[,-1]

ocorrencias.atlantico<-ocorrencias[ocorrencias[,"Lat"]>-56,]
ocorrencias.atlantico[,3]<-paste(ocorrencias.atlantico[,3],ocorrencias.atlantico[,4],sep="_")
ocorrencias.atlantico<-ocorrencias.atlantico[,-4]
ocorrencias.rear<-ocorrencias.atlantico[,c(3,2,1)]
ocorrencias.xy<-ocorrencias.rear[,-1]


ocorrencias.var<-as.data.frame(extract(environmental,ocorrencias.xy,cellnumber=T))
ocorrencias.var<-cbind(ocorrencias,ocorrencias.var)
ocorrencias.var.sna<-na.omit(ocorrencias.var)
ocorrencias.na<-ocorrencias.var[(is.na(ocorrencias.var[,5])),]
row.names(ocorrencias.na)<-NULL
write.xlsx(ocorrencias.na,"Ocorrencias_Poliquetas_Fora_Mascara.xls",row.names=F)


#Arredar ocorrencias para dentro da mascara
ocorrencias.arredado<-ocorrencias.na
mask.tabela<-rasterToPoints(environmental)
d1<-NULL

mask.tabela.x<-mask.tabela[,1]
mask.tabela.y<-mask.tabela[,2]

ocorrencias.na.x<-ocorrencias.na[,2]
ocorrencias.na.y<-ocorrencias.na[,3]

for (a in 1:nrow(ocorrencias.na)){
  print(ocorrencias.arredado[a,])
  plot(environmental)
  points(ocorrencias.arredado[a,c(2,3)])
  for (b in 1:nrow(mask.tabela)){
    d1[b]<-sqrt((mask.tabela.x[b]-ocorrencias.na.x[a])^2+(mask.tabela.y[b]-ocorrencias.na.y[a])^2)
  }
  posicao<-match(min(d1),d1)
  menor.dist<-mask.tabela[posicao,]
  #tabela.distancia<-cbind(mask.tabela,d1)
  #menor.dist<-tabela.distancia[tabela.distancia[,"d1"] == min(tabela.distancia[,"d1"]),]
  ocorrencias.arredado[a,c(2,3)]<-menor.dist[c(1,2)]
  plot(environmental)
  points(ocorrencias.arredado[a,c(2,3)])
  print(ocorrencias.arredado[a,])
  print(paste("Linha_",a,"....OK",sep=""))
  d1<-NULL
}
ocorrencias.xy.arredado<-ocorrencias.arredado[,c(2,3)]


ocorrencias<-read.table(file.choose(),h=T,sep=",")

ocorrencias.var<-extract(environmental,ocorrencias.xy.arredado,cellnumber=T)
ocorrencias.var.arredado<-cbind(ocorrencias.arredado[,c(1:3)],ocorrencias.var)

ocorrencias.finais<-rbind(ocorrencias.var.sna,ocorrencias.var.arredado)
ocorrencias.finais<-ocorrencias.finais[,c(1:3)]
write.xlsx(ocorrencias.finais,"Ocorrencias_Poliquetas_Arredado.xls",row.names=F)

plot(environmental)
points(ocorrencias.arredado[,c(2,3)])
    #ocorrencias.na.sp<-distancias.celula[distancias.celula[,3]==min(distancias.celula[,3]),c(1,2)]
    #ocorrencias.arredadas<-rbind(ocorrencias.arredadas,ocorrencias.na.sp)


#Cria um vetor com as esp?cies presentes na planilha
especies<-unique(ocorrencias.rear[,1])
length(especies)

#Remover presen?as duplicadas


plot(environmental)
points(ocorrencias.na[,c(2,3)])

ocorrencias.unicas<-NULL

for (b in 1:length(especies.na)){
  ocorrencias.especie<-ocorrencias.na[ocorrencias.na[,1]==especies.na[b],]
  duplicadas<-which(duplicated(ocorrencias.especie[,"cells"])==T)
  ocorrencias.unicas<-rbind(ocorrencias.unicas,ocorrencias.especie[-duplicadas,])
}
rownames(ocorrencias.unicas)<-NULL
write.xlsx(ocorrencias.unicas[,c(1:3)],"OcorrenciasUnicasPoliquetas.xls",col.names=T,row.names=F)

-10.48770  -70.869200

#3.Arquivo final para SAM
composicao<-read.table(file.choose(),h=T,sep="\t")
riqueza<-NULL

riqueza.celula<-composicao[,c(1,2,26)]
poli.sam.completo<-cbind(tabela.env,riqueza.celula)
poli.sam.dados<-poli.sam.completo[poli.sam.completo[,"riqueza.celula"]!=0,]

#4.Calcular distancia a Centros Especificos
distancia.USP<-NULL
distancia.UFPB<-NULL
distancia.MarPlata<-NULL

USP<-c(-23.83,-45.42)
UFPB<-c(-7.14,-34.85)
MarPlata<-c(-38.0,-57.6)

for (c in 1:nrow(poli.sam.dados)){
  distancia.USP[c]<-sqrt((USP[2]-poli.sam.dados[c,1])^2+(USP[1]-poli.sam.dados[c,2])^2)
  distancia.UFPB[c]<-sqrt((UFPB[2]-poli.sam.dados[c,1])^2+(UFPB[1]-poli.sam.dados[c,2])^2)
  distancia.MarPlata[c]<-sqrt((MarPlata[2]-poli.sam.dados[c,1])^2+(MarPlata[1]-poli.sam.dados[c,2])^2)
}

distancia.minima<-NULL
for (c in 1:length(distancia.USP)){
  distancia.minima[c]<-min(distancia.USP[c],distancia.UFPB[c],distancia.MarPlata[c])  
}

poli.sam.dados<-cbind(poli.sam.dados,distancia.minima)

write.table(poli.sam.dados,"Dados_Poliquetas_Sam.txt",sep="\t",row.names=F)

#4.Calculo do esforco amostral

ocorrencias<-read.table(file.choose(),h=T,sep=",")
ocorrencias.xy<-ocorrencias[,-1]
ocorrencias.var<-as.data.frame(extract(environmental,ocorrencias.xy,cellnumber=T))
ocorrencias.var<-cbind(ocorrencias,ocorrencias.var)
ocorrencias.crescente<-ocorrencias.var[order(ocorrencias.var[,"cells"], decreasing = FALSE), ]
ocorrencias.var<-ocorrencias.crescente
celula.numero<-unique(ocorrencias.var[,"cells"])
celula.xy<-unique(ocorrencias.var[,c("Long","Lat")])
esforco.celula<-NULL

for (d in 1:length(celula.numero)){
  esforco.celula[d]<-nrow(ocorrencias.var[ocorrencias.var[,"cells"]==celula.numero[d],])
}

esforco<-matrix(NA,nrow=length(esforco.celula),ncol=2)
esforco[,1]<-celula.numero
esforco[,2]<-esforco.celula
esforco<-as.data.frame(esforco)
colnames(esforco)<-c("cells","EsforcoAmostral")
esforco<-cbind(riqueza.celula,esforco)
esforco<-esforco[,c(1,2,5)]

write.table(esforco,"Esforco_Amostral_Poliquetas.txt",sep="\t",row.names=F)
