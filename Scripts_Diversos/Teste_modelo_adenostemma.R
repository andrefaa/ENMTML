setwd("C:/Users/Andre/Desktop/Andre/Mestrado/Ocorrencias/Ocorrencias_CSV")
library(XML)

Adenostemma<-read.table("Adenostemma_bras_mod.csv", sep=",", header=TRUE)
coord<-data.frame(Adenostemma$lon,Adenostemma$lat)


setwd("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma")
var<-list.files(pattern=".bil")
rel<-list.files("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/Teste Adenostemma/mosaic_nodata")
r<-raster("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/Teste Adenostemma/mosaic_nodata/w001001.adf")

bio5_34<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio5_34.bil")
bio5_44<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio5_44.bil")
bio5<-merge(bio5_34,bio5_44)

bio6_34<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio6_34.bil")
bio6_44<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio6_44.bil")
bio6<-merge(bio6_34,bio6_44)

bio13_34<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio13_34.bil")
bio13_44<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio13_44.bil")
bio13<-merge(bio13_34,bio13_44)

bio14_34<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio14_34.bil")
bio14_44<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio14_44.bil")
bio14<-merge(bio14_34,bio14_44)

bio15_34<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio15_34.bil")
bio15_44<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio15_44.bil")
bio15<-merge(bio15_34,bio15_44)

bio18_34<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio18_34.bil")
bio18_44<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio18_44.bil")
bio18<-merge(bio18_34,bio18_44)

bio19_34<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio19_34.bil")
bio19_44<- raster("C:/Andre_Macroeco/climatic layers/Modelo Adenostemma/bio19_44.bil")
bio19<-merge(bio19_34,bio19_44)


biocl<-stack(bio5,bio6,bio13,bio14,bio15,bio18,bio19)
bioclcr<- crop(biocl,extentras)
relcr<- crop (r,e)

resamp<- raster(extentras, nrows=19200, ncols=14400, crs=relcr@crs)

bioclim<-resample(bioclcr, resamp, method="ngb")

variaveis<- addLayer(bioclim,relcr)

#Dividir a amostra em 5 partes, 4 para calibrar, 1 para validar
fold<- kfold(coord, k=5) 
cbind(coord,fold)

#Selecionar todas as Adenostemmaenadas, menos as que sejam do grupo 1, 80%(diferentes de = "!")
occtrain<- coord[fold !=1,]

#Selecionar 20% das ocorrencias para teste (que sejam igual a = "==")
occtest<- coord[fold ==1,]

#Extent das vari?veis
e<- extent(-60, -30, -50, -10)
extentras<- alignExtent(e, relcr,snap="near")

#rJava funciona com r em versao 32-bits (mudar em Tools->Global Options->Rversion)
library(rJava)

#Rodar o modelo maxent,da como resposta a contribuicao das variáveis
max_aden<-maxent(variaveis,occtrain)

library(maps)
#Criar o mapa da previsão, predict
map_aden_maxent<-predict(max_aden,variaveis)
plot(map_aden_maxent,col=mypallet(100), colNA= "light blue", axes=T, box=F,xlab="Longitude", ylab="Latidude")
map(add=T,fill=F,xlim=c(-60,-30),ylim=c(-50,-10))
#title("Myrmecophaga tridactyla Present MaxEnt")
points(Adenostemma$lon,Adenostemma$lat, col="black",pch=19,cex=0.7)
savePlot("Myrmecophaga_maxmodelo_atual.tiff","tiff")

#Para ver os gráficos da resposta de cada variável
response(max_aden)

#Para alguma coisa..huehuehue
str(max_aden)

#Para ver a tabela de resultados do maxent
max_aden@results

#Para ver o valor de AUC de treino
training.AUC<-max_aden@results[5]

#Threshold com menos erros
thr<-max_aden@results[38]

#Threshold para espécies raras, presença com o menor valor
thr2<-max_aden@results[60]

#Reclassificar a partir do threshold escolhido, no caso utilizando o thr(pode ser inserido o valor desejado)
pa<-reclassify (map_aden_maxent, c(-Inf, thr, 0, thr, 1,1))
plot(pa,col=mypallet2(100),colNA= "light blue", axes=T, box=F,xlab="Longitude", ylab="Latidude")
map(add=T,fill=F,xlim=c(-90,-20),ylim=c(-50,15))
#title("Myrmecophaga tridactyla Present MaxEnt")
points(Adenostemma$lon,Adenostemma$lat, col="black",pch=19,cex=0.7)
savePlot("Myrmecophaga_maxent_atual2.tiff","tiff")

#Testar com dados externos
pseudoausencias<-randomPoints(variaveis,100)
e1<-evaluate(max_aden, p=occtest, a=pseudoausencias,x=variaveis)
