library(raster)

mat1<-read.table(file.choose(),h=T,sep='\t') #Nova matriz de composicao de especies
mat2<-read.table(file.choose(),h=T,sep=',') #Matriz com os dados de riqueza

div<-rep(NA,nrow(mat1))
mat1<-cbind(mat1,div)

lin<-which(is.na(mat1[,238])==F) #A coluna 238 aqui e a coluna de qualquer especie
div<-mat2[,240] #A coluna 240 e a coluna com os valores de diversidade

mat1[lin,239]<-div
mat.gr<-mat1[,c(1,2,239)]
mat.gr[mat.gr == -9999] <- NA #Volta os -9999 a valer NA
gridded(mat.gr)<-~x+y
mat.gr<-raster(mat.gr)
plot(mat.gr)

setwd("C:\\Users\\decoa\\Desktop\\Livia")
writeRaster(mat.gr,"Diversidade_Presente.asc",format="ascii")
