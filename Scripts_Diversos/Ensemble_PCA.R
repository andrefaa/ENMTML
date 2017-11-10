Ensemble_PCA<- function(diretorio){
  #diretorio:diretorio dos modelos de adequabilidade
  
  library(raster)
  library(SDMTools)
  setwd(diretorio)

#1.Ler os mapas de adequabilidade e cria um stack

adeq_prep<-list.files(pattern='.asc')
asc_base<-read.asc(adeq_prep[1])
adequabilidade<-raster(asc_base)
camadas<-adeq_prep[-1]

for (a in 1:length(camadas)){
  asc<-read.asc(camadas[a])
  raster<-raster(asc)
  adequabilidade<-stack(adequabilidade,raster)
}
names(adequabilidade)<-substr(adeq_prep,1,nchar(adeq_prep)-4)

#Z-transform nos mapas de adequabilidade
adequabilidade<-(adequabilidade-minValue(adequabilidade))/(maxValue(adequabilidade)-minValue(adequabilidade))

#2.Realiza a PCA desses mapas

data.frame<-rasterToPoints(adequabilidade)
data.frame<-na.omit(data.frame)
pca.raw<-data.frame[,-c(1:2)]

# Realizar a PCA 
data.pca <- prcomp(pca.raw,retx=TRUE,scale=T)

#Salvar os coeficientes
coeficientes<-data.pca$rotation
write.table(coeficientes,"result_coeff.txt",row.names=T,sep="\t")

#Salvar os Biplots
tiff(file = "BiplotPCA.tiff", width = 3200, height = 3200, units = "px", res = 800)
plot(coeficientes[,1],coeficientes[,2],
     xlim=c((min(coeficientes[,1])-0.5),(max(coeficientes[,1])+0.5)),
     ylim=c((min(coeficientes[,2])-0.5),(max(coeficientes[,2])+0.5)),
     xlab="PC1",ylab="PC2")
text(coeficientes[,1],coeficientes[,2],labels=rownames(coeficientes),pos=3)
arrows(0,0,coeficientes[,1],coeficientes[,2],col="red")
dev.off()

tiff(file = "BiplotPCA_1x3.tiff", width = 3200, height = 3200, units = "px", res = 800)
plot(coeficientes[,1],coeficientes[,3],
     xlim=c((min(coeficientes[,1])-0.5),(max(coeficientes[,1])+0.5)),
     ylim=c((min(coeficientes[,3])-0.5),(max(coeficientes[,3])+0.5)),
     xlab="PC1",ylab="PC3")
text(coeficientes[,1],coeficientes[,3],labels=rownames(coeficientes),pos=3)
arrows(0,0,coeficientes[,1],coeficientes[,3],col="red")
dev.off()

#Salvar e plotar os eixos
eigenvectors<-data.pca$x
eigenvectors.xy<-cbind(data.frame[,c(1:2)],eigenvectors)
pc1<-eigenvectors.xy[,1:3]
pc2<-eigenvectors.xy[,c(1:2,4)]
pc1<-data.frame(pc1)
pc2<-data.frame(pc2)
gridded(pc1)<-~x+y
gridded(pc2)<-~x+y
pc1.raster<-raster(pc1)
pc2.raster<-raster(pc2)

writeRaster(pc1.raster,"PC1.asc",format="ascii")
writeRaster(pc2.raster,"PC2.asc",format="ascii")


#Porcentagem da variacaoo explicada
n.eixos<-length(summary(data.pca)$importance[3,])
cumulat.var<-summary(data.pca)$importance[3,]
variacao.explicada<-data.frame(cumulat.var)
write.table(variacao.explicada,'result_latent.txt',row.names=T,sep='\t')

}

Ensemble_PCA("C:\\Poliquetas\\Poliquetas_Pacifico\\Resultados\\TodasOcorrencias\\PCA")


