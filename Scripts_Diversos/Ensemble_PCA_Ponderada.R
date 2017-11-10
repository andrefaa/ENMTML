Ensemble_PCA<- function(diretorio,dir_output){
  #diretorio:diretorio dos modelos de adequabilidade
  #dir_output: diretorio onde serao salvos os resultados
  
  library(raster)
  library(SDMTools)

#1.Importa tabela com as especies(Linhas) e os TSS dos diferentes algoritmos(Colunas)
  tss<-read.table(file.choose(),header=T,sep="\t")
  ssp<-unique(tss[,1])

#1.Ler os mapas de adequabilidade e cria um stack(por espécie)
for (a in 1:length(ssp)){
  setwd(diretorio)
  adeq_prep<-list.files(pattern= paste(ssp[a]))
  adeq_prep<-grep(".asc",adeq_prep,value=T)
  asc_base<-read.asc(adeq_prep[1])
  adequabilidade<-raster(asc_base)
  #Z-transform nos mapas de adequabilidade
  adequabilidade<-(adequabilidade-minValue(adequabilidade))/(maxValue(adequabilidade)-minValue(adequabilidade))
  camadas<-adeq_prep[-1]

  for (b in 1:length(camadas)){
    asc<-read.asc(camadas[b])
    raster<-raster(asc)
    
    #Z-transform nos mapas de adequabilidade
    raster<-(raster-minValue(raster))/(maxValue(raster)-minValue(raster))
    
    adequabilidade<-stack(adequabilidade,raster)
  }
  names(adequabilidade)<-substr(adeq_prep,1,nchar(adeq_prep)-4)


#Ponderar a PCA pelo valor de TSS de cada modelo
  adeq.pond<-stack()
  models<-tss[,-1]
  for (c in 1:ncol(models)){
    adeq.pond.temp<-adequabilidade[[c]]*models[a,c]
    adeq.pond<-stack(adeq.pond,adeq.pond.temp)
  }
  adequabilidade<-adeq.pond

#2.Realiza a PCA desses mapas

  data.frame<-rasterToPoints(adequabilidade)
  data.frame<-na.omit(data.frame)
  pca.raw<-data.frame[,-c(1:2)]

# Realizar a PCA 
  data.pca <- prcomp(pca.raw,retx=TRUE,scale=T)

#Salvar os coeficientes
  setwd(dir_output)
  coeficientes<-data.pca$rotation
  write.table(coeficientes,paste(ssp[a],"result_coeff.txt",sep="_"),row.names=T,sep="\t")

#Salvar os Biplots
  tiff(file = paste(ssp[a],"BiplotPCA.tiff",sep="_"), width = 3200, height = 3200, units = "px", res = 800)
  plot(coeficientes[,1],coeficientes[,2],
      xlim=c((min(coeficientes[,1])-0.5),(max(coeficientes[,1])+0.5)),
      ylim=c((min(coeficientes[,2])-0.5),(max(coeficientes[,2])+0.5)),
      xlab="PC1",ylab="PC2")
  text(coeficientes[,1],coeficientes[,2],labels=rownames(coeficientes),pos=3)
  arrows(0,0,coeficientes[,1],coeficientes[,2],col="red")
  dev.off()

  tiff(file = paste(ssp[a],"BiplotPCA_1x3.tiff",sep="_"), width = 3200, height = 3200, units = "px", res = 800)
  plot(coeficientes[,1],coeficientes[,3],
      xlim=c((min(coeficientes[,1])-0.5),(max(coeficientes[,1])+0.5)),
      ylim=c((min(coeficientes[,3])-0.5),(max(coeficientes[,3])+0.5)),
      xlab="PC1",ylab="PC3")
  text(coeficientes[,1],coeficientes[,3],labels=rownames(coeficientes),pos=3)
  arrows(0,0,coeficientes[,1],coeficientes[,3],col="red")
  dev.off()

#Salvar e plotar os eixos
  eigenvectors<-NULL
  for (d in 1:ncol(coeficientes)){
    if (sum(coeficientes[,d])<0){
      eigenvectors<-cbind(eigenvectors,(data.pca$x[,d])*-1)
    }else{
      eigenvectors<-cbind(eigenvectors,data.pca$x[,d])
    }
  }
  colnames(eigenvectors)<-colnames(coeficientes)
  eigenvectors.xy<-cbind(data.frame[,c(1:2)],eigenvectors)
  pc1<-eigenvectors.xy[,1:3]
  pc2<-eigenvectors.xy[,c(1:2,4)]
  pc1<-data.frame(pc1)
  pc2<-data.frame(pc2)
  gridded(pc1)<-~x+y
  gridded(pc2)<-~x+y
  pc1.raster<-raster(pc1)
  pc2.raster<-raster(pc2)

  writeRaster(pc1.raster,paste("ENS",ssp[a],"PC1.asc",sep="_"),format="ascii")
  writeRaster(pc2.raster,paste("ENS",ssp[a],"PC2.asc",sep="_"),format="ascii")


#Porcentagem da variacaoo explicada
  n.eixos<-length(summary(data.pca)$importance[3,])
  cumulat.var<-summary(data.pca)$importance[3,]
  variacao.explicada<-data.frame(cumulat.var)
  write.table(variacao.explicada,paste(ssp[a],'result_latent.txt',sep="_"),row.names=T,sep='\t')
}
}

Ensemble_PCA("C:\\Users\\decoa\\Desktop\\Daniel_Ens",
             "C:\\Users\\decoa\\Desktop\\Daniel_Ens")
