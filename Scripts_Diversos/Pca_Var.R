#Funcao para realizar a PCA nas variaveis ambientais
#Parametros
  #dir_input: Diretorio com as variaveis ambientais(.asc)
  #dir_out: Diretorio onde ser?o salvos os eixos da PCA (.asc)
#Outputs:
  #Eixos da PCA(.asc)
  #Correlacao das variaveis com os eixos(result_coef.txt)
  #Variacao explicada por cada eixo (result_latent.txt)

Pca_var<-function(dir_input,dir_out){

  library(raster)
  library(SDMTools)
  
setwd(dir_input)
  
env_prep<-list.files(pattern='.asc')
asc_base<-read.asc(env_prep[1])
environmental<-raster(asc_base)
camadas<-env_prep[-1]

for (a in 1:length(camadas)){
  asc<-read.asc(camadas[a])
  raster<-raster(asc)
  environmental<-stack(environmental,raster)
}
names(environmental)<-substr(env_prep,1,nchar(env_prep)-4)


#1.2.Realizar a PCA nas camadas ambientais(utilizar PCA como input)

  setwd(dir_out)

  data.frame<-rasterToPoints(environmental)
  pca.raw<- na.omit(data.frame)
  pca.raw.xy<-pca.raw[,1:2]
  pca.raw<-pca.raw[,-c(1:2)]
  
  #Scale transform 
  data.scaled <- data.frame(apply(pca.raw,2,scale))
  
  # Realizar a PCA 
  data.pca <- prcomp(data.scaled,retx=TRUE)
  
  #Salvar os coeficientes
  coeficientes<-data.pca$rotation
  write.table(coeficientes,'result_coef.txt',row.names=T,sep='\t')

  #Porcentagem da variacaoo explicada
  n.eixos<-length(summary(data.pca)$importance[3,])
  cumulat.var<-summary(data.pca)$importance[3,]
  variacao.explicada<-data.frame(cumulat.var)
  write.table(variacao.explicada,'result_latent.txt',row.names=T,sep='\t')

  
  #Recuperar os loadings e transformar em um data.frame
  eixos<-as.data.frame(data.pca$x)
  eixos.xy<-cbind(pca.raw.xy,eixos)
  gridded(eixos.xy)<- ~x+y
  environmental<-stack(eixos.xy)
  writeRaster(environmental,"PCA.asc",bylayer=T,suffix='numbers',format="ascii")
  
}  

Pca_var('C:\\Users\\decoa\\Desktop\\Disciplinas 2016_2\\Macroecologia\\Projeto\\Env\\Neo','C:\\Users\\decoa\\Desktop\\Disciplinas 2016_2\\Macroecologia\\Projeto\\Env\\PCA')
