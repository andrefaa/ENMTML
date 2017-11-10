# M_SDM a priori
# Feito por Poliana Mendes e Santiago Argentino
# Objetivo: gerar variáveis asc que representam o espaço para serem utilizadas na modelagem.

#Método 1: X_Y, cria um asc de longitude e outro de latitude
library(raster)
# abra um raster de uma outra variável ambiental que vai ser utilizada na modelagem
var<-raster(file.choose())
#escolha um diretorio para salvar as variáveis
setwd(choose.dir())
#Encontrando o tamanho do seu raster
row<-nrow(var)
col<-ncol(var)
#criando uma matrix com mesmo tamanho de linhas e colunas do seu raster
lat<-matrix(0,row,col)
long<-matrix(0,row,col)
#calculando a latitude do centroide de cada pixel do raster
for (l in 1:row){
  lat[l,]=(row-l+0.5)*yres(var)+ymin(var)
  lat2<-var
  lat2[!is.na(lat2[])]<-lat
  }
writeRaster(lat2,'latitude.asc',format="ascii",overwrite=TRUE)
#calculando a longitude do centroide de cada pixel do raster
for(c in 1:col){
  long[,c]=(c-0.5)*xres(var)+xmin(var)
  long2<-var
  long2[!is.na(long2[])]<-long
}
writeRaster(long2,'longitude.asc',format="ascii",overwrite=TRUE)

#Metodo 2- Allouche Método das distâncias mínimas, cria um asc de distâncias até o ponto de ocorrência mais próximo para todas as espécies
require(flexclust)
library(raster)
#abra a tabela com pontos de ocorrencia das especies, que está em txt delimitado por tabulação, 
#essa tabela deve conter três colunas, Species, long, lat 
Species<-read.table(file.choose(),head=T)
#encontre o nome das espécies
namesp<-levels(Species[,1])
#abra um raster ambiental que será utilizado na modelagem
var<-raster(file.choose())
#escolha o diretório para salvar seu asc
setwd(choose.dir())
for (b in 1:length(levels(Species[,1]))){
  spi<-as(var,'SpatialPixels')@coords
  r1<-rasterize(Species[Species[,1]==namesp[b],2:3],var,field=1)
  r1<-as(r1,'SpatialPixels')@coords
  #calcula a distancia de todos os centroides dos pixels aos pontos de ocorrencia
  distr<-dist2(spi,r1,method='euclidean',p=2)
  #dá o valor minimo de distancia
  distr2<-apply(distr,1,min)
  spdist<-var
  spdist[!is.na(spdist[])]<-distr2
  writeRaster(spdist,paste('mindist_',namesp[b],'.asc',sep=""),format="ascii")
}

#Metodo 3: Allouche distancia cumulativa
require(flexclust)
library(raster)
#tabela com pontos de ocorrencia das especies
Species<-read.table(file.choose(),head=T)
namesp<-levels(Species[,1])
#abra um raster ambiental que será utilizado na modelagem
var<-raster(file.choose())
#escolha o diretório para salvar seu asc
setwd(choose.dir())
for (b in 1:length(levels(Species[,1]))){
  spi<-as(var,'SpatialPixels')@coords
  r1<-rasterize(Species[Species[,1]==namesp[b],2:3],var,field=1)
  r1<-as(r1,'SpatialPixels')@coords
  #calcula a distancia de todos os centroides dos pixels aos pontos de ocorrencia
  distr<-dist2(spi,r1,method='euclidean',p=2)
  #dá o valor cumulativo de distancia
  distr2<-1/(distr^2)
  distr2<-apply(distr2,1,sum)
  spdist<-var
  spdist[!is.na(spdist[])]<-distr2
  writeRaster(spdist,paste('cumudist_',namesp[b],'.asc',sep=""),format="ascii")
  
}

#Metodo 4 Allouche  kernel-gaussiano
require(flexclust)
library(raster)
#tabela com pontos de ocorrencia das especies
Species<-read.table(file.choose(),head=T)
namesp<-levels(Species[,1])
#abra um raster ambiental que será utilizado na modelagem
var<-raster(file.choose())
#escolha o diretório para salvar seu asc
setwd(choose.dir())
for (b in 1:length(levels(Species[,1]))){
  spi<-as(var,'SpatialPixels')@coords
  r1<-rasterize(Species[Species[,1]==namesp[b],2:3],var,field=1)
  r1<-as(r1,'SpatialPixels')@coords
  #calcula a distancia de todos os centroides dos pixels aos pontos de ocorrencia
  distr<-dist2(spi,r1,method='euclidean',p=2)
  #sd_gaus<- fazer distancia maxima entre distancias minimas entre pontos
  distp<-dist2(r1,r1,method='euclidean',p=2)
  distp1<-matrix(0,nrow(distp),1)
  for (c in 1:nrow(distp)) {
    vec<-distp[c,]
    distp1[c]<-min(vec[vec!=min(vec)])
  }
  sd_graus<-max(distp1)
  #dá o valor medio de distancia
  distr2<-distr
  distr2<-(1/sqrt(2*pi*sd_graus)*exp(-1*(distr/(2*sd_graus^2))))
  distr3<-apply(distr2,1,sum)
  spdist<-var
  spdist[!is.na(spdist[])]<-distr3
  writeRaster(spdist,paste('gauss_',namesp[b],'.asc',sep=""),format="ascii")
  }
