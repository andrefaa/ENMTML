# M_SDM a priori
# Feito por Poliana Mendes, Santiago Argentino e Andr? Mineiro
# Objetivo: gerar vari?veis asc que representam o espa?o para serem utilizadas na modelagem.
install.packages("raster")
install.packages("flexclust")
install.packages("SDMTools")

#M?todo 1: X_Y, cria um asc de longitude e outro de latitude

#definir pasta para salvar arquivos longitude e latitude
setwd("C:\\Users\\PONDS01\\Dropbox\\pos_doc\\FASE1\\SPS_SIMU_TESTE\\SDMs\\M-SDM\\layer\\XY")
#carregar pacotes necess?rios
library(SDMTools)
library(raster)
print("escolha um asc ambiental para servir como m?scara")
#abrir um asc ambiental para servir como m?scara
asc<-read.asc(file.choose())
#lst<-list.files(pattern=".asc")
#asc<-read.asc(lst[1])
ras<-raster(asc)
rsd<-as.data.frame(ras,xy=T,centroids=T)
long<-rsd
lat<-rsd
lins<-which(is.na(long[,3])==F)

long[lins,3]<-long[lins,1]
gridded(long)<-~x+y
long<-raster(long)
writeRaster(long,"Longitude.asc",format="ascii")

lat[lins,3]<-lat[lins,2]
gridded(lat)<-~x+y
lat<-raster(lat)
writeRaster(lat,"Latitude.asc",format="ascii")


#library(raster)
# abra um raster de uma outra vari?vel ambiental que vai ser utilizada na modelagem
#var<-raster(file.choose())
#escolha um diretorio para salvar as vari?veis
#setwd(choose.dir())
#Encontrando o tamanho do seu raster
#row<-nrow(var)
#col<-ncol(var)
#criando uma matrix com mesmo tamanho de linhas e colunas do seu raster
#lat<-matrix(0,row,col)
#long<-matrix(0,row,col)
#calculando a latitude do centroide de cada pixel do raster
#for (l in 1:row){
#  lat[l,]=(row-l+0.5)*yres(var)+ymin(var)
#  lat2<-var
#  lat2[!is.na(lat2[])]<-lat
#  }
#writeRaster(lat2,'latitude.asc',format="ascii",overwrite=TRUE)
#calculando a longitude do centroide de cada pixel do raster
#for(c in 1:col){
 # long[,c]=(c-0.5)*xres(var)+xmin(var)
  #long2<-var
  #long2[!is.na(long2[])]<-long
#}
#writeRaster(long2,'longitude.asc',format="ascii",overwrite=TRUE)

#Metodo 2- Allouche M?todo das dist?ncias m?nimas, cria um asc de dist?ncias at? o ponto de ocorr?ncia mais pr?ximo para todas as esp?cies
require(flexclust)
library(raster)
#abra a tabela com pontos de ocorrencia das especies, que est? em txt delimitado por tabula??o, 
#essa tabela deve conter tr?s colunas, Species, long, lat 
Species<-read.table(file.choose(),head=T)
#encontre o nome das esp?cies
namesp<-levels(Species[,1])
#abra um raster ambiental que ser? utilizado na modelagem
#var<-ras
var<-raster(file.choose())
#escolha o diret?rio para salvar seu asc
setwd(choose.dir())
#setwd('C:\\Users\\PONDS01\\Dropbox\\pos_doc\\FASE1\\SPS_SIMU_TESTE\\SDMs\\M-SDM\\layer\\MIN')
for (b in 1:length(levels(Species[,1]))){
  spi<-as(var,'SpatialPixels')@coords
  r1<-rasterize(Species[Species[,1]==namesp[b],2:3],var,field=1)
  r1<-as(r1,'SpatialPixels')@coords
  #calcula a distancia de todos os centroides dos pixels aos pontos de ocorrencia
  distr<-dist2(spi,r1,method='euclidean',p=2)
  #d? o valor minimo de distancia
  distr2<-apply(distr,1,min)
  distr3<-(distr2-min(distr2))/(max(distr2)-min(distr2))
  spdist<-var
  spdist[!is.na(spdist[])]<-distr3
  writeRaster(spdist,paste('mindist_',namesp[b],'.asc',sep=""),format="ascii")
}

#Metodo 3: Allouche distancia cumulativa
require(flexclust)
library(raster)
#tabela com pontos de ocorrencia das especies
Species<-read.table(file.choose(),head=T)
namesp<-levels(Species[,1])
#abra um raster ambiental que ser? utilizado na modelagem
var<-raster(file.choose())
#escolha o diret?rio para salvar seu asc
setwd(choose.dir())
#setwd('C:\\Users\\PONDS01\\Dropbox\\pos_doc\\FASE1\\SPS_SIMU_TESTE\\SDMs\\M-SDM\\layer\\CUMULATIVE')
spi<-as(var,'SpatialPixels')@coords
for (b in 1:length(levels(Species[,1]))){
  r1<-rasterize(Species[Species[,1]==namesp[b],2:3],var,field=1)
  r1<-as(r1,'SpatialPixels')@coords
  #calcula a distancia de todos os centroides dos pixels aos pontos de ocorrencia
  distr<-dist2(spi,r1,method='euclidean',p=2)
  distr2<-distr
  distr2<-distr+1
  #d? o valor cumulativo de distancia
  distr3<-1/(1/distr2^2)
  distr4<-apply(distr2,1,sum)
  distr5<-(distr4-min(distr4))/(max(distr4)-min(distr4))
  spdist<-var
  spdist[!is.na(spdist[])]<-distr5
  writeRaster(spdist,paste(namesp[b],'.asc',sep=""),format="ascii",overwrite=T)
}


#Metodo 4 Allouche  kernel-gaussiano
require(flexclust)
library(raster)
#tabela com pontos de ocorrencia das especies
Species<-read.table(file.choose(),head=T)
namesp<-levels(Species[,1])
#abra um raster ambiental que ser? utilizado na modelagem
var<-raster(file.choose())
print('escolha o diret?rio com as adequabilidades')
setwd(choose.dir())
setwd('C:\\Users\\PONDS01\\Dropbox\\pos_doc\\FASE1\\SPS_SIMU_TESTE\\SDMs\\M-SDM\\layer\\GAUSSIAN')
spi<-as(var,'SpatialPixels')@coords
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
  #d? o valor medio de distancia
  distr2<-distr
  distr2<-(1/sqrt(2*pi*sd_graus)*exp(-1*(distr/(2*sd_graus^2))))
  distr3<-apply(distr2,1,sum)
  spdist<-var
  spdist[!is.na(spdist[])]<-distr3
  writeRaster(spdist,paste('gauss_',namesp[b],'.asc',sep=""),format="ascii")
  }


#M-SDM a posteriori

#carregando pacotes
library(raster)
library(flexclust)
library(SDMTools)
library(GISTools)
library(maptools)
library(sp)
library(rgdal)
library(rgeos)

#preparando os dados 

  #essa tabela deve conter tr?s colunas, Species, long, lat 
  Species<-read.table(file.choose(),head=T)
  #encontre o nome das esp?cies
  namesp<-levels(Species[,1])
  #encontre a tabela dos limiares de corte do seu modelo
  ths<-read.table(file.choose(),head=T)
  print('escolha pasta com asc de adequabilidade')
  setwd(choose.dir())
  #abra os arquivos com a adequabilidade
  list<-list.files(pattern='.asc')
  adeqs<-stack(list.files(pattern='.asc'))

#M?todo 1- pixel 

  print('escolha pasta para salvar')
  setwd(choose.dir())
  for (b in 1:length(levels(Species[,1]))){
  th<-ths[b,6]
  #tentar fazer isso autom?tico
  #print('abrir adeq para cada esp?cie')
  adeq<-subset(adeqs,b)
  adeq3<-reclassify(adeq,c(-Inf,th,NA,th,1,1))
  spraster<-rasterize(Species[Species[,1]==namesp[b],2:3],adeq3,field=1)
  sps<-as(spraster,'SpatialPixels')@coords
  spi<-as(adeq3,'SpatialPixels')@coords
  dist<-dist2(spi,sps,method='euclidean',p=2)
  distmin<-apply(dist,1,min)
  distmini<--1*distmin
  dist3<-(distmini-min(distmini))/(max(distmini)-min(distmini))
  spdist<-adeq3
  spdist[!is.na(spdist[])]<-dist3
  #criando os cen?rios
    cen25<-reclassify(spdist,c(-Inf,0.25,NA,0.25,1,1))
    cen50<-reclassify(spdist,c(-Inf,0.50,NA,0.50,1,1))
    cen75<-reclassify(spdist,c(-Inf,0.75,NA,0.75,1,1))
    cen95<-reclassify(spdist,c(-Inf,0.95,NA,0.95,1,1))
    writeRaster(cen25,paste('pixel_cen25_',namesp[b],'.asc',sep=""),format="ascii",overwrite=TRUE) 
    writeRaster(cen50,paste('pixel_cen50_',namesp[b],'.asc',sep=""),format="ascii",overwrite=TRUE) 
    writeRaster(cen75,paste('pixel_cen75_',namesp[b],'.asc',sep=""),format="ascii",overwrite=TRUE) 
    writeRaster(cen95,paste('pixel_cen95_',namesp[b],'.asc',sep=""),format="ascii",overwrite=TRUE) 
}

#Metodo 2: Paisagem


print('escolha pasta para salvar')
setwd(choose.dir())
for (b in 1:length(levels(Species[,1]))){
  th<-ths[b,6]
  adeq<-subset(adeqs,b)
  adeq3<-reclassify(adeq,c(-Inf,th,NA,th,1,1))
  #adeq4<-projectRaster(adeq3,crs= "+proj=merc +lat_0=0 +lon_0=-80")
  pts1<-SpatialPointsDataFrame(Species[Species[,1]==namesp[b],2:3],Species[Species[,1]==namesp[b],],proj4string=CRS("+proj=longlat +datum=WGS84"), bbox = NULL)
  #pts2<-spTransform(pts1,CRS("+proj=merc +lat_0=0 +lon_0=-80"))
  polygon<-rasterToPolygons(adeq3, fun=NULL, n=8, na.rm=TRUE, digits=12, dissolve=TRUE)
  polygon2 <- disaggregate(polygon)
  #plot(polygon2)
  #plot(pts1, add=T)
  polypoint<-intersect(polygon2,pts1)
  polypoint<-aggregate(polypoint)
  #plot(polypoint)
  a<-gDistance(polygon2,polypoint, byid=TRUE)
  distmini<--1*a
  dist3<-(distmini-min(distmini))/(max(distmini)-min(distmini))
  polygon2$NewAT<-dist3[1:length(dist3)]
  fist<-rasterize(polygon2,adeq3,field=polygon2$NewAT)
  cen25<-reclassify(fist,c(-Inf,0.25,NA,0.25,1,1))
  cen50<-reclassify(fist,c(-Inf,0.50,NA,0.50,1,1))
  cen75<-reclassify(fist,c(-Inf,0.75,NA,0.75,1,1))
  cen95<-reclassify(fist,c(-Inf,0.95,NA,0.95,1,1))
  writeRaster(cen25,paste('paisagem_cen25_',namesp[b],'.asc',sep=""),format="ascii",overwrite=TRUE) 
  writeRaster(cen50,paste('paisagem_cen50_',namesp[b],'.asc',sep=""),format="ascii",overwrite=TRUE) 
  writeRaster(cen75,paste('paisagem_cen75_',namesp[b],'.asc',sep=""),format="ascii",overwrite=TRUE) 
  writeRaster(cen95,paste('paisagem_cen95_',namesp[b],'.asc',sep=""),format="ascii",overwrite=TRUE)
  }

# Metodo grafos (Saura 2006)
print('escolha pasta para salvar')
setwd(choose.dir())
for (b in 1:length(levels(Species[,1]))){
  th<-ths[b,6]
  adeq<-subset(adeqs,b)
  adeq3<-reclassify(adeq,c(-Inf,th,NA,th,1,1))
  pts1<-SpatialPointsDataFrame(Species[Species[,1]==namesp[b],2:3],Species[Species[,1]==namesp[b],],proj4string=CRS("+proj=longlat +datum=WGS84"), bbox = NULL)
  polygon<-rasterToPolygons(adeq3, fun=NULL, n=8, na.rm=TRUE, digits=12, dissolve=TRUE)
  polygon2 <- disaggregate(polygon)
  #area cada patch
  areai<-gArea(polygon2,byid=TRUE)
  #areatotal
  areaT<-gArea(polygon2,byid=FALSE)
  #desagregando o j para calcular a ?rea
  polypoint<-intersect(polygon2,pts1)
  polypoint<-aggregate(polypoint)
  polypoint2<-disaggregate(polypoint)
  areaj<-gArea(polypoint2,byid=TRUE)
  polygon2$NewAT<-1:length(areai)
  polypoint2$NewAT<-1:length(areaj)
  res<-matrix(data=NA,nrow=length(areai),ncol=length(areaj))
  for (c in 1:length(areai)){
  for (d in 1:length(areaj)){
    polygon3<-polygon2
    polygon3<-polygon3[polygon3$NewAT==c,]
    polypoint3<-polypoint2
    polypoint3<-polypoint2[polypoint3$NewAT==d,]
    dij=gDistance(polygon3,polypoint3,byid=TRUE)
    pc=(areai[c]*areaj[d]*dij*-1)
    res[c,d]=pc
  }}
  resT<-apply(res,1,sum)
  resti<-(resT-min(resT))/(max(resT)-min(resT))
  polygon2$NewAD<-resti[1:length(resti)]
  fist<-rasterize(polygon2,adeq3,field=polygon2$NewAT)
  cen25<-reclassify(fist,c(-Inf,0.25,NA,0.25,1,1))
  cen50<-reclassify(fist,c(-Inf,0.50,NA,0.50,1,1))
  cen75<-reclassify(fist,c(-Inf,0.75,NA,0.75,1,1))
  cen95<-reclassify(fist,c(-Inf,0.95,NA,0.95,1,1))
  writeRaster(cen25,paste('grafos_cen25_',namesp[b],'.asc',sep=""),format="ascii",overwrite=TRUE) 
  writeRaster(cen50,paste('grafos_cen50_',namesp[b],'.asc',sep=""),format="ascii",overwrite=TRUE) 
  writeRaster(cen75,paste('grafos_cen75_',namesp[b],'.asc',sep=""),format="ascii",overwrite=TRUE) 
  writeRaster(cen95,paste('grafos_cen95_',namesp[b],'.asc',sep=""),format="ascii",overwrite=TRUE)
}

######



