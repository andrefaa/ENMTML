setwd("C:\\Users\\decoa\\Dropbox\\Temperatura")
library(raster)

lista<-list.files(pattern=".tif")
occ<-read.table(file.choose(),header=T,sep="\t")
occ.xy<-occ[,3:2]

resu<-occ.xy
for (x in 1:length(lista)){
  print(x)
  ras<-raster(lista[x])
  temp<-extract(ras,occ.xy)
  resu<-cbind(resu,temp)
}
resu<-resu[,-c(1,2)]
colnames(resu)<-substr(lista,1,nchar(lista)-4)
resu<-cbind(occ[,2:3],resu)
write.table(resu,"Temp_Mensal_Localidades.txt",row.names=F,sep="\t")
