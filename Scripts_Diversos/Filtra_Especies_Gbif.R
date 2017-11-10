library(dismo)
library(maps)


spBR<-read.table(file.choose(),sep="\t",header=T) #Ocorrencias DB.Karla
spSL<-read.table(file.choose(),sep="\t",header=T)
spo1<-read.table(file.choose(),sep="\t",header=T)
spo2<-read.table(file.choose(),sep="\t",header=T)
spBR<-spBR[,c(1,3,2)]
uBR<-unique(spBR[,1])

poli.w<-NULL
for (x in 1:length(uBR)){
  print(as.character(uBR[x]))
  poli<-try(gbif(as.character(uBR[x])),silent=F)
  if ('try-error' %in% class(poli)) next
  else poli<-poli[,which(names(poli) %in% c("species","lon","lat"))]
  poli.w<-rbind(poli.w,poli)
  print(x/length(uBR))
}
spWR<-na.omit(poli.w)
spWR<-spWr[,c(3,2,1)]
colnames(spWR)<-c("Species","Long","Lat")
colnames(spBR)<-c("Species","Long","Lat")

spWR<-spWR[spWR[,1] %in% spBR[,1],] #Ocorrencias GBIF
spSL<-spSL[spSL[,1] %in% spBR[,1],] #Ocorrencias SpeciesLink
spo1<-spo1[spo1[,1] %in% spBR[,1],] #Ocorrencias Obis1
spo2<-spo2[spo2[,1] %in% spBR[,1],] #Ocorrencias Obis2
spOB<-spOB[spOB[,1] %in% spBR[,1],]

spOB<-rbind(spo1,spo2)              #Ocorrencias OBIS

sps<-unique(spWR[,1])
spsBR<-as.character(unique(spBR[,1]))

setwd("C:\\OneDrive\\Poliquetas\\Art1-Lacunas\\Imagens_Dist")

for (y in 1:length(spsBR)){
  tiff(file = paste(spsBR[y],".tiff",sep=""), width = 800, height = 800, units = "px", res = 200)
  map()
  points(spBR[spBR==spsBR[y],2:3],pch=15,col="green",cex=1)
  points(spWR[spWR==spsBR[y],2:3],pch=16,col="red",cex=1)
  points(spSL[spSL==spsBR[y],2:3],pch=17,col="blue",cex=1)
  points(spOB[spOB==spsBR[y],2:3],pch=18,col="orange",cex=1)
  title(spsBR[y])
  dev.off()
}

spBR<-cbind(spBR,rep("NB_e_Nosso",nrow(spBR)))
colnames(spBR)<-c("Species","Long","Lat","Source")
spOB<-cbind(spOB,rep("OBIS",nrow(spOB)))
colnames(spOB)<-colnames(spBR)
spWR<-cbind(spWR,rep("WorldClim",nrow(spWR)))
colnames(spWR)<-colnames(spBR)
spSL<-cbind(spSL,rep("SpeciesLink",nrow(spSL)))
colnames(spSL)<-colnames(spBR)
sps.all<-rbind(spBR,spWR,spOB,spSL)
sp.all<-sps.all[order(sps.all[,1]),]

write.table(sp.all,"Occorrencias_Filtradas.txt",row.names=F,sep="\t")
