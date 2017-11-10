#ENFA
library(adehabitatHS)

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

ocorrencias.r<-read.table(file.choose(),h=T,sep='\t')
ocorrencias<-ocorrencias.r[,which(names(ocorrencias.r) %in% c("Species","Long","Lat"))]
ocorrencias.xy<-ocorrencias[,-1]

############## MATRIZ DUALIDADE & ENFA PARA DADOS NOSSOS ########################
env.data.frame<-as.data.frame(rasterToPoints(environmental))
env.data.frame.var<-as.data.frame(env.data.frame[,-c(1,2)])
env.data.frame.xy<-as.data.frame(env.data.frame[,c(1,2)])

spat.pixels<-SpatialPixelsDataFrame(env.data.frame.xy,data=env.data.frame.xy)
spat.points<-SpatialPointsDataFrame(ocorrencias.xy,data=ocorrencias.xy)

cp <- count.points(spat.points,spat.pixels) # Conta o numero de pontos por celula

pr <- slot(cp, "data")[,1]
histniche(env.data.frame.var, pr)

pc <- dudi.pca(env.data.frame.var, scannf=FALSE) # Cria a matriz de dualidade

en3 <- enfa(pc, pr, scan=FALSE,nf=6) # Realiza o ENFA
barplot(en3$s) # Barplot dos autovalores, decréscimo brusco nos eixos significa um grau de especializacao alto
scatter(en3) # Scatterplot do ENFA mostrando a relaco da espécie com o espaço ambiental
coeficientes<-en3$s #Autovalores de cada eixo(utilizados para projecao do resultado)
marginalidade.eixos<-en3$tab # Aqui o programa gera a tabela com os valores de marginalidade para cada variavel por celula

#Selecionar a marginalidade apenas para as presenças
pr.df<-as.data.frame(pr)
presenca<-which(pr.df!=0)
marginalidade.especie<-marginalidade.eixos[presenca,]

#Calculo da Mediana e da Distancia a Mediana para Divisao em classes
medianas<-NULL
for (u in 1:ncol(marginalidade.especie)){
  medianas[u]<-median(marginalidade.especie[,u])
}

dist.mediana<-marginalidade.especie
for (k in 1:length(medianas)){
  for (m in 1:nrow(marginalidade.especie)){
  dist.mediana[m,k]<-marginalidade.especie[m,k]-medianas[k]
  }
}

summary(dist.mediana)
cuts<-NULL
for (n in 1:ncol(dist.mediana)){
  cuts[[n]]<-list(seq(round(min(dist.mediana[,n])-0.5),round(max(dist.mediana[,n])+0.5), by = 0.2))
}

#Distancia de Cada celula para a mediana

dist.global<-NULL

for (o in 1:length(medianas)){
  dist.global.vect<-marginalidade.eixos[,o]
  dist.global.vect<-dist.global.vect-medianas[o]
  dist.global<-cbind(dist.global,dist.global.vect)
}
colnames(dist.global)<-colnames(marginalidade.eixos)

#Definir a classe na qual cada celula do mapa cai(TESTAR TAMBEM ENFA ALTERNATIVA 2 A PARTIR DAQUI!)
classe.cell<-NULL
classe.cell.mat<-NULL

for (w in 1:length(cuts)){
    cut<-unlist(cuts[[w]])
    dist.global.vect<-dist.global[,w]
  for (z in 1:length(dist.global.vect)){
     if(dist.global.vect[z]>max(cut) | dist.global.vect[z]<min(cut)){
       classe.cell[z]<-NA
     }else{
       classe.cell[z]<-cut[which(abs(cut-dist.global[z,w])==min(abs(cut-dist.global[z,w])))]
     }
  }
  classe.cell.mat<-cbind(classe.cell.mat,classe.cell)
}
colnames(classe.cell.mat)<-colnames(dist.global)

dist.med.matrix<-NULL
dist.n.cell<-NULL
for (g in 1:ncol(classe.cell.mat)){
  classe.cell.mat.vect<-classe.cell.mat[,g]
  classe.cell.mat.vect.na<-na.omit(classe.cell.mat.vect)
  for (h in 1:length(classe.cell.mat.vect)){
    if (is.na(classe.cell.mat.vect[h])==FALSE & classe.cell.mat.vect[h]>=0){
      classe.cell.mat.vect[h]<-sum(classe.cell.mat.vect.na>=classe.cell.mat.vect[h])
    }
    if (is.na(classe.cell.mat.vect[h])==FALSE & classe.cell.mat.vect[h]<0){
      classe.cell.mat.vect[h]<-sum(classe.cell.mat.vect.na<=classe.cell.mat.vect[h])
    }
    if (is.na(classe.cell.mat.vect[h])==TRUE){
      classe.cell.mat.vect[h]<-0
    }
  }
  print(paste("Coluna_",g,"_OK",sep=""))
  dist.med.matrix<-cbind(dist.med.matrix,classe.cell.mat.vect)
}

######################## FUNCIONANDO ATÉ AQUI(PARA CIMA) #################

gn <- gnesfa(pc, Focus = pr)
scatterniche(gn$li, pr, pts=TRUE)

en3 <- enfa(pc, pr, scan=FALSE,nf=7)
en1$s
barplot(en1$s)
scatter(en1)

RMD1 <- data.frame(RMD1=(en1$li[,1]^2)/max(en1$li[,1]^2))

coordinates(RMD1) <- coordinates(map)
gridded(RMD1) <- TRUE
image(RMD1)


li<-en1$li
RMD1 <- data.frame(RMD1=(en1$li))
coordinates(RMD1)<-coordinates(map)
gridded(RMD1)<-TRUE
par(mfrow=c(1,1))
for (a in 2:ncol(RMD1)){
  image(RMD1[,a])
}

contagem<-NULL
for (x in 1:nrow(li)){
  contagem[x]<-sum(li[,2]<=li[x,2])
}
median(li[,2])
hist(li[,1])
median(li[,1])

tabela<-as.data.frame(RMD1$RMD1)
tabela<-cbind(RMD1@coords,tabela)
gridded(tabela)~x+y
colnames(tabela[3])<-"z"
plot(tabela)
raster(tabela)
RMD1@proj4string
eigenv<-en1$s
eimar<-en1$mar

mad <- madifa(pc, pr, scan=FALSE)
li.m<-mad$li
?enfa
lw<-as.data.frame(en1$lw)


############ENiRG

lista_env<-list.files(pattern='.asc')
initGRASS("/usr/bin/grass-6.4.0", home=tempdir())
initGRASS("C:/GRASS", home=tempdir())
esse<-import.egvs(lista_env[19], "tann")


