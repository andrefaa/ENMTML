rm(list=ls())
redacao<-read.table(file.choose(),header=F,sep="\t")
nrow(redacao)/5

for (a in 1:(nrow(redacao)/5)){
  if (a==1){
    redacao[a,2]<-redacao[2,1]
  }
  if(a>1){
  redacao[a,2]<-redacao[(2+5*(a-1)),1]
  }
}
  
for (b in 1:(nrow(redacao)/5)){
  if (b==1){
    redacao[b,3]<-redacao[3,1]
  }
  if(b>1){
  redacao[b,3]<-redacao[(3+5*(b-1)),1]
  }
}

for (c in 1:(nrow(redacao)/5)){
  if (c==1){
    redacao[c,4]<-redacao[4,1]
  }
  if(c>1){
  redacao[c,4]<-redacao[(4+5*(c-1)),1]
  }
}

for (d in 1:(nrow(redacao)/5)){
  if (d==1){
    redacao[d,5]<-redacao[5,1]
  }
  if(d>1){
  redacao[d,5]<-redacao[(5+5*(d-1)),1]
  }
}

for (e in 1:(nrow(redacao)/5)){
  if (e==1){
    redacao[e,6]<-redacao[1,1]
  }
  if(e>1){
    redacao[e,6]<-redacao[(1+5*(e-1)),1]  
  }
}

redacao.na<-na.omit(redacao)
redacao.na[,1]<-redacao.na[,6]
redacao.na<-redacao.na[,-6]
redacao.na.pcd<-redacao.na[redacao.na[,5]!="Reprovado - PcD",]
redacao.na.pcd<-redacao.na.pcd[redacao.na.pcd[,5]!="Aprovado - PcD",]
redacao.na.alfab<-redacao.na.pcd[order(redacao.na.pcd$V2),]

row.names(redacao.na.aprovados)<-NULL
redacao.classificacao<-redacao.na.aprovados[order(redacao.na.aprovados$V4,decreasing=TRUE),]
redacao.classificacao[,6]<-1:nrow(redacao.classificacao)
redacao.classificacao$V4<-apply(apply(redacao.classificacao[,4], 2, gsub, patt=",", replace="."), 2, as.numeric)
redacao.classificacao.P2<-redacao.classificacao[,c(3,4)]
redacao.classificacao.P2<-apply(apply(redacao.classificacao.P2, 2, gsub, patt=",", replace="."), 2, as.numeric)
redacao.classificacao[,4]<-redacao.classificacao.P2[,2]

for (x in 1:nrow(redacao.classificacao)){
  redacao.classificacao[x,4]<-redacao.classificacao[x,4]*2
}

objetivas<-read.table(file.choose(),header=T,sep="\t")
objetivas<-redacao[,-c(9,10)]
objetivas<-objetivas[objetivas[,8]=="SIM",]
objetivas.dot<-objetivas[,c(3:6)]
objetivas.dot<-apply(apply(objetivas.dot, 2, gsub, patt=",", replace="."), 2, as.numeric)
objetivas[,c(3:6)]<-objetivas.dot

redacao.classificacao.alfab<-redacao.classificacao[order(redacao.classificacao$V2),]
objetivas.alfab<-objetivas[order(objetivas$Nome),]
objetivas.alfab<-objetivas.alfab[,c(1:3,5)]
redacao.na.alfab<-redacao.na.alfab[,c(1:4)]
redacao.na.alfab<-redacao.na.alfab[,-3]
redacao.na.alfab<-redacao.na.alfab[-c(307),]
classificacao<-cbind(redacao.na.alfab,objetivas.alfab)
classificacao<-classificacao[,-c(4:5)]
classificacao.dot<-apply(apply(classificacao, 2, gsub, patt=",", replace="."), 2, as.numeric)7
classificacao[,c(3:5)]<-classificacao.dot[,c(3:5)]
classificacao[,6]<-0
classificacao[,6]<-((classificacao[,3]*2)+(classificacao[,4]*1)+(classificacao[,5]*3))/6
classificacao.decres<-classificacao[order(classificacao$V6,decreasing=TRUE),]
classificacao.decres[,7]<-1:nrow(classificacao.decres)
