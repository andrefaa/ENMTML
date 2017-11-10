setwd("C:/Users/Andre/Google Drive/Mestrado/Mestrado UFG/Modelos de Distribuicao-MatheusRibeiro/Projeto_RolaBosta")
rm(list=ls())

ocorrencias<-read.table(file.choose(),h=T,sep=",")
especies<-unique(ocorrencias[,1])
occ.treino<-NULL
occ.teste<-NULL



for (x in 1:4){  
  for (i in 1:length(especies)){
    
    ocorrencias.especie<-ocorrencias[ocorrencias[,1]==especies[i],]
    occ.grupo<-kfold(ocorrencias.especie,4)
    
    occ.treino.sp<-ocorrencias.especie[(occ.grupo!=x),]
    occ.teste.sp<-ocorrencias.especie[(occ.grupo==x),]
    occ.treino<-rbind(occ.treino,occ.treino.sp)
    occ.teste<-rbind(occ.teste,occ.teste.sp)
  }
  nome.treino<-paste("Ocorrencias_Treino_",x,".csv",sep="")
  nome.teste<-paste("Ocorrencias_Teste_",x,".csv",sep="")
  write.table(occ.treino,nome.treino,sep=",",row.names=F)
  write.table(occ.teste,nome.teste,sep=",",row.names=F)
  occ.treino<-NULL
  occ.teste<-NULL
}


