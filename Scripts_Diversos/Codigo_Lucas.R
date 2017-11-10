#Funcao para selecionar apenas os registros unicos observados, baseado em um limiar pre-definido:
#Parametros
  #limit:tempo minimo entre registros (EM MINUTOS!)

Registros.Unicos <- function(limit) {

tab<-read.table(file.choose(),header = 1,sep="\t")
tab$HORA<-as.character(tab$HORA)
tab$HORA<-as.POSIXct(tab$HORA, format="%H:%M:%S")
tab<-tab[-c(which(is.na(tab$HORA))),]
ssp<-unique(tab[,8])
tab.dia.un<-NULL
count<-NULL

for (a in ssp){
  tab.sp<-tab[tab[,8]==a,]
  dia<-unique(tab.sp[,14])
  
  for (b in dia){
    tab.dia<-tab.sp[tab.sp[,14]==b,]
    tab.dia<-tab.dia[order(tab.dia$HORA),]
    difer<-diff(tab.dia$HORA)
    keep <- c(T, rep(F, length(difer)))
    acc <- 0
    
    if (nrow(tab.dia)==0){
      keep<-T
    }else{
      for (i in seq_along(difer)) {
        acc <- acc + difer[i]
        if (acc > limit) {
          keep[i+1] <- T
          acc <- 0
        }
      }
    }
    tab.dia.un<-rbind(tab.dia.un,tab.dia[keep,])
  }
  count<-c(count,nrow(tab.dia.un[tab.dia.un[,8]==a,]))
}
print("Selecione o Diretorio para salvar os Resultados:")
dir.out<-choose.dir()
tab.count<-cbind(ssp,as.data.frame(count))
write.table(tab.dia.un,paste(dir.out,"Registros_Unicos.txt",sep="\\"),sep="\t",row.names=F)
write.table(tab.count,paste(dir.out,"Numero_Registros.txt",sep="\\"),sep="\t",row.names=F)
}

