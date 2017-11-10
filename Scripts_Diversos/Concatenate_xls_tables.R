#Definir o diretório onde estao os outputs do Maxent
setwd ("C:\\OneDrive\\Doutorado\\Bloco1 - ENM\\Cap2 - Virtual Species ENM\\Modelos\\Narrow\\Distribuicao\\2EIXOS\\NaoEquilibrio\\Occ_50")
library(xlsx)
library(plyr)
#Listar as pastas(uma para cada espécie trabalhada, por exemplo)
filelist<-list.files(pattern=".xls")
nomes <- strsplit(filelist, "_")
nomes <- rapply(nomes,function(x) x[2])

#Criar tabela base
tabela<-data.frame("Species"=character(1),"Long"=numeric(1),"Lat"=numeric(1),stringsAsFactors=FALSE)

Automatic.Concatenate.Xls.Tables<-function(lista_de_pastas){
  for (x in 1:length(filelist)){
    tabela_especie<-data.frame(read.xlsx(filelist[x],sheetIndex="Planilha1"),stringsAsFactors=FALSE)
    tabela_especie[,1] <- nomes[x]
    colnames(tabela_especie) <- colnames(tabela)
    tabela<-rbind.data.frame(tabela,tabela_especie)
  }
  tabela <- tabela[-1,]
  write.table(tabela,paste("Ocorrencias",as.character(nrow(tabela)/100),".txt",sep=""),sep="\t",row.names=F)
}

Automatic.Concatenate.Xls.Tables(getwd())
