
setwd ("C:/Users/Andre/Desktop/Andre/Mestrado/Ocorrencias")

######################## A PARTIR DE TABELA CSV ###################################

tabela<-read.table("C:/Users/Andre/Desktop/Andre/Mestrado/Ocorrencias/Myrmecophaga_tridactyla.csv", sep=",", header=TRUE)

dataframe<- data.frame(tabela$scientificname,tabela$latitude, tabela$longitude)

Myrmecophaga_tridactyla<-na.exclude(dataframe)

colnames(Myrmecophaga_tridactyla)<- c("species","lat","lon")

write.table(Myrmecophaga_tridactyla, "C:/Users/Andre/Desktop/Andre/Mestrado/Ocorrencias/Ocorrencias CSV/Myrmecophaga_tridactyla.csv", row.names=F)


######################## A PARTIR DE TABELA XLSX ###################################

library(XLConnect)

theData<- readWorksheet(loadWorkbook("C:/Users/Andre/Desktop/Andre/Mestrado/Ocorrencias/Adenostemma_brasilianum.xlsx"),sheet=1)

write.table(theData, "provisorio.csv",sep=",",row.names=F)

tabela2<-read.table("C:/Users/Andre/Desktop/Andre/Mestrado/Ocorrencias/provisorio.csv", sep=",", header=TRUE)

names(tabela2)

dataframe2<- data.frame(tabela2$scientificname,tabela2$latitude, tabela2$longitude)

Adenostemma_brasilianum<- na.exclude(dataframe2)

colnames(Adenostemma_brasilianum)<- c("species","lat","lon")

write.table(Adenostemma_brasilianum, "C:/Users/Andre/Desktop/Andre/Mestrado/Ocorrencias/Ocorrencias CSV/Adenostemma_brasilianum.csv", row.names=F)
