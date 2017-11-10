Split_Xls <- function(dir_output,teste){

  #Parametros
    #dir_output: diretorio que serao salvos os arquivos xls por especie
  
library(xlsx)
  print("Select the training dataset .xls file:")
  input.tr<-read.xlsx(file.choose(),sheetIndex=1)

  if(teste==1){
    print("Select the testing dataset .xls file:")
    input.ts<-read.xlsx(file.choose(),sheetIndex=1)
  }

  ssp<-unique(input.tr[,1])

  setwd(dir_output)

  for (x in 1:length(ssp)){
    sp.treino<-input.tr[input.tr[,1]==ssp[x],]
    name.tr<-paste(ssp[x],"_par.xls",sep="")
    write.xlsx(sp.treino,name.tr,row.names=F)
    if(teste==1){
      sp.teste<-input.ts[input.ts[,1]==ssp[x],]
      name.ts<-paste(ssp[x],"_impar.xls",sep="")
      write.xlsx(sp.teste,name.ts,row.names=F)
    }
  }
}

Split_Xls("C:/Users/Andre/Google Drive/Mestrado/Capitulo 1(Transferabilidade)/Ocorrencias_20/Aleatorio/Replicas/Split Ocorrencias/Rep10",1)
