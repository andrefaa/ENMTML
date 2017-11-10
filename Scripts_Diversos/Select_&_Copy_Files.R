setwd ("C:/Users/Andre/Google Drive/Mestrado/Capitulo 1(Transferabilidade)/Teste_Caribe/Resultados Maxent_Caribe/2Q/Agaricia fragilis/Impar-Caribe")


filelist<-list.files(pattern="Projet.asc")
file.copy(filelist, "C:/Users/Andre/Google Drive/Mestrado/Capitulo 1(Transferabilidade)/Teste_Caribe/Validacao/2Q/Agaricia fragilis/Par-Caribe")

##########################PARA TRABALHAR EM VÁRIOS DIRETÓRIOS############################

#Definir o diretório onde estao os outputs do Maxent
setwd ("C:/Users/Andre/Google Drive/Mestrado/Capitulo 1(Transferabilidade)/Teste_Caribe/Resultados Maxent_Caribe/10Q/")

#Listar as pastas(uma para cada espécie trabalhada, por exemplo)
filelist<-list.files()

Automatic.Copy.to.Evaluate.Maxent.Models<-function(lista_de_especies){
  for (x in 1:length(filelist)){
    especie<-list.files(filelist[x])
    for (y in 2:length(especie)){
      projecoes<-list.files(paste(getwd(),"/",sep="",filelist[x],"/",especie[y]),pattern="Projet.asc")
      input<-paste("C:/Users/Andre/Google Drive/Mestrado/Capitulo 1(Transferabilidade)/Teste_Caribe/Resultados Maxent_Caribe/10Q","/",filelist[x],"/",especie[y],sep="")
      output<- paste("C:/Users/Andre/Google Drive/Mestrado/Capitulo 1(Transferabilidade)/Teste_Caribe/Validacao/10Q/Par")
      input_projecoes<-paste(input,"/",projecoes,sep="")
      #for (z in 1:length(input_projecoes)){
      file.copy(input_projecoes,output,overwrite=FALSE)
      #}
    }
  }
}

Automatic.Copy.to.Evaluate.Maxent.Models(filelist)

