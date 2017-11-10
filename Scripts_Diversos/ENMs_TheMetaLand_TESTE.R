ENMs_TheMetaLand<-function (diretorio,background,psa){
  
  #Essa ? a fun?ao padronizada para cria?ao de modelos de distribui?ao do laborat?rio
  #TheMetaLand
  #Par?metro iniciais:
    #diretorio:O diretorio com as ocorrencias(.csv) e as camadas ambientais
    #env: o formato das camadas ambientais (raster,asc ou txt)
    #PCA:Deseja realizar uma PCA nas var.ambientais e utilizar os eixos?(S/N)
    #background:n?mero de pontos utilizados para criar o background
    #psa:n?mero de pontos utilizados como pseudo-ausencias
    #k.fold: N?mero de grupos para dividir as ocorrencias(k-fold)
    #Algoritmos: os algoritmos utilizados para criacao dos modelos (GLM,SVM,Mahalanobis,Domains,MaxEnt,Random Forest)
 #OBSERVACOES IMPORTANTES:
  #O formato das ocorrencias deve ser (species,x(long),y(lat))
  
  #Carregar os pacotes
  library(raster)
  library(sp)
  library(dismo)
  library(rJava)#MAXENT
  library(kernlab)#SVM
  library(maps)
  library(SDMTools)#Ler ASCs
  library(adehabitatHS)#ENFA
  
  setwd(diretorio)
  
  rm(list=ls())
  
#1.VARI?VEIS AMBIENTAIS

  #if (env=='raster'){
  #  environmental<-stack(list.files(pattern=".bil"))
    
  #}
  
  #if(env=='asc'){
    env_prep<-list.files(pattern=".asc")
    asc_base<-read.asc(env_prep[1])
    environmental<-raster(asc_base)
    camadas<-env_prep[-1]
    
    for (a in 1:length(camadas)){
      asc<-read.asc(camadas[a])
      raster<-raster(asc)
      environmental<-stack(environmental,raster)
    }
    
    names(environmental)<-env_prep
  #}
  
  #if (env=='txt'){
   #env_prep<-read.table(list.files(pattern=".txt"),h=T)
  #  gridded(env_prep)<- ~x+y
   # environmental<-stack(env_prep)
  #}

  #1.1.Recortar as vari?veis ambientais para regi?o de interesse(Definir background)
        #e: xmin,xmax,ymin,ymax da regiao de interesse
          #Comumente utilizadas:
            #Neotr?pico: (-90,-30,-60,15)
            #Corais.AMS: (-81,-28,-43,14)
            #Corais.Caribe:(-100,-28,-43,33)

  #extent<- c(-90,-30,-60,15)
  #environmental<- crop(environmental,extent)

  
  #1.2.Realizar a PCA nas camadas ambientais(utilizar PCA como input)

  #if (PCA="S"){
   # data.frame<-rasterToPoints(environmental)
    #pca.raw<-data.frame[,-c(1:2)]
    #pca.raw<- na.omit(pca.raw)
    
    # Scale transform 
    #data.scaled <- data.frame(apply(pca.raw,2,scale))
    #nrow(data.scaled)

    # Realizar a PCA 
    #data.pca <- prcomp(data.scaled,retx=TRUE)
    #str(data.pca)

  
    #Porcentagem da varia??o explicada
    #str(summary(data.pca))
    #n.eixos<-length(summary(data.pca)$importance[3,])
    #cumulat.var<-summary(data.pca)$importance[3,]
    
    #Eixos que explicam 95% da varia??o
    #var.95<-cumulat.var<=0.95
    
    #Recuperar os loadings e transformar em um data.frame
    #eixos<-as.data.frame(data.pca$x)
    #eixos.95<-eixos[,var.95]==T
    #eixos.95.var<-eixos[,1:ncol(eixos.95)]
    #eixos.xy<-cbind(data.frame[,(1:2)],eixos.95.var)
    #gridded(eixos.xy)<- ~x+y
    #environmental<-stack(eixos.xy)
  #}

  #plot(environmental$PC1)
  #points(ocorrencias.xy)

  #1.3.Projecoes para outros per?odos de tempo ou outra regiao

##### P/INCLUIR(#1) #####

  #1.4.Definir Background
  background.env<-rasterToPoints(environmental)
  id.back<-sample(1:nrow(background.env),10000)
  background.env<-background.env[id.back,]
  

#2.ARQUIVOS DE OCORR?NCIA

  ocorrencias<-read.table(list.files(pattern=".txt"),h=T,sep="\t")
  ocorrencias.xy<-ocorrencias[,-1]

  #Cria um vetor com as esp?cies presentes na planilha
  especies<-unique(ocorrencias[,1])
  
#Remover presen?as duplicadas
  ocorrencias.var<-extract(environmental,ocorrencias.xy,cellnumber=T)
  ocorrencias.var<-cbind(ocorrencias,ocorrencias.var)
  ocorrencias.var<-na.omit(ocorrencias.var)

  ocorrencias.unicas<-NULL

  for (b in 1:length(especies)){
    ocorrencias.especie<-ocorrencias.var[ocorrencias.var[,1]==especies[b],]
    duplicadas<-which(duplicated(ocorrencias.especie[,"cells"])==T)
    ocorrencias.unicas<-rbind(ocorrencias.unicas,ocorrencias.especie[-duplicadas,])
  }
  rownames(ocorrencias.unicas)<-NULL
  
  #Caso queira trabalhar apenas com long e lat sem valores ambientais nas ocorrencias
  #ocorrencias<- cbind(cbind(ocorrencias,ocorrencias.var[,"cells"]))
  #colnames(ocorrencias)<-c("species","long","lat","cells")
  #ocorrencias<-na.omit(ocorrencias)
  #duplicadas<-which(duplicated(ocorrencias[,"cells"])==T)
  #ocorrencias.unicas<-ocorrencias.var[-duplicadas,]
  #rownames(ocorrencias.unicas)<-NULL

  #2.1.Separar as ocorrencias em grupos(k-fold)
  grupos<-NULL

  #Definir o n?mero de grupos que se deseja criar(k-fold)
  k<-2

  #Grupo de esp?cies que nao pode ser dividida(devido ao baixo N) e ser? exclu?da
  #ocorrencias.excluidas<-NULL
  
  #Separar em grupos
  for (c in 1:length(especies)){
    ocorrencias.unicas.especie<-ocorrencias.unicas[ocorrencias.unicas[,1]==especies[c],]
    
    if (nrow(ocorrencias.unicas.especie)<k){
      print(paste(especies[c],":N?mero de ocorrencias insuficiente para realizar k-fold",sep=" "))
      #ocorrencias.excluidas<-rbind(ocorrencias.unicas[-ocorrencias.unicas[,1]==especies[c],])
      #ocorrencias.excluidas<-na.omit(ocorrencias.excluidas)
      
    } else{
    occ.grupo<-kfold(ocorrencias.unicas.especie,k)
    grupos<-c(grupos,occ.grupo)
    
    }
  }
  nrow(ocorrencias.unicas)  
  length(grupos)
  #ocorrencias.unicas<-ocorrencias.unicas[-ocorrencias.excluidas,]
  ocorrencias.especie.grupos<-cbind(ocorrencias.unicas,grupos)


  #Criar listas onde cada elemento ser? utilizado para treino/teste nos modelos
  
  lista.treino<-NULL
  lista.teste<-NULL
  
  for (d in 1:k){
    occ.treino<-ocorrencias.especie.grupos[(ocorrencias.especie.grupos[,"grupos"]!=d),]
    occ.teste<-ocorrencias.especie.grupos[(ocorrencias.especie.grupos==d),]
    lista.treino[[d]]<-list(occ.treino)
    lista.teste[[d]]<-list(occ.teste)   
  }


  #2.2.Definir pseudo-ausencias
  pseudo.ausencias<-rasterToPoints(environmental)
  especies.unicas<-unique(ocorrencias.especie.grupos[,1])
  pseudo.ausencias.base<-NULL
  list.pseudo<-NULL

  for (e in 1:length(lista.treino)){
    ocorrencia.unicas.especie.psa<-lista.treino[[e]]

    for (f in 1:length(especies.unicas)){
      ocorrencias.unicas.especie.psa<-ocorrencias.especie.grupos[ocorrencias.especie.grupos[,1]==especies[f],]
      ocorrencias.unicas.especie.psa<-ocorrencias.unicas.especie.psa[ocorrencias.unicas.especie.psa[,"grupos"]!=1,]
      id.psa<-sample(1:nrow(pseudo.ausencias),nrow(ocorrencias.unicas.especie.psa))
      pseudo.ausencias.select<-pseudo.ausencias[id.psa,]
      pseudo.ausencias.select<-as.data.frame(pseudo.ausencias.select)
      pseudo.ausencias.sp<-cbind(ocorrencias.unicas.especie.psa[,1],pseudo.ausencias.select)
      pseudo.ausencias.base<-rbind(pseudo.ausencias.base,pseudo.ausencias.sp)
    }
  
    list.pseudo[[e]]<-pseudo.ausencias
  }  
  


ocorrencias.treino<-as.data.frame(lista.treino[[1]])
ocorrencias.treino<-ocorrencias.treino[,-11]
especies.treino<-unique(ocorrencias.treino[,1])

#return(list(ocorrencias.treino,pseudo.ausencias,background.env,environmental))

#BIOCLIM
Bioclim<-stack()

for (g in 1:length(especies.treino)){
  ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[g],]
  Bioclim.modelo<-bioclim(ocorrencias.treino.sp[,-c(1:4)])
  Bioclim.sp<-predict(Bioclim.modelo,environmental)
  Bioclim<-stack(Bioclim,Bioclim.sp)
}

names(Bioclim)<-especies.treino
#return(Bioclim)

#MAXENT

Maxent<-stack()
Response<-NULL
Sys.setenv(NOAWT=TRUE)

for (h in 1:length(especies.treino)){
  ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[1],]
  Maxent.modelo<-maxent(x=environmental,p=ocorrencias.treino.sp[,c(2:3)])
  Maxent.sp<-predict(Maxent.modelo,environmental)
  Maxent<-stack(Maxent,Maxent.sp)
}

names(Maxent)<-especies.treino
return(Maxent)


#pb.1<-rep(1,nrow(ocorrencias.treino))
#ocorrencias.treino.max<-cbind(ocorrencias.treino,pb.1)
#pb.0<-rep(0,nrow(background.env))
#background.env.max<-cbind(background.env,pb.0)

} # fecha a funcao

#PARA INCLUIR NA FUNCAO:
  #1.Trabalhar com per?odos de tempo diferente

output2<-ENMs_TheMetaLand("C:/Users/Andre/Desktop/Andre/Teste_codigo",background=10000,psa=10000)
output

#Para corrigir erros
#ocorrencias.treino<-as.data.frame(output[[1]])
#environmental<-output[[4]]
