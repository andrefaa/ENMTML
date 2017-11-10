ENMs_TheMetaLand<-function (diretorio,env='',pca='',background,k.fold,algoritmo=''){
  
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
  library(xlsx)#Criar arquivos xls com ocorrencias treino/teste
  
  setwd(diretorio)
  
#1.VARIAVEIS AMBIENTAIS
  
  if(env == 'raster'){
    environmental<-stack(list.files(pattern='.bil'))
    
  }
  
  if(env == 'asc'){
    env_prep<-list.files(pattern='.asc')
    asc_base<-read.asc(env_prep[1])
    environmental<-raster(asc_base)
    camadas<-env_prep[-1]
    
    for (a in 1:length(camadas)){
      asc<-read.asc(camadas[a])
      raster<-raster(asc)
      environmental<-stack(environmental,raster)
    }
    names(environmental)<-env_prep
  }
  
  if(env == 'txt'){
   env_prep<-read.table(list.files(pattern='.txt'),h=T)
    gridded(env_prep)<- ~x+y
    environmental<-stack(env_prep)
  }

  #1.1.Definir Background
  background.env<-rasterToPoints(environmental)
  id.back<-sample(1:nrow(background.env),background)
  background.env<-background.env[id.back,]

  #1.2.Recortar as vari?veis ambientais para regi?o de interesse(Definir background)
    #e: xmin,xmax,ymin,ymax da regiao de interesse
    #Comumente utilizadas:
    #Neotr?pico: (-90,-30,-60,15)
    #Corais.AMS: (-81,-28,-43,14)
    #Corais.Caribe:(-100,-28,-43,33)
    
    #extent<- c(-90,-30,-60,15)
    #environmental<- crop(environmental,extent)
  
  
  1.3.Realizar a PCA nas camadas ambientais(utilizar PCA como input)
  
    if (PCA="S"){
     data.frame<-rasterToPoints(environmental)
    pca.raw<-data.frame[,-c(1:2)]
    pca.raw<- na.omit(pca.raw)
    
    #Scale transform 
    data.scaled <- data.frame(apply(pca.raw,2,scale))
    nrow(data.scaled)
    
    # Realizar a PCA 
    data.pca <- prcomp(data.scaled,retx=TRUE)
    str(data.pca)
        
    #Porcentagem da variacaoo explicada
    str(summary(data.pca))
    n.eixos<-length(summary(data.pca)$importance[3,])
    cumulat.var<-summary(data.pca)$importance[3,]
    
    #Eixos que explicam 95% da varia??o
    var.95<-cumulat.var<=0.95
    
    #Recuperar os loadings e transformar em um data.frame
    eixos<-as.data.frame(data.pca$x)
    eixos.95<-eixos[,var.95]==T
    eixos.95.var<-eixos[,1:ncol(eixos.95)]
    eixos.xy<-cbind(data.frame[,(1:2)],eixos.95.var)
    gridded(eixos.xy)<- ~x+y
    environmental<-stack(eixos.xy)
    }

#2.ARQUIVOS DE OCORRENCIA

  ocorrencias<-read.table(list.files(pattern='.txt'),h=T,sep='\t')
  ocorrencias.xy<-ocorrencias[,-1]

  #Cria um vetor com as esp?cies presentes na planilha
  especies<-unique(ocorrencias[,1])
  
  #Remover presencas duplicadas
  ocorrencias.var<-extract(environmental,ocorrencias.xy,cellnumber=T)
  ocorrencias.var<-cbind(ocorrencias,ocorrencias.var)
  ocorrencias.var<-na.omit(ocorrencias.var)

  ocorrencias.unicas<-NULL

  for (b in 1:length(especies)){
    ocorrencias.especie<-ocorrencias.var[ocorrencias.var[,1]==especies[b],]
    duplicadas<-which(duplicated(ocorrencias.especie[,'cells'])==T)
    ocorrencias.unicas<-rbind(ocorrencias.unicas,ocorrencias.especie[-duplicadas,])
  }
  rownames(ocorrencias.unicas)<-NULL
  
  #2.1.Separar as ocorrencias em grupos(k-fold)
  grupos<-NULL

  #Grupo de especies que nao pode ser dividida(devido ao baixo N) e ser? exclu?da
  ocorrencias.excluidas<-NULL
  
  #Separar em grupos
  for (c in 1:length(especies)){
    ocorrencias.unicas.especie<-ocorrencias.unicas[ocorrencias.unicas[,1]==especies[c],]
    
    if (nrow(ocorrencias.unicas.especie)<k.fold){
      print(paste(especies[c],':N?mero de ocorrencias insuficiente para realizar k-fold',sep=' '))
      ocorrencias.excluidas<-rbind(ocorrencias.unicas[-ocorrencias.unicas[,1]==especies[c],])
      ocorrencias.excluidas<-na.omit(ocorrencias.excluidas)
      
    } else{
    occ.grupo<-kfold(ocorrencias.unicas.especie,k.fold)
    grupos<-c(grupos,occ.grupo)
    
    }
  }
  nrow(ocorrencias.unicas)  
  length(grupos)
  #ocorrencias.unicas<-ocorrencias.unicas[-ocorrencias.excluidas,]
  ocorrencias.especie.grupos<-cbind(ocorrencias.unicas,grupos)


  #Criar listas onde cada elemento sera utilizado para treino/teste nos modelos
  
  lista.treino<-NULL
  lista.teste<-NULL
  
  for (d in 1:k.fold){
    if(k.fold == 1){
      occ.treino<-ocorrencias.especie.grupos[(ocorrencias.especie.grupos[,'grupos']==d),]
      lista.treino[[d]]<-list(occ.treino)
    }else{
    occ.treino<-ocorrencias.especie.grupos[(ocorrencias.especie.grupos[,'grupos']!=d),]
    occ.teste<-ocorrencias.especie.grupos[(ocorrencias.especie.grupos[,'grupos']==d),]
    lista.treino[[d]]<-list(occ.treino)
    lista.teste[[d]]<-list(occ.teste)   
    }
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

#3. Criar os modelos:

for (l in 1:length(lista.treino)){
  
  #Escolhe cada elemento da lista para realizar o modelo
  ocorrencias.treino<-as.data.frame(lista.treino[[l]])
  ocorrencias.treino<-ocorrencias.treino[,-11]
  especies.treino<-unique(ocorrencias.treino[,1])
  write.xlsx(ocorrencias.treino[,c(1:3)],paste("Ocorrencias_treino",l,".xls",sep=""),row.names=F)
  
  if (k.fold != 1){
    ocorrencias.teste<-as.data.frame(lista.teste[[l]])
    ocorrencias.teste<-ocorrencias.teste[,-11]
    write.xlsx(ocorrencias.teste[,c(1:3)],paste("Ocorrencias_teste",l,".xls",sep=""),row.names=F)  
  }
   
  #BIOCLIM
  if (any(algoritmo=='bioclim')){
  
    for (g in 1:length(especies.treino)){
      ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[g],]
      Bioclim.modelo<-bioclim(ocorrencias.treino.sp[,-c(1:4)])
      Bioclim.sp<-predict(Bioclim.modelo,environmental)
      writeRaster(Bioclim.sp,paste(especies.treino[g],l,'bioclim.asc',sep='_'),format='ascii')
    }
  }
  
  #MAXENT
  if (any(algoritmo=='maxent')){
  Sys.setenv(NOAWT=TRUE)
  
    for (h in 1:length(especies.treino)){
      ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[h],]
      Maxent.modelo<-maxent(x=environmental,p=ocorrencias.treino.sp[,c(2:3)])
      Maxent.sp<-predict(Maxent.modelo,environmental)
      writeRaster(Maxent.sp,paste(especies.treino[h],l,'maxent.asc',sep='_'),format='ascii')
    }
  }
  
  #SVM
  if (any(algoritmo=='svm')){
  
    for (i in 1:length(especies.treino)){
      ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[i],]
        
      #Agrupar background e presencas
      pb.1<-rep(1,nrow(ocorrencias.treino.sp))
      ocorrencias.treino.max<-cbind(ocorrencias.treino.sp,pb.1)
      ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
      species.pb<-rep(especies.treino[i],nrow(background.env))
      pb.0<-rep(0,nrow(background.env))
      background.env.max<-as.data.frame(cbind(background.env,pb.0))
      background.env.max<-cbind(species.pb,background.env.max)
      colnames(background.env.max)<-colnames(ocorrencias.treino.max)
      input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
      
      #Rodar o modelo SVM
      svm.modelo<-ksvm(pb.1 ~ pca1.asc+pca2.asc+pca3.asc+pca4.asc+pca5.asc+pca6.asc,data=input.occ.back)
      
      SVM.sp<-predict(environmental,svm.modelo)
      writeRaster(SVM.sp,paste(especies.treino[i],l,'svm.asc',sep='_'),format='ascii')
    }
  }
  
  #GLM
  if (any(algoritmo=='glm')){
  
    for (j in 1:length(especies.treino)){
      ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[j],]
      
      #Agrupar background e presencas
      pb.1<-rep(1,nrow(ocorrencias.treino.sp))
      ocorrencias.treino.max<-cbind(ocorrencias.treino.sp,pb.1)
      ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
      species.pb<-rep(especies.treino[j],nrow(background.env))
      pb.0<-rep(0,nrow(background.env))
      background.env.max<-as.data.frame(cbind(background.env,pb.0))
      background.env.max<-cbind(species.pb,background.env.max)
      colnames(background.env.max)<-colnames(ocorrencias.treino.max)
      input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
      return(input.occ.back)
      
      #Rodar o modelo GLM
      glm.modelo<-glm(pb.1 ~ pca1.asc+pca2.asc+pca3.asc+pca4.asc+pca5.asc+pca6.asc,data=input.occ.back,family=binomial(link='logit'))
      
      GLM.sp<-predict(environmental,glm.modelo)
      writeRaster(GLM.sp,paste(especies.treino[j],l,'glm.asc',sep='_'),format='ascii')
    }
  }
  
  #ENFA
  #if (any(algoritmo=="enfa")){
    #for (j in 1:length(especies.treino)){
      #ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[j],]
    #}
  #}
}

} # fecha a funcao

#PARA INCLUIR NA FUNCAO:
  #1.Trabalhar com per?odos de tempo diferente

ENMs_TheMetaLand('C:/Users/Andre/Desktop/Andre/Teste_codigo','asc',pca='N',10000,1,algoritmo='glm')

rm(list=ls())
