ENMs_TheMetaLand<-function(diretorio,env='',pca='',background,particao,algoritmo=''){
  
  #Essa ? a fun?ao padronizada para cria?ao de modelos de distribui?ao do laborat?rio
  #TheMetaLand
  #Parametro iniciais:
    #diretorio: Diretorio com as camadas ambientais
    #env: Formato das camadas ambientais (raster/asc/txt)
    #PCA: Deseja realizar uma PCA nas var.ambientais e utilizar os eixos?(S/N)
    #background: Numero de pontos utilizados para criar o background
    #particao: Tipo de particao de dados (k.fold/bootstrap)
    #Algoritmos: os algoritmos utilizados para criacao dos modelos (GLM,SVM,MaxEnt,Bioclim)
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
    env_prep<-list.files(pattern='.bil')
    environmental<-stack(list.files(pattern='.bil'))
    names(environmental)<-substr(env_prep,1,nchar(env_prep)-4)
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
    names(environmental)<-substr(env_prep,1,nchar(env_prep)-4)
  }

  if(env == 'txt'){
    environmental<-read.table(list.files(pattern='.txt'),h=T)
    env_prep<-colnames(environmental)
    gridded(environmental)<- ~x+y
    environmental<-stack(environmental)
    names(environmental)<-substr(env_prep,1,nchar(env_prep)-4)
  }

  #1.1.Definir Background
  background.env<-rasterToPoints(environmental)
  id.back<-sample(1:nrow(background.env),background)
  background.env<-background.env[id.back,]


  #1.2.Realizar a PCA nas camadas ambientais(utilizar PCA como input)
  
    if (pca=="S"){
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

  print("Select the occurrences file(.txt):")
  ocorrencias.r<-read.table(file.choose(),h=T,sep='\t')
  ocorrencias<-ocorrencias.r[,which(names(ocorrencias.r) %in% c("Species","Long","Lat"))]
  ocorrencias.xy<-ocorrencias[,-1]

  #Cria um vetor com as esp?cies presentes na planilha
  especies<-unique(ocorrencias[,1])
  
  #Remover presencas duplicadas
  ocorrencias.var<-extract(environmental,ocorrencias.xy,cellnumber=T)
  ocorrencias.var<-cbind(ocorrencias,ocorrencias.var)
  ocorrencias.var<-na.omit(ocorrencias.var)

  ocorrencias.unicas<-NULL
  registros.ssp<-NULL
  registros.unicos<-NULL

  for (b in 1:length(especies)){
    
    ocorrencias.especie<-ocorrencias.var[ocorrencias.var[,1]==especies[b],]
    duplicadas<-which(duplicated(ocorrencias.especie[,'cells'])==T)
    registros.ssp<-nrow(ocorrencias.especie[-duplicadas])
    registros.unicos<-rbind(registros.unicos,registros.ssp)
    
    if (length(duplicadas)!=0){
      ocorrencias.unicas<-rbind(ocorrencias.unicas,ocorrencias.especie[-duplicadas,])
    }else{
      ocorrencias.unicas<-rbind(ocorrencias.unicas,ocorrencias.especie)
    }    
  }

  rownames(ocorrencias.unicas)<-NULL
  registros.unicos<-cbind(as.character(especies),registros.unicos)
  colnames(registros.unicos)<-c("Especie","Registros_Unicos")
  write.table(registros.unicos,"Registros_Unicos.txt",row.names=F,sep="\t")
 

  #2.1.Separar as ocorrencias em grupos(k-fold)
    
  if (particao=="k.fold"){
    
    cat("Select the number of groups to be created (>=1): ")
    k.fold <- as.integer(readLines(n = 1))
    grupos<-NULL
  
    #Grupo de especies que nao pode ser dividida(devido ao baixo N) e ser? exclu?da
    ocorrencias.excluidas<-NULL
    
    #Separar em grupos
    for (c in 1:length(especies)){
      ocorrencias.unicas.especie<-ocorrencias.unicas[ocorrencias.unicas[,1]==especies[c],]
      
      if (nrow(ocorrencias.unicas.especie)<k.fold){
        print(paste(especies[c],':Numero de ocorrencias insuficiente para realizar k-fold',sep=' '))
        ocorrencias.excluidas<-rbind(ocorrencias.excluidas,ocorrencias.unicas.especie)        
      } else{
      occ.grupo<-kfold(ocorrencias.unicas.especie,k.fold)
      grupos<-c(grupos,occ.grupo)
      
      }
    }
    
    ocorrencias.unicas<-ocorrencias.unicas[-c(as.integer(row.names(ocorrencias.excluidas))),]
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
  }

#3. Criar os modelos:

if (particao =="bootstrap"){
  cat("Select the number of replicates (>=1): ")
  replicas <- as.integer(readLines(n = 1))
  ocorrencias.treino<-NULL
  ocorrencias.teste<-NULL
  for (y in 1:replicas){
    for (x in 1:length(especies)){
      ocorrencias.unicas.especie<-ocorrencias.unicas[ocorrencias.unicas[,1]==especies[x],]
      id.ocor<-sample(1:nrow(ocorrencias.unicas.especie),round(1*nrow(ocorrencias.unicas.especie)))
      occ.treino.sp<-ocorrencias.unicas.especie[id.ocor,]
      #occ.teste.sp<-ocorrencias.unicas.especie[-id.ocor,]
      ocorrencias.treino<-rbind(ocorrencias.treino,occ.treino.sp)
      #ocorrencias.teste<-rbind(ocorrencias.teste,occ.teste.sp)
    }
    write.xlsx(ocorrencias.treino[,c(1:3)],paste("Ocorrencias_treino",y,".xls",sep=""),row.names=F)
    #write.xlsx(ocorrencias.teste[,c(1:3)],paste("Ocorrencias_teste",y,".xls",sep=""),row.names=F)
    
    especies.treino<-unique(ocorrencias.treino[,1])
      
      #BIOCLIM
      if (any(algoritmo=='bioclim')){
        
        for (g in 1:length(especies.treino)){
          ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[g],]
          Bioclim.modelo<-bioclim(ocorrencias.treino.sp[,-c(1:4)])
          Bioclim.sp<-predict(Bioclim.modelo,environmental)
          writeRaster(Bioclim.sp,paste(especies.treino[g],y,'bioclim.asc',sep='_'),format='ascii')
        }
      }
      
      #MAXENT
      if (any(algoritmo=='maxent')){
        Sys.setenv(NOAWT=TRUE)
        
        for (h in 1:length(especies.treino)){
          ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[h],]
          Maxent.modelo<-maxent(x=environmental,p=ocorrencias.treino.sp[,c(2:3)],args=c(
                                "randomseed","noproduct","nothreshold","nohinge","noautofeature","nodoclamp","maximumiterations=1000"))
          Maxent.sp<-predict(Maxent.modelo,environmental)
          writeRaster(Maxent.sp,paste(especies.treino[h],y,'maxent.asc',sep='_'),format='ascii')
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
          input.occ.back<-input.occ.back[,-c(1:3)]

          #Rodar o modelo SVM
          svm.modelo<-ksvm(pb.1~.,data=input.occ.back)
          
          SVM.sp<-predict(environmental,svm.modelo)
          writeRaster(SVM.sp,paste(especies.treino[i],y,'svm.asc',sep='_'),format='ascii')
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
          input.occ.back<-input.occ.back[,-c(1:3)]
        
          #Rodar o modelo GLM
          glm.modelo<-glm(pb.1~.,data=input.occ.back,family=binomial(link='logit'))
          
          GLM.sp<-predict(environmental,glm.modelo)
          writeRaster(GLM.sp,paste(especies.treino[j],y,'glm.asc',sep='_'),format='ascii')
        }
      }
      
      #ENFA
      #if (any(algoritmo=="enfa")){
      #for (j in 1:length(especies.treino)){
      #ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[j],]
      #}
      #}
  }
}
    
if (particao=="k.fold"){
  for (l in 1:length(lista.treino)){
    
    #Escolhe cada elemento da lista para realizar o modelo
    ocorrencias.treino<-as.data.frame(lista.treino[[l]])
    ocorrencias.treino<-ocorrencias.treino[,-which(names(ocorrencias.treino) %in% c("grupos"))]
    write.xlsx(ocorrencias.treino[,c(1:3)],paste("Ocorrencias_treino",l,".xls",sep=""),row.names=F)
    
    if (k.fold != 1){
      ocorrencias.teste<-as.data.frame(lista.teste[[l]])
      ocorrencias.teste<-ocorrencias.teste[,-which(names(ocorrencias.teste) %in% c("grupos"))]
      write.xlsx(ocorrencias.teste[,c(1:3)],paste("Ocorrencias_teste",l,".xls",sep=""),row.names=F)  
    }

    especies.treino<-unique(ocorrencias.treino[,1])
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
      Maxent.modelo<-maxent(x=environmental,p=ocorrencias.treino.sp[,c(2:3)],args=c(
                            "randomseed","nowarnings","notooltips","noaskoverwrite","nowriteclampgrid","nowritemess",
                            "noproduct","nothreshold","nohinge","writeplotdata","noautofeature","nodoclamp","maximumiterations=1000"))
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
      input.occ.back<-input.occ.back[,-c(1:3)]
      
      #Rodar o modelo SVM
      svm.modelo<-ksvm(pb.1~.,data=input.occ.back)
      
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
      input.occ.back<-input.occ.back[,-c(1:3)]

      #Rodar o modelo GLM
      glm.modelo<-glm(pb.1~.,data=input.occ.back,family=binomial(link='logit'))
      
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
}
} # fecha a funcao

#PARA INCLUIR NA FUNCAO:
  #1.Trabalhar com per?odos de tempo diferente
  #2.A funcao response("resultado do modelo") mostra as curvas de resposa para cada variável

ENMs_TheMetaLand('C:/Users/Andre/Google Drive/Mestrado/PeixesOsseos/GLM_49','asc',pca='N',10000,"bootstrap",algoritmo='glm')

rm(list=ls())
