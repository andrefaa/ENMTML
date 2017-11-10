ENMs_TheMetaLand<-function(diretorio,env='',pca='',project='',background,particao,algoritmo=''){
  
  #Funcao padronizada para criar modelos de distribuicao do laboratorio TheMetaLand
  
  #Parametro iniciais:
    #diretorio: Diretorio com as camadas ambientais
    #env: Formato das camadas ambientais (raster/asc/txt)
    #PCA: Deseja realizar uma PCA nas var.ambientais e utilizar os eixos?(S/N)
    #project : Deseja projetar o modelo em outra area? (S/N)
    #background: Numero de pontos utilizados para criar o background
    #particao: Tipo de particao de dados 
      #k.fold : Separa as ocorrencias em grupos aleatorios utilizando k-1 para criar os modelos e o outro para testar
      #bootstrap : Particao aleatoria dos dados baseada em uma porcentagem
      #leave1out : Aconselhado para poucos dados(<5);cria os modelos com n-1 ocorrencia e testa na ocorrencia restante
    #Algoritmos: lista de algoritmos atualmente suportados
      #bioclim : Bioclim
      #maxentS : Maxent Simples
      #maxentD : Maxent Default
      #svm : Support Vector Machine(SVM)
      #glm : Generalized Linear Model (GLM)
      #RF : Random Forest
      #MDA : Mixture Discriminant Analysis (MDA)
      #mahal : Mahalanobis Distance
  
 ## OBSERVACOES IMPORTANTES!!!!!!!!!!!!!!!:
  ######## O formato das ocorrencias deve ser (Species,Long,Lat) ########
  ######## Espécies com número de ocorrencias <5 serão desconsideradas no momento de criar os modelos ########

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
  library(randomForest)#Random Forest
  library(mda)#MDA
  #library(AICcmodavg)#Calculo de AIC
  
  
  setwd(diretorio)
  
#1.VARIAVEIS AMBIENTAIS----
  
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

  #1.1.Realizar a PCA nas camadas ambientais(utilizar PCA como input)
  
    if (pca=="S"){
      data.frame<-rasterToPoints(environmental)
      data.frame<-na.omit(data.frame)
      pca.raw<-data.frame[,-c(1:2)]
      
      means<-colMeans(pca.raw)
      stds<-apply(pca.raw,2,sd)

      #Scale transform 
      data.scaled <- data.frame(apply(pca.raw,2,scale))
      
      # Realizar a PCA 
      data.pca <- prcomp(data.scaled,retx=TRUE)
      
      #Salvar os coeficientes
      coeficientes<-data.pca$rotation
      write.table(coeficientes,'result_coef.txt',row.names=T,sep='\t')
                
      #Porcentagem da variacaoo explicada
      n.eixos<-length(summary(data.pca)$importance[3,])
      cumulat.var<-summary(data.pca)$importance[3,]
      variacao.explicada<-data.frame(cumulat.var)
      write.table(variacao.explicada,'result_latent.txt',row.names=T,sep='\t')      
      
      #Eixos que explicam 95% da varia??o
      var.95<-cumulat.var<=0.96
      
      #Recuperar os loadings e transformar em um data.frame
      eixos<-as.data.frame(data.pca$x)
      eixos.95<-eixos[,var.95]==T
      eixos.95.var<-eixos[,1:ncol(eixos.95)]
      eixos.xy<-cbind(data.frame[,(1:2)],eixos.95.var)
      gridded(eixos.xy)<- ~x+y
      environmental<-stack(eixos.xy)
    }

    #1.2.Definir Background
    background.env<-rasterToPoints(environmental)
    id.back<-sample(1:nrow(background.env),background)
    background.env.sample<-background.env[id.back,]
    
    #1.3.Área a ser projetada
    if (project=='S'){
      print("Select the folder with ASC variables for projection:")
      dir.project<-choose.dir()
      setwd(dir.project)
      proj_prep<-list.files(pattern='.asc')
      asc_proj<-read.asc(proj_prep[1])
      project<-raster(asc_proj)
      camadas.proj<-proj_prep[-1]
      
      for (a in 1:length(camadas.proj)){
        asc.p<-read.asc(camadas.proj[a])
        raster.p<-raster(asc.p)
        project<-stack(project,raster.p)
      }
      names(project)<-substr(proj_prep,1,nchar(proj_prep)-4)
      project.env<-rasterToPoints(project)
      project.env.raw<-project.env[,-c(1:2)]
      
      if (pca =="S"){
      project<-stack()
      nomes.PC<-NULL
      scale<-sweep(project.env.raw,2,means)
      scale<-scale %*% diag(1/stds)
      
        for (x in 1:ncol(project.env.raw)){
          coeficientes.pc<-coeficientes[,x]
          names(coeficientes.pc)<-NULL
          PC<-project.env.raw %*% diag(coeficientes.pc)
          PC<-rowSums(PC)
          PC<-cbind(project.env[,1:2],PC)
          colnames(PC)<-c("x","y","PC")
          PC<-as.data.frame(PC)
          gridded(PC)<- ~x+y
          PC<-raster(PC)
          project<-stack(project,PC)
          nomes.PC[x]<-paste("PC",x,sep="")
        }
      names(project)<-nomes.PC
      project<-project[[1:nlayers(environmental)]]
      project.env<-rasterToPoints(project)   
      }
    }else{
      project<-environmental
      project.env<-background.env
    }


#2.ARQUIVOS DE OCORRENCIA----
  
  print("Select a folder to save occurrence data for model evaluation:")
  dir2<-choose.dir()
  setwd(dir2)

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
 

  #2.1.Separar as ocorrencias em grupos(k-fold)----
    
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
  
    #2.2.Definir pseudo-ausencias----
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

#3. Criar os modelos:----

  #3.1.Selecionar as pastas onde salvar os outputs

  if (any(algoritmo=='bioclim')){
    print("Select the folder to save Bioclim results")
    dir.b<-choose.dir()
  }
  
  if (any(algoritmo=='mahal')){
    print("Select the folder to save Bioclim results")
    dir.mahal<-choose.dir()
  }
  
  if (any(algoritmo=='maxentS')){
    print("Select the folder to save Maxent Simplified results")
    dir.maxS<-choose.dir()
  }

  if (any(algoritmo=='maxentD')){
    print("Select the folder to save Maxent Default results")
    dir.maxD<-choose.dir()
  }

  if (any(algoritmo=='svm')){
    print("Select the folder to save SVM results")
    dir.s<-choose.dir()
  }

  if (any(algoritmo=='glm')){
    print("Select the folder to save GLM results")
    dir.g<-choose.dir()
  }

  if (any(algoritmo=='RF')){
    print("Select the folder to save Random Forest results")
    dir.rf<-choose.dir()
  }

  if (any(algoritmo=='MDA')){
    print("Select the folder to save MDA results")
    dir.mda<-choose.dir()
  }

#3.2.Particao----

if (particao =="bootstrap"){
  cat("Select the number of replicates (>=1): ")
  replicas <- as.integer(readLines(n = 1))
  cat("Select the percentage of occurrences used for training (0-1): ")
  porcentagem<-readLines(n = 1)
  porcentagem<-as.numeric(porcentagem)
  
  for (y in 1:replicas){
    
    ocorrencias.treino<-NULL
    ocorrencias.teste<-NULL
    
    if (porcentagem==1){
      for (x in 1:length(especies)){
        ocorrencias.unicas.especie<-ocorrencias.unicas[ocorrencias.unicas[,1]==especies[x],]
        
        #Retirar espécies com numero de ocorrencias < x
        if (nrow(ocorrencias.unicas.especie)<5){
          print(paste("Ocorrencias de ",ocorrencias.unicas.especie[,1]," insuficientes para criar um modelo",sep=""))
        }else{
        id.ocor<-sample(1:nrow(ocorrencias.unicas.especie),round(porcentagem*nrow(ocorrencias.unicas.especie)))
        occ.treino.sp<-ocorrencias.unicas.especie[id.ocor,]
        ocorrencias.treino<-rbind(ocorrencias.treino,occ.treino.sp)
        }
      }
      setwd(dir2)
      ocorrencias.treino[,1]<-paste(ocorrencias.treino[,1],y,sep="")
      write.xlsx(ocorrencias.treino[,c(1:3)],paste("Ocorrencias_treino",y,".xls",sep=""),row.names=F)  
      especies.treino<-unique(ocorrencias.treino[,1])
      
    }else{
    
    for (x in 1:length(especies)){
      ocorrencias.unicas.especie<-ocorrencias.unicas[ocorrencias.unicas[,1]==especies[x],]
      
      #Retirar espécies com numero de ocorrencias < x
      if (nrow(ocorrencias.unicas.especie)<5){
        print(paste("Ocorrencias de ",ocorrencias.unicas.especie[,1]," insuficientes para criar um modelo",sep=""))
      }else{
      id.ocor<-sample(1:nrow(ocorrencias.unicas.especie),round(porcentagem*nrow(ocorrencias.unicas.especie)))
      occ.treino.sp<-ocorrencias.unicas.especie[id.ocor,]
      occ.teste.sp<-ocorrencias.unicas.especie[-id.ocor,]
      ocorrencias.treino<-rbind(ocorrencias.treino,occ.treino.sp)
      ocorrencias.teste<-rbind(ocorrencias.teste,occ.teste.sp)
      }
    }
    setwd(dir2)
    ocorrencias.treino[,1]<-paste(ocorrencias.treino[,1],y,sep="")
    ocorrencias.teste[,1]<-paste(ocorrencias.teste[,1],y,sep="")
    write.xlsx(ocorrencias.treino[,c(1:3)],paste("Ocorrencias_treino",y,".xls",sep=""),row.names=F)
    write.xlsx(ocorrencias.teste[,c(1:3)],paste("Ocorrencias_teste",y,".xls",sep=""),row.names=F)
    
    especies.treino<-unique(ocorrencias.treino[,1])
    }  
      #BIOCLIM
      if (any(algoritmo=='bioclim')){
        
        setwd(dir.b)
        
        for (g in 1:length(especies.treino)){
          ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[g],]
          Bioclim.modelo<-bioclim(ocorrencias.treino.sp[,-c(1:4)])
          Bioclim.sp<-predict(Bioclim.modelo,project)
          writeRaster(Bioclim.sp,paste(especies.treino[g],'.asc',sep=''),format='ascii')
          print(paste("Bioclim ",especies.treino[g],"......OK",sep=""))
        }
      }
      
      #MAHALANOBIS
      if (any(algoritmo=='mahal')){
      
       setwd(dir.mahal)
      
        for (g in 1:length(especies.treino)){
          ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[g],]
          Mahal.modelo<-mahal(ocorrencias.treino.sp[,-c(1:4)])
          Mahal.sp<-predict(Mahal.modelo,project)
          writeRaster(Mahal.sp,paste(especies.treino[g],'.asc',sep=''),format='ascii')
          print(paste("Mahalanobis ",especies.treino[g],"......OK",sep=""))
        }
      }
      
      #MAXENT SIMPLES
      if (any(algoritmo=='maxentS')){
        Sys.setenv(NOAWT=TRUE)

        setwd(dir.maxS)
        
        for (h in 1:length(especies.treino)){
          ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[h],]
          Maxent.modelo<-maxent(x=environmental,p=ocorrencias.treino.sp[,c(2:3)],args=c(
                                "randomseed","noproduct","nothreshold","nohinge","noautofeature","nodoclamp","maximumiterations=1000"))
          Maxent.sp<-raster::predict(object=project,model=Maxent.modelo)
          writeRaster(Maxent.sp,paste(especies.treino[h],'.asc',sep=''),format='ascii')
          print(paste("Maxent Simples ",especies.treino[h],"......OK",sep=""))
        }
      }
    
      #MAXENT DEFAULT
      if (any(algoritmo=='maxentD')){
        Sys.setenv(NOAWT=TRUE)
      
        setwd(dir.maxD)
      
        for (h in 1:length(especies.treino)){
          ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[h],]
          Maxent.default<-maxent(x=environmental,p=ocorrencias.treino.sp[,c(2:3)])
          Maxent.sp.default<-raster::predict(object=project,model=Maxent.default)
          writeRaster(Maxent.sp.default,paste(especies.treino[h],'.asc',sep=''),format='ascii')
          print(paste("Maxent Default ",especies.treino[h],"......OK",sep=""))
        }
      }
      
      #SVM
      if (any(algoritmo=='svm')){
        
        setwd(dir.s)
        
        for (i in 1:length(especies.treino)){
          ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[i],]
          
          #Agrupar background e presencas
          pb.1<-rep(1,nrow(ocorrencias.treino.sp))
          ocorrencias.treino.max<-cbind(ocorrencias.treino.sp,pb.1)
          ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
          species.pb<-rep(especies.treino[i],nrow(ocorrencias.treino.max))
          pb.0<-rep(0,nrow(ocorrencias.treino.max))
          background.sample<-sample(1:nrow(background.env),length(pb.0))
          background.sample<-background.env[background.sample,]
          background.env.max<-as.data.frame(cbind(background.sample,pb.0))
          background.env.max<-cbind(species.pb,background.env.max)
          colnames(background.env.max)<-colnames(ocorrencias.treino.max)
          input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
          input.occ.back<-input.occ.back[,-c(1:3)]
          row.names(input.occ.back)<-NULL

          #Rodar o modelo SVM
          svm.modelo<-ksvm(pb.1~.,data=input.occ.back)
          
          SVM.sp<-predict(project,svm.modelo)
          writeRaster(SVM.sp,paste(especies.treino[i],'.asc',sep=''),format='ascii')
          print(paste("SVM ",especies.treino[i],"......OK",sep=""))
        }
      }
      
      #GLM
      if (any(algoritmo=='glm')){
        
        setwd(dir.g)
        
        for (j in 1:length(especies.treino)){
          ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[j],]
          
          #Agrupar background e presencas
          pb.1<-rep(1,nrow(ocorrencias.treino.sp))
          ocorrencias.treino.max<-cbind(ocorrencias.treino.sp,pb.1)
          ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
          species.pb<-rep(especies.treino[j],nrow(ocorrencias.treino.max))
          pb.0<-rep(0,nrow(ocorrencias.treino.max))
          background.sample<-sample(1:nrow(background.env),length(pb.0))
          background.sample<-background.env[background.sample,]
          background.env.max<-as.data.frame(cbind(background.sample,pb.0))
          background.env.max<-cbind(species.pb,background.env.max)
          colnames(background.env.max)<-colnames(ocorrencias.treino.max)
          input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
          input.occ.back<-input.occ.back[,-c(1:3)]
          row.names(input.occ.back)<-NULL
          input.occ.back.quad<-data.frame(input.occ.back^2)
        
          
          #Rodar o modelo GLM
          glm.modelo<-glm(pb.1~.,data=input.occ.back,family=binomial(link='logit'))#Modelo Linear
          
          glm.modelo.quad<-glm(pb.1~.,data=input.occ.back.quad,family=binomial(link='logit'))#Modelo Quadrático
          
          #AIC.GLM<-AICc(glm.modelo.quad, return.K = FALSE, second.ord = TRUE,nobs = NULL, c.hat = 1)
          
          GLM.sp<-predict(project,glm.modelo)
          writeRaster(GLM.sp,paste(especies.treino[j],'.asc',sep=''),format='ascii')
          print(paste("GLM ",especies.treino[j],"......OK",sep=""))
        }
      }
    
    #Random Forest
    if (any(algoritmo=='RF')){
      
      setwd(dir.rf)
      
      for (j in 1:length(especies.treino)){
        ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[j],]
        
        #Agrupar background e presencas
        pb.1<-rep(1,nrow(ocorrencias.treino.sp))
        ocorrencias.treino.max<-cbind(ocorrencias.treino.sp,pb.1)
        ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
        species.pb<-rep(especies.treino[j],nrow(ocorrencias.treino.max))
        pb.0<-rep(0,nrow(ocorrencias.treino.max))
        background.sample<-sample(1:nrow(background.env),length(pb.0))
        background.sample<-background.env[background.sample,]
        background.env.max<-as.data.frame(cbind(background.sample,pb.0))
        background.env.max<-cbind(species.pb,background.env.max)
        colnames(background.env.max)<-colnames(ocorrencias.treino.max)
        input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
        input.occ.back<-input.occ.back[,-c(1:3)]
        row.names(input.occ.back)<-NULL
        
        
        #Rodar o modelo Random Forest
        RF.modelo<-tuneRF(input.occ.back[,-ncol(input.occ.back)],input.occ.back$pb.1,nTreeTry=2000,stepFactor=2,doBest=TRUE,plot=FALSE,trace=FALSE)
        
        RF.sp<-predict(object=RF.modelo,newdata=project.env,type="response")
        RF.xy<-data.frame(cbind(project.env[,1:2],RF.sp))
        gridded(RF.xy)<-~x+y
        RF.raster<-raster(RF.xy)
        writeRaster(RF.raster,paste(especies.treino[j],'.asc',sep=''),format='ascii')
        print(paste("Random Forest ",especies.treino[j],"......OK",sep=""))
      }
    }
    
    #Mixture Discriminant Analysis (MDA)
    if (any(algoritmo=='MDA')){
      
      setwd(dir.mda)
      
      for (j in 1:length(especies.treino)){
        ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[j],]
        
        #Agrupar background e presencas
        pb.1<-rep(1,nrow(ocorrencias.treino.sp))
        ocorrencias.treino.max<-cbind(ocorrencias.treino.sp,pb.1)
        rm(pb.1)
        ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
        species.pb<-rep(especies.treino[j],nrow(ocorrencias.treino.max))
        pb.0<-rep(0,nrow(ocorrencias.treino.max))
        background.sample<-sample(1:nrow(background.env),length(pb.0))
        background.sample<-background.env[background.sample,]
        background.env.max<-as.data.frame(cbind(background.sample,pb.0))
        background.env.max<-cbind(species.pb,background.env.max)
        colnames(background.env.max)<-colnames(ocorrencias.treino.max)
        input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
        input.occ.back<-input.occ.back[,-c(1:3)]
        row.names(input.occ.back)<-NULL
        
        
        #Rodar o modelo MDA
        MDA.modelo<-mda(pb.1~.,data=input.occ.back)
        
        MDA.sp<-predict(project,MDA.modelo)
        writeRaster(MDA.sp,paste(especies.treino[j],'.asc',sep=''),format='ascii')
        print(paste("MDA ",especies.treino[j],"......OK",sep=""))
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
      
      setwd(dir.b)
      
      for (g in 1:length(especies.treino)){
        ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[g],]
        Bioclim.modelo<-bioclim(ocorrencias.treino.sp[,-c(1:4)])
        Bioclim.sp<-predict(Bioclim.modelo,project)
        writeRaster(Bioclim.sp,paste(especies.treino[g],'.asc',sep=''),format='ascii')
      }
    }
    
    #MAHALANOBIS
    if (any(algoritmo=='mahal')){
      
      setwd(dir.mahal)
      
      for (g in 1:length(especies.treino)){
        ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[g],]
        Mahal.modelo<-mahal(ocorrencias.treino.sp[,-c(1:4)])
        Mahal.sp<-predict(Mahal.modelo,project)
        writeRaster(Mahal.sp,paste(especies.treino[g],'.asc',sep=''),format='ascii')
      }
    }
    
    #MAXENT SIMPLES
    if (any(algoritmo=='maxentS')){
      Sys.setenv(NOAWT=TRUE)
      
      setwd(dir.maxS)
      
      for (h in 1:length(especies.treino)){
        ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[h],]
        Maxent.modelo<-maxent(x=environmental,p=ocorrencias.treino.sp[,c(2:3)],args=c(
          "randomseed","noproduct","nothreshold","nohinge","noautofeature","nodoclamp","maximumiterations=1000"))
        Maxent.sp<-predict(Maxent.modelo,project)
        writeRaster(Maxent.sp,paste(especies.treino[h],'.asc',sep=''),format='ascii')
      }
    }
    
    #MAXENT DEFAULT
    if (any(algoritmo=='maxentD')){
      Sys.setenv(NOAWT=TRUE)
      
      setwd(dir.maxD)
      
      for (h in 1:length(especies.treino)){
        ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[h],]
        Maxent.default<-maxent(x=environmental,p=ocorrencias.treino.sp[,c(2:3)])
        Maxent.sp.default<-predict(Maxent.default,project)
        writeRaster(Maxent.sp.default,paste(especies.treino[h],'.asc',sep=''),format='ascii')
      }
    }
    
    #SVM
    if (any(algoritmo=='svm')){
      
      setwd(dir.s)
      
      for (i in 1:length(especies.treino)){
        ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[i],]
        
        #Agrupar background e presencas
        pb.1<-rep(1,nrow(ocorrencias.treino.sp))
        ocorrencias.treino.max<-cbind(ocorrencias.treino.sp,pb.1)
        ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
        species.pb<-rep(especies.treino[i],nrow(ocorrencias.treino.max))
        pb.0<-rep(0,nrow(ocorrencias.treino.max))
        background.sample<-sample(1:nrow(background.env),length(pb.0))
        background.sample<-background.env[background.sample,]
        background.env.max<-as.data.frame(cbind(background.sample,pb.0))
        background.env.max<-cbind(species.pb,background.env.max)
        colnames(background.env.max)<-colnames(ocorrencias.treino.max)
        input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
        input.occ.back<-input.occ.back[,-c(1:3)]
        row.names(input.occ.back)<-NULL
        
        #Rodar o modelo SVM
        svm.modelo<-ksvm(pb.1~.,data=input.occ.back)
        
        SVM.sp<-predict(project,svm.modelo)
        writeRaster(SVM.sp,paste(especies.treino[i],'.asc',sep=''),format='ascii')
      }
    }
    
    #GLM
    if (any(algoritmo=='glm')){
      
      setwd(dir.g)
      
      for (j in 1:length(especies.treino)){
        ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[j],]
        
        #Agrupar background e presencas
        pb.1<-rep(1,nrow(ocorrencias.treino.sp))
        ocorrencias.treino.max<-cbind(ocorrencias.treino.sp,pb.1)
        ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
        species.pb<-rep(especies.treino[j],nrow(ocorrencias.treino.max))
        pb.0<-rep(0,nrow(ocorrencias.treino.max))
        background.sample<-sample(1:nrow(background.env),length(pb.0))
        background.sample<-background.env[background.sample,]
        background.env.max<-as.data.frame(cbind(background.sample,pb.0))
        background.env.max<-cbind(species.pb,background.env.max)
        colnames(background.env.max)<-colnames(ocorrencias.treino.max)
        input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
        input.occ.back<-input.occ.back[,-c(1:3)]
        row.names(input.occ.back)<-NULL
        input.occ.back.quad<-data.frame(input.occ.back^2)
        
        
        #Rodar o modelo GLM
        glm.modelo<-glm(pb.1~.,data=input.occ.back,family=binomial(link='logit'))#Modelo Linear
        
        glm.modelo.quad<-glm(pb.1~.,data=input.occ.back.quad,family=binomial(link='logit'))#Modelo Quadrático
        
        AIC.GLM<-AICc(glm.modelo.quad, return.K = FALSE, second.ord = TRUE,nobs = NULL, c.hat = 1)
        
        GLM.sp<-predict(project,glm.modelo)
        writeRaster(GLM.sp,paste(especies.treino[j],'.asc',sep=''),format='ascii')
      }
    }
    
    #Random Forest
    if (any(algoritmo=='RF')){
      
      setwd(dir.rf)
      
      for (j in 1:length(especies.treino)){
        ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[j],]
        
        #Agrupar background e presencas
        pb.1<-rep(1,nrow(ocorrencias.treino.sp))
        ocorrencias.treino.max<-cbind(ocorrencias.treino.sp,pb.1)
        ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
        species.pb<-rep(especies.treino[j],nrow(ocorrencias.treino.max))
        pb.0<-rep(0,nrow(ocorrencias.treino.max))
        background.sample<-sample(1:nrow(background.env),length(pb.0))
        background.sample<-background.env[background.sample,]
        background.env.max<-as.data.frame(cbind(background.sample,pb.0))
        background.env.max<-cbind(species.pb,background.env.max)
        colnames(background.env.max)<-colnames(ocorrencias.treino.max)
        input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
        input.occ.back<-input.occ.back[,-c(1:3)]
        row.names(input.occ.back)<-NULL
        
        
        #Rodar o modelo Random Forest
        RF.modelo<-tuneRF(input.occ.back[,-ncol(input.occ.back)],input.occ.back$pb.1,nTreeTry=2000,stepFactor=2,doBest=TRUE,plot=FALSE,trace=FALSE)
        
        RF.sp<-predict(object=RF.modelo,newdata=project.env,type="response")
        RF.xy<-data.frame(cbind(project.env[,1:2],RF.sp))
        gridded(RF.xy)<-~x+y
        RF.raster<-raster(RF.xy)
        writeRaster(RF.raster,paste(especies.treino[j],'.asc',sep=''),format='ascii')
      }
    }
    
    #Mixture Discriminant Analysis (MDA)
    if (any(algoritmo=='MDA')){
      
      setwd(dir.mda)
      
      for (j in 1:length(especies.treino)){
        ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[j],]
        
        #Agrupar background e presencas
        pb.1<-rep(1,nrow(ocorrencias.treino.sp))
        ocorrencias.treino.max<-cbind(ocorrencias.treino.sp,pb.1)
        ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
        species.pb<-rep(especies.treino[j],nrow(ocorrencias.treino.max))
        pb.0<-rep(0,nrow(ocorrencias.treino.max))
        background.sample<-sample(1:nrow(background.env),length(pb.0))
        background.sample<-background.env[background.sample,]
        background.env.max<-as.data.frame(cbind(background.sample,pb.0))
        background.env.max<-cbind(species.pb,background.env.max)
        colnames(background.env.max)<-colnames(ocorrencias.treino.max)
        input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
        input.occ.back<-input.occ.back[,-c(1:3)]
        row.names(input.occ.back)<-NULL
        
        
        #Rodar o modelo MDA
        MDA.modelo<-mda(pb.1~.,data=input.occ.back)
        
        MDA.sp<-predict(project,MDA.modelo)
        
        writeRaster(RF.raster,paste(especies.treino[j],'.asc',sep=''),format='ascii')
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
  
if (particao =="leave1out"){
  
  ocorrencias.treino<-NULL
  ocorrencias.teste<-NULL
  
  for (o in 1:length(especies)){
    
    ocorrencias.unicas.sp<-ocorrencias.unicas[ocorrencias.unicas[,1]==especies[o],]

    for (y in 1:nrow(ocorrencias.unicas.sp)){
      
      occ.treino.sp<-ocorrencias.unicas.sp[-y,]
      occ.teste.sp<-ocorrencias.unicas.sp[y,]
      occ.treino.sp[,1]<-paste(occ.treino.sp[,1],y,sep="")
      occ.teste.sp[,1]<-paste(occ.teste.sp[,1],y,sep="")
      ocorrencias.treino<-rbind(ocorrencias.treino,occ.treino.sp)
      ocorrencias.teste<-rbind(ocorrencias.teste,occ.teste.sp)
        
      ocorrencias.treino.y<-occ.treino.sp
      especies.treino<-unique(ocorrencias.treino.y[,1])
        
      #BIOCLIM
      if (any(algoritmo=='bioclim')){
        
        setwd(dir.b)
        
          Bioclim.modelo<-bioclim(ocorrencias.treino.y[,-c(1:4)])
          Bioclim.sp<-predict(Bioclim.modelo,project)
          writeRaster(Bioclim.sp,paste(especies.treino,'.asc',sep=''),format='ascii')
      }
      
      #MAHALANOBIS
      if (any(algoritmo=='mahal')){
        
        setwd(dir.mahal)
          Mahal.modelo<-mahal(ocorrencias.treino.y[,-c(1:4)])
          Mahal.sp<-predict(Mahal.modelo,project)
          writeRaster(Mahal.sp,paste(especies.treino,'.asc',sep=''),format='ascii')
      }
      
      #MAXENT SIMPLES
      if (any(algoritmo=='maxentS')){
        Sys.setenv(NOAWT=TRUE)
        
        setwd(dir.maxS)
        
          Maxent.modelo<-maxent(x=environmental,p=ocorrencias.treino.y[,c(2:3)],args=c(
            "randomseed","noproduct","nothreshold","nohinge","noautofeature","nodoclamp","maximumiterations=1000"))
          Maxent.sp<-predict(Maxent.modelo,project)
          writeRaster(Maxent.sp,paste(especies.treino,'Simples.asc',sep=''),format='ascii')
      }
      
      #MAXENT DEFAULT
      if (any(algoritmo=='maxentD')){
        Sys.setenv(NOAWT=TRUE)
        
        setwd(dir.maxD)
        
          Maxent.default<-maxent(x=environmental,p=ocorrencias.treino.y[,c(2:3)])
          Maxent.sp.default<-predict(Maxent.default,project)
          writeRaster(Maxent.sp.default,paste(especies.treino,'Default.asc',sep=''),format='ascii')
      }
      
      #SVM
      if (any(algoritmo=='svm')){
        
        setwd(dir.s)
        
          #Agrupar background e presencas
          pb.1<-rep(1,nrow(ocorrencias.treino.y))
          ocorrencias.treino.max<-cbind(ocorrencias.treino.y,pb.1)
          ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
          species.pb<-rep(especies.treino,nrow(ocorrencias.treino.max))
          pb.0<-rep(0,nrow(ocorrencias.treino.max))
          background.sample<-sample(1:nrow(background.env),length(pb.0))
          background.sample<-background.env[background.sample,]
          background.env.max<-as.data.frame(cbind(background.sample,pb.0))
          background.env.max<-cbind(species.pb,background.env.max)
          colnames(background.env.max)<-colnames(ocorrencias.treino.max)
          input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
          input.occ.back<-input.occ.back[,-c(1:3)]
          row.names(input.occ.back)<-NULL
          
          #Rodar o modelo SVM
          svm.modelo<-ksvm(pb.1~.,data=input.occ.back)
          
          SVM.sp<-predict(project,svm.modelo)
          writeRaster(SVM.sp,paste(especies.treino,'.asc',sep=''),format='ascii')
      }
      
      #GLM
      if (any(algoritmo=='glm')){
        
        setwd(dir.g)
        
          #Agrupar background e presencas
          pb.1<-rep(1,nrow(ocorrencias.treino.y))
          ocorrencias.treino.max<-cbind(ocorrencias.treino.y,pb.1)
          ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
          species.pb<-rep(especies.treino,nrow(ocorrencias.treino.max))
          pb.0<-rep(0,nrow(ocorrencias.treino.max))
          background.sample<-sample(1:nrow(background.env),length(pb.0))
          background.sample<-background.env[background.sample,]
          background.env.max<-as.data.frame(cbind(background.sample,pb.0))
          background.env.max<-cbind(species.pb,background.env.max)
          colnames(background.env.max)<-colnames(ocorrencias.treino.max)
          input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
          input.occ.back<-input.occ.back[,-c(1:3)]
          row.names(input.occ.back)<-NULL
          input.occ.back.quad<-data.frame(input.occ.back^2)
          
          
          #Rodar o modelo GLM
          glm.modelo<-glm(pb.1~.,data=input.occ.back,family=binomial(link='logit'))#Modelo Linear
          
          glm.modelo.quad<-glm(pb.1~.,data=input.occ.back.quad,family=binomial(link='logit'))#Modelo Quadrático
          
          AIC.GLM<-AICc(glm.modelo.quad, return.K = FALSE, second.ord = TRUE,nobs = NULL, c.hat = 1)
          
          GLM.sp<-predict(project,glm.modelo)
          writeRaster(GLM.sp,paste(especies.treino,'.asc',sep=''),format='ascii')
      }
      
      #Random Forest
      if (any(algoritmo=='RF')){
        
        setwd(dir.rf)
        
          #Agrupar background e presencas
          pb.1<-rep(1,nrow(ocorrencias.treino.y))
          ocorrencias.treino.max<-cbind(ocorrencias.treino.y,pb.1)
          ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
          species.pb<-rep(especies.treino,nrow(ocorrencias.treino.max))
          pb.0<-rep(0,nrow(ocorrencias.treino.max))
          background.sample<-sample(1:nrow(background.env),length(pb.0))
          background.sample<-background.env[background.sample,]
          background.env.max<-as.data.frame(cbind(background.sample,pb.0))
          background.env.max<-cbind(species.pb,background.env.max)
          colnames(background.env.max)<-colnames(ocorrencias.treino.max)
          input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
          input.occ.back<-input.occ.back[,-c(1:3)]
          row.names(input.occ.back)<-NULL
          
          
          #Rodar o modelo Random Forest
          RF.modelo<-tuneRF(input.occ.back[,-ncol(input.occ.back)],input.occ.back$pb.1,nTreeTry=2000,stepFactor=2,doBest=TRUE,plot=FALSE,trace=FALSE)
          
          RF.sp<-predict(object=RF.modelo,newdata=project.env,type="response")
          RF.xy<-data.frame(cbind(project.env[,1:2],RF.sp))
          gridded(RF.xy)<-~x+y
          RF.raster<-raster(RF.xy)
          writeRaster(RF.raster,paste(especies.treino,'.asc',sep=''),format='ascii')
      }
      
      #Mixture Discriminant Analysis (MDA)
      if (any(algoritmo=='MDA')){
        
        setwd(dir.mda)

          #Agrupar background e presencas
          pb.1<-rep(1,nrow(ocorrencias.treino.y))
          ocorrencias.treino.max<-cbind(ocorrencias.treino.y,pb.1)
          ocorrencias.treino.max<-ocorrencias.treino.max[,-4]
          species.pb<-rep(especies.treino,nrow(ocorrencias.treino.max))
          pb.0<-rep(0,nrow(ocorrencias.treino.max))
          background.sample<-sample(1:nrow(background.env),length(pb.0))
          background.sample<-background.env[background.sample,]
          background.env.max<-as.data.frame(cbind(background.sample,pb.0))
          background.env.max<-cbind(species.pb,background.env.max)
          colnames(background.env.max)<-colnames(ocorrencias.treino.max)
          input.occ.back<-rbind(ocorrencias.treino.max,background.env.max)
          input.occ.back<-input.occ.back[,-c(1:3)]
          row.names(input.occ.back)<-NULL
          
          
          #Rodar o modelo MDA
          MDA.modelo<-mda(pb.1~.,data=input.occ.back)
          
          MDA.sp<-predict(project,MDA.modelo)
          
          writeRaster(RF.raster,paste(especies.treino,'.asc',sep=''),format='ascii')
      }
      
      #ENFA
      #if (any(algoritmo=="enfa")){
      #for (j in 1:length(especies.treino)){
      #ocorrencias.treino.sp<-ocorrencias.treino[ocorrencias.treino[,1]==especies.treino[j],]
      #}
      #}
    }
  }
  
  setwd(dir2)
  write.table(ocorrencias.treino[,c(1:3)],paste("Ocorrencias_treino",".txt",sep=""),row.names=F,sep="\t")
  write.table(ocorrencias.teste[,c(1:3)],paste("Ocorrencias_teste",".txt",sep=""),row.names=F,sep="\t")
  
  }
} # fecha a funcao

#PARA INCLUIR NA FUNCAO:
  #2.A funcao response("resultado do modelo") mostra as curvas de resposa para cada variável

ENMs_TheMetaLand('D:\\Users\\Andre\\Google Drive\\Mestrado\\Layers\\Bioclim\\Presente_10km','raster',pca='S',project='N',10000,"bootstrap",algoritmo=c('maxentS','svm','RF'))

diretorio<-'C:\\R\\Teste_codigo\\Var_Originais'
rm(list=ls())
gc()
