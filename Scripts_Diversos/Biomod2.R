library(biomod2)
?biomod2

#Importar as ocorrencias
DataSpecies <- read.csv("C:/Users/Andre/Desktop/Andre/Mestrado/Ocorrencias/Corais/GBIF+Literatura/Mussismilia hispida/M_hispida_GBIF+REV.csv",header=TRUE)
DataSpecies <- DataSpecies[,c(1:3)]
colnames(DataSpecies) <- c("MussismiliaHispida","latitude","longitude")
head(DataSpecies)

  #Dividir a Amostra para 5-fold cross-validation
  #library(dismo)
  #fold<- kfold(DataSpecies, k=5)
  #DataSpecies<-cbind(DataSpecies,fold)

  #Selecionar todas as coordenadas, menos as que sejam do grupo 1, 80%(diferentes de = "!")
  #occtrain<- DataSpecies[fold !=1,]

  #Selecionar 20% das ocorrencias para teste (que sejam igual a = "==")
  #occtest<- DataSpecies[fold ==1,]

#Ocorrencias treino para inserçao no biomod
myRespName <- 'MussismiliaHispida'
myResp <- as.numeric(DataSpecies[,myRespName])
myRespXY <- DataSpecies[,c("longitude","latitude")]

  #Ocorrencias teste para inserçao no biomod
  #evalResp <- as.numeric(occtest[,myRespName])
  #evalRespXY <- occtest[,c("latitude","longitude")]

#Importar/stack das variáveis ambientais
myExpl = stack("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca1.asc",
               "C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca2.asc",
               "C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca3.asc",
               "C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca4.asc",
               "C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca5.asc",
               "C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca6.asc",
               "C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca7.asc",
               "C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca8.asc",
               "C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca9.asc")

#Formatar os dados de entrada para o biomod
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep=1,
                                     PA.nb.absences=1000,
                                     PA.strategy="random",
                                     na.rm=TRUE)

#Checar o arquivo de entrada
myBiomodData

#Definir os parametros para o MAXENT
myBiomodOptions <- BIOMOD_ModelingOptions()

#Check parameters
myBiomodOptions

#Run model
myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
                                     models=c("MAXENT"),
                                     models.options=myBiomodOptions,
                                     NbRunEval=1,
                                     DataSplit=80,
                                     Yweights=NULL,
                                     Prevalence=NULL,
                                     VarImport=0,
                                     models.eval.meth=c("ROC","TSS"),
                                     SaveObj=TRUE,
                                     rescal.all.models=FALSE,
                                     do.full.models=FALSE,
                                     modeling.id="teste")

#Summary of modeling stuff
myBiomodModelOut

#Projetar o output dos modelos
myBiomodProjection <- BIOMOD_Projection(modeling.output=myBiomodModelOut,
                                        new.env=myExpl,
                                        proj.name="current",
                                        selected.models="all",
                                        binary.meth="TSS",
                                        compress=FALSE,
                                        build.clamping.mask=TRUE)


# print summary and plot projections
myBiomodProjection
windows()
plot(myBiomodProjection)
