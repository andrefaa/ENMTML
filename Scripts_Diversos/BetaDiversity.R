library(vegan)
library(betapart)
setwd(choose.dir())

#########################################################################################
############## Copie essas funcoeses (ate a proxima linha de "#") no seu R ################
#########################################################################################

beta.sor<-function(x){
  ## x must be a data frame
  matr<-as.matrix(x);
  nr<-dim(matr)[1];
  result<-matrix(nrow=nr,ncol=nr);
  rownames(result)<-rownames(matr);
  colnames(result)<-rownames(matr);
  for(i in 1:nr) {
    for(j in i:nr) {
      bij<-sum(matr[i,]&(!matr[j,]));
      cij<-sum(matr[j,]&(!matr[i,]));
      aij<-sum(matr[i,]&(matr[j,]));
      Sorensen.ij<-(bij+cij)/((2*aij)+bij+cij);
      result[i,j]<-Sorensen.ij;
      result[j,i]<-Sorensen.ij;
    }
  }
  d<-as.dist(result);
  d
}

beta.sim<-function(x){
  ## x must be a data frame
  matr<-as.matrix(x);
  nr<-dim(matr)[1];
  result<-matrix(nrow=nr,ncol=nr);
  rownames(result)<-rownames(matr);
  colnames(result)<-rownames(matr);
  for(i in 1:nr) {
    for(j in i:nr) {
      bij<-sum(matr[i,]&(!matr[j,]));
      cij<-sum(matr[j,]&(!matr[i,]));
      aij<-sum(matr[i,]&(matr[j,]));
      Simpson.ij<-min(bij,cij)/(min(bij,cij)+aij);
      result[i,j]<-Simpson.ij;
      result[j,i]<-Simpson.ij;
    }
  }
  d<-as.dist(result);
  d
}

beta.nes<-function(x){
  ## x must be a data frame
  matr<-as.matrix(x);
  nr<-dim(matr)[1];
  result<-matrix(nrow=nr,ncol=nr);
  rownames(result)<-rownames(matr);
  colnames(result)<-rownames(matr);
  for(i in 1:nr) {
    for(j in i:nr) {
      bij<-sum(matr[i,]&(!matr[j,]));
      cij<-sum(matr[j,]&(!matr[i,]));
      aij<-sum(matr[i,]&(matr[j,]));
      nestedness.ij<-((max(bij,cij)-min(bij,cij))/((2*aij)+bij+cij))*(aij/(min(bij,cij)+aij));
      result[i,j]<-nestedness.ij;
      result[j,i]<-nestedness.ij;
    }
  }
  d<-as.dist(result);
  d
}

#########################################################################################
################################# So Ate Aqui!!!!! ######################################
#########################################################################################

# Carregue os dados, que voce ja salvou como texto(MSDOS), nao serve outro formato!!!!!
# A sua planilha tem o nome das especies (sp1, spn) mas nao tem rotulo de linha, identificador 
# dos locais, caso tenha, nao vai funcionar!!!!!!!

dados<-read.table(file=choose.files(), head=T,sep="\t")
dados.xy <- dados[,c(1:2)]
dados<-dados[,-c(1:2)]


## Calcula a diversidade por sorense
dir()

#dados <-read.table("Zygo_All_teste1.txt", header=TRUE)

a<-beta.sor(dados)
b<-beta.sim(dados)
c<-beta.nes(dados)


## Criando um vetor com as medias das diversidade
sor<-as.vector(rowSums(as.matrix(a))/(nrow(as.matrix(a))-1))
sim<-as.vector(rowSums(as.matrix(b))/(nrow(as.matrix(b))-1))
nes<-as.vector(rowSums(as.matrix(c))/(nrow(as.matrix(c))-1))

# Juntando tudo e salvando as coisas como um planilha para que vc possa fazer analises no Excell (no Statistica e melhor!!!)!!!

res<-cbind(sor,sim,nes)
res<-cbind(dados[,1:2],res)
write.table(res,"BetaAMZ.txt",sep="\t",row.names=F)

res.a<-cbind(dados.xy,sor)
res.a2 <- res.a
gridded(res.a) <- ~ x+y
res.a <- raster(res.a)
plot(res.a)
writeRaster(res.a,"SorensenBeta.asc",format="ascii")
#########################################################################################
################################# TESTE DE HIPOTESE #####################################
#########################################################################################

dados <- as.matrix(dados)

beta_res <- beta.multi(dados, index.family="sorensen")
beta_res

beta_null <- oecosimu(decostand(dados,"pa"),
                      nestedbetasor,
                      method="r1", nsimul = 999)


beta_null
