##
## written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
## 2nd of November 2010. University of Lausanne. Department of ecology and evolution, Switzerland
##
## DESCRIPTION
## This script investigates ecological niche overlap from occurrence and spatial environmental data
## 
## INSTALLATION
## 
## Put script.R, niche.overlap.functions.R and occ.prep.functions.R in a folder with a R shortcut.
## Use this folder as workspace by setting the path in the proprieties of the shortcut (right click)
## In this folder also put your 2 datasets of occurences data (column names should be x,y) and
## 2 datesets of points representing the study areas with environmental values (column names should be x,y,X1,X2,...,Xn)
## (these datasets are usually derived from a sampling of the study areas in a GIS) 
## The main script script.R allows setting the analyses and sources the functions necessary for the calculations 
## (niche.overlap.functions.R and occ.prep.functions.R).
##


#################################################################################################
############################## load functions and packages ######################################
#################################################################################################

source("niche.overlap.functions.R")
source("occ.prep.functions.R")
source("plot_niche_blue.R")

#install.packages("BIOMOD",repos="http://R-Forge.R-project.org")
library(biomod2)
library(ade4)
library(adehabitat)
library(sp)
library(gam)
library(MASS)
library(mvtnorm)
library(gbm)
library(dismo)

#################################################################################################
############################## preparation of datasets ##########################################
#################################################################################################

# load climate variable for all site of the study area 1 (column names should be x,y,X1,X2,...,Xn)
clim1<-na.exclude(read.delim("Env_area_1.txt",h=T,sep="\t"))

# load climate variable for all site of the study area 2 (column names should be x,y,X1,X2,...,Xn)
clim2<-na.exclude(read.delim("Env_area_2.txt",h=T,sep="\t")) 

# global climate for both ranges
clim12<-rbind(clim1,clim2)

# loading occurrence sites for the species (column names should be x,y)
occ.sp.aggr<-na.exclude(read.delim("Ocorrencias_Mhispida_xy.txt",h=T,sep="\t"))

# remove occurrences closer than a minimum distance to each other (remove aggregation). Setting min.dist=0 will remove no occurrence.

occ.sp<-occ.desaggragation(df=occ.sp.aggr,colxy=1:2,min.dist=0.16666,plot=F)

# create sp occurrence dataset by adding climate variables from the global climate datasets
# resolution should be the resolution of the climate data grid

occ.sp1<-na.exclude(sample.sp.globvar(dfsp=occ.sp,colspxy=1:2,colspkept=NULL,dfvar=clim1,colvarxy=1:2,colvar="all",resolution=0.16666))
occ.sp2<-na.exclude(sample.sp.globvar(dfsp=occ.sp,colspxy=1:2,colspkept=NULL,dfvar=clim2,colvarxy=1:2,colvar="all",resolution=0.16666))

# create presence/absence datasets (used in ENFA and SDMs)
row.pa1<-sample.sp.globvar(dfsp=clim1,colspxy=1:2,colspkept=NULL,dfvar=occ.sp1,colvarxy=1:2,colvar=3,resolution=0) 
#find rows of clim1 where the species is present
pa<-data.frame((!is.na(row.pa1))*1);names(pa)<-"pa" #create 01 column
pa1<-cbind(clim1,pa)

row.pa2<-sample.sp.globvar(dfsp=clim2,colspxy=1:2,colspkept=NULL,dfvar=occ.sp2,colvarxy=1:2,colvar=3,resolution=0) 
#find rows of clim1 where the species is present
pa<-data.frame((!is.na(row.pa2))*1);names(pa)<-"pa" #create 01 column
pa2<-cbind(clim2,pa)

#################################################################################################
############################## ANALYSIS - selection of parameters ###############################
#################################################################################################

# selection of the type of analysis.
# If PROJ =F, the models are calibrated on both ranges.
# If PROJ =T, the models are calibrated on species 1 range only and projected to range 2. 
# Analyses where both ranges are needed (ex: LDA) are not done
PROJ = T

# selection of variables to include in the analyses
names(clim12)
Xvar<-c(3:11)
nvar<-length(Xvar)

#number of interation for the tests of equivalency and similarity
iterations<-100

#resolution of the gridding of the climate space
R=1000

#################################################################################################
################### row weigthing and grouping factors for ade4 functions  ######################
#################################################################################################

# if PROJ = F
row.w.1.occ<-1-(nrow(occ.sp1)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ1
row.w.2.occ<-1-(nrow(occ.sp2)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ2
row.w.occ<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(row.w.1.occ, nrow(occ.sp1)),rep(row.w.2.occ, nrow(occ.sp2)))

row.w.1.env<-1-(nrow(clim1)/nrow(clim12))  # prevalence of clim1
row.w.2.env<-1-(nrow(clim2)/nrow(clim12))  # prevalence of clim2
row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))

fac<-as.factor(c(rep(1, nrow(clim1)),rep(2, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(2, nrow(occ.sp2))))

# if PROJ = T

row.w.occ.PROJT<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
row.w.env.PROJT<-c(rep(1, nrow(clim1)),rep(0, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))

# global dataset for the analysis and rows for each sub dataset
data.env.occ<-rbind(clim1,clim2,occ.sp1,occ.sp2)[Xvar]
row.clim1<-1:nrow(clim1)
row.clim2<-(nrow(clim1)+1):(nrow(clim1)+nrow(clim2))
row.clim12<-1:(nrow(clim1)+nrow(clim2))
row.sp1<-(nrow(clim1)+nrow(clim2)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1))
row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+nrow(occ.sp2))


#################################################################################################
#################################### one-variable ###############################################
#################################################################################################

      # measures niche overlap along one chosen variable

onevar<-5 #choose the variable to analyse here
nameonevar<-names(clim12)[onevar]

			# predict the scores on the axes
scores.clim12<- data.frame(clim12[,onevar]);names(scores.clim12)<-nameonevar
scores.clim2<- data.frame(clim2[,onevar]);names(scores.clim2)<-nameonevar
scores.clim1<- data.frame(clim1[,onevar]);names(scores.clim1)<-nameonevar
scores.sp2<- data.frame(occ.sp2[,onevar]);names(scores.sp2)<-nameonevar
scores.sp1<- data.frame(occ.sp1[,onevar]);names(scores.sp1)<-nameonevar

			# calculation of occurence density and test of niche equivalency and similarity 
z1<- grid.clim(scores.clim12,scores.clim1,scores.sp1,R)
z2<- grid.clim(scores.clim12,scores.clim2,scores.sp2,R)
a<-niche.equivalency.test(z1,z2,rep=100)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
plot.niche(z1,title=paste(nameonevar,"- EU niche"),name.axis1="PC1",name.axis2="PC2")
plot.niche(z2,title=paste(nameonevar,"- NA niche"),name.axis1="PC1",name.axis2="PC2")
plot.new()
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")

#################################################################################################
#################################### PCA-occ ####################################################
#################################################################################################

      # measures niche overlap along the two first axes of a PCA calibrated on occurence data
			
if(PROJ == F){	#fit of the analyse using occurences from both ranges
	pca.cal <-dudi.pca(data.env.occ,row.w = row.w.occ, center = T, scale = T, scannf = F, nf = 2)
}
if(PROJ == T) {   #fit of the analyse using occurences from range 1
	pca.cal <-dudi.pca(data.env.occ,row.w = row.w.occ.PROJT, center = T, scale = T, scannf = F, nf = 2)
}
			# predict the scores on the axes
scores.clim12<- pca.cal$li[row.clim12,]
scores.clim1<- pca.cal$li[row.clim1,]
scores.clim2<- pca.cal$li[row.clim2,]
scores.sp1<- pca.cal$li[row.sp1,]
scores.sp2<- pca.cal$li[row.sp2,]


			# calculation of occurence density and test of niche equivalency and similarity 
z1<- grid.clim(scores.clim12,scores.clim1,scores.sp1,R)
z2<- grid.clim(scores.clim12,scores.clim2,scores.sp2,R)
a<-niche.equivalency.test(z1,z2,rep=1)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
plot.niche(z1,title="PCA - EU niche",name.axis1="PC1",name.axis2="PC2")
plot.niche(z2,title="PCA - EU niche",name.axis1="PC1",name.axis2="PC2")
plot.contrib(pca.cal$co,pca.cal$eig)
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")

#################################################################################################
#################################### PCA-ENV ####################################################
#################################################################################################

      # measures niche overlap along the two first axes of a PCA calibrated on all the pixels of the study areas
      			
if(PROJ == F){	#fit of the analyse using occurences from both ranges		
	pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)
}
if(PROJ == T){	#fit of the analyse using occurences from range 1		
	pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env.PROJT, center = T, scale = T, scannf = F, nf = 2)
}
			# predict the scores on the axes
scores.clim12<- pca.cal$li[row.clim12,]
scores.clim1<- pca.cal$li[row.clim1,]
scores.clim2<- pca.cal$li[row.clim2,]
scores.sp1<- pca.cal$li[row.sp1,]
scores.sp2<- pca.cal$li[row.sp2,]

			# calculation of occurence density and test of niche equivalency and similarity 
z1<- grid.clim(scores.clim12,scores.clim1,scores.sp1,R)
z2<- grid.clim(scores.clim12,scores.clim2,scores.sp2,R)
a<-niche.equivalency.test(z1,z2,rep=1)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
x11()
plot.niche.blue(z1,title="PCA-env - EU niche",name.axis1="PC1",name.axis2="PC2")
par(new=T)
plot.niche(z2,title="PCA-env - NA niche",name.axis1="PC1",name.axis2="PC2")
plot.contrib(pca.cal$co,pca.cal$eig)
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")

#################################################################################################
#################################### WITHIN-GROUP ###############################################
#################################################################################################

      # measures niche overlap along the two first axes of a within analysis on occurence data with with ranges as a priori groups  			
if(PROJ == F){	
			#fit of the analyse using occurences from both ranges
pca.cal <-dudi.pca(data.env.occ,row.w = row.w.occ, center = T, scale = T, scannf = F, nf = nvar)
within.cal<-within(pca.cal, fac, scannf = F, nf = 2)

scores.clim12<- within.cal$li[row.clim12,]
scores.clim1<- within.cal$li[row.clim1,]
scores.clim2<- within.cal$li[row.clim2,]
scores.sp1<- within.cal$li[row.sp1,]
scores.sp2<- within.cal$li[row.sp2,]

			# calculation of occurence density and test of niche equivalency and similarity 
z1<- grid.clim(scores.clim12,scores.clim1,scores.sp1,R)
z2<- grid.clim(scores.clim12,scores.clim2,scores.sp2,R)
a<-niche.equivalency.test(z1,z2,rep=1)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
plot.niche(z1,title="WITHIN - EU niche",name.axis1="PC1",name.axis2="PC2")
plot.niche(z2,title="WITHIN - NA niche",name.axis1="PC1",name.axis2="PC2")
plot.contrib(within.cal$co,within.cal$eig)
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")
}

#################################################################################################
#################################### WITHIN-ENVIRONMENT #########################################
#################################################################################################

      # measures niche overlap along the two first axes of a within analysis based on all the sites of the study areas 
      # with with ranges as a priori groups
      
if(PROJ == F){	
			#fit of the analyse using global environments from both ranges
pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = nvar)
within.cal<-within(pca.cal, fac, scannf = F, nf = 2)

scores.clim12<- within.cal$li[row.clim12,]
scores.clim1<- within.cal$li[row.clim1,]
scores.clim2<- within.cal$li[row.clim2,]
scores.sp1<- within.cal$li[row.sp1,]
scores.sp2<- within.cal$li[row.sp2,]

			# calculation of occurence density and test of niche equivalency and similarity 
z1<- grid.clim(scores.clim12,scores.clim1,scores.sp1,R)
z2<- grid.clim(scores.clim12,scores.clim2,scores.sp2,R)
a<-niche.equivalency.test(z1,z2,rep=1)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
plot.niche(z1,title="WITHIN-env - EU niche",name.axis1="PC1",name.axis2="PC2")
plot.niche(z2,title="WITHIN-env - NA niche",name.axis1="PC1",name.axis2="PC2")
plot.contrib(within.cal$co,within.cal$eig)
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")
}
#################################################################################################
#################################### BETWEEN-GROUP ##############################################
#################################################################################################

      # measures niche overlap along the axis of a between analysis based on the occurence data 
      # with with ranges as a priori groups
      
if(PROJ == F){	
		#fit of the analyse using global environments from both ranges

pca.cal <-dudi.pca(data.env.occ,row.w = row.w.occ, center = T, scale = T, scannf = F, nf = nvar)
between.cal<-between(pca.cal, fac, scannf = F, nf = nvar)

scores.clim12<- data.frame(between.cal$ls[row.clim12,])
scores.clim1<- data.frame(between.cal$ls[row.clim1,])
scores.clim2<- data.frame(between.cal$ls[row.clim2,])
scores.sp1<- data.frame(between.cal$ls[row.sp1,])
scores.sp2<- data.frame(between.cal$ls[row.sp2,])

			# calculation of occurence density and test of niche equivalency and similarity 
z1<- grid.clim(scores.clim12,scores.clim1,scores.sp1,R)
z2<- grid.clim(scores.clim12,scores.clim2,scores.sp2,R)
a<-niche.equivalency.test(z1,z2,rep=100)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
plot.niche(z1,title="BETWEEN - EU niche",name.axis1="PC1",name.axis2="PC2")
plot.niche(z2,title="BETWEEN - NA niche",name.axis1="PC1",name.axis2="PC2")
plot.contrib(between.cal$co)
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")
}
#################################################################################################
#################################### LINEAR DISCRIMINANT ########################################
#################################################################################################

      # measures niche overlap along the axis of a discriminant analysis based on the occurence data 
      # with with ranges as a priori groups
      
if(PROJ == F){	
			#fit of the analyse using occurences from both ranges
fac<-as.factor(c(rep(1, nrow(occ.sp1)),rep(2, nrow(occ.sp2))))
lda.cal<-lda(rbind(occ.sp1,occ.sp2)[Xvar],fac)

			# predict the scores on the axes
scores.clim12<- predict(lda.cal,clim12[Xvar])$x
scores.clim2<- predict(lda.cal,clim2[Xvar])$x
scores.clim1<- predict(lda.cal,clim1[Xvar])$x
scores.sp2<- predict(lda.cal,occ.sp2[Xvar])$x
scores.sp1<- predict(lda.cal,occ.sp1[Xvar])$x

			# calculation of occurence density and test of niche equivalency and similarity 
z1<- grid.clim(scores.clim12,scores.clim1,scores.sp1,R)
z2<- grid.clim(scores.clim12,scores.clim2,scores.sp2,R)
a<-niche.equivalency.test(z1,z2,rep=100)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
plot.niche(z1,title="LDA - EU niche",name.axis1="PC1",name.axis2="PC2")
plot.niche(z2,title="LDA - NA niche",name.axis1="PC1",name.axis2="PC2")
plot.contrib(lda.cal$scaling)
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")
}
#################################################################################################
#################################### ENFA #######################################################
#################################################################################################

      # measures niche overlap along the two first axes of an Ecological Factor Analysis
      
if(PROJ == F){
  row.w.2<-1-(nrow(clim2)/nrow(rbind(clim1,clim2)))
  row.w.1<-1-(nrow(clim1)/nrow(rbind(clim1,clim2)))    
  row.seq<-c(rep(row.w.1,nrow(clim1)),rep(row.w.2,nrow(clim2)))
}

if(PROJ == T){  
  row.seq<-c(rep(1,nrow(clim1)),rep(0,nrow(clim2)))
}

pc<- dudi.pca(rbind(clim1,clim2)[Xvar],row.w=row.seq,scannf=F,nf=nvar)
enfa<-enfa(pc,c(pa1$pa,pa2$pa),scannf=F,nf=nvar)

            # predict the scores on the axes
scores.clim12<- enfa$li[,1:2]
scores.clim2<- enfa$li[(1+length(pa1$pa)):(nrow(clim2)+length(pa1$pa)),1:2]
scores.clim1<- enfa$li[1:nrow(clim1),1:2]
scores.sp2<- enfa$li[which(pa2$pa==1)+length(pa1),1:2]
scores.sp1<- enfa$li[which(pa1$pa==1),1:2]

			# calculation of occurence density and test of niche equivalency and similarity 
z1<- grid.clim(scores.clim12,scores.clim1,scores.sp1,R)
z2<- grid.clim(scores.clim12,scores.clim2,scores.sp2,R)
a<-niche.equivalency.test(z1,z2,rep=1)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
plot.niche(z1,title="ENFA - EU niche",name.axis1="marginality",name.axis2="tolerance")
plot.niche(z2,title="ENFA - NA niche",name.axis1="marginality",name.axis2="tolerance")
plot.contrib(enfa$co,enfa$eig)
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")


#################################################################################################
#################################### MDS ########################################################
#################################################################################################


      # measures niche overlap along the two first axes of a Multi Dimensional Scaling analysis calibrated on occurence data
      
if(PROJ == F){	
data<-rbind(occ.sp1[Xvar],occ.sp2[Xvar])
fac<-as.factor(c(rep(1, nrow(occ.sp1)),rep(2, nrow(occ.sp2))))
mds<-cmdscale(dist(data),ncol(data),eig = T)
}

if(PROJ == T){	
data<-occ.sp1[Xvar]
mds<-cmdscale(dist(data),ncol(data),eig = T)
}

# GAMs are used here to interpolate the relationship found in the mds model to new data because no predict function is not available!!

gam1<-gam(mds$point[,1]~s(aetpet,3)+s(gdd,3)+s(p,3)+s(pet,3)+s(ppi,3)+s(stdp,3)+s(t,3)+s(tmax,3)+s(tmin,3),data=data)
gam2<-gam(mds$point[,2]~s(aetpet,3)+s(gdd,3)+s(p,3)+s(pet,3)+s(ppi,3)+s(stdp,3)+s(t,3)+s(tmax,3)+s(tmin,3),data=data)

scores.clim12<-data.frame(unlist(cbind(predict(gam1,newdata=data.frame(clim12[,Xvar]),type="response"),predict(gam2,type="response",newdata=data.frame(clim12[,Xvar])))))
scores.clim1<-data.frame(unlist(cbind(predict(gam1,type="response",newdata=clim1[Xvar]),predict(gam2,type="response",newdata=clim1[Xvar]))))
scores.clim2<-data.frame(unlist(cbind(predict(gam1,type="response",newdata=clim2[Xvar]),predict(gam2,type="response",newdata=clim2[Xvar]))))
scores.sp1<- data.frame(unlist(cbind(predict(gam1,type="response",newdata=occ.sp1[Xvar]),predict(gam2,type="response",newdata=occ.sp1[Xvar]))))
scores.sp2<- data.frame(unlist(cbind(predict(gam1,type="response",newdata=occ.sp2[Xvar]),predict(gam2,type="response",newdata=occ.sp2[Xvar]))))

			# calculation of occurence density and test of niche equivalency and similarity 
z1<- grid.clim(scores.clim12,scores.clim1,scores.sp1,R)
z2<- grid.clim(scores.clim12,scores.clim2,scores.sp2,R)
a<-niche.equivalency.test(z1,z2,rep=100)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
plot.niche(z1,title="MDS - EU niche",name.axis1="PC1",name.axis2="PC2")
plot.niche(z2,title="MDS - NA niche",name.axis1="PC1",name.axis2="PC2")
contrib<-t(cor(mds$points,data))[,1:2]
plot.contrib(contrib,mds$eig)
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")

#################################################################################################
############ SDM (using BIOMOD for GLM,GBM and RF and using dismo for MaxEnt) ###################
#################################################################################################

      # measures niche overlap along the response gradient (probabilities) of a set of SDMs (i.e., GLM, GBM, RF and MaxEnt)
      
if(PROJ == F){data<-rbind(pa1,pa2)}
if(PROJ == T){data<-pa1}

# convert climate data to raster objects for MaxEnt modeling  
xyz <- data[,c(1:3)]
rasts <- rasterFromXYZ(xyz)
for(r in 4:(ncol(data)-1)){
	xyz <- data[,c(1,2,r)]
	nrast <- rasterFromXYZ(xyz)
	rasts <- stack(rasts, nrast)
}
rasts@layernames <- names(data)[3:(ncol(data)-1)]

#BIOMOD modeling
Response<-data.frame(data[,ncol(data)]);names(Response)<-"sp"
nbpresence<-sum(data["pa"]==1)
Initial.State(Response = Response, Explanatory = data[,Xvar])	
Models(GLM=T,TypeGLM ="quad",GBM=T,GAM=T,RF=T, 
   NbRunEval = 1, DataSplit =80, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=F, TSS=F, VarImport=10,
   NbRepPA=1, strategy="random", nb.absences=nbpresence)

Projection(Proj = clim12[,Xvar],Proj.name="clim12",GLM=T,GBM=T,GAM=F,RF=T,ANN=F,CTA=F,MARS=F,FDA=F,SRE=F)
load("proj.clim12/Proj_clim12_sp")
scores.clim12.GLM<-data.frame(Proj_clim12_sp[,"GLM",1,])/1000
scores.clim12.GBM<-data.frame(Proj_clim12_sp[,"GBM",1,])/1000
scores.clim12.RF<-data.frame(Proj_clim12_sp[,"RF",1,])/1000

Projection(Proj = clim1[,Xvar],Proj.name="clim1",GLM=T,GBM=T,GAM=F,RF=T,ANN=F,CTA=F,MARS=F,FDA=F,SRE=F)
load("proj.clim1/Proj_clim1_sp")
scores.clim1.GLM<-data.frame(Proj_clim1_sp[1:100,"GLM",1,])/1000
scores.clim1.GBM<-data.frame(Proj_clim1_sp[1:100,"GBM",1,])/1000
scores.clim1.RF<-data.frame(Proj_clim1_sp[1:100,"RF",1,])/1000

Projection(Proj = clim2[,Xvar],Proj.name="clim2",GLM=T,GBM=T,GAM=F,RF=T,ANN=F,CTA=F,MARS=F,FDA=F,SRE=F)
load("proj.clim2/Proj_clim2_sp")
scores.clim2.GLM<-data.frame(Proj_clim2_sp[1:100,"GLM",1,])/1000
scores.clim2.GBM<-data.frame(Proj_clim2_sp[1:100,"GBM",1,])/1000
scores.clim2.RF<-data.frame(Proj_clim2_sp[1:100,"RF",1,])/1000

Projection(Proj = occ.sp1[,Xvar],Proj.name="sp1",GLM=T,GBM=T,GAM=F,RF=T,ANN=F,CTA=F,MARS=F,FDA=F,SRE=F)
load("proj.sp1/Proj_sp1_sp")
scores.sp1.GLM<-data.frame(Proj_sp1_sp[1:100,"GLM",1,])/1000
scores.sp1.GBM<-data.frame(Proj_sp1_sp[1:100,"GBM",1,])/1000
scores.sp1.RF<-data.frame(Proj_sp1_sp[1:100,"RF",1,])/1000

Projection(Proj = occ.sp2[,Xvar],Proj.name="sp2",GLM=T,GBM=T,GAM=F,RF=T,ANN=F,CTA=F,MARS=F,FDA=F,SRE=F)
load("proj.sp2/Proj_sp2_sp")
scores.sp2.GLM<-data.frame(Proj_sp2_sp[1:100,"GLM",1,])/1000
scores.sp2.GBM<-data.frame(Proj_sp2_sp[1:100,"GBM",1,])/1000
scores.sp2.RF<-data.frame(Proj_sp2_sp[1:100,"RF",1,])/1000

#MaxEnt modeling
occur <- data[which(data$pa==1),1:2]
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='') #check if maxent.jar is located in the right folder
me <- maxent(rasts, occur, args=c("-X", 20)) # fit maxent model, with 20% reserved for testing

scores.clim12.MAXENT <- data.frame(predict(me, clim12[,Xvar], progress="text"))
scores.clim1.MAXENT <- data.frame(predict(me, clim1[,Xvar], progress="text"))
scores.clim2.MAXENT <- data.frame(predict(me, clim2[,Xvar], progress="text"))
scores.sp1.MAXENT <- data.frame(predict(me, occ.sp1[,Xvar], progress="text"))
scores.sp2.MAXENT <- data.frame(predict(me, occ.sp2[,Xvar], progress="text"))

############### measures and plot overlap along prediction gradient

#GLM
z1<- grid.clim(scores.clim12.GLM,scores.clim1.GLM,scores.sp1.GLM,R)
z2<- grid.clim(scores.clim12.GLM,scores.clim2.GLM,scores.sp2.GLM,R)
a<-niche.equivalency.test(z1,z2,rep=100)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
plot.niche(z1,title="GLM - EU niche",name.axis1="Probability of occurence",name.axis2="PC2")
plot.niche(z2,title="GLM - NA niche",name.axis1="Probability of occurence",name.axis2="PC2")
plot.contrib(t(VarImportance$sp[3,]))
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")

#GBM
z1<- grid.clim(scores.clim12.GBM,scores.clim1.GBM,scores.sp1.GBM,R)
z2<- grid.clim(scores.clim12.GBM,scores.clim2.GBM,scores.sp2.GBM,R)
a<-niche.equivalency.test(z1,z2,rep=100)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
plot.niche(z1,title="GBM - EU niche",name.axis1="Probability of occurence",name.axis2="PC2")
plot.niche(z2,title="GBM - NA niche",name.axis1="Probability of occurence",name.axis2="PC2")
plot.contrib(t(VarImportance$sp[2,]))
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")

#RF
z1<- grid.clim(scores.clim12.RF,scores.clim1.RF,scores.sp1.RF,R)
z2<- grid.clim(scores.clim12.RF,scores.clim2.RF,scores.sp2.RF,R)
a<-niche.equivalency.test(z1,z2,rep=100)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
plot.niche(z1,title="RF - EU niche",name.axis1="Probability of occurence",name.axis2="PC2")
plot.niche(z2,title="RF - NA niche",name.axis1="Probability of occurence",name.axis2="PC2")
plot.contrib(t(VarImportance$sp[4,]))
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")

#MaxEnt
z1<- grid.clim(scores.clim12.MAXENT,scores.clim1.MAXENT,scores.sp1.MAXENT,R)
z2<- grid.clim(scores.clim12.MAXENT,scores.clim2.MAXENT,scores.sp2.MAXENT,R)
a<-niche.equivalency.test(z1,z2,rep=100)# test of niche equivalency and similarity according to Warren et al. 2008
b<-niche.similarity.test(z1,z2,rep=100)
b2<-niche.similarity.test(z2,z1,rep=100)

			#plot			
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
plot.niche(z1,title="Maxent - EU niche",name.axis1="Probability of occurence",name.axis2="PC2")
plot.niche(z2,title="Maxent - NA niche",name.axis1="Probability of occurence",name.axis2="PC2")
plot(me)
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)))
plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")
