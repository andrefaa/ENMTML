
##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
## remove occurences in a dataframe that are closer to each other than a specified distance threshold
##
##ARGUMENTS
##df: dataframe with x, y, and variables
##colxy: the range of columns for x and y in df
##colvar: the range of columns for variables in df	
##min.dist: minimun distance threshold in the sub-dataframe
##

occ.desaggragation <-function(df,colxy,colvar=NULL,min.dist,plot=T){

initial<-df
train<-initial
xx<-colxy[1]
yy<-colxy[2]
kept<-0 ;out<-0; keep<-c()
x11(2,2,pointsize = 12); par(mar=c(0,0,0,0)); plot.new()

while(nrow(train)>0){
		
	i<-sample(1:nrow(train),1)
		
	if(sum(sqrt((train[,xx]-train[i,xx])^2 + (train[,yy]-train[i,yy])^2)<=min.dist)>1) {
		out<-out+1
		plot.new(); text(0.5,0.8,paste("# initial:",nrow(initial))); text(0.5,0.5,paste("# kept: ",kept)); text(0.5,0.2,paste("# out: ",out))
	}
	else {
		keep<-c(keep,row.names(train[i,]))
		kept<-kept+1
		plot.new(); text(0.5,0.8,paste("# initial:",nrow(initial))); text(0.5,0.5,paste("# kept: ",kept)); text(0.5,0.2,paste("# out: ",out))
	}

	train<-train[-i,]
}
keep.row<-rep(F,nrow(initial))

for(k in 1:nrow(initial)){
	if( sum(row.names(initial)[k]==keep)==1) keep.row[k]<-T
}
dev.off()

if(is.null(colvar))final<-initial[keep.row,colxy]
if(ncol(df)==2)final<-initial[keep.row,colxy]
if(!is.null(colvar)&ncol(df)>2)final<-initial[keep.row,c(colxy,colvar)]

if(plot==T){
	x11()
	plot(initial[,colxy],main="distribution of occurences",sub=paste("# initial (black):",nrow(initial)," | # kept (red): ",kept),pch=19,col="black",cex=0.2)
	points(final[,colxy],pch=19,col="red",cex=0.2)
}
return(final)
}

##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
##
## add environmental values to a species dataframe.
## the xy (lat/long) coordinates of the species occurrences are compared to those of the environment dataframe
## and the value of the closest pixel is added to the species dataframe. 
## when the closest environment pixel is more distant than resolution, NA is added instead of the value.
## (similar to sample() in ArcGIS)

##ARGUMENTS
##dfsp: species dataframe with x, y and optional other variables
##colspxy: the range of columns for x and y in dfsp
##colspkept: the columns of dfsp that should be kept in the final dataframe (by default: xy )
##dfvar: environmental dataframe with x, y and environmental variables
##colvarxy: the range of columns for x and y in dfvar
##colvar: the range of enviromental variables columns in dfvar. (by default: all exept xy )
##resolution: distance between x,y of species and environmental datafreme after which values shouldn't be added 
##(typically, the resolution of the data in dfvar)

sample.sp.globvar <-function(dfsp,colspxy,colspkept="xy",dfvar,colvarxy,colvar="all",resolution){

if(sum(colspkept=="xy")==1)colspkept<-colspxy
if(sum(colvar=="all")==1) {
	if(!is.null(colspkept)) colvar<-(1:ncol(dfvar))[-colvarxy]
	if(is.null(colspkept))	colvar<-(1:ncol(dfvar))
}
colspx<-colspxy[1];colspy<-colspxy[2];colvarx<-colvarxy[1];colvary<-colvarxy[2]

x<-dfsp[,colspx]
X<-dfvar[,colvarx]
y<-dfsp[,colspy]
Y<-dfvar[,colvary]

train<-data.frame(matrix(nrow=nrow(dfsp),ncol=length(colvar)))
names(train)<-names(dfvar)[colvar]

x11(2,2,pointsize = 12); par(mar=c(0,0,0,0));
for (i in 1:nrow(dfsp)){
	dist<-sqrt((X-x[i])^2 + (Y-y[i])^2)
	min<-min(dist)
	if(min<=resolution){
		if(length(colvar)>1)train[i,]<-dfvar[dist==min,colvar][1,]
		if(length(colvar)==1) train[i,]<-dfvar[dist==min,colvar][1]
	}
	plot.new(); text(0.5,0.5,paste(paste("sampling:","\n","runs to go: ",nrow(dfsp)-i))); 
}
dev.off()

if(!is.null(colspkept))final<-cbind(dfsp[,colspkept],train)
if(is.null(colspkept))final<-train

return(final)
}

##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
##Investigate spatial autocorrelation by drawing a mantel Correlogram (autocorrelation vs distance)
##
##ARGUMENTS
##df: dataframe with x, y, and variables
##colxy: the range of columns for x and y in df
##colvar: the range of columns for variables in df
##n: number of random occurences used for the test (computation time increase tremendiously when using more than 500occ.) 	
##max: maximum distance to be computed in the correlogram
##nclass: number of class of distance to be computed in the correlogram
##nperm: number of permutation in the randomization process


mantel.correlogram <- function(df,colxy,n,colvar,max,nclass,nperm){

library(ecodist)

envnorm<-data.frame(t((t(df[,colvar])-mean(df[,colvar]))/sd(df[,colvar])))
row.rand<-sample(1:nrow(df),n,replace=T)
envdist<-dist(envnorm[row.rand,])
geodist<-dist(df[row.rand,colxy])
b<- seq(from = min(geodist), to = max, length.out = nclass)
crlg<-mgram(envdist,geodist,breaks=b,nperm=nperm)
plot(crlg)
abline(h=0)
}

##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
##randomly sample pseudoabsences from an environmental dataframe covering the study area
##A minimum distance from presences can be set.
##ARGUMENTS
##nbabsences: number of pseudoabsences desired 
##glob: environmental dataframe covering the study area to sample, with x,y 
##colxyglob: the range of columns for x and y in glob
##colvar: the range of columns for x and y in glob. colvar="all" keeps all the variables in glob in the final dataframe. colvar=NULL keeps only x and y
##presence: occurence dataframe 
##colxypresence: the range of columns for x and y in presence
##mindist: minimum distance from prensences closer to wich pseudoabsences shouldn't be drawn (buffer distance around presences)

rand.pseudoabsences<-function(nbabsences, glob, colxyglob,colvar="all", presence, colxypresence, mindist){

colxglob<-colxyglob[1]
colyglob<-colxyglob[2]
colxpresence<-colxypresence[1]
colypresence<-colxypresence[2]

keep<-c()

no.i<-1
while(no.i <= nbabsences){
	ki<-sample(1:nrow(glob),1)
	if(sum(((glob[ki,colxglob]- presence[,colxpresence])^2 + (glob[ki,colyglob]- presence[,colypresence])^2) <= mindist^2)==0) {
		keep[no.i]<-ki
 		no.i<-no.i+1
	}
}
if(sum(colvar=="all")==1) colvar<-(1:ncol(glob))[-colvarxy]
if(!is.null(colvar))pseudoabs<-glob[keep,c(colxyglob,colvar)]
if(is.null(colvar))pseudoabs<-glob[keep,colxyglob]

return(pseudoabs)
}
