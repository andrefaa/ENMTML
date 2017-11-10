Moran_for_Quadrants_Pair<- function(occurrence,pc1){
#Moran_for_Quadrants <- function(occurrence,pc1,m){
  #Parametros:
    #Occurrence:Occurrence data already with the quadrant it belongs
    #pc1:First axis of a PCA
    #m:Percentage of points used for Moran Calculation
  
    library(ape)
  
    #Read the environmental asc file from where env_data will be extracted
    asc.base<-pc1

    #Full dataset
    occs<-occurrence
    occs.xy<-occs[,c(2,3)]
    occs.env<-extract(asc.base,occs.xy)
    occs<-cbind(occs,occs.env)
    
    #Calculate Z values for each observation
    occs.z<-occs$occs.env-mean(occs$occs.env)
    occs<-cbind(occs,occs.z)

    #Moran Correction

    sp.dists <- as.matrix(dist(cbind(occs$Long, occs$Lat)))
    sp.dists.inv <- 1/sp.dists
    diag(sp.dists.inv) <- 0
    sp.dists.inv[sp.dists.inv ==Inf] <- 0
    sp.dists.inv.t<-sp.dists.inv
    for (x in 1:max(occs[,4])){
       sp.dists.inv.t[which(occs[,4]==x),]<-pmax((sp.dists.inv[which(occs[,4]==x),]*(sp.dists.inv[which(occs[,4]==x),]==max(sp.dists.inv[which(occs[,4]==x),which((occs[,4]==x+1))]))),(sp.dists.inv[which(occs[,4]==x),]*(sp.dists.inv[which(occs[,4]==x),]==max(sp.dists.inv[which(occs[,4]==x),which((occs[,4]==x-1))]))))
    }
    sp.dists.t<-sp.dists.inv.t
    Moran<-as.numeric(Moran.I(occs$occs.env,sp.dists.t,scaled=T)[1])

    return(Moran)
}