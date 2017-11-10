Moran_for_Quadrants<- function(occurrence,pc1){
#Moran_for_Quadrants <- function(occurrence,pc1,m){
  #Parametros:
    #Occurrence:Occurrence data already with the quadrant it belongs
    #pc1:First axis of a PCA
    #m:Percentage of points used for Moran Calculation
  
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
    sel.odd<-which((occs[,4] %% 2 != 0))
    sel.even<-which((occs[,4] %% 2 == 0))
    sp.dists.odd<-sp.dists.inv
    sp.dists.even<-sp.dists.inv
    
    sp.dists.odd[sel.even,]<-0
    sp.dists.even[sel.odd,]<-0
    
    if (length(sel.odd)==1){
      sp.dists.odd[sel.odd,]<-0
      sp.dists.odd[,sel.odd]<-0
    }else{
      sp.dists.odd[sel.odd,sel.odd]<-0  
    }
    
    if (length(sel.even)==1){
      sp.dists.even[sel.even,]<-0
      sp.dists.even[,sel.even]<-0
    }else{
      sp.dists.even[sel.even,sel.even]<-0  
    }
    
    maximo.colE<-max.col(sp.dists.even,ties.method="first")
    maximo.colO<-max.col(sp.dists.odd,ties.method="first")
    
    for (a in 1:nrow(sp.dists.even)){
      sp.dists.even[a,]<-ifelse(sp.dists.even[a,]==sp.dists.even[a,maximo.colE[a]],sp.dists.even[a,maximo.colE[a]],0)
    }
    
    for (a in 1:nrow(sp.dists.odd)){
      sp.dists.odd[a,]<-ifelse(sp.dists.odd[a,]==sp.dists.odd[a,maximo.colO[a]],sp.dists.odd[a,maximo.colO[a]],0)
    }
    
    sp.dists.t<-pmax(sp.dists.odd,sp.dists.even)
    Moran<-as.numeric(Moran.I(occs$occs.env,sp.dists.t,scaled=T)[1])

    return(Moran)
}