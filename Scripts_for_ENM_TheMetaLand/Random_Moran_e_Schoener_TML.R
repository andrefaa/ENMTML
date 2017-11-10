### This is the R code to conduct ENMs evaluation by the Odds-and-Evens Framework
# In the moment, the code is configured to create the splits only by its latitudinal axis
# Read carefully the Initialization and function parameters and the messages that appear on the console
# Developers: Paulo De Marco Junior & Andre Andrade

Random_Moran_e_Schoener_TML <- function(env,occ.tr,occ.tes){
  
  #Function Parameters:
  #env:environmental variables stack
  #occ: occurrence data

  #Important Observations:
  #Occurence File MUST have the columns named Species (containing species name),Long (containing longitude),Lat (containing latitude)
  
  #Development
    occ.tr<-cbind(occ.tr[,1:3],rep(1,nrow(occ.tr)))
    colnames(occ.tr)<-c("Species","Long","Lat","Subset")
    occ.tes<-cbind(occ.tes[,1:3],rep(2,nrow(occ.tes)))
    colnames(occ.tes)<-colnames(occ.tr)
    occ.ssp<-rbind(occ.tr,occ.tes)
  
  #Calculate Moran's I
    Moran<-Moran_for_Quadrants(occ.ssp,env[[1]])
    print(Moran)

    #Calculate Schoener's D
    
    pcs<-env

    clim<-rasterToPoints(pcs)
    colnames(clim)<-c("long","lat","pc1","pc2")
    
    occ.sp1<-na.exclude(cbind(occ.tr[,2:3],extract(pcs,occ.tr[,2:3])))
    colnames(occ.sp1)<-colnames(clim)
    occ.sp2<-na.exclude(cbind(occ.tes[,2:3],extract(pcs,occ.tes[,2:3])))
    colnames(occ.sp2)<-colnames(clim)
    
    if (nrow(occ.sp1)<5 || nrow(occ.sp2)<5){
      D<-NA
      
    }else{
    
    PROJ = T
    
    # selection of variables to include in the analyses
    colnames(clim)
    Xvar<-c(3:4)
    nvar<-length(Xvar)
    
    #number of interation for the tests of equivalency and similarity
    iterations<-100
    
    #resolution of the gridding of the climate space
    R=1000
    
    # if PROJ = T
    row.w.occ.PROJT<-c(rep(0, nrow(clim)),rep(1, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
    row.w.env.PROJT<-c(rep(1, nrow(clim)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
    
    # global dataset for the analysis and rows for each sub dataset
    data.env.occ<-rbind(clim,occ.sp1,occ.sp2)[Xvar]
    row.clim<-1:nrow(clim)
    row.sp1<-(nrow(clim)+1):(nrow(clim)+nrow(occ.sp1))
    row.sp2<-(nrow(clim)+nrow(occ.sp1)+1):(nrow(clim)+nrow(occ.sp1)+nrow(occ.sp2))
    
    #################################################################################################
    #################################### PCA-ENV ####################################################
    #################################################################################################
    
    # measures niche overlap along the two first axes of a PCA calibrated on all the pixels of the study areas
    
    if(PROJ == T) {   #fit of the analyse using occurences from range 1
      pca.cal <-dudi.pca(data.env.occ,row.w = row.w.occ.PROJT, center = T, scale = T, scannf = F, nf = 2)
    }
    
    # predict the scores on the axes
    scores.clim<- pca.cal$li[row.clim,]
    scores.sp1<- pca.cal$li[row.sp1,]
    scores.sp2<- pca.cal$li[row.sp2,]
    
    # calculation of occurence density and test of niche equivalency and similarity 
    z1<- grid.clim(scores.clim,scores.clim,scores.sp1,R)
    z2<- grid.clim(scores.clim,scores.clim,scores.sp2,R)
    D<-round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)
    }
    print(Moran)
    print(D)
    res<-cbind(Moran,D)
    return(res)
}