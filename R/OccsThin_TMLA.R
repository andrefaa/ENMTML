OccsThin <- function(occ,
                     envT,
                     ThinMethod,
                     VarColin,
                     DirR){
  #Function to thin occurrence data for ENM_TMLA
  #Parameters:
    #occ: Species list of occurrence data
    #envT: Predictors
    #ThinMethod: Methods chosen by user to thin occurrences
    #VarColin: Method chosen to deal with Variables Colinearity
    #DirR: Directory to save TXT with thinned occuurences
  
  #Convert from decimals to km
  spN <- names(occ)
  occDF <- lapply(occ, function(x) cbind(latlong2grid(x[,1:2]),x[,4]))
  
  if(ThinMethod==1){
  #1.Defined by variogram----
    #Check if there is a PC
    if(VarColin!="PCA" && names(envT)[1]!="PC1"){
      pc1 <- PCA_env_TMLA(env=envT,Dir=pred_dir)[[1]]
    }else{
      pc1 <- envT[[1]]
    }

    #Optimal distance for each species
    ocsD <- lapply(occDF, function(x) dist(x[,1:2]))
    maxD <- lapply(ocsD, function(x) max(x))
    breaksD <- lapply(maxD, function(x) seq(0,x,l=10))
    v1 <- vector("list", length = length(breaksD))
    for(i in 1:length(breaksD)){
      v1[[i]] <- variog(coords=occDF[[i]][,1:2],data=occDF[[i]][,3],uvec=breaksD[[i]])
      v1[[i]] <- v1[[i]]$u[which(v1[[i]]$v==min(v1[[i]]$v[-1]))]
    }
    
    #Data Frame for thining
    occDF <- ldply(occDF,data.frame)
    
    #Thinning
    occPOS <- vector("list", length = length(breaksD))
    for(i in 1:length(v1)){
      occT <- thin(occDF[occDF$.id==spN[i],],lat.col = "y", long.col="x",spec.col=".id", thin.par = v1[[i]],reps=20,
                  write.files=F,locs.thinned.list.return=T,write.log.file=F)
      occT <- occT[[which(sapply(occT, function(x) nrow(x))==max(sapply(occT, function(x) nrow(x))))[1]]]
      occPOS[[i]] <- as.integer(row.names(occT))
    }
    
    #Select Thinned Occurrences
    for(i in 1:length(occPOS)){
      occ[[i]] <- occ[[i]][occPOS[[i]],]
    }
    
    #Number of occurrences after Thining
    uni <- data.frame(Species=spN,UniqueOcc=sapply(occ,function(x) nrow(x)))
    write.table(uni,file.path(DirR,"N_Occ_Thinned.txt"),sep="\t",row.names=F)
    return(occ)

  } else if (ThinMethod==2){
  #2.Defined by user----
    cat("Select distance for thining(in km):")
    D <- as.integer(readLines(n=1))
    
    #Data Frame for thining
    occDF <- ldply(occDF,data.frame)
    
    #Thinning
    occPOS <- vector("list", length = length(occ))
    for(i in 1:length(occPOS)){
      occT <- thin(occDF[occDF$.id==spN[i],],lat.col = "y", long.col="x",spec.col=".id", thin.par = D,reps=20,
                   write.files=F,locs.thinned.list.return=T,write.log.file=F)
      occT <- occT[[which(sapply(occT, function(x) nrow(x))==max(sapply(occT, function(x) nrow(x))))[1]]]
      occPOS[[i]] <- as.integer(row.names(occT))
    }
    
    #Select Thinned Occurrences
    for(i in 1:length(occPOS)){
      occ[[i]] <- occ[[i]][occPOS[[i]],]
    }
    
    #Number of occurrences after Thining
    uni <- data.frame(Species=spN,UniqueOcc=sapply(occ,function(x) nrow(x)))
    write.table(uni,file.path(DirR,"N_Occ_Thinned.txt"),sep="\t",row.names=F)
    
    return(occ)
    
  } else if (ThinMethod==3){
  #3.Based on cellsize----
    #Haversine Transformation
    D <- rasterToPoints(envT[[1]])[,-3]
    D <- haversine(D[1,], D[2,], R = 6371.0)*2
    
    #Data Frame for thining
    occDF <- ldply(occDF,data.frame)
    
    #Thinning
    occPOS <- vector("list", length = length(occ))
    for(i in 1:length(occPOS)){
      occT <- thin(occDF[occDF$.id==spN[i],],lat.col = "y", long.col="x",spec.col=".id", thin.par = D,reps=20,
                   write.files=F,locs.thinned.list.return=T,write.log.file=F)
      occT <- occT[[which(sapply(occT, function(x) nrow(x))==max(sapply(occT, function(x) nrow(x))))[1]]]
      occPOS[[i]] <- as.integer(row.names(occT))
    }
    
    #Select Thinned Occurrences
    for(i in 1:length(occPOS)){
      occ[[i]] <- occ[[i]][occPOS[[i]],]
    }
    
    #Number of occurrences after Thining
    uni <- data.frame(Species=spN,UniqueOcc=sapply(occ,function(x) nrow(x)))
    write.table(uni,file.path(DirR,"N_Occ_Thinned.txt"),sep="\t",row.names=F)
    
    return(occ)
  }
}