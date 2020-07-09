OccsThin <- function(occ,
                     envT,
                     ThinMethod,
                     VarColin,
                     DirR,
                     pred_dir,
                     distance = NULL) {
  #Function to thin occurrence data for ENM_TMLA
  #Parameters:
  #occ: Species list of occurrence data
  #envT: Predictors
  #ThinMethod: Methods chosen by user to thin occurrences
  #VarColin: Method chosen to deal with Variables Colinearity
  #DirR: Directory to save TXT with thinned occuurences
  
  #Convert from decimals to km
  spN <- names(occ)
  occDF <-
    lapply(occ, function(x)
      cbind(lat2grd(x[, 1:2]), x[, 4]))
  
  if (ThinMethod == "MORAN") {
    #1.Defined by variogram----
    #Check if there is a PC
    if (!is.null(VarColin)) {
      if (VarColin != "PCA" && names(envT)[1] != "PC1") {
        pc1 <- PCA_env_TMLA(env = envT, Dir = pred_dir)[[1]]
      } else{
        pc1 <- envT[[1]]
      }
    }else{
      if (names(envT)[1] != "PC1") {
        pc1 <- PCA_env_TMLA(env = envT, Dir = pred_dir)[[1]]
      }
    }
    
    #Optimal distance for each species
    # ocsD <- lapply(occDF, function(x)
    #   dist(x[, 1:2]))
    # maxD <- lapply(ocsD, function(x)
    #   max(x))
    # breaksD <- lapply(maxD, function(x)
    #   seq(0, x, l = 10))
    options(warn = -1)
    v1 <- vector("list", length = length(occDF))
    names(v1) <- spN
    for (i in 1:length(occDF)) {
      v1[[i]] <- tryCatch({ 
        v1[[i]] <- pgirmess::correlog(coords=occDF[[i]][, 1:2], z=occDF[[i]][, 3], method="Moran", nbclass=10)
      }, error = function(err) { 
        v1[[i]] <- pgirmess::correlog(coords=occDF[[i]][, 1:2], z=occDF[[i]][, 3], method="Moran", nbclass=NULL)
      }, error = function(e){
        v1[[i]] <- NA
      })
      
      #Select valid distances
      if(!is.na(v1[[i]])){
        if (any(abs(v1[[i]][,2])<=0.2)){
          v1[[i]] <- subset(v1[[i]],abs(v1[[i]][,2])<=0.2)
          v1[[i]] <- t(data.frame(v1[[i]][which(abs(v1[[i]][,1])==min(abs(v1[[i]][,1]))),]))
        }else{
          v1[[i]] <- t(data.frame(v1[[i]][which(abs(v1[[i]][,2])==min(abs(v1[[i]][,2]))),]))
        }
      }else{
          v1[[i]] <- 0
      }
    }
    options(warn = 1)
      # v1[[i]] <-
      #   geoR::variog(coords = occDF[[i]][, 1:2],
      #                data = occDF[[i]][, 3],
      #                uvec = breaksD[[i]])
      # v1[[i]] <- v1[[i]]$u[which(v1[[i]]$v == min(v1[[i]]$v[-length(v1[[i]]$v)]))]
    
    #Data Frame for thining
    occDF <- plyr::ldply(occ, data.frame,.id="sp")
    Moran <- plyr::ldply(v1, data.frame,.id="sp")
    
    #Thinning
    occPOS <- vector("list", length = length(occ))
    for (i in 1:length(occ)) {
      invisible(utils::capture.output(
        occT <-
          spThin::thin(
            loc.data =  occDF[occDF$sp == spN[i],],
            lat.col = "y",
            long.col = "x",
            spec.col = "sp",
            thin.par = v1[[i]][1],
            reps = 20,
            write.files = F,
            locs.thinned.list.return = T,
            write.log.file = F
          )
      ))
      occT <-
        occT[[which(sapply(occT, function(x)
          nrow(x)) == max(sapply(occT, function(x)
            nrow(x))))[1]]]
      occPOS[[i]] <- as.integer(row.names(occT))
    }
    
    #Select Thinned Occurrences
    for (i in 1:length(occPOS)) {
      occ[[i]] <- occ[[i]][occPOS[[i]], ]
    }
    
    #Number of occurrences after Thining
    uni <-
      data.frame(Species = spN,
                 UniqueOcc = sapply(occ, function(x)
                   nrow(x)))
    utils::write.table(
      uni,
      file.path(DirR, "N_Occ_Thinned.txt"),
      sep = "\t",
      row.names = F
    )
    
    #Record Thinning Distance
    Moran <- Moran[,1:4]
    colnames(Moran) <- c("sp","Distancia_km","I_Moran","p_value")
    utils::write.table(Moran,file.path(DirR,"ThinningDistance_Moran.txt"),sep="\t",row.names=F)
    
    return(occ)
    
  } else if (ThinMethod == "USER-DEFINED") {
    #2.Defined by user----
    # cat("Select distance for thining(in km):")
    # distance <- as.integer(readLines(n=1))
    
    #Data Frame for thining
    occDF <- plyr::ldply(occDF, data.frame)
    
    #Thinning
    occPOS <- vector("list", length = length(occ))
    for (i in 1:length(occPOS)) {
      invisible(utils::capture.output(
        occT <-
          spThin::thin(
            occDF[occDF$.id == spN[i], ],
            lat.col = "y",
            long.col = "x",
            spec.col = ".id",
            thin.par = distance,
            reps = 20,
            write.files = F,
            locs.thinned.list.return = T,
            write.log.file = F
          )
      ))
      occT <-
        occT[[which(sapply(occT, function(x)
          nrow(x)) == max(sapply(occT, function(x)
            nrow(x))))[1]]]
      occPOS[[i]] <- as.integer(row.names(occT))
    }
    
    #Select Thinned Occurrences
    for (i in 1:length(occPOS)) {
      occ[[i]] <- occ[[i]][occPOS[[i]], ]
    }
    
    #Number of occurrences after Thining
    uni <-
      data.frame(Species = spN,
                 UniqueOcc = sapply(occ, function(x)
                   nrow(x)))
    utils::write.table(
      uni,
      file.path(DirR, "N_Occ_Thinned.txt"),
      sep = "\t",
      row.names = F
    )
    
    return(occ)
    
  } else if (ThinMethod == "CELLSIZE") {
    #3.Based on cellsize----
    #Haversine Transformation
    distance <-
      raster::xyFromCell(envT[[1]], 1:raster::ncell(envT[[1]]))
    df <-
      data.frame(x = c(distance[1, c(2, 1)]), y = c(distance[2, c(2, 1)]))
    distance <- pracma::haversine(df$x, df$y) * 2
    
    #Data Frame for thining
    occDF <- plyr::ldply(occDF, data.frame)
    
    #Thinning
    occPOS <- vector("list", length = length(occ))
    for (i in 1:length(occPOS)) {
      invisible(utils::capture.output(
        occT <-
          spThin::thin(
            occDF[occDF$.id == spN[i],],
            lat.col = "y",
            long.col = "x",
            spec.col = ".id",
            thin.par = distance,
            reps = 20,
            write.files = F,
            locs.thinned.list.return = T,
            write.log.file = F
          )
      ))
      occT <-
        occT[[which(sapply(occT, function(x)
          nrow(x)) == max(sapply(occT, function(x)
            nrow(x))))[1]]]
      occPOS[[i]] <- as.integer(row.names(occT))
    }
    
    #Select Thinned Occurrences
    for (i in 1:length(occPOS)) {
      occ[[i]] <- occ[[i]][occPOS[[i]], ]
    }
    
    #Number of occurrences after Thining
    uni <-
      data.frame(Species = spN,
                 UniqueOcc = sapply(occ, function(x)
                   nrow(x)))
    utils::write.table(
      uni,
      file.path(DirR, "N_Occ_Thinned.txt"),
      sep = "\t",
      row.names = F
    )
    
    return(occ)
  }
}
