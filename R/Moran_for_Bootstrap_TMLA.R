Moran_for_Bootstrap_TMLA <- function(pc1,
                                     occTR) {
  #Moran_for_Quadrants <- function(occurrence,pc1,m){
  #Parametros:
  #Occurrence:Occurrence data already with the quadrant it belongs
  #pc1:First axis of a PCA
  #m:Percentage of points used for Moran Calculation
  
  #Full dataset
  occs.xy <- occTR[, c(2, 3)]
  occs.env <- raster::extract(pc1, occs.xy)
  occ <- cbind(occTR, occs.env)
  spN <- unique(occ[, 1])
  
  #Calculate Z values for each observation
  occs.z <- occ$occs.env - mean(occ$occs.env)
  occ <- cbind(occ, occs.z)
  
  #Calculate Moran
  
  Moran <- list()
  for (i in spN) {
    occSP <- occ[occ$sp == i, ]
    
    #Moran Correction
    sp.dists <- 1 / (as.matrix(stats::dist(cbind(
      occSP$x, occSP$y
    ))))
    diag(sp.dists) <- 0
    sp.dists[sp.dists == Inf] <- 0
    
    sel.odd <- which((occSP$Partition == 1))
    sel.even <- which((occSP$Partition == 2))
    sp.dists.odd <- sp.dists
    sp.dists.even <- sp.dists
    
    sp.dists.odd[sel.even, ] <- 0
    sp.dists.even[sel.odd, ] <- 0
    
    if (length(sel.odd) == 1) {
      sp.dists.odd[sel.odd, ] <- 0
      sp.dists.odd[, sel.odd] <- 0
    } else{
      sp.dists.odd[sel.odd, sel.odd] <- 0
    }
    
    if (length(sel.even) == 1) {
      sp.dists.even[sel.even, ] <- 0
      sp.dists.even[, sel.even] <- 0
    } else{
      sp.dists.even[sel.even, sel.even] <- 0
    }
    
    maximo.colE <- max.col(sp.dists.even, ties.method = "first")
    maximo.colO <- max.col(sp.dists.odd, ties.method = "first")
    
    for (a in 1:nrow(sp.dists.even)) {
      sp.dists.even[, a] <-
        ifelse(sp.dists.even[, a] == sp.dists.even[maximo.colO[a], a], sp.dists.even[maximo.colO[a], a], 0)
    }
    
    for (a in 1:nrow(sp.dists.odd)) {
      sp.dists.odd[, a] <-
        ifelse(sp.dists.odd[, a] == sp.dists.odd[maximo.colE[a], a], sp.dists.odd[maximo.colO[a], a], 0)
    }
    
    sp.dists <-
      pmax(sp.dists.odd, sp.dists.even)
    
    Moran[[i]] <-
      as.numeric(ape::Moran.I(occSP$occs.z, sp.dists, scaled = T)[1])
  }
  Moran <- ldply(Moran, data.frame, .id = "sp")
  colnames(Moran)[2] <- "Moran's I"
  return(Moran)
}
