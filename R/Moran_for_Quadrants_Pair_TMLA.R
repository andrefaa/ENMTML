Moran_for_Quadrants_Pair_TMLA <- function(occ, pc1, quad, type = "") {
  #Moran_for_Quadrants <- function(occurrence,pc1,m){
  #Parametros:
  #Occurrence:Occurrence data already with the quadrant it belongs
  #pc1:First axis of a PCA
  #m:Percentage of points used for Moran Calculation

  #Full dataset
  occ <- cbind(occ, raster::extract(pc1, occ[, 1:2]))
  colnames(occ) <- c("x", "y", "Seg", "Partition", "Env")

  #Calculate Z values for each observation
  occ$Env <- occ$Env - mean(occ$Env)

  #Moran Correction

  sp.dists <- 1 / (as.matrix(stats::dist(cbind(occ$x, occ$y))))
  diag(sp.dists) <- 0
  sp.dists[sp.dists == Inf] <- 0

  if (type == "nearest") {
    for (y in 1:quad) {
      sp.dists[which(occ$Seg == y), ] <-
        pmax((sp.dists[which(occ$Seg == y), ] * (sp.dists[which(occ$Seg == y), ] ==
                                                   max(sp.dists[which(occ$Seg == y), which((occ$Seg == y + 1))]))),
             (sp.dists[which(occ$Seg == y), ] * (sp.dists[which(occ$Seg == y), ] ==
                                                   max(sp.dists[which(occ$Seg == y), which((occ$Seg == y - 1))]))))
    }
  }

  if (type == "all") {
    odd <- which((occ$Partition == 1))
    even <- which((occ$Partition == 2))
    sp.dists[odd, odd] <- 0
    sp.dists[even, even] <- 0
    mins <- apply(sp.dists, 2, max)
    for (i in 1:length(mins)) {
      sp.dists[, i] <- ifelse(sp.dists[, i] == mins[i], mins[i], 0)
    }
  }

  if (sum(sp.dists) == 0) {
    Moran <- NA
  } else{
    Moran <- abs(as.numeric(ape::Moran.I(occ$Env, sp.dists, scaled = T)$observed))
  }

  return(Moran)
}
