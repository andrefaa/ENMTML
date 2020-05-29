utils::globalVariables("j")
MOP <- function(Variables,
                RecordsData,
                Algorithm,
                VarCol,
                DirProj,
                DirMask) {
  #Function to calculate MOP metrics
  #Variables: Raster stack in which MESS\MOP will be calculated (from upper function)
  #RecordsData: Species occurrence data\absence (from upper function)
  #VarCol: Name of columns with the predictors information (from upper function)
  #DirProj: Output folder
  
  mop <- function(g_raster, m_occ, percent = 10) {
    m0 <- stats::na.omit(raster::rasterToPoints(g_raster))
    m2 <- g <- m0[, -c(1:2)] # G region
    m0 <- m0[, c(1:2)]
    m1 <- m_occ # data (pres-abs-backgrou) from M region
    set <- c(seq(1, nrow(m2), round(nrow(m2) / 2000)), nrow(m2) + 1)
    mop1 <- lapply(seq_len((length(set) - 1)), function(x) {
      seq_rdist <- set[x]:(set[x + 1] - 1)
      auclidean <- flexclust::dist2(m2[seq_rdist,], m1)
      meanq <- apply(auclidean, 1, function(x) {
        d <- stats::quantile(x, probs = percent / 100,
                      na.rm = TRUE)
        return(mean(x[which(x <= d)]))
      })
      return(meanq)
    })
    mop_vals <- unlist(mop1)
    mop_vals <- 1 - (mop_vals / max(mop_vals))
    mopr <- data.frame(m0, mop_vals)
    suppressWarnings({
      sp::gridded(mopr) <- ~ x + y
    })
    mopr <- raster::raster(mopr)
    return(mopr)
  }
  
  spN <- unique(RecordsData$sp)
  
  #Initialisation
  for (i in 1:length(DirProj)) {
    foreach (
      j = 1:length(spN),
      .packages = c("raster", "dismo"),
      .export = "mop"
    ) %dopar% {
      #Check for existence
      if (!file.exists(file.path(DirProj[i], paste0(spN[j], "_MOP.tif")))) {
        spOccS <- RecordsData[RecordsData$sp == spN[j], ]
        if (is.null(DirMask)) {
          MM <- Variables[[i]][1]
        } else{
          MM <- raster::raster(file.path(DirMask, paste0(spN[j], ".tif")))
        }
        MP <- raster::rasterize(spOccS[, 2:3], MM)
        MM[!is.na(MP[])] <- NA
        if (raster::cellStats(MM, sum) < 10000) {
          occM <-
            raster::sampleRandom(MM, size = raster::cellStats(MM, sum), xy = T)[, 1:2]
          occM <- raster::extract(Variables[[i]], occM)
        } else{
          occM <- raster::sampleRandom(MM, size = 10000, xy = T)[, 1:2]
          occM <- raster::extract(Variables[[i]], occM)
        }
        rm(MP)
        rm(MM)
        occM <- rbind(spOccS[VarCol], occM[VarCol])
        occM <- stats::na.omit(occM)
        MOPr <- mop(g_raster = Variables[[i]], m_occ = occM)
        plot(MOPr < 0.8)
        raster::writeRaster(
          MOPr,
          file.path(DirProj[i], paste0(spN[j], "_MOP.tif")),
          format = "GTiff",
          NAflag = -9999,
          overwrite = TRUE
        )
      }
    }
  }
}
