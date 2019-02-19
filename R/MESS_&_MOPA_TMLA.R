MESS_and_MOP <- function(Variables,
                         RecordsData, 
                         VarCol,
                         ProjDir,
                         Mehtods = c("MESS", "MOP")) {
  #Function to calculate MESS and MOPA metrics
  #Parameters:
  #Variables: Raster stack in which the model will be projected (other region/time period)
  #SpOccE: Species data with environmental information, created by ENM_TMLA
  #VarCol: Name of columns with the predictors information in species data
  #FoldF: Directory of envF variables (same structure as ENM_TMLA)
  #Background: Is the algorithm background-based
  
  # MOP function
  mop <- function(g_raster, m_occ, percent = 10) {
    m0 <- na.omit(rasterToPoints(g_raster))
    m2 <- g <- m0[,-c(1:2)] # G region
    m0 <- m0[, c(1:2)]
    m1 <-
      m <- m_occ # data (pres-abs-backgrou) from M region
    set <- c(seq(1, nrow(m2), round(nrow(m2) / 100)), nrow(m2) + 1)
    mop1 <- lapply(seq_len((length(set) - 1)), function(x) {
      seq_rdist <- set[x]:(set[x + 1] - 1)
      auclidean <- flexclust::dist2(m2[seq_rdist, ], m1)
      meanq <- apply(auclidean, 1, function(x) {
        d <- quantile(x, probs = percent / 100,
                      na.rm = TRUE)
        return(mean(x[which(x <= d)]))
      })
      return(meanq)
    })
    mop_vals <- unlist(mop1)
    mop_vals <- 1 - (mop_vals / max(mop_vals))
    mopr <- data.frame(m0, mop_vals)
    suppressWarnings({
      gridded(mopr) <- ~ x + y
    })
    mopr <- raster(mopr)
    return(mopr)
  }
  
  if (Mehtods == "MESS") {
    spN <- unique(RecordsData$sp)
    #Initialisation
    for (i in 1:length(ProjDir)) {
      if (!file.exists(file.path(ProjDir[i], "MESS"))) {
        dir.create(file.path(ProjDir[i], "MESS"))
      }
      for (j in spN) {
        if (any(algorithm %in% c("BIO", "MAH", "DOM"))) {
          spOccS <- RecordsData[RecordsData$sp == j, ]
          MESS <- mess(Variables[[i]], spOccS[spOccS$PresAbse == 1, VarCol])
          MESS <- stand_mess(MESS)
          writeRaster(
            MESS,
            file.path(ProjDir[i], "MESS", paste0(j, "_MESS_Presence.tif")),
            format = "GTiff",
            NAflag = -9999
          )
        }
        if (any(algorithm %in% c("GLM", "GAM", "SVM", "BRT", "RDF", "MLK", "GAU"))) {
          spOccS <- RecordsData[RecordsData$sp == j, ]
          MESS <- mess(Variables[[i]], spOccS[VarCol])
          MESS <- stand_mess(MESS)
          writeRaster(
            MESS,
            file.path(ProjDir[i], "MESS", paste0(j, "_MESS_PresAbse.tif")),
            format = "GTiff",
            NAflag = -9999
          )
        }
        if (any(algorithm %in% c("MXS", "MXD", "ENF"))) {
          spOccS <- RecordsDataM[RecordsDataM$sp == j, ]
          MESS <- mess(Variables[[i]], RecordsDataM[VarCol])
          MESS <- stand_mess(MESS)
          writeRaster(
            MESS,
            file.path(ProjDir[i], "MESS", paste0(j, "_MESS_Background.tif")),
            format = "GTiff",
            NAflag = -9999
          )
        }
      }
    }
  }
  
  
  if (Mehtods == "MOP") {
    spN <- unique(RecordsData$sp)
    #Initialisation
    for (i in 1:length(ProjDir)) {
      if (!file.exists(file.path(ProjDir[i], "MOP"))) {
        dir.create(file.path(ProjDir[i], "MOP"))
      }
      for (j in spN) {
        if (any(algorithm %in% c("BIO", "MAH", "DOM"))) {
          spOccS <- RecordsData[RecordsData$sp == j,]
          MOP <-
            mop(Variables[[i]], spOccS[spOccS$PresAbse == 1, VarCol])
          writeRaster(
            MOP,
            file.path(ProjDir[i], "MOP", paste0(j, "_MOP_Presence.tif")),
            format = "GTiff",
            NAflag = -9999
          )
        }
        if (any(algorithm %in% c("GLM", "GAM", "SVM", "BRT", "RDF", "MLK", "GAU"))) {
          spOccS <- RecordsData[RecordsData$sp == j,]
          MOP <- mop(Variables[[i]], spOccS[VarCol])
          writeRaster(
            MOP,
            file.path(ProjDir[i], "MOP", paste0(j, "_MOP_PresAbse.tif")),
            format = "GTiff",
            NAflag = -9999
          )
        }
        if (any(algorithm %in% c("MXS", "MXD", "ENF"))) {
          spOccS <- RecordsDataM[RecordsDataM$sp == j,]
          MOP <- mop(Variables[[i]], RecordsDataM[VarCol])
          writeRaster(
            MOP,
            file.path(ProjDir[i], "MOP", paste0(j, "_MOP_Background.tif")),
            format = "GTiff",
            NAflag = -9999
          )
        }
      }
    }
  }
}
