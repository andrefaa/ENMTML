MESS_and_MOP <- function(Variables,
                         RecordsData,
                         RecordsDataM,
                         algorithm,
                         VarCol,
                         DirProj,
                         Methods = c("MESS", "MOP")) {
  #Function to calculate MESS and MOPA metrics
  #Parameters:
  #Variables: Raster stack in which MESS\MOP will be calculated (from upper function)
  #RecordsData: Species occurrence data\absence (from upper function)
  #VarCol: Name of columns with the predictors information (from upper function)
  #DirProj: Output folder
  #Mehtods = Calculate which metric? (MESS\MOP)
  
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
  
  if ("MESS" %in% Methods) {
    spN <- unique(RecordsData$sp)
    #Initialisation
    for (i in 1:length(DirProj)) {
      foreach (j =1: length(spN),.packages=c("raster","dismo"))%dopar% {
        if (any(algorithm %in% c("BIO", "MAH", "DOM"))) {
          spOccS <- RecordsData[RecordsData$sp == spN[j], ]
          MESS <- mess(Variables[[i]], spOccS[spOccS$PresAbse == 1, VarCol])
          MESS[MESS==Inf] <- NA
          MESS[!is.na(MESS[])] <- 1+(na.omit(MESS[]) - max(na.omit(MESS[])))/(max(na.omit(MESS[])) - min(na.omit(MESS[])))
          # MESS[!is.na(MESS[])] <- 1 - (na.omit(MESS[]) / min(na.omit(MESS[])))
          writeRaster(
            MESS,
            file.path(DirProj[i], paste0(spN[j], "_MESS_Presence.tif")),
            format = "GTiff",
            NAflag = -9999,
            overwrite=TRUE
          )
        }
        if (any(algorithm %in% c("GLM", "GAM", "SVM", "BRT", "RDF", "GAU"))) {
          spOccS <- RecordsData[RecordsData$sp == spN[j], ]
          MESS <- mess(Variables[[i]], spOccS[VarCol])
          MESS[MESS==Inf] <- NA
          MESS[!is.na(MESS[])] <- 1+(na.omit(MESS[]) - max(na.omit(MESS[])))/(max(na.omit(MESS[])) - min(na.omit(MESS[])))
          # MESS[!is.na(MESS[])] <- 1 - (na.omit(MESS[]) / min(na.omit(MESS[])))
          writeRaster(
            MESS,
            file.path(DirProj[i], paste0(spN[j], "_MESS_PresAbse.tif")),
            format = "GTiff",
            NAflag = -9999,
            overwrite=TRUE
          )
        }
        if (any(algorithm %in% c("MXS", "MXD", "ENF", "MLK"))) {
          spOccS <- RecordsDataM[RecordsDataM$sp == spN[j], ]
          MESS <- mess(Variables[[i]], spOccS[VarCol])
          MESS[MESS==Inf] <- NA
          MESS[!is.na(MESS[])] <- 1+(na.omit(MESS[]) - max(na.omit(MESS[])))/(max(na.omit(MESS[])) - min(na.omit(MESS[])))
          # MESS[!is.na(MESS[])] <- 1 - (na.omit(MESS[]) / min(na.omit(MESS[])))
          writeRaster(
            MESS,
            file.path(DirProj[i], paste0(spN[j], "_MESS_Background.tif")),
            format = "GTiff",
            NAflag = -9999,
            overwrite=TRUE
          )
        }
      }
    }
  }
  
  
  if ("MOP" %in% Methods) {
    spN <- unique(RecordsData$sp)
    #Initialisation
    for (i in 1:length(DirProj)) {
      foreach (j =1: length(spN),.packages=c("raster","dismo"),.export="mop")%dopar% {
        if (any(algorithm %in% c("BIO", "MAH", "DOM"))) {
          spOccS <- RecordsData[RecordsData$sp == spN[j],]
          MOP <-
            mop(Variables[[i]], spOccS[spOccS$PresAbse == 1, VarCol])
          writeRaster(
            MOP,
            file.path(DirProj[i], paste0(spN[j], "_MOP_Presence.tif")),
            format = "GTiff",
            NAflag = -9999,
            overwrite=TRUE
          )
        }
        if (any(algorithm %in% c("GLM", "GAM", "SVM", "BRT", "RDF", "GAU"))) {
          spOccS <- RecordsData[RecordsData$sp == spN[j],]
          MOP <- mop(Variables[[i]], spOccS[VarCol])
          writeRaster(
            MOP,
            file.path(DirProj[i], paste0(spN[j], "_MOP_PresAbse.tif")),
            format = "GTiff",
            NAflag = -9999,
            overwrite=TRUE
          )
        }
        if (any(algorithm %in% c("MXS", "MXD", "MLK","ENF"))) {
          spOccS <- RecordsDataM[RecordsDataM$sp == spN[j],]
          MOP <- mop(Variables[[i]], spOccS[VarCol])
          writeRaster(
            MOP,
            file.path(DirProj[i], paste0(spN[j], "_MOP_Background.tif")),
            format = "GTiff",
            NAflag = -9999,
            overwrite=TRUE
          )
        }
      }
    }
  }
}
