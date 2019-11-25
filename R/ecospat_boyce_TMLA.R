#### boyce index

#### functions calculating Boyce index (Hirzel et al. 2006) by Blaise Petitpierre & Frank Breiner(28.06.2013)
# fit: A vector or Raster-Layer containing the predicted suitability values
#obs: A vector containing the predicted suitability values or xy-coordinates (if fit is a Raster-Layer) of the validation points (presence records)

#### internal function calculating predicted-to-expected ratio for each class-interval

boycei <- function(interval, obs, fit) {
  fit.bin <- fit
  obs.bin <- obs
  fit.bin[fit[] >= interval[1] & fit[] <= interval[2]] <- "i"
  fit.bin[fit.bin != "i"] <- 0
  obs.bin[obs[] >= interval[1] & obs[] <= interval[2]] <- "i"
  obs.bin[obs.bin != "i"] <- 0
  pi <- length(which(obs.bin == "i")) / length(obs)
  ei <- length(which(fit.bin == "i")) / length(fit.bin)
  fi <- pi / ei
  return(fi)
}

#### Calculating Boyce index as in Hirzel et al. 2006
# fit: A vector or Raster-Layer containing the predicted suitability values
# obs: A vector containing the predicted suitability values or xy-coordinates (if fit is a Raster-Layer) of the validation points (presence records)
# nclass : number of classes or vector with classes threshold. If nclass=0, Boyce index is calculated with a moving window (see next parameters)
# windows.w : width of the moving window (by default 1/10 of the suitability range)
# res : resolution of the moving window (by default 100 focals)
# PEplot : if True, plot the predicted to expected ratio along the suitability class

ecospat.boyce <-
  function(fit,
           obs,
           nclass = 0,
           window.w = "default",
           res = 100,
           PEplot = TRUE) {
    if (class(fit) == "RasterLayer") {
      #    if (class(obs) == "data.frame") {
      if (class(obs) == "data.frame" | class(obs) == "matrix") {
        obs <- raster::extract(fit, obs)
      }
      fit <- raster::getValues(fit)
      fit <- fit[!is.na(fit)]
    }

    if (window.w == "default") {
      window.w <- (max(fit) - min(fit)) / 10
    }

    interval <- c(min(fit), max(fit))
    mini <- interval[1]
    maxi <- interval[2]

    if (nclass == 0) {
      vec.mov <-
        seq(
          from = mini,
          to = maxi - window.w,
          by = (maxi - mini - window.w) / res
        )

      vec.mov[res + 1] <-
        vec.mov[res + 1] + 1  #Trick to avoid error with closed interval in R

      interval <- cbind(vec.mov, vec.mov + window.w)
    } else if (length(nclass) > 1) {
      vec.mov <- c(mini, nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    } else if (nclass > 0 & length(nclass) < 2) {
      vec.mov <- seq(from = mini,
                     to = maxi,
                     by = (maxi - mini) / nclass)
    }

    f <- apply(interval, 1, boycei, obs, fit)
    to.keep <- which(f != "NaN")  # index to keep no NaN data
    f <- f[to.keep]

    if (length(f) < 2) {
      b <- NA  #at least two points are necessary to draw a correlation
    } else {
      r <-
        c(1:length(f))[f != c(f[-1], FALSE)]  #index to remove successive duplicates
      b <-
        cor(f[r], vec.mov[to.keep][r], method = "spearman")  # calculation of the spearman correlation (i.e. Boyce index) after removing successive duplicated values
    }

    HS <-
      apply(interval, 1, sum) / 2  # mean habitat suitability in the moving window
    HS[length(HS)] <-
      HS[length(HS)] - 1  #Correction of the 'trick' to deal with closed interval
    HS <- HS[to.keep]  #exlude the NaN

    if (PEplot == TRUE) {
      plot(
        HS,
        f,
        xlab = "Habitat suitability",
        ylab = "Predicted/Expected ratio",
        col = "grey",
        cex = 0.75
      )

      points(HS[r], f[r], pch = 19, cex = 0.75)

    }

    results <- list(F.ratio = f,
                    Spearman.cor = round(b, 3),
                    HS = HS)
    return(results)
  }
