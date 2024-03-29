# Written by Andre Andrade
utils::globalVariables("x")
BandsPartition_TMLA <- function(evnVariables = NULL,
                                RecordsData = NULL,
                                N = NULL,
                                pseudoabsencesMethod = NULL,
                                PrAbRatio = NULL,
                                DirSave = NULL,
                                DirM = NULL,
                                MRst = NULL,
                                type = NULL,
                                Geo_Buf = NULL,
                                cores = NULL) {
  # Parameters
  # evnVariables: Predictors
  # RecordsData: Occurrence List
  # N: Longitudinal(1) or Latitudinal(2) bands

  # Development
  # Start Cluster
  # if (Sys.getenv("RSTUDIO") == "1" &&
  #   !nzchar(Sys.getenv("RSTUDIO_TERM")) &&
  #   Sys.info()["sysname"] == "Darwin" &&
  #   as.numeric(gsub("[.]", "", getRversion())) >= 360) {
  #   cl <- parallel::makeCluster(cores, outfile = "", setup_strategy = "sequential")
  # } else {
  #   cl <- parallel::makeCluster(cores, outfile = "")
  # }
  # doParallel::registerDoParallel(cl)


  # Separate data by groups
  RecordsData <-
    lapply(RecordsData, function(x) {
      cbind(x, rep(0, nrow(x)))
    })
  colnames <- c("x", "y", "Seg")
  RecordsData <-
    lapply(RecordsData, function(x) {
      colnames(x) <- colnames
      return(x)
    })
  RecordsData <- lapply(RecordsData, function(x) {
    round(x, digits = 5)
  })

  grid <- seq(2, 20, 2)

  # Start species loop----
  results <- list()
  #   foreach(
  #   x = 1:length(RecordsData),
  #   .packages = c("raster", "dismo", "plyr", "rgdal"),
  #   .export = c(
  #     "Moran_for_Quadrants_Pair_TMLA", "inv_bio", "inv_geo",
  #     "KM", "OptimRandomPoints"
  #   )
  # ) %dopar% 
    for (x in  1:length(RecordsData)){

    # Prevent Auxiliary files from rgdal
    # rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")

    opt <- NULL
    # print(names(RecordsData)[x])
    RecordsData.s <- RecordsData[[x]]

    for (quad in seq(2, 20, 2)) {
      RecordsData.st <- RecordsData.s
      # print(paste(quad, "Quadrants", sep = " "))

      # Bands----
      for (segm in 1:quad) {
        axfin <-
          min(RecordsData.st[, N]) + ((max(RecordsData.st[, N]) - min(RecordsData.st[, N])) /
            quad) * segm
        axin <-
          min(RecordsData.st[, N]) + ((max(RecordsData.st[, N]) - min(RecordsData.st[, N])) /
            quad) * (segm - 1)

        RecordsData.st[RecordsData.st[, N] >= axin &
          RecordsData.st[, N] <= axfin, "Seg"] <- segm
      }

      RecordsData.st <-
        cbind(RecordsData.st, ifelse((RecordsData.st$Seg / 2) %% 1, 1, 2))
      colnames(RecordsData.st) <- c("x", "y", "Seg", "Partition")

      # Moran's I----
      pc1 <- one_layer_pca(env_layer = evnVariables)
      Moran <-
        Moran_for_Quadrants_Pair_TMLA(
          occ = RecordsData.st,
          pc1 = pc1,
          quad = quad,
          type = type
        )
      Moran <- data.frame(cbind(quad, Moran), stringsAsFactors = F)

      # Euclidean Distance----
      RecordsData_e <-
        cbind(RecordsData.st, extract(evnVariables, RecordsData.st[, 1:2]))
      RecordsData_e <-
        split(RecordsData_e, f = RecordsData_e$Partition)

      eucl <-
        flexclust::dist2(RecordsData_e[[1]][, -c(1:4)], RecordsData_e[[2]][, -c(1:4)], method = "euclidean")
      eucl <- mean(eucl, na.rm = TRUE)

      # SD of number of records per Band----
      pa <-
        RecordsData.st$Partition # Vector with presences and absences
      Sd <- stats::sd(table(pa)) /
        mean(table(pa))

      # Combine calculations in a data frame
      res.t <- data.frame(cbind(Moran, eucl, Sd), stringsAsFactors = F)
      colnames(res.t) <- c("Partition", "Moran", "Euclidean", "SD")
      opt <- rbind(opt, res.t)
    }
    names(opt) <- names(res.t)

    # SELLECTION OF THE BEST NUMBER OF BANDS----
    Opt2 <- opt
    Dup <- !duplicated(Opt2[c("Moran", "Euclidean", "SD")])
    Opt2 <- Opt2[Dup, ]

    while (nrow(Opt2) > 1) {
      # I MORAN
      if (nrow(Opt2) == 1) {
        break
      }
      Opt2 <- Opt2[which(Opt2$Moran <= summary(Opt2$Moran)[2]), ]
      if (nrow(Opt2) == 1) {
        break
      }
      # MESS
      Opt2 <- Opt2[which(Opt2$Euclidean >= summary(Opt2$Euclidean)[5]), ]
      if (nrow(Opt2) == 1) {
        break
      }
      # SD
      Opt2 <- Opt2[which(Opt2$Sd <= summary(Opt2$Sd)[2]), ]

      if (nrow(Opt2) == 2) {
        break
      }

      if (unique(Opt2$Imoran.Grid.P) &&
        unique(Opt2$Mess.Grid.P) && unique(Opt2$Sd.Grid.P)) {
        Opt2 <- Opt2[nrow(Opt2), ]
      }
    }

    if (nrow(Opt2) > 1) {
      Opt2 <- Opt2[nrow(Opt2), ]
    }

    resOpt <- data.frame(Sp = names(RecordsData)[x], Opt2)

    # Create Bands Mask

    msk <- evnVariables[[1]]
    msk[!is.na(msk[, ])] <- 0

    quad <- Opt2$Partition

    for (segm in 1:quad) {
      # axfin <- min(RecordsData.s[,N])+((max(RecordsData.s[,N])-min(RecordsData.s[,N]))/quad)*segm
      # axin <- min(RecordsData.s[,N])+((max(RecordsData.s[,N])-min(RecordsData.s[,N]))/quad)*(segm-1)
      # RecordsData.s[RecordsData.s[,N]>=axin & RecordsData.s[,N]<=axfin, "Seg"] <- segm
      if (N %in% 1) {
        ex <- raster::extent(msk)
        axfin <- ex[1] + ((ex[2] - ex[1]) / quad) * segm
        axin <- ex[1] + ((ex[2] - ex[1]) / quad) * (segm - 1)
        RecordsData.s[RecordsData.s[, N] >= axin &
          RecordsData.s[, N] <= axfin, "Seg"] <- segm
        if (((segm / 2) - round(segm / 2)) != 0) {
          y1 <- ncol(msk) - floor((axfin - msk@extent[1]) / xres(msk))
          y0 <- ncol(msk) - floor((axin - msk@extent[1]) / xres(msk))
          for (y in y1:y0) {
            msk[, y] <- 1
          }
          msk[][is.na(evnVariables[[1]][])] <- NA
        } else {
          y1 <- ncol(msk) - floor((axfin - msk@extent[1]) / xres(msk))
          y0 <- ncol(msk) - floor((axin - msk@extent[1]) / xres(msk))

          for (y in y1:y0) {
            msk[, y] <- 2
          }
          msk[][is.na(evnVariables[[1]][])] <- NA
        }
      }

      if (N %in% 2) {
        ex <- raster::extent(msk)
        axfin <- ex[3] + ((ex[4] - ex[3]) / quad) * segm
        axin <- ex[3] + ((ex[4] - ex[3]) / quad) * (segm - 1)
        RecordsData.s[RecordsData.s[, N] >= axin &
          RecordsData.s[, N] <= axfin, "Seg"] <- segm
        if (((segm / 2) - round(segm / 2)) != 0) {
          y1 <- nrow(msk) - floor((axfin - msk@extent[3]) / yres(msk))
          y0 <- nrow(msk) - floor((axin - msk@extent[3]) / yres(msk))

          for (y in y1:y0) {
            msk[y, ] <- 1
          }
          msk[][is.na(evnVariables[[1]][])] <- NA
        } else {
          y1 <- nrow(msk) - floor((axfin - msk@extent[3]) / yres(msk))
          y0 <- nrow(msk) - floor((axin - msk@extent[3]) / yres(msk))

          for (y in y1:y0) {
            msk[y, ] <- 2
          }
          msk[][is.na(evnVariables[[1]][])] <- NA
        }
      }
    }

    RecordsData.s <-
      cbind(RecordsData.s, ifelse((RecordsData.s$Seg / 2) %% 1, 2, 1))
    RecordsData.s <-
      cbind(rep(names(RecordsData)[x], nrow(RecordsData.s)), RecordsData.s)
    colnames(RecordsData.s) <- c("sp", "x", "y", "Seg", "Partition")
    RecordsData.s <- RecordsData.s[, c("sp", "x", "y", "Partition")]
    RecordsData.s <-
      cbind(RecordsData.s, rep(1, nrow(RecordsData.s)))
    colnames(RecordsData.s) <-
      c("sp", "x", "y", "Partition", "PresAbse")
    msk[msk[] == 0] <- NA # ASC with the ODD-EVEN quadrants
    writeRaster(
      msk,
      paste(DirSave, paste(names(RecordsData)[x], ".tif", sep = ""), sep = "\\"),
      format = "GTiff",
      NAflag = -9999,
      overwrite = T
    )

    ############################################################
    #                                                          #
    #              Pseudoabsences allocation-----              #####
    #                                                          #
    ############################################################

    # Pseudo-Absences with Random allocation-----
    if (pseudoabsencesMethod == "RND") {
      # Clip the mask raster to generate rando pseudoabsences
      pseudo.mask <- msk
      pseudo.mask2 <- list()

      for (i in 1:2) {
        mask3 <- pseudo.mask
        mask3[!mask3[] == i] <- NA
        pseudo.mask2[[i]] <- mask3
      }

      pseudo.mask <- brick(pseudo.mask2)
      rm(pseudo.mask2)
      # Random allocation of pseudoabsences
      absences <- list()
      for (i in 1:2) {
        set.seed(x)
        if (!is.null(MRst)) {
          SpMask <-
            raster::raster(file.path(DirM, paste0(names(
              RecordsData
            )[x], ".tif")))
          pseudo.mask[[i]] <- pseudo.mask[[i]] * SpMask
          if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(RecordsData[[x]]))) {
            warning(
              "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
            )
            stop(
              "Please try again with another restriction type or without restricting the extent"
            )
          }
        }
        # absences.0 <-
        #   dismo::randomPoints(
        #     pseudo.mask[[i]],
        #     (1 / PrAbRatio) * sum(RecordsData.s[, "Partition"] == i),
        #     ext = extent(pseudo.mask[[i]]),
        #     prob = FALSE
        #   )

        absences.0 <- OptimRandomPoints(r = pseudo.mask[[i]], n = (1 / PrAbRatio) * sum(RecordsData.s[, "Partition"] == i), p = RecordsData.s[RecordsData.s[, "Partition"] == i, 2:3])

        colnames(absences.0) <- c(x, y)
        absences[[i]] <-
          data.frame(absences.0, stringsAsFactors = F)
        rm(absences.0)
      }
      for (i in 1:length(absences)) {
        absences[[i]] <- cbind(absences[[i]], rep(i, nrow(absences[[i]])))
      }
      absences <-
        lapply(absences, function(x) {
          cbind(x, rep(0, nrow(x)))
        })
      absences <-
        lapply(absences, function(y) {
          cbind(rep(names(RecordsData)[x], nrow(y)), y)
        })
      absences <- plyr::ldply(absences, data.frame)
      colnames(absences) <- colnames(RecordsData.s)
    }

    # Pseudo-Absences allocation with Environmental constrain ----
    if (pseudoabsencesMethod == "ENV_CONST") {
      pseudo.mask <-
        inv_bio(e = evnVariables, p = RecordsData.s[, c("x", "y")])

      # Split the raster of environmental layer with grids
      pseudo.mask <- pseudo.mask * msk
      pseudo.mask2 <- list()

      for (i in 1:2) {
        mask3 <- pseudo.mask
        mask3[!msk[] == i] <- NA
        mask3[is.na(msk[])] <- NA
        pseudo.mask2[[i]] <- mask3
      }
      pseudo.mask <- raster::brick(pseudo.mask2)
      rm(pseudo.mask2)

      absences <- list()
      for (i in 1:2) {
        set.seed(x)
        if (!is.null(MRst)) {
          SpMask <-
            raster::raster(file.path(DirM, paste0(names(
              RecordsData
            )[x], ".tif")))
          pseudo.mask[[i]] <- pseudo.mask[[i]] * SpMask
          if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(RecordsData[[x]]))) {
            warning(
              "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
            )
            stop(
              "Please try again with another restriction type or without restricting the extent"
            )
          }
        }
        # absences.0 <-
        #   dismo::randomPoints(
        #     pseudo.mask[[i]],
        #     (1 / PrAbRatio) * sum(RecordsData.s[, "Partition"] == i),
        #     ext = raster::extent(pseudo.mask[[i]]),
        #     prob = FALSE
        #   )

        absences.0 <- OptimRandomPoints(r = pseudo.mask[[i]], n = (1 / PrAbRatio) * sum(RecordsData.s[, "Partition"] == i), p = RecordsData.s[RecordsData.s[, "Partition"] == i, 2:3])


        colnames(absences.0) <- c(x, y)
        absences[[i]] <-
          data.frame(absences.0, stringsAsFactors = F)
        rm(absences.0)
      }
      for (i in 1:length(absences)) {
        absences[[i]] <- cbind(absences[[i]], rep(i, nrow(absences[[i]])))
      }
      absences <-
        lapply(absences, function(x) {
          cbind(x, rep(0, nrow(x)))
        })
      absences <-
        lapply(absences, function(y) {
          cbind(rep(names(RecordsData)[x], nrow(y)), y)
        })
      absences <- plyr::ldply(absences, data.frame)
      colnames(absences) <- colnames(RecordsData.s)
    }

    # Pseudo-Absences allocation with Geographical constrain-----
    if (pseudoabsencesMethod == "GEO_CONST") {
      pseudo.mask <-
        inv_geo(e = evnVariables[[1]], p = RecordsData.s[, c("x", "y")], d = Geo_Buf)

      # Split the raster of environmental layer with grids
      pseudo.mask <- pseudo.mask * msk
      pseudo.mask2 <- list()

      for (i in 1:2) {
        mask3 <- pseudo.mask
        mask3[!msk[] == i] <- NA
        mask3[is.na(msk[])] <- NA
        pseudo.mask2[[i]] <- mask3
      }
      pseudo.mask <- raster::brick(pseudo.mask2)
      rm(pseudo.mask2)

      absences <- list()
      for (i in 1:2) {
        set.seed(x)
        if (!is.null(MRst)) {
          SpMask <-
            raster(file.path(DirM, paste0(names(
              RecordsData
            )[x], ".tif")))
          pseudo.mask[[i]] <- pseudo.mask[[i]] * SpMask
          if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(RecordsData[[x]]))) {
            warning(
              "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
            )
            stop(
              "Please try again with a smaller geographical buffer or without restricting the accessible area"
            )
          }
        }
        # absences.0 <-
        #   dismo::randomPoints(
        #     pseudo.mask[[i]],
        #     (1 / PrAbRatio) * sum(RecordsData.s[, "Partition"] == i),
        #     ext = extent(pseudo.mask[[i]]),
        #     prob = FALSE
        #   )
        absences.0 <- OptimRandomPoints(r = pseudo.mask[[i]], n = (1 / PrAbRatio) * sum(RecordsData.s[, "Partition"] == i), p = RecordsData.s[RecordsData.s[, "Partition"] == i, 2:3])

        colnames(absences.0) <- c(x, y)
        absences[[i]] <-
          data.frame(absences.0, stringsAsFactors = F)
        rm(absences.0)
      }
      for (i in 1:length(absences)) {
        absences[[i]] <- cbind(absences[[i]], rep(i, nrow(absences[[i]])))
      }
      absences <-
        lapply(absences, function(x) {
          cbind(x, rep(0, nrow(x)))
        })
      absences <-
        lapply(absences, function(y) {
          cbind(rep(names(RecordsData)[x], nrow(y)), y)
        })
      absences <- plyr::ldply(absences, data.frame)
      colnames(absences) <- colnames(RecordsData.s)
    }

    # Pseudo-Absences allocation with Environmentla and Geographical  constrain-----
    if (pseudoabsencesMethod == "GEO_ENV_CONST") {
      pseudo.mask <-
        inv_geo(e = evnVariables[[1]], p = RecordsData.s[, c("x", "y")], d = Geo_Buf)
      pseudo.mask2 <-
        inv_bio(e = evnVariables, p = RecordsData.s[, c("x", "y")])
      pseudo.mask <- pseudo.mask * pseudo.mask2

      # Split the raster of environmental layer with grids
      pseudo.mask <- pseudo.mask * msk
      pseudo.mask2 <- list()

      for (i in 1:2) {
        mask3 <- pseudo.mask
        mask3[!msk[] == i] <- NA
        mask3[is.na(msk[])] <- NA
        pseudo.mask2[[i]] <- mask3
      }
      pseudo.mask <- raster::brick(pseudo.mask2)
      rm(pseudo.mask2)

      absences <- list()
      for (i in 1:2) {
        set.seed(x)
        if (!is.null(MRst)) {
          SpMask <-
            raster(file.path(DirM, paste0(names(
              RecordsData
            )[x], ".tif")))
          pseudo.mask[[i]] <- pseudo.mask[[i]] * SpMask
          if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(RecordsData[[x]]))) {
            warning(
              "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
            )
            stop(
              "Please try again with a smaller geographical buffer or without restricting the accessible area"
            )
          }
        }
        # absences.0 <-
        #   dismo::randomPoints(
        #     pseudo.mask[[i]],
        #     (1 / PrAbRatio) * sum(RecordsData.s[, "Partition"] == i),
        #     ext = extent(pseudo.mask[[i]]),
        #     prob = FALSE
        #   )
        absences.0 <- OptimRandomPoints(r = pseudo.mask[[i]], n = (1 / PrAbRatio) * sum(RecordsData.s[, "Partition"] == i), p = RecordsData.s[RecordsData.s[, "Partition"] == i, 2:3])

        colnames(absences.0) <- c(x, y)
        absences[[i]] <-
          data.frame(absences.0, stringsAsFactors = F)
        rm(absences.0)
      }
      for (i in 1:length(absences)) {
        absences[[i]] <- cbind(absences[[i]], rep(i, nrow(absences[[i]])))
      }
      absences <-
        lapply(absences, function(x) {
          cbind(x, rep(0, nrow(x)))
        })
      absences <-
        lapply(absences, function(y) {
          cbind(rep(names(RecordsData)[x], nrow(y)), y)
        })
      absences <- plyr::ldply(absences, data.frame)
      colnames(absences) <- colnames(RecordsData.s)
    }

    # Pseudo-Absences allocation with Environmentla and Geographical and k-mean constrain-----
    if (pseudoabsencesMethod == "GEO_ENV_KM_CONST") {
      pseudo.mask <- inv_geo(
        e = evnVariables[[1]],
        p = RecordsData.s[, c("x", "y")],
        d = Geo_Buf
      )

      pseudo.mask2 <- inv_bio(
        e = evnVariables,
        p = RecordsData.s[, c("x", "y")]
      )

      pseudo.mask <- pseudo.mask * pseudo.mask2

      # Split the raster of environmental layer with grids
      pseudo.mask <- pseudo.mask * msk
      pseudo.mask2 <- list()

      for (i in 1:2) {
        mask3 <- pseudo.mask
        mask3[!msk[] == i] <- NA
        mask3[is.na(msk[])] <- NA
        pseudo.mask2[[i]] <- mask3
      }
      pseudo.mask <- raster::brick(pseudo.mask2)
      rm(pseudo.mask2)

      absences <- list()
      for (i in 1:2) {
        set.seed(x)
        if (!is.null(MRst)) {
          SpMask <-
            raster::raster(file.path(DirM, paste0(names(
              RecordsData
            )[x], ".tif")))
          pseudo.mask[[i]] <- pseudo.mask[[i]] * SpMask
          if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(RecordsData[[x]]))) {
            warning(
              "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
            )
            stop(
              "Please try again with a smaller geographical buffer or without restricting the accessible area"
            )
          }
        }

        absences.0 <- KM(
          rasterToPoints(pseudo.mask[[i]])[, -3],
          raster::mask(evnVariables, pseudo.mask[[i]]),
          (1 / PrAbRatio) * sum(RecordsData.s[, "Partition"] ==
            i)
        )

        absences[[i]] <-
          data.frame(absences.0, stringsAsFactors = F)
        rm(absences.0)
      }

      for (i in 1:length(absences)) {
        absences[[i]] <- cbind(absences[[i]], rep(i, nrow(absences[[i]])))
      }

      absences <-
        lapply(absences, function(x) {
          cbind(x, rep(0, nrow(x)))
        })
      absences <-
        lapply(absences, function(y) {
          cbind(rep(names(RecordsData)[x], nrow(y)), y)
        })
      absences <- plyr::ldply(absences, data.frame)
      colnames(absences) <- colnames(RecordsData.s)
    }

    RecordsData.s <- rbind(RecordsData.s, absences)

    # Final Data Frame Results
    out <- list(
      ResultList = RecordsData.s,
      ResOptm = resOpt
    )
    # return(out)
    results[[x]] <- out
  }

  FinalResult <- dplyr::bind_rows(lapply(results, function(x) x[[1]]))
  FinalInfoGrid <- dplyr::bind_rows(lapply(results, function(x) x[[2]]))

  utils::write.table(
    FinalInfoGrid,
    paste(DirSave, "BestPartitions.txt", sep = "\\"),
    sep = "\t",
    row.names = F
  )
  utils::write.table(
    FinalResult,
    paste(DirSave, "OccBands.txt", sep = "\\"),
    sep = "\t",
    row.names = F
  )
  # parallel::stopCluster(cl)
  return(FinalResult)
}
