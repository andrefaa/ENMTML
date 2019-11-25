ENS_Posterior <- function(RecordsData,
                          Algorithm,
                          PredictType,
                          Threshold,
                          DirAlg,
                          DirSave) {
  #Ensemble directory
  #Binary ensemble directories
  ensFCat <- file.path(DirSave, Threshold)
  for (i in 1:length(ensFCat)) {
    dir.create(ensFCat[i])
    assign(paste("Dir", PredictType[PredictType != "N"][i], "Cat", sep =
                   ""), ensFCat[i])
  }

  # Lists for validation-----
  for (i in 1:length(PredictType)) {
    assign(paste("Validation_", PredictType[i], sep = ""), list())
  }

  # Create Lists to Read Algorithms MSDMs
  ListRaster <- as.list(Algorithm)
  names(ListRaster) <- Algorithm
  spN <- as.character(unique(RecordsData$sp))
  Validation_Res <- as.list(spN)

  for (i in 1:length(spN)) {
    for (k in 1:length(ListRaster)) {
      ListRaster[[k]] <-
        raster::brick(raster::raster(file.path(
          DirAlg[k], list.files(DirAlg[k], pattern = ".tif")
        )[i]))
    }

    SpData <- RecordsData[RecordsData[, "sp"] == spN[i],]
    PAtest <- SpData[SpData[, "Partition"] == 2,]

    #Evaluate MSDM Models
    Eval <- list()
    ValidALL <- list()
    for (k in 1:length(ListRaster)) {
      PredPoint <- raster::extract(ListRaster[[k]], SpData[, c("x", "y")])
      PredPoint <-
        data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval[[k]] <-
        dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                        PredPoint[PredPoint$PresAbse == 0, 2])
      Thr <- dismo::threshold(Eval[[k]])[Threshold]
      names(Thr) <- Threshold
      Validation <- SUMMRES(Eval, 1, Thr)
      Validation <-
        cbind(rep(Algorithm[k], nrow(Validation)), Validation)
      colnames(Validation)[1] <- "Algorithm"
      ValidALL[[k]] <- Validation
    }
    ValidALL <- plyr::ldply(ValidALL, data.frame)
    ValidALL <- cbind(rep(spN[i], nrow(Validation)), ValidALL)
    colnames(ValidALL)[1] <- "sp"
    Validation_Res[[i]] <- ValidALL

    # With Mean Ensemble
    if (any(PredictType == "Mean")) {
      Final <- raster::brick(ListRaster)
      names(Final) <- Algorithm
      Final <- STANDAR(Final)
      Final <- round(mean(Final), 4)

      # Threshold
      PredPoint <- raster::extract(Final, SpData[, c("x", "y")])
      PredPoint <-
        data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Thr <- dismo::threshold(Eval)[Threshold]
      names(Thr) <- Threshold
      Validation <- SUMMRES(list(Eval), N = 1, Thr)
      Validation_Mean[[i]] <-
        data.frame(sp = spN[i], Algorithm = "MEA", Validation)

      raster::writeRaster(
        Final,
        paste(DirMean, '/', spN[i], ".tif", sep = ""),
        format = 'GTiff',
        overwrite = TRUE
      )
      TYPE <- Validation_Mean[[i]]["TYPE"]
      for (h in 1:length(Thr)) {
        raster::writeRaster(
          Final >= as.numeric(Thr)[h],
          paste(
            DirMeanCat,
            '/',
            spN[i],
            "_",
            as.character(TYPE[h, ]),
            ".tif",
            sep = ""
          ),
          format = 'GTiff',
          overwrite = TRUE
        )
      }
    }

    # With Over the Mean(Superior) Ensemble
    if (any(PredictType == 'Sup')) {
      SpValidation <- ValidALL
      Nom <- NULL
      if ("no_omission" %in% Threshold) {
        Nom <- c(Nom, "LPT")
      }
      if ("spec_sens" %in% Threshold) {
        Nom <- c(Nom, "MAX")
      }
      Validation <- NULL
      for (h in 1:length(Nom)) {
        SpValidationT <- SpValidation[SpValidation$TYPE == Nom[h], ]
        SpValidationT$Algorithm <-
          as.character(SpValidationT$Algorithm)

        Best <- which(SpValidationT$TSS >= mean(SpValidationT$TSS))
        Best <- SpValidationT$Algorithm[Best]

        W <- names(ListRaster) %in% Best
        Final <- raster::brick(ListRaster[W])
        Final <- STANDAR(Final)
        Final <- round(mean(Final), 4)

        # Threshold
        PredPoint <- raster::extract(Final, SpData[, c("x", "y")])
        PredPoint <-
          data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
        Eval <-
          dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                          PredPoint[PredPoint$PresAbse == 0, 2])
        Thr <- dismo::threshold(Eval)[Threshold[h]]
        names(Thr) <- Threshold[h]
        Validation <- rbind(Validation, SUMMRES(list(Eval), N = 1, Thr))

        raster::writeRaster(
          Final,
          paste(DirSup, '/', paste(spN[i], Nom[h], sep = "_"), ".tif", sep =
                  ""),
          format = 'GTiff',
          overwrite = TRUE
        )
        raster::writeRaster(
          Final >= as.numeric(Thr),
          paste(DirSupCat, '/', spN[i], "_", Nom[h], ".tif", sep =
                  ""),
          format = 'GTiff',
          overwrite = TRUE
        )
      }
      Validation_Sup[[i]] <-
        data.frame(sp = spN[i], Algorithm = "SUP", Validation)
    }

    # With PCA ------
    if (any(PredictType == 'PCA')) {
      # Selection of best algorithms based on TSS
      Final <- raster::brick(ListRaster)

      #PCA
      Final <- PCA_ENS_TMLA(Final)

      # Threshold
      PredPoint <- raster::extract(Final, SpData[, c("x", "y")])
      PredPoint <-
        data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Thr <- dismo::threshold(Eval)[Threshold]
      names(Thr) <- Threshold
      Validation <- SUMMRES(list(Eval), N = 1, Thr)
      Validation_PCA[[i]] <-
        data.frame(sp = spN[i], Algorithm = "PCA", Validation)

      raster::writeRaster(
        Final,
        paste(DirPCA, '/', spN[i], ".tif", sep = ""),
        format = 'GTiff',
        overwrite = TRUE
      )
      TYPE <- Validation_PCA[[i]]["TYPE"]
      for (h in 1:length(Thr)) {
        raster::writeRaster(
          Final >= as.numeric(Thr)[h],
          paste(
            DirPCACat,
            '/',
            spN[i],
            "_",
            as.character(TYPE[h, ]),
            ".tif",
            sep = ""
          ),
          format = 'GTiff',
          overwrite = TRUE
        )
      }
    }

    # With PCA over the Mean(Superior) Ensemble
    if (any(PredictType == 'PCASup')) {
      # Selection of best algorithms based on TSS
      SpValidation <- ValidALL

      Nom <- NULL
      if ("no_omission" %in% Threshold) {
        Nom <- c(Nom, "LPT")
      }
      if ("spec_sens" %in% Threshold) {
        Nom <- c(Nom, "MAX")
      }

      Validation <- NULL
      for (h in 1:length(Nom)) {
        SpValidationT <- SpValidation[SpValidation$TYPE == Nom[h], ]
        SpValidationT$Algorithm <-
          as.character(SpValidationT$Algorithm)

        Best <- which(SpValidationT$TSS >= mean(SpValidationT$TSS))
        Best <- SpValidationT$Algorithm[Best]

        W <- names(ListRaster) %in% Best
        Final <- brick(ListRaster[W])

        #PCA
        Final <- PCA_ENS_TMLA(Final)

        # Threshold
        PredPoint <- raster::extract(Final, SpData[, c("x", "y")])
        PredPoint <-
          data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
        Eval <-
          dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                          PredPoint[PredPoint$PresAbse == 0, 2])
        Thr <- dismo::threshold(Eval)[Threshold[h]]
        names(Thr) <- Threshold[h]
        Validation <- rbind(Validation, SUMMRES(list(Eval), N = 1, Thr))

        raster::writeRaster(
          Final,
          paste(DirPCASup, '/', spN[i], "_", Nom[h], ".tif", sep =
                  ""),
          format = 'GTiff',
          overwrite = TRUE
        )
        raster::writeRaster(
          Final >= as.numeric(Thr),
          paste(DirPCASupCat, '/', spN[i], "_", Nom[h], ".tif", sep =
                  ""),
          format = 'GTiff',
          overwrite = TRUE
        )
      }
      Validation_PCASup[[i]] <-
        data.frame(sp = spN[i], Algorithm = "PCS", Validation)
    }

    #With PCA over the threshold Ensemble
    if (any(PredictType == 'PCAThr')) {
      SpValidation <- ValidALL
      TYPE <- NULL
      if ("no_omission" %in% Threshold) {
        TYPE <- c(TYPE, "LPT")
      }
      if ("spec_sens" %in% Threshold) {
        TYPE <- c(TYPE, "MAX")
      }

      for (n in 1:length(TYPE)) {
        Final <- raster::brick(ListRaster)
        SpValType <- SpValidation[SpValidation$TYPE == TYPE[n], ]

        #Select only values above the Threshold
        for (k in Algorithm) {
          FinalSp <- Final[[k]]
          FinalSp[FinalSp < SpValType[SpValType$Algorithm == k, "THR"]] <-
            0
          Final[[k]] <- FinalSp
        }

        #PCA
        Final <- PCA_ENS_TMLA(Final)

        # Threshold
        PredPoint <- raster::extract(Final, SpData[, c("x", "y")])
        PredPoint <-
          data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
        Eval <-
          dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                          PredPoint[PredPoint$PresAbse == 0, 2])
        Thr <- c(dismo::threshold(Eval))[Threshold[n]]
        names(Thr) <- Threshold[n]
        Validation <- SUMMRES_ENS(Eval, Thr)
        Validation_PCAThr[[i]] <-
          data.frame(sp = spN[i], Algorithm = "PCT", Validation)

        raster::writeRaster(
          Final,
          paste(
            DirPCAThr,
            '/',
            paste(spN[i], TYPE[n], sep = "_"),
            ".tif",
            sep = ""
          ),
          format = 'GTiff',
          overwrite = TRUE
        )
        raster::writeRaster(
          Final >= unlist(Thr),
          paste(DirPCAThrCat, '/', spN[i], "_", TYPE[n], ".tif", sep =
                  ""),
          format = 'GTiff',
          overwrite = TRUE
        )
      }
    }
  }

  # Save .txt with the models performance----
  Obj <- ls(pattern = 'Validation_')
  result <- list()
  for (i in 1:length(Obj)) {
    result[[i]] <- ldply(get(Obj[i]))
  }
  result <- ldply(result)
  write.table(
    result,
    paste(DirSave, "Validation_MSDM.txt", sep = '/'),
    sep = "\t",
    col.names = T,
    row.names = F
  )
}
