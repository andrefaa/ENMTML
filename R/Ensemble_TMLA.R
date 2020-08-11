utils::globalVariables("s")
Ensemble_TMLA <- function(DirR,
                          ValTable,
                          ThrTable,
                          PredictType,
                          RecordsData,
                          spN = spN,
                          Threshold,
                          sensV,
                          Proj,
                          cores,
                          ensemble_metric = NULL) {
  #Start Cluster----
  if (Sys.getenv("RSTUDIO") == "1" &&
      !nzchar(Sys.getenv("RSTUDIO_TERM")) &&
      Sys.info()["sysname"] == "Darwin" &&
      as.numeric(gsub('[.]', '', getRversion())) >= 360) {
        cl <- parallel::makeCluster(cores,outfile="", setup_strategy = "sequential")
      }else{
        cl <- parallel::makeCluster(cores,outfile="")
      }
  doParallel::registerDoParallel(cl)

  #Create Folders----
  DirENS <- paste(DirR, "Ensemble", sep = "/")
  dir.create(DirENS)

  ensF <- paste(DirENS, PredictType[PredictType != "N"], sep = "/")
  for (i in 1:length(ensF)) {
    dir.create(ensF[i])
  }

  #Binary ensemble directories
  ensFCat <-
    file.path(sort(rep(ensF, length(Threshold))), Threshold)
  for (i in 1:length(ensFCat)) {
    dir.create(ensFCat[i])
  }

  #Projection directories
  if (!is.null(Proj)) {
    ProjN <- list.files(file.path(DirR, "Projection"), full.names = T)
    sapply(ProjN, function(x)
      dir.create(file.path(x, "Ensemble"), recursive = T))
    ModFut <- file.path(ProjN, "Ensemble")
    for (i in 1:length(ModFut)) {
      DirENS <- file.path(ModFut[i], PredictType)
      sapply(DirENS, function(x)
        dir.create(x, recursive = T))
      #Binary ensemble projection
      DirENSFCat <-
        file.path(sort(rep(DirENS, length(Threshold))), Threshold)
      sapply(DirENSFCat, function(x)
        dir.create(x, recursive = T))
    }
    DirENSFCat <-
      file.path(sort(rep(ModFut, length(PredictType))), PredictType, Threshold)
  }

  #Load Rasters----
  Folders <- list.files(file.path(DirR, "Algorithm"), full.names = T)
  # spN <-
  #   substr(list.files(Folders[1], pattern = "\\.tif$"),
  #          1,
  #          nchar(list.files(Folders[1], pattern = "\\.tif$")) - 4)
  if (!is.null(Proj)) {
    ProjN <- as.list(ProjN)
    names(ProjN) <- list.files(file.path(DirR, "Projection"))
    ProjN <-
      lapply(ProjN, function(x)
        grep(
          list.files(x, full.names = T),
          pattern = 'Ensemble',
          invert = T,
          value = T
        ))
  }

  #Thresholds List----
  ListSummary <- as.list(PredictType)
  names(ListSummary) <- PredictType

  #Species Loop----
  result <- foreach(
    s = 1:length(spN),
    .packages = c("raster", "dismo",
                  "plyr", "RStoolbox"),
    .export = c(
      "PCA_ENS_TMLA",
      "STANDAR",
      "STANDAR_FUT",
      "Eval_Jac_Sor_TMLA",
      "Thresholds_TMLA"
    )
  ) %dopar% {
    #Species Occurrence Data & Evaluation Table----
    SpData <- RecordsData[RecordsData["sp"] == spN[s],]
    SpVal <- ValTable[ValTable["Sp"] == spN[s],]
    SpVal <-
      SpVal[SpVal$Algorithm %in% c("BIO",
                                   "MAH",
                                   "DOM",
                                   "ENF",
                                   "GLM",
                                   "GAM",
                                   "SVM",
                                   "BRT",
                                   "RDF",
                                   "MXS",
                                   "MXD",
                                   "MLK",
                                   "GAU"),]
    SpThr <- ThrTable[ThrTable["Sp"] == spN[s],]

    #Raster Stack----
    List <- file.path(Folders, paste0(spN[s], ".tif"))
    List <- List[file.exists(List)]
    ListRaster <- raster::stack(List)
    names(ListRaster) <-
      list.files(file.path(DirR, "Algorithm"))[file.exists(file.path(Folders, paste0(spN[s], ".tif")))]

    if (!is.null(Proj)) {
      ListF <-
        lapply(ProjN, function(x)
          file.path(x, paste0(spN[s], ".tif")))
      ListF <- lapply(ListF, function(x)
        x[file.exists(x)])
      ListFut <- lapply(ListF, function(x)
        raster::stack(x))
      for (i in 1:length(ListFut)) {
        FilesF <- file.exists(file.path(ProjN[[i]], paste0(spN[s], ".tif")))
        names(ListFut[[i]]) <-
          list.files(file.path(DirR, "Algorithm"))[FilesF]
      }
      names(ListFut) <- list.files(file.path(DirR, "Projection"))
    }

    ##%######################################################%##
    #                                                          #
    ####                  Ensemble Methods                  ####
    #                                                          #
    ##%######################################################%##

    # Mean Ensemble----
    if (any(PredictType == "MEAN")) {
      Final <- raster::brick(raster::stack(ListRaster))
      FinalT <- raster::calc(Final, mean)
      Final <- STANDAR(FinalT)
      PredPoint <- raster::extract(Final, SpData[, c("x", "y")])
      PredPoint <-
        data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Eval_JS <-
        Eval_Jac_Sor_TMLA(p = PredPoint[PredPoint$PresAbse == 1, 2],
                          a = PredPoint[PredPoint$PresAbse == 0, 2])

      #Final Thresholds
      Thr <- Thresholds_TMLA(Eval, Eval_JS, sensV)
      ListSummary[["MEAN"]] <-
        data.frame(Sp = spN[s], Algorithm = "MEA", Thr)
      DirMEAN <- grep(pattern = "Ensemble/MEAN",
                         x = ensF,
                         value = T)
      raster::writeRaster(
        Final,
        paste(DirMEAN, '/', spN[s], ".tif", sep = ""),
        format = 'GTiff',
        overwrite = TRUE
      )
      Thr_Alg <- Thr[Thr$THR %in% Threshold, 2]
      DirMEANCat <- grep(pattern = "/MEAN/",
                         x = ensFCat,
                         value = T)
      for (t in 1:length(Thr_Alg)) {
        raster::writeRaster(
          Final >= Thr_Alg[t],
          paste(DirMEANCat[t], '/', spN[s], ".tif", sep = ""),
          format = 'GTiff',
          overwrite = TRUE
        )

        #Future Projection
        if (!is.null(Proj)) {
          for (p in 1:length(ListFut)) {
            Final <- raster::brick(raster::stack(ListFut[[p]]))
            Final <- raster::calc(Final, mean)
            Final <- STANDAR_FUT(Final, FinalT)

            raster::writeRaster(
              Final,
              file.path(ModFut[p], "MEAN", spN[s]),
              format = 'GTiff',
              overwrite = TRUE
            )
            DirMEANFCat <-
              grep(
                pattern = paste0(names(ListFut)[p], "/Ensemble/MEAN/", Threshold),
                x = DirENSFCat,
                value = T
              )
            for (t in 1:length(DirMEANFCat)) {
              raster::writeRaster(
                Final >= Thr_Alg[t],
                file.path(DirMEANFCat, paste(spN[s], sep = "_")),
                format = 'GTiff',
                overwrite = TRUE
              )
            }
          }
        }
      }
    }

    # Weighted Mean Ensemble----
    if (any(PredictType == "W_MEAN")) {
      ThResW <- unlist(SpVal[ensemble_metric])
      Final <- raster::brick(raster::stack(ListRaster))
      Final <- raster::calc(Final, function(x)
        x * ThResW)
      FinalT <- raster::calc(Final, mean)
      Final <- STANDAR(FinalT)

      PredPoint <- raster::extract(Final, SpData[, c("x", "y")])
      PredPoint <-
        data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Eval_JS <-
        Eval_Jac_Sor_TMLA(p = PredPoint[PredPoint$PresAbse == 1, 2],
                          a = PredPoint[PredPoint$PresAbse == 0, 2])

      #Final Thresholds
      Thr <- Thresholds_TMLA(Eval, Eval_JS, sensV)
      ListSummary[["W_MEAN"]] <-
        data.frame(Sp = spN[s], Algorithm = "WMEA", Thr)

      #Save Maps
      DirW_MEAN <- grep(pattern = "Ensemble/W_MEAN",
                      x = ensF,
                      value = T)
      raster::writeRaster(
        Final,
        paste(DirW_MEAN, '/', spN[s], ".tif", sep = ""),
        format = 'GTiff',
        overwrite = TRUE
      )
      Thr_Alg <- Thr[Thr$THR %in% Threshold, 2]
      DirMEANCat <- grep(pattern = "/W_MEAN/",
                         x = ensFCat,
                         value = T)
      for (t in 1:length(Thr_Alg)) {
        raster::writeRaster(
          Final >= Thr_Alg[t],
          paste(DirMEANCat[t], '/', spN[s], ".tif", sep = ""),
          format = 'GTiff',
          overwrite = TRUE
        )
      }

      #Future Projection
      if (!is.null(Proj)) {
        for (p in 1:length(ListFut)) {
          Final <- raster::brick(raster::stack(ListFut[[p]]))
          Final <- raster::calc(Final, function(x)
            x * ThResW)
          Final <- raster::calc(Final, mean)
          Final <- STANDAR_FUT(Final, FinalT)

          raster::writeRaster(
            Final,
            file.path(ModFut[p], "W_MEAN", spN[s]),
            format = 'GTiff',
            overwrite = TRUE
          )
          DirMEANFCat <-
            grep(
              pattern = paste0(names(ListFut)[p], "/Ensemble/W_MEAN/", Threshold),
              x = DirENSFCat,
              value = T
            )
          for (t in 1:length(Thr_Alg)) {
            raster::writeRaster(
              Final >= Thr_Alg[t],
              file.path(DirMEANFCat, paste0(spN[s], ".tif")),
              format = 'GTiff',
              overwrite = TRUE
            )
          }
        }
      }
    }

    # Superior Ensemble----
    if (any(PredictType == 'SUP')) {
      Best <-
        as.character(SpVal[which(unlist(SpVal[ensemble_metric]) >= mean(unlist(SpVal[ensemble_metric]))), "Algorithm"])
      W <- names(ListRaster)[names(ListRaster) %in% Best]
      nom <- names(ListRaster)

      ListRaster2 <- raster::brick(raster::stack(ListRaster))
      names(ListRaster2) <- nom
      Final <- raster::subset(ListRaster2, subset = W)
      FinalT <- raster::calc(Final, mean)
      Final <- STANDAR(FinalT)
      PredPoint <- raster::extract(Final, SpData[, c("x", "y")])
      PredPoint <-
        data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Eval_JS <-
        Eval_Jac_Sor_TMLA(p = PredPoint[PredPoint$PresAbse == 1, 2],
                          a = PredPoint[PredPoint$PresAbse == 0, 2])

      #Final Thresholds
      Thr <- Thresholds_TMLA(Eval, Eval_JS, sensV)
      ListSummary[["SUP"]] <-
        data.frame(Sp = spN[s], Algorithm = "SUP", Thr)

      #Save Final Maps
      DirSUP <- grep(pattern = "Ensemble/SUP",
                        x = ensF,
                        value = T)
      raster::writeRaster(
        Final,
        paste(DirSUP, '/', spN[s], ".tif", sep = ""),
        format = 'GTiff',
        overwrite = TRUE
      )
      Thr_Alg <- Thr[Thr$THR %in% Threshold, 2]
      DirMEANCat <- grep(pattern = "/SUP/",
                         x = ensFCat,
                         value = T)
      for (t in 1:length(Thr_Alg)) {
        raster::writeRaster(
          Final >= Thr_Alg[t],
          paste(DirMEANCat[t], '/', spN[s], ".tif", sep = ""),
          format = 'GTiff',
          overwrite = TRUE
        )
      }

      #Projection
      if (!is.null(Proj)) {
        for (p in 1:length(ListFut)) {
          ListFut[[p]] <- raster::brick(raster::stack(ListFut[[p]]))
          names(ListFut[[p]]) <- nom
          Final <- raster::subset(ListFut[[p]], subset = W)
          Final <- raster::calc(Final, mean)
          Final <- STANDAR_FUT(Final, FinalT)

          raster::writeRaster(
            Final,
            file.path(ModFut[p], "SUP", paste(spN[s], sep = "_")),
            format = 'GTiff',
            overwrite = TRUE
          )
          DirMEANFCat <-
            grep(
              pattern = paste0(names(ListFut)[p], "/Ensemble/SUP/", Threshold),
              x = DirENSFCat,
              value = T
            )
          for (t in 1:length(Thr_Alg)) {
            raster::writeRaster(
              Final >= Thr_Alg[t],
              file.path(DirMEANFCat, paste(spN[s], sep = "_")),
              format = 'GTiff',
              overwrite = TRUE
            )
          }
        }
      }
    }

    # With PCA ------
    if (any(PredictType == 'PCA')) {
      Final <- raster::brick(raster::stack(ListRaster))
      Final <- PCA_ENS_TMLA(Final)
      PredPoint <- raster::extract(Final, SpData[, c("x", "y")])
      PredPoint <-
        data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Eval_JS <-
        Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                          PredPoint[PredPoint$PresAbse == 0, 2])
      #Final Thresholds
      Thr <- Thresholds_TMLA(Eval, Eval_JS, sensV)
      ListSummary[["PCA"]] <-
        data.frame(Sp = spN[s], Algorithm = "PCA", Thr)

      #Save Final Maps
      DirPCA <- grep(pattern = "Ensemble/PCA",
                        x = ensF,
                        value = T)
      raster::writeRaster(
        Final,
        paste(DirPCA, '/', spN[s], ".tif", sep = ""),
        format = 'GTiff',
        overwrite = TRUE
      )
      Thr_Alg <- Thr[Thr$THR %in% Threshold, 2]
      DirMEANCat <- grep(pattern = "/PCA/",
                         x = ensFCat,
                         value = T)
      for (t in 1:length(Thr_Alg)) {
        raster::writeRaster(
          Final >= Thr_Alg[t],
          paste(DirMEANCat[t], '/', spN[s], ".tif", sep = ""),
          format = 'GTiff',
          overwrite = TRUE
        )
      }

      #Future Projection
      if (!is.null(Proj)) {
        for (p in 1:length(ListFut)) {
          Final <- raster::brick(raster::stack(ListFut[[p]]))
          # Final <- ListFut[[p]]
          Final <- PCA_ENS_TMLA(Final)

          raster::writeRaster(
            Final,
            file.path(ModFut[p], "PCA", spN[s]),
            format = 'GTiff',
            overwrite = TRUE
          )
          DirMEANFCat <-
            grep(
              pattern = paste0(names(ListFut)[p], "/Ensemble/PCA/", Threshold),
              x = DirENSFCat,
              value = T
            )
          for (t in 1:length(Thr_Alg)) {
            raster::writeRaster(
              Final >= Thr_Alg[t],
              file.path(DirMEANFCat, paste(spN[s], sep = "_")),
              format = 'GTiff',
              overwrite = TRUE
            )
          }
        }
      }
    }

    # With PCA over the Mean(Superior) Ensemble----
    if (any(PredictType == 'PCA_SUP')) {
      Best <-
        SpVal[which(unlist(SpVal[ensemble_metric]) >= mean(unlist(SpVal[ensemble_metric]))), "Algorithm"]
      W <- names(ListRaster)[names(ListRaster) %in% Best]
      nom <- names(ListRaster)

      ListRaster2 <- raster::brick(raster::stack(ListRaster))
      names(ListRaster2) <- nom
      Final <- raster::subset(ListRaster2, subset = W)
      Final <- PCA_ENS_TMLA(Final)
      PredPoint <- raster::extract(Final, SpData[, c("x", "y")])
      PredPoint <-
        data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Eval_JS <-
        Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                          PredPoint[PredPoint$PresAbse == 0, 2])
      #Final Thresholds
      Thr <- Thresholds_TMLA(Eval, Eval_JS, sensV)
      ListSummary[["PCS"]] <-
        data.frame(Sp = spN[s], Algorithm = "PCS", Thr)

      #Save Final Maps
      DirPCA_SUP <- grep(pattern = "Ensemble/PCA_SUP",
                        x = ensF,
                        value = T)
      raster::writeRaster(
        Final,
        paste(DirPCA_SUP, '/', spN[s], ".tif", sep = ""),
        format = 'GTiff',
        overwrite = TRUE
      )
      Thr_Alg <- Thr[Thr$THR %in% Threshold, 2]
      DirMEANCat <- grep(pattern = "/PCA_SUP/",
                         x = ensFCat,
                         value = T)
      for (t in 1:length(Thr_Alg)) {
        raster::writeRaster(
          Final >= Thr_Alg[t],
          paste(DirMEANCat[t], '/', spN[s], ".tif", sep = ""),
          format = 'GTiff',
          overwrite = TRUE
        )
      }

      #Future Projection
      if (!is.null(Proj)) {
        for (p in 1:length(ListFut)) {
          ListFut[[p]] <- raster::brick(raster::stack(ListFut[[p]]))
          names(ListFut[[p]]) <- nom
          Final <- raster::subset(ListFut[[p]], subset = W)
          Final <- PCA_ENS_TMLA(Final)

          raster::writeRaster(
            Final,
            file.path(ModFut[p], "PCA_SUP", paste(spN[s], sep =
                                                    "_")),
            format = 'GTiff',
            overwrite = TRUE
          )
          DirMEANFCat <-
            grep(
              pattern = paste0(names(ListFut)[p], "/Ensemble/PCA_SUP/", Threshold),
              x = DirENSFCat,
              value = T
            )
          for (t in 1:length(Thr_Alg)) {
            raster::writeRaster(
              Final >= Thr_Alg[t],
              file.path(DirMEANFCat, paste(spN[s], sep = "_")),
              format = 'GTiff',
              overwrite = TRUE
            )
          }
        }
      }
    }

    #With PCA over the threshold Ensemble----
    if (any(PredictType == 'PCA_THR')) {
      Final <- raster::brick(raster::stack(ListRaster))
      ValidTHR <-
        SpThr[grepl(Threshold, SpThr[, "THR"], ignore.case = T), "THR_VALUE"]
      for (k in 1:raster::nlayers(Final)) {
        FinalSp <- Final[[k]]
        FinalSp[FinalSp < ValidTHR[k]] <- 0
        Final[[k]] <- FinalSp
      }
      Final <- PCA_ENS_TMLA(Final)
      PredPoint <- raster::extract(Final, SpData[, c("x", "y")])
      PredPoint <-
        data.frame(PresAbse = SpData[, "PresAbse"], PredPoint)
      Eval <- dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                              PredPoint[PredPoint$PresAbse == 0, 2])
      Eval_JS <-
        Eval_Jac_Sor_TMLA(PredPoint[PredPoint$PresAbse == 1, 2],
                          PredPoint[PredPoint$PresAbse == 0, 2])

      #Final Thresholds
      Thr <- Thresholds_TMLA(Eval, Eval_JS, sensV)
      ListSummary[["PCT"]] <-
        data.frame(Sp = spN[s],
                   Algorithm = "PCT",
                   Threshold = Thr)
      #Save Final Maps
      DirPCA_THR <- grep(pattern = "Ensemble/PCA_THR",
                        x = ensF,
                        value = T)
      raster::writeRaster(
        Final,
        paste(DirPCA_THR, '/', paste(spN[s], sep = "_"), ".tif", sep =
                ""),
        format = 'GTiff',
        overwrite = TRUE
      )
      Thr_Alg <- Thr[Thr$THR %in% Threshold, 2]
      DirMEANCat <- grep(pattern = "/PCA_THR/",
                         x = ensFCat,
                         value = T)
      for (t in 1:length(Thr_Alg)) {
        raster::writeRaster(
          Final >= Thr_Alg[t],
          paste(DirMEANCat[t], '/', spN[s], ".tif", sep = ""),
          format = 'GTiff',
          overwrite = TRUE
        )
      }

      #Future Projection
      if (!is.null(Proj)) {
        for (p in 1:length(ListFut)) {
          Final <- raster::brick(raster::stack(ListFut[[p]]))

          #Select only values above the Threshold
          for (k in 1:raster::nlayers(ListRaster)) {
            FinalSp <- Final[[k]]
            FinalSp[FinalSp < ValidTHR[k]] <- 0
            if (all(parallel::makeCluster(FinalSp[]) == 0)) {
              Final[Final[[k]]] <- NULL
            } else{
              Final[[k]] <- FinalSp
            }
          }
          Final <- PCA_ENS_TMLA(Final)

          raster::writeRaster(
            Final,
            file.path(ModFut[p], "PCA_THR", paste(spN[s], sep =
                                                    "_")),
            format = 'GTiff',
            overwrite = TRUE
          )
          DirMEANFCat <-
            grep(
              pattern = paste0(names(ListFut)[p], "/Ensemble/PCA_THR/", Threshold),
              x = DirENSFCat,
              value = T
            )
          for (t in 1:length(Thr_Alg)) {
            raster::writeRaster(
              Final >= Thr_Alg[t],
              file.path(DirMEANFCat, paste(spN[s], sep = "_")),
              format = 'GTiff',
              overwrite = TRUE
            )
          }
        }
      }
    }

    #Final Data Frame Results
    result <- ldply(ListSummary, data.frame, .id = NULL)
    return(result)
  }
  # Save .txt with the models performance----
  FinalSummary <- do.call("rbind", result)
  utils::write.table(
    FinalSummary,
    paste(DirR, "Thresholds_Ensemble.txt", sep = '/'),
    sep = "\t",
    col.names = T,
    row.names = F
  )
  parallel::stopCluster(cl)
}
