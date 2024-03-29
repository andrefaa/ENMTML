## Written by Santiago Velazco
utils::globalVariables("s")
BlockPartition_TMLA <- function(evnVariables = NULL,
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
  # RecordsData: matrix or data frame with presences records
  # N: 2 (dafault). interger  Number of group for data  paritioning
  # pseudoabsences: logical, TRUE (dafault).
  # mask: Raster object. Preferible on wiht the same resolution and extent than
  #       variable used in the future
  # pseudoabsencesMethod: a character string indicating which pseudo-absences
  #                       method is to be computed
  # PrAbRatio: numeric. value of PrAbRatio to be computed.
  # evnVariables: Raster object. Variable set to be used in pseusoabsences
  # cellSize: numeric vector. a vector of values with different cell grid sizes

  #Cellsize
  cellSize = seq(res(evnVariables[[1]])[1]*2, res(evnVariables[[1]])[1]*100, length.out = 30)
  
  
  # Mask
  mask <- evnVariables[[1]]
  if (class(mask) != "brick") {
    mask <- raster::brick(mask)
  }
  names(mask) <- "Group"
  mask[!is.na(mask[, ])] <- 1
  
  #Extent
  e <- raster::extent(mask)
  
  # Loop for each species-----
  SpNames <- names(RecordsData)
  
  #Start Cluster
  # if (Sys.getenv("RSTUDIO") == "1" &&
  #     !nzchar(Sys.getenv("RSTUDIO_TERM")) &&
  #     Sys.info()["sysname"] == "Darwin" &&
  #     as.numeric(gsub('[.]', '', getRversion())) >= 360) {
  #   cl <- parallel::makeCluster(cores,outfile="", setup_strategy = "sequential")
  # }else{
  #   cl <- parallel::makeCluster(cores,outfile="")
  # }
  # doParallel::registerDoParallel(cl)
  
  # LOOP----
  results <- list()
    # foreach(
    #   s = 1:length(RecordsData),
    #   .packages = c("raster", "dismo",'rgdal'),
    #   .export = c("inv_bio", "inv_geo", "KM_BLOCK", "OptimRandomPoints")
    # ) %dopar% 
      for(s in 1:length(RecordsData)){
      
      #Prevent Auxiliary files from rgdal
        suppressWarnings(rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE"))
      
      # print(paste(s, SpNames[s]))
      # Extract coordinates----
      presences <- RecordsData[[s]]
      mask2 <- mask
      mask2[] <- 0
      presences2 <- data.frame(pa = c(rep(1, nrow(presences))), presences)
      
      # Transform the presences points in a DataFrameSpatialPoints
      sp::coordinates(presences2) = presences2[, c("x", "y")]
      raster::crs(presences2) <- raster::projection(mask)
      
      #### Data partitioning using a grid approach ####
      
      # raster resolution
      DIM <-
        matrix(0, length(cellSize), 2) # the number of rows and columns of each grid
      colnames(DIM) <- c("R", "C")
      DIM
      part <- data.frame(matrix(0, nrow(presences2@data), nrow(DIM)))
      part2 <- list()
      
      for (i in 1:length(cellSize)) {
        mask3 <- mask2
        res(mask3) <- cellSize[i]
        DIM[i,] <- dim(mask3)[1:2]
        raster::values(mask3) <- 1 # Add values to cells /
        NAS <-
          c(raster::extract(mask3, presences2)) # Extract values to test if exist NAs
        if (any(is.na(NAS))) {
          while (any(is.na(NAS))) {
            raster::extent(mask3) <- raster::extent(mask3) + cellSize[i]
            raster::res(mask3) <- cellSize[i] # Give to cells a size
            DIM[i,] <- dim(mask3)[1:2]
            raster::values(mask3) <- 1
            NAS <- raster::extract(mask3, presences2)
          }
        }
        
        # In this section is assigned the group of each cell
        group <- rep(c(rep(1:N, DIM[i, 2])[1:DIM[i, 2]],
                       rep(c((N / 2 + 1):N, 1:(N / 2)
                       ), DIM[i, 2])[1:DIM[i, 2]])
                     , ncell(mask3) / DIM[i, 2])[1:ncell(mask3)]
        
        mask3[] <-
          data.frame(group, expand.grid(C = 1:DIM[i, 2], R = 1:DIM[i, 1]))$group
        
        # Matrix within each columns represent the partitions of points
        # for each grid resolution
        
        part[, i] <- raster::extract(mask3, presences2)
        part2[[i]] <- data.frame(raster::extract(mask3, presences2,), presences)
        
      }#
      
      # Here will be deleted grids that assigned partitions less than the number
      # of groups
      pp <- sapply(part[1:nrow(presences), ], function(x)
        length(unique(range(x))))
      pp <- ifelse(pp == N, TRUE, FALSE)
      # Elimination of those partition that have one record in some group
      pf <- sapply(part[1:nrow(presences), ], table)
      if (is.list(pf) == TRUE) {
        pf <- which(sapply(pf, min) == 1)
      } else{
        pf <- which(apply(pf, 2, min) == 1)
      }
      pp[pf] <- FALSE
      # grid <- grid[pp]
      part <- data.frame(part[, pp])
      names(part) <- names(which(pp == T))
      part2 <- part2[pp]
      
      # Performace of cells ----
      # SD of number of records per cell size-----
      pa <- presences2@data[, 1] # Vector wiht presences and absences
      Sd.Grid.P <- rep(NA, ncol(part))
      for (i in 1:ncol(part)) {
        Sd.Grid.P[i] <- stats::sd(table(part[pa == 1, i])) /
          mean(table(part[pa == 1, i]))
      }
      
      # Euclidean distance -----
      Eucl.Grid.P <- rep(NA, ncol(part))
      Env.P <- raster::extract(evnVariables, presences)
      for (i in 1:ncol(part)) {
        Env.P1 <- cbind(part[i], Env.P)
        Env.P2 <- split(Env.P1[, -1], Env.P1[, 1])
        Eucl.Grid.P[i] <- mean(flexclust::dist2(Env.P2[[1]], Env.P2[[2]]), method = "euclidean")
        rm(Env.P1)
      }
      
      # Imoran-----
      if(grepl("PC",names(evnVariables)[1])){
        pc1 <- evnVariables[[1]]
      } else{
        pc1 <- one_layer_pca(env_layer = evnVariables)
      }
      Imoran.Grid.P <- rep(NA, length(part2))
      for (p in 1:length(part2)) {
        part3 <- part2[[p]]
        # rownames(part3) <-
        #   paste(part3$group, part3$C, part3$R, part3$lon, part3$lat)
        if (type == "nearest") {
          mineucli <- list()
          for (r in 1:DIM[p, 'R']) {
            part4 <- part3[part3$R == r,]
            if (nrow(part4) >= 2) {
              for (c in 1:DIM[p, 'C']) {
                A <- part4[part4$C == c, ]
                B <- part4[part4$C == (c + 1), ]
                euclilist <- matrix(0, (nrow(A) * nrow(B)), 3)
                if ((nrow(A) >= 1 & nrow(B) >= 1)) {
                  coord <- data.frame(
                    expand.grid(rownames(A), rownames(B)),
                    cbind(expand.grid(A[, 5], B[, 5]),
                          expand.grid(A[, 4], B[, 4]))
                  )
                  colnames(coord) <- c("A", "B", "xa", "xb", "ya", "yb")
                  coord$eucli <-
                    sqrt((coord$xa - coord$xb) ^ 2 + (coord$ya - coord$yb) ^ 2)
                  mineucli[[length(mineucli) + 1]] <-
                    coord[which.min(coord$eucli),]
                }
              }
            }
          }
          for (r in 1:DIM[p, 'C']) {
            part4 <- part3[part3$C == r,]
            if (nrow(part4) >= 2) {
              for (c in 1:DIM[p, 'R']) {
                A <- part4[part4$R == c, ]
                B <- part4[part4$R == (c + 1), ]
                if ((nrow(A) >= 1 & nrow(B) >= 1)) {
                  coord <- data.frame(
                    expand.grid(rownames(A), rownames(B)),
                    cbind(expand.grid(A[, 5], B[, 5]),
                          expand.grid(A[, 4], B[, 4]))
                  )
                  colnames(coord) <- c("A", "B", "xa", "xb", "ya", "yb")
                  coord$eucli <-
                    sqrt((coord$xa - coord$xb) ^ 2 + (coord$ya - coord$yb) ^ 2)
                  mineucli[[length(mineucli) + 1]] <-
                    coord[which.min(coord$eucli),]
                }
              }
            }
          }
          
          euclida <- plyr::ldply(mineucli, data.frame)
          euclida$A <- as.character(euclida$A)
          euclida$B <- as.character(euclida$B)
          species2 <- presences[c(euclida$A, euclida$B), ]
          dist <- as.matrix(dist(species2))
          dist <- 1 / dist
          diag(dist) <- 0
          dist[which(dist == Inf)] <- 0
          species2$pc1 <- raster::extract(evnVariables[[1]], species2)
        }
        
        if (type == "all") {
          odd <- which((part3$Group == 1))
          even <- which((part3$Group == 2))
          dist <- as.matrix(dist(presences))
          dist <- 1 / dist
          diag(dist) <- 0
          dist[which(dist == Inf)] <- 0
          dist[odd, odd] <- 0
          dist[even, even] <- 0
          mins <- apply(dist, 2, max)
          for (i in 1:length(mins)) {
            dist[, i] <- ifelse(dist[, i] == mins[i], mins[i], 0)
          }
          species2 <-
            cbind(presences, pc1 = raster::extract(pc1, presences))
        }
        
        if (nrow(species2) < 3) {
          Imoran.Grid.P[p] <- NA
        } else{
          Imoran.Grid.P[p] <-
            Moran.I(species2$pc1,
                    dist,
                    na.rm = T,
                    scaled = T)$observed
        }
      }
      
      Imoran.Grid.P <-
        abs(Imoran.Grid.P) 
      N.grid <- 1:length(cellSize[pp])
      
      Opt <-
        data.frame(N.grid, cellSize[pp], round(data.frame(
          Imoran.Grid.P, Eucl.Grid.P, Sd.Grid.P
        ), 3))
      # Cleaning those variances based in data divided in a number of partition less than
      # the number of groups
      
      # SELLECTION OF THE BEST CELL SIZE----
      Opt2 <- Opt
      Dup <-
        !duplicated(Opt2[c("Imoran.Grid.P", "Eucl.Grid.P", "Sd.Grid.P")])
      Opt2 <- Opt2[Dup, ]
      
      while (nrow(Opt2) > 1) {
        # I MORAN
        if (nrow(Opt2) == 1)
          break
        Opt2 <-
          Opt2[which(Opt2$Imoran.Grid.P <= summary(Opt2$Imoran.Grid.P)[2]), ]
        if (nrow(Opt2) == 1)
          break
        # Euclidean distance
        Opt2 <-
          Opt2[which(Opt2$Eucl.Grid.P >= summary(Opt2$Eucl.Grid.P)[5]), ]
        if (nrow(Opt2) == 1)
          break
        # SD
        Opt2 <-
          Opt2[which(Opt2$Sd.Grid.P <= summary(Opt2$Sd.Grid.P)[2]), ]
        if (nrow(Opt2) == 2)
          break
        
        if (unique(Opt2$Imoran.Grid.P) &&
            unique(Opt2$Eucl.Grid.P) && unique(Opt2$Sd.Grid.P)) {
          Opt2 <- Opt2[nrow(Opt2), ]
        }
      }
      
      if (nrow(Opt2) > 1) {
        Opt2 <- Opt2[nrow(Opt2), ]
      }
      
      # Optimum size for presences
      # print(Opt2)
      presences <- data.frame(Partition = as.numeric(part[, Opt2$N.grid]), presences)
      
      #Save blocks raster
      mask3 <- mask2
      res(mask3) <- Opt2$cellSize.pp.
      raster::values(mask3) <- 1
      
      NAS <-
        c(raster::extract(mask3, presences2)) # Extract values to test if exist NAs
      if (any(is.na(NAS))) {
        while (any(is.na(NAS))) {
          raster::extent(mask3) <- raster::extent(mask3) + Opt2$cellSize.pp.
          raster::res(mask3) <- Opt2$cellSize.pp. # Give to cells a size
          raster::values(mask3) <- 1
          NAS <- raster::extract(mask3, presences2)
        }
      }
      
      DIM <- dim(mask3)[1:2]
      
      group <- rep(c(rep(1:N, DIM[2])[1:DIM[2]],
                     rep(c((N / 2 + 1):N, 1:(N / 2)
                     ), DIM[2])[1:DIM[2]])
                   , ncell(mask3) / DIM[2])[1:ncell(mask3)]
      mask3[] <-
        data.frame(group, expand.grid(C = 1:DIM[2], R = 1:DIM[1]))$group
      
      pseudo.mask <- mask
      pseudo.mask2 <- list()
      RtoP <- data.frame(raster::rasterToPoints(mask)[, -3])
      sp::coordinates(RtoP) = c("x", "y")
      raster::crs(RtoP) <- raster::projection(mask)
      FILTER <- raster::extract(mask3,RtoP)
      pseudo.mask[which(pseudo.mask[] == 1)] <- as.matrix(FILTER)
      # writeRaster(pseudo.mask, paste(DirSave, paste(SpNames[s],'.tif',sep=""),sep='/'),
      #             format = 'GTiff', NAflag = -9999, overwrite = TRUE)
      #
      for (i in 1:N) {
        mask3 <- pseudo.mask
        mask3[!mask3[] == i] <- 0
        pseudo.mask2[[i]] <- mask3
      }
      #
      pseudo.mask <- raster::brick(pseudo.mask2)
      # rm(pseudo.mask2)
      
      ##%######################################################%##
      #                                                          #
      ####             Pseudoabsences allocation              ####
      #                                                          #
      ##%######################################################%##
      
      
      # Pseudo-Absences with Random allocation-----
      if (pseudoabsencesMethod == "RND") {
        pseudo.mask_p <- pseudo.mask
        pseudo.mask <- sum(pseudo.mask_p)
        pseudo.mask_p[pseudo.mask_p==0] <- NA
        
        raster::writeRaster(
          pseudo.mask,
          paste(DirSave, paste(SpNames[s], '.tif', sep = ""), sep = '/'),
          format = 'GTiff',
          NAflag = -9999,
          overwrite = TRUE
        )
        
        # Random allocation of Pseudo-Absences
        absences <- list()
        for (i in 1:N) {
          set.seed(s)
          if (!is.null(MRst)) {
            SpMask <- raster::raster(file.path(DirM, paste0(SpNames[s], ".tif")))
            pseudo.mask_p[[i]] <- pseudo.mask_p[[i]] * SpMask
            # if(sum(is.na(SpMask[])==F)<(PrAbRatio*nrow(RecordsData[[s]]))){
            # warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences")
            # stop("Please try again with another restriction type or without restricting the extent")
            # }
          }
          # absences.0 <-
          #   dismo::randomPoints(
          #     pseudo.mask_p[[i]],
          #     (1 / PrAbRatio) * sum(presences[, 1] == i),
          #     ext = e,
          #     prob = FALSE
          #   )
          
          absences.0 <- OptimRandomPoints(r=pseudo.mask_p[[i]], n=(1 / PrAbRatio)*sum((presences[, 1] == i)),p=presences[presences[, 1] == i, 2:3] )
          colnames(absences.0) <- c("lon", "lat")
          absences[[i]] <- as.data.frame(absences.0)
        }
        names(absences) <- 1:N
      }
      
      # Pseudo-Absences allocation with Environmental constrain ----
      if (pseudoabsencesMethod == "ENV_CONST") {
        pseudo.mask_p <- inv_bio(evnVariables, presences[, -1])
        
        # Split the raster of environmental layer with grids
        pseudo.mask_p <- raster::mask(pseudo.mask, pseudo.mask_p)
        pseudo.mask <- sum(pseudo.mask_p)
        pseudo.mask_p[pseudo.mask_p==0] <- NA
        
        raster::writeRaster(
          pseudo.mask,
          paste(DirSave, paste(SpNames[s], '.tif', sep = ""), sep = '/'),
          format = 'GTiff',
          NAflag = -9999,
          overwrite = TRUE
        )
        
        absences <- list()
        for (i in 1:N) {
          set.seed(s)
          if (!is.null(MRst)) {
            SpMask <-
              raster::raster(file.path(DirM, paste0(SpNames[s], ".tif")))
            pseudo.mask_p[[i]] <- pseudo.mask_p[[i]] * SpMask
            if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(RecordsData[[s]]))) {
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
          #     pseudo.mask_p[[i]],
          #     (1 / PrAbRatio) * sum(presences[, 1] == i),
          #     ext = e,
          #     prob = FALSE
          #   )
          absences.0 <- OptimRandomPoints(r=pseudo.mask_p[[i]], n=(1 / PrAbRatio)*sum((presences[, 1] == i)),p=presences[presences[, 1] == i, 2:3] )
          
          colnames(absences.0) <- c("lon", "lat")
          absences[[i]] <- as.data.frame(absences.0)
        }
        
        names(absences) <- 1:N
      }
      
      # Pseudo-Absences allocation with Geographical constrain-----
      if (pseudoabsencesMethod == "GEO_CONST") {
        pseudo.mask_p <-
          inv_geo(e = evnVariables, p = presences[, -1], d = Geo_Buf)
        
        # Split the raster of environmental layer with grids
        pseudo.mask_p <- raster::mask(pseudo.mask, pseudo.mask_p)
        pseudo.mask <- sum(pseudo.mask_p)
        pseudo.mask_p[pseudo.mask_p==0] <- NA
        
        raster::writeRaster(
          pseudo.mask,
          paste(DirSave, paste(SpNames[s], '.tif', sep = ""), sep = '/'),
          format = 'GTiff',
          NAflag = -9999,
          overwrite = TRUE
        )
        
        absences <- list()
        for (i in 1:N) {
          set.seed(s)
          if (!is.null(MRst)) {
            SpMask <-
              raster::raster(file.path(DirM, paste0(SpNames[s], ".tif")))
            pseudo.mask_p[[i]] <- pseudo.mask_p[[i]] * SpMask
            if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(RecordsData[[s]]))) {
              warning(
                "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
              )
              stop(
                "Please try again with a smaller geographical buffer or without restricting the accessible area"
              )
            }
          }
          
          absences.0 <-
            OptimRandomPoints(r = pseudo.mask_p[[i]],
                              n = (1 / PrAbRatio) * sum((presences[, 1] == i)),
                              p = presences[presences[, 1] == i, 2:3])
          
          colnames(absences.0) <- c("lon", "lat")
          absences[[i]] <- as.data.frame(absences.0)
        }
        
        names(absences) <- 1:N
      }
      
      # Pseudo-Absences allocation with Environmental and Geographical  constrain-----
      if (pseudoabsencesMethod == "GEO_ENV_CONST") {
        pseudo.mask_p <- inv_bio(evnVariables, presences[, -1])
        pseudo.mask_pg <-
          inv_geo(e = evnVariables, p = presences[, -1], d = Geo_Buf)
        pseudo.mask_p <- pseudo.mask_p * pseudo.mask_pg
        
        # Split the raster of environmental layer with grids
        pseudo.mask_p <- raster::mask(pseudo.mask, pseudo.mask_p)
        pseudo.mask <- sum(pseudo.mask_p)
        pseudo.mask_p[pseudo.mask_p==0] <- NA
        
        raster::writeRaster(
          pseudo.mask,
          paste(DirSave, paste(SpNames[s], '.tif', sep = ""), sep = '/'),
          format = 'GTiff',
          NAflag = -9999,
          overwrite = TRUE
        )
        
        absences <- list()
        for (i in 1:N) {
          set.seed(s)
          if (!is.null(MRst)) {
            SpMask <-
              raster::raster(file.path(DirM, paste0(SpNames[s], ".tif")))
            pseudo.mask_p[[i]] <- pseudo.mask_p[[i]] * SpMask
            if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(RecordsData[[s]]))) {
              warning(
                "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
              )
              stop(
                "Please try again with another restriction type or without restricting the extent"
              )
            }
          }
          
          absences.0 <- OptimRandomPoints(r=pseudo.mask_p[[i]], n=(1 / PrAbRatio)*sum((presences[, 1] == i)),p=presences[presences[, 1] == i, 2:3] )
          colnames(absences.0) <- c("lon", "lat")
          absences[[i]] <- as.data.frame(absences.0)
        }
        
        names(absences) <- 1:N
      }
      
      # Pseudo-Absences allocation with Environmental and Geographical and k-mean constrain-----
      if (pseudoabsencesMethod == "GEO_ENV_KM_CONST") {
        pseudo.mask_p <- inv_bio(evnVariables, presences[, -1])
        pseudo.mask_pg <-
          inv_geo(e = evnVariables, p = presences[, -1], d = Geo_Buf)
        pseudo.mask_p <- pseudo.mask_p * pseudo.mask_pg
        
        # Split the raster of environmental layer with grids
        pseudo.mask_p <- raster::mask(pseudo.mask, pseudo.mask_p)
        pseudo.mask <- sum(pseudo.mask_p)
        pseudo.mask_p[pseudo.mask_p==0] <- NA
        
        raster::writeRaster(
          pseudo.mask,
          paste(DirSave, paste(SpNames[s], '.tif', sep = ""), sep = '/'),
          format = 'GTiff',
          NAflag = -9999,
          overwrite = TRUE
        )
        
        absences <- list()
        for (i in 1:N) {
          set.seed(s)
          if (!is.null(MRst)) {
            SpMask <- raster::raster(file.path(DirM, paste0(SpNames[s], ".tif")))
            pseudo.mask_p[[i]] <- pseudo.mask_p[[i]] * SpMask
            if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(RecordsData[[s]]))) {
              warning(
                "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
              )
              stop(
                "Please try again with another restriction type or without restricting the extent"
              )
            }
          }
          
          absences.0 <-
            KM_BLOCK(
              raster::rasterToPoints(pseudo.mask_p[[i]])[, -3],
              raster::mask(evnVariables, pseudo.mask_p[[i]]),
              (1 / PrAbRatio) * sum(presences[, 1] == i)
            )
          colnames(absences.0) <- c("lon", "lat")
          
          absences[[i]] <- as.data.frame(absences.0)
        }
        names(absences) <- 1:N
      }
      
      absences <- plyr::ldply(absences, data.frame)
      names(absences) <- c("Partition", "x", "y")
      absences[, c("x", "y")] <- round(absences[, c("x", "y")], 4)
      colnames(absences) <- colnames(presences)
      # Final data.frame result----
      PresAbse <-
        rep(c(1, 0), sapply(list(presences, absences), nrow))
      result <-
        data.frame(
          Sp = SpNames[s],
          PresAbse,
          rbind(presences, absences),
          stringsAsFactors = F
        )
      result <- result[, c("Sp", "x", "y", "Partition", "PresAbse")]
      
      Opt2 <- data.frame(Sp = SpNames[s], Opt2)
      
      # Final data.frame result2----
      out <- list(ResultList = result,
                  BestGridList = Opt2)
      # utils::write.table(result,paste(DirSave, paste0(SpNames[s],".txt"), sep="\\"), sep="\t",row.names=F)
      # return(out)
      results[[s]] <- out
    }
  
  # parallel::stopCluster(cl)
  FinalResult <- dplyr::bind_rows(lapply(results, function(x) x[[1]]))
  FinalInfoGrid <- dplyr::bind_rows(lapply(results, function(x) x[[2]]))
  
  colnames(FinalResult) <- c("sp", "x", "y", "Partition", "PresAbse")
  utils::write.table(
    FinalResult,
    paste(DirSave, "OccBlocks.txt", sep = "\\"),
    sep = "\t",
    row.names = F
  )
  utils::write.table(
    FinalInfoGrid,
    paste(DirSave, "BestPartitions.txt", sep = '/'),
    sep = "\t",
    col.names = T,
    row.names = F
  )
  
  return(FinalResult)
}
