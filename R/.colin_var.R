#' Collinearity reduction of predictor variables
#'
#' @param env SpatRaster An object of class SpatRaster containing the predictors.
#' This function does not allow categorical variables
#'
#' @param colin character. Collinearity reduction method. It is necessary to
#' provide a vector for this argument. The next methods are implemented:
#' \itemize{
#'   \item pearson: Highlights correlated variables according to Pearson correlation. A threshold of maximum correlation
#'   must be specified. Otherwise, a threshold of 0.7 is defined as default.
#'   Usage method = c('pearson', th='0.7').
#'   \item vif: Select variables by Variance Inflation Factor, a threshold can be specified by
#'   user. Otherwise, a threshold of 10 is defined as default.Usage method = c('vif', th = '10').
#'   \item pca: Perform a Principal Component Analysis and use the principal components as the
#'   new predictors. The selected components account for 95\% of the whole variation in the system.
#'   Usage method = c('pca').
#'   \item fa: Perform a Factorial Analysis and select, from the original predictors, the number of factors is defined by Broken-Stick and variables with the highest correlation to the factors are selected.  Usage method = c('fa').
#' }
#' @param proj character. Only used for pca method. Path to a folder that contains sub-folders for the different projection
#' scenarios. Variables names must have the same names as in the raster used in env argument. Usage proj = "C:/User/Desktop/Projections" (see in Details more about the use of this argument)
#' @param maxcell numeric. Number of raster cells to be randomly sampled. Taking a sample could be
#' useful to reduce memory usage for large rasters. If NULL, the function will use all
#' raster cells. Default NULL. Usage maxcell = 50000.
#'
#' @return
#' #' If 'pearson', the following elements are created:
#' \itemize{
#' \item cor_table: a matrix object with pairwise correlation values of the environmental variables
#' \item cor_variables: a list object with the same length of the number of environmental values containing the pairwise relations that exceeded the correlation threshold for each one of the environmental variables
#' }
#'
#' If 'vif' method, the following elements are created:
#' \itemize{
#' \item env: a SpatRaster object with selected environmental variables
#' \item removed_variables: a character vector with removed environmental variables
#' \item vif_table: a data frame with VIF values for all environmental variables
#' }
#'
#' If 'pca' method, the following elements are created:
#' \itemize{
#' \item env: SpatRaster with scores of selected principal component (PC) that sum up 95\% of the
#' whole variation or original environmental variables
#' \item coefficients: a matrix with the coefficient of principal component (PC) for predictors
#' \item cumulative_variance: a tibble with the cumulative variance explained in selected principal component (PC)
#' }
#'
#' If 'fa' method, returns a list with the following elements:
#' \itemize{
#' \item env: SpatRaster with scores of selected variables due to correlation to factors.
#' \item number_factors: number of factors selected according to the Broken-Stick criteria,
#' \item removed_variables: removed variables,
#' \item uniqueness: uniqueness of each environmental variable according to the factorial analysis,
#' \item loadings: environmental variables loadings in each of the chosen factors
#' }
#'
#' @details In the case of having environmental variables for the current conditions and other time
#' periods (future or present), it is recommended to perform the PCA analysis with the current
#' environmental condition and project the PCA for the other time periods. To do so, it is necessary
#' to use “proj” argument. Path to a folder (e.g., projections) that contains sub-folders for the
#' different projection scenarios (e.g., years and emissions). Within each sub-folder must be stored
#' single or multiband rasters with the environmental variables.
#'
#' For example:
#'
#' C:/Users/my_pc/projections/ \cr
#'     ├── MRIESM_2050_ssp126 \cr
#'     │   └── var1.tif\cr
#'     │   └── var2.tif\cr
#'     │   └── var3.tif\cr
#'     ├── MRIESM_2080_ssp585\cr
#'     │   └── var1.tif\cr
#'     │   └── var2.tif\cr
#'     │   └── var3.tif\cr
#'     ├── UKESM_2050_ssp370\cr
#'     │   └── var1.tif\cr
#'     │   └── var2.tif\cr
#'     │   └── var3.tif
#'
#' If pca method is run with time projections, colin_var function will create the
#' Projection_PCA (the exact path is in the path object returned by the function) with the same
#' system of sub-folders and multiband raster with the principal components (pcs.tif)
#'
#' C:/Users/my_pc/Projection_PCA/\cr
#'       ├── MRIESM_2050_ssp126\cr
#'       │   └── pcs.tif           # a multiband tif with principal components\cr
#'       ├── MRIESM_2080_ssp585\cr
#'       │   └── pcs.tif\cr
#'       ├── UKESM_2050_ssp370\cr
#'       │   └── pcs.tif
#'
#' @export
#' @importFrom dplyr tibble
#' @importFrom stats na.omit cor lm prcomp factanal
#' @importFrom terra rast as.data.frame spatSample subset predict scale writeRaster
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#'
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' # Perform pearson collinearity control
#' var <- .colin_var(env = somevar, colin = c("pearson", th = "0.7"))
#' var$cor_table
#' var$cor_variables
#'
#' # For all .colin_var methods it is possible to take a sample or raster to reduce memory
#' var <- .colin_var(env = somevar, colin = c("pearson", th = "0.7"), maxcell = 10000)
#' var$cor_table
#' var$cor_variables
#'
#' # Perform vif collinearity control
#' var <- .colin_var(env = somevar, colin = c("vif", th = "8"))
#' var$env
#' var$removed_variables
#' var$vif_table
#'
#' # Perform pca collinearity control
#' var <- .colin_var(env = somevar, colin = c("pca"))
#' plot(var$env)
#' var$env
#' var$coefficients
#' var$cumulative_variance
#'
#'
#' # Perform pca collinearity control with different projections
#' ## Below will be created a set of folders to simulate the structure of the directory where
#' ## environmental variables are stored for different scenarios
#' dir_sc <- file.path(tempdir(), "projections")
#' dir.create(dir_sc)
#' dir_sc <- file.path(dir_sc, c('scenario_1', 'scenario_2'))
#' sapply(dir_sc, dir.create)
#'
#' somevar <-
#'   system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' terra::writeRaster(somevar, file.path(dir_sc[1], "somevar.tif"), overwrite=TRUE)
#' terra::writeRaster(somevar, file.path(dir_sc[2], "somevar.tif"), overwrite=TRUE)
#'
#' ## Perform pca with projections
#' dir_w_proj <- dirname(dir_sc[1])
#' dir_w_proj
#' var <- .colin_var(env = somevar, colin = "pca", proj = dir_w_proj)
#' var$env
#' var$coefficients
#' var$cumulative_variance
#' var$proj
#'
#'
#' # Perform fa colinearity control
#' var <- .colin_var(env = somevar, colin = c("fa"))
#' var$env
#' var$number_factors
#' var$removed_variables
#' var$uniqueness
#' var$loadings
#' }
#'
.colin_var <- function(env,
                       colin,
                       proj = NULL,
                       maxcell = NULL,
                       pred_dir) {
  . <- NULL
  
  #Create sub-folder for collinearity results
  dir.create(file.path(pred_dir,'./Collinearity/Tables'),recursive=T)
 
  #1.Pearson Correlation----
  if (any(colin %in% "pearson")) {
    if (is.na(colin["th"])) {
      th <- 0.7
    } else {
      th <- as.numeric(colin["th"])
    }
    
    if(is.null(maxcell)){
      h <- terra::as.data.frame(env,na.rm=T)
    } else {
      # Raster random sample
      set.seed(10)
      h <- env %>%
        terra::spatSample(., size = maxcell, method="random", na.rm=TRUE, as.df=TRUE)
    }
    
    h2 <- base::abs(stats::cor(h, method = "pearson"))
    diag(h2) <- 0
    #Save Correlation Table
    utils::write.csv(h2,file.path(pred_dir,'/Collinearity/Tables/Pearson_Correlation_table.csv'),row.names = T)
    
    #Select only variables below threshold
    LOOP <- T
    while(any(h2>th)){
      h <- h[,-which.max(rowSums(h2))[1]]
      h2 <- base::abs(stats::cor(h, method = "pearson"))
      diag(h2) <- 0
      if(ncol(h)==2){
        h <- h[,1]
        break
      }
    }
    
    utils::write.csv(h2,file.path(pred_dir,'/Collinearity/Tables/Pearson_Correlation_table_selected_variables.csv'),row.names = T)
    
    #SpatRaster subset
    env <- terra::subset(env,subset=names(h))
    terra::writeRaster(env,
                       file.path(pred_dir,'/Collinearity',paste0(names(env),'.tif')),
                       filetype='GTiff',
                       overwrite = T)
    
    
    result <- list(
      env = env,
      sel_variables = names(env)
    )
  }
  
  #2.VIF----
  if (any(colin %in% "vif")) {
    if (is.null(colin["th"])) {
      th <- 10
    } else {
      th <- as.numeric(colin["th"])
    }
    
    if(is.null(maxcell)){
      h <- terra::as.data.frame(env)
    } else {
      # Raster random sample
      set.seed(10)
      h <- env %>%
        terra::spatSample(., size = maxcell, method="random", na.rm=TRUE, as.df=TRUE) %>%
        stats::na.omit()
    }
    
    LOOP <- TRUE
    if (nrow(h) > 200000) {
      h <- h[sample(1:nrow(h), 200000), ]
    }
    n <- list()
    n$variables <- colnames(h)
    exc <- c()
    
    while (LOOP) {
      v <- rep(NA, ncol(h))
      names(v) <- colnames(h)
      for (i in 1:ncol(h)) {
        v[i] <- 1 / (1 - summary(lm(h[, i] ~ ., data = h[,-i]))$r.squared)
      }
      if(length(v)==terra::nlyr(env)){
        vif_tab <- data.frame(Variables = names(v), VIF = as.vector(v))
        utils::write.csv(vif_tab,file.path(pred_dir,'/Collinearity/Tables/VIF_table.csv'),row.names = F)
      }
      
      if (v[which.max(v)] >= th) {
        ex <- names(v[which.max(v)])
        exc <- c(exc, ex)
        h <- h[, -which(colnames(h) == ex)]
      } else {
        LOOP <- FALSE
      }
    }
    if (length(exc) > 0) {
      n$excluded <- exc
    }
    
    v <- rep(NA, ncol(h))
    names(v) <- colnames(h)
    for (i in 1:ncol(h)) {
      v[i] <- 1 / (1 - summary(stats::lm(h[, i] ~ ., data = h[-i]))$r.squared)
    }
    
    # n$corMatrix <- stats::cor(x, method = "pearson")
    vif_tab <- data.frame(Variables = names(v), VIF = as.vector(v))
    utils::write.csv(vif_tab,file.path(pred_dir,'/Collinearity/Tables/VIF_table_selected_variables.csv'),row.names = F)
    
    # diag(n$corMatrix) <- 0
    env <-
      terra::subset(env, subset = n$results$Variables)
    
    result <- list(
      env = env,
      sel_variables = names(v)
      )
  }
  
  #3.PCA----
  if (any(colin %in% "pca")) {
    
    if(is.null(maxcell)){
      h <- terra::as.data.frame(env, xy = FALSE, na.rm = TRUE)
    } else {
      # Raster random sample
      set.seed(10)
      h <- env %>%
        terra::spatSample(., size = maxcell, method="random", na.rm=TRUE, as.df=TRUE) %>%
        stats::na.omit()
    }
    
    #Means and SD for Projection
    means <- colMeans(h)
    stds <- apply(h, 2, stats::sd)
    
    # Standardize values
    vnmes <- names(means)
    for(nl in 1:length(vnmes)){
      env[[vnmes[nl]]] <- (env[[vnmes[nl]]]-means[vnmes[nl]])/stds[vnmes[nl]]
    }
    
    #PCA
    h <- stats::prcomp(h,
                       retx = TRUE,
                       scale. = FALSE,
                       center = FALSE
    )
    cof <- h$rotation
    coefficients <- data.frame(cof) %>% dplyr::tibble(variable = rownames(.), .)
    utils::write.csv(coefficients,file.path(pred_dir,'/Collinearity/Tables/PCA_Coefficients.csv'),row.names = F)
    
    cvar <- summary(h)$importance["Cumulative Proportion", ]
    naxis <- Position(function(x) {
      x >= 0.95
    }, cvar)
    cvar <- data.frame(Cumulative_Variance = cvar) %>% dplyr::tibble(Axis = rownames(.), .)
    utils::write.csv(cvar,file.path(pred_dir,'/Collinearity/Tables/PCA_Cumulative_Variance.csv'),row.names = F)
    
    env <- terra::predict(env, h, index = 1:naxis)
    names(env) <- cvar$Axis[1:naxis]
    terra::writeRaster(
      env,
      file.path(pred_dir,'/Collinearity', "pcs.tif"),
      filetype='GTiff',
      overwrite=TRUE
    )
    
    result <- list(
      env = env
    )
    
    if (!is.null(proj)) {
      dpca <- file.path(dirname(proj), "Projection_PCA")
      dir.create(dpca)
      subfold <- list.files(proj)
      subfold <- as.list(file.path(dpca, subfold))
      sapply(subfold, function(x) {
        dir.create(x)
      })
      
      proj <- list.files(proj, full.names = TRUE)
      for (i in 1:length(proj)) {
        scen <- terra::rast(list.files(proj[i], full.names = TRUE))
        scen <- scen[[names(means)]]
        scen <- terra::scale(scen, center = means, scale = stds)
        scen <- terra::predict(scen, p)
        terra::writeRaster(
          scen,
          file.path(subfold[[i]], "pcs.tif"),
          filetype='GTiff',
          overwrite=TRUE
        )
      }
      
      result$proj <- dpca
    }
  }
  
  #4.Factor Analysis----
  if (any(colin %in% "fa")) {
    p <- terra::scale(env, center = TRUE, scale = TRUE)
    
    if(is.null(maxcell)){
      p <- terra::as.data.frame(p, xy = FALSE, na.rm = TRUE)
    } else {
      # Raster random sample
      set.seed(10)
      p <- p %>%
        terra::spatSample(., size = maxcell, method="random", na.rm=TRUE, as.df=TRUE) %>%
        stats::na.omit()
    }
    
    if (nrow(p) > 200000) {
      p <- p[sample(1:nrow(p), 200000), ]
    }
    
    e <- eigen(stats::cor(p))
    len <- length(e$values)
    a <- NULL
    r <- NULL
    
    for (j in 1:len) {
      a[j] <- 1 / len * sum(1 / (j:len))
      r[j] <- e$values[j] / (sum(e$values))
    }
    
    ns <- length(which(r > a))
    
    fit <-
      tryCatch(
        stats::factanal(
          x = p,
          factors = ns,
          rotation = "varimax",
          lower = 0.001
        ),
        error = function(e) {
          stats::factanal(
            x = p,
            factors = ns - 1,
            rotation = "varimax",
            lower = 0.001
          )
        }
      )
    
    sel <-
      row.names(fit$loadings)[apply(fit$loadings, 2, which.max)]
    rem <-
      row.names(fit$loadings)[!row.names(fit$loadings) %in% sel]
    
    env <- terra::subset(env, sel)
    
    h <- fit$loadings %>%
      matrix() %>%
      data.frame()
    colnames(h) <- paste("Factor", 1:ncol(h), sep = "_")
    
    result <- list(
      env = env,
      number_factors = fit$factors,
      removed_variables = rem,
      uniqueness = fit$uniquenesses,
      loadings = dplyr::tibble(Variable = rownames(fit$loadings), h)
    )
  }
  
  return(result)
}