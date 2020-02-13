##%######################################################%##
#                                                          #
####                    Partial_ROC                     ####
#                                                          #
##%######################################################%##

#' Modified from the original function from 
#' ellipseenm (https://github.com/marlonecobos/ellipsenm/blob/master/R/partial_roc.R)

#PartialROC
partial_roc <- function(prediction, test_data, longitude, latitude, error = 5,
                        
                        iterations = 500, percentage = 50) {
  
  
  calc_aucDF <- function(big_classpixels, fractional_area, test_data, n_data,
                         
                         n_samp, error_sens) {
    
    rowsID <- sample(x = n_data, size = n_samp, replace = TRUE)
    
    test_data1 <- test_data[rowsID]
    
    omssion_matrix <- big_classpixels > test_data1
    
    sensibility <- 1 - colSums(omssion_matrix) / n_samp
    
    xyTable <- data.frame(fractional_area, sensibility)
    
    less_ID <- which(xyTable$sensibility <= error_sens)
    
    xyTable <- xyTable[-less_ID, ]
    
    xyTable <- xyTable[order(xyTable$fractional_area, decreasing = F), ]
    
    auc_pmodel <- trap_roc(x=xyTable$fractional_area, y=xyTable$sensibility)
    
    auc_prand <- trap_roc(x=xyTable$fractional_area, y=xyTable$fractional_area)
    
    auc_ratio <- auc_pmodel / auc_prand
    
    auc_table <- data.frame(auc_pmodel, auc_prand, auc_ratio = auc_ratio )
    
    return(auc_table)
    
  }
  
  
  
  # -----------
  
  # preparing data
  c_pred <- class(prediction)[1]
  c_tdat <- class(test_data)[1]
  min_pred <- ifelse(c_pred == "numeric", min(prediction, na.rm = TRUE),
                     
                     prediction@data@min)
  
  max_pred <- ifelse(c_pred == "numeric", max(prediction, na.rm = TRUE),
                     
                     prediction@data@max)
  
  
  
  prediction <- round((prediction / max_pred) * 1000)
  
  if (c_pred == "RasterLayer") {
    
    if (c_tdat != "numeric") {
      
      test_data <- na.omit(raster::extract(prediction,
                                           
                                           test_data[, c(longitude, latitude)]))
      
    } else {
      
      test_data <- round((test_data / max_pred) * 1000)
      
    }
    
    classpixels <- data.frame(raster::freq(prediction, useNA = "no"))
    
  } else {
    
    test_data <- round((test_data / max_pred) * 1000)
    
    vals <- na.omit(unique(prediction))
    
    classpixels <- data.frame(value = vals, count = c(table(prediction)),
                              
                              row.names = 1:length(vals))
    
  }
  
  
  
  # -----------
  
  # analysis
  
  if(min_pred == max_pred){
    
    warning("\nprediction has no variability, pROC will return NA.\n")
    
    
    
    p_roc <- rep(NA, 2)
    
    names(p_roc) <- c(paste0("Mean_AUC_ratio_at_", error, "%"), "pval_pROC")
    
    
    
    auc_ratios <- rep(NA, 3)
    
    names(auc_ratios) <- c("Prediction_partial_AUC", "Random_curve_partial_AUC",
                           
                           "AUC_ratio")
    
    
    
    p_roc_res <- list(pROC_summary = p_roc, pROC_results = auc_ratios)
    
    
    
    return(p_roc_res)
    
  }else {
    
    classpixels <- classpixels %>%
      
      dplyr::mutate_(value = ~rev(value),
                     
                     count = ~rev(count),
                     
                     totpixperclass = ~cumsum(count),
                     
                     percentpixels = ~totpixperclass/sum(count)) %>%
      
      dplyr::arrange(value)
    
    
    
    
    
    error_sens <- 1 - (error / 100)
    
    prediction_errors <- classpixels[, "value"]
    
    fractional_area <- classpixels[, "percentpixels"]
    
    n_data <- length(test_data)
    
    n_samp <- ceiling((percentage / 100) * n_data)
    
    
    
    big_classpixels <- matrix(rep(prediction_errors, each = n_samp),
                              
                              ncol = length(prediction_errors))
    
    
    
    partial_AUC <- 1:iterations %>%
      
      purrr::map_df(~calc_aucDF(big_classpixels, fractional_area, test_data,
                                
                                n_data, n_samp, error_sens))
    
    
    
    naID <- !is.na(partial_AUC$auc_ratio)
    
    nona_valproc <- partial_AUC$auc_ratio[naID]
    
    mauc <- mean(nona_valproc)
    
    proc <- sum(nona_valproc <= 1) / length(nona_valproc)
    
    valid_iter <- length(nona_valproc)
    
    
    
    p_roc <- c(mauc, proc, valid_iter)
    
    names(p_roc) <- c(paste0("Mean_AUC_ratio_at_", error, "%"), "pval_pROC",
                      
                      "Valid_iterations")
    
    
    
    auc_ratios <- partial_AUC
    
    names(auc_ratios) <- c("Prediction_partial_AUC", "Random_curve_partial_AUC",
                           
                           "AUC_ratio")
    
    
    
    p_roc_res <- list(pROC_summary = p_roc, pROC_results = auc_ratios)
    
    
    
    return(p_roc_res)
    
  }
  
}


#Trap_roc
trap_roc <- function(x, y) {
  
  x_s = length(x)
  
  y_s = length(y)
  
  if(x_s != y_s){
    
    stop("x  and y must have the same length!")}
  
  auc = 0
  
  for(i in 2:x_s) {
    
    auc <- auc+0.5*(y[i-1] + y[i])*(x[i]-x[i-1])
    
  }
  
  return (auc)
  
}