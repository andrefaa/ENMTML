
start_cluster <- function(cores) {
  cl <- NULL
  if (cores > 1) {
    if (Sys.getenv("RSTUDIO") == "1" &&
        !nzchar(Sys.getenv("RSTUDIO_TERM")) &&
        Sys.info()["sysname"] == "Darwin" &&
        as.numeric(gsub('[.]', '', getRversion())) >= 360) {
      cl <- parallel::makeCluster(cores,outfile="", setup_strategy = "sequential")
    }else{
      cl <- parallel::makeCluster(cores,outfile="")
    }
    doParallel::registerDoParallel(cl)
  }
  return(cl)
}

stop_cluster <- function(cluster) {
  if (!is.null(cluster)) {
    parallel::stopCluster(cluster)
  }
}
