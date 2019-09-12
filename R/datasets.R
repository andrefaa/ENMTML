#' Virtual species occurrences.
#'
#' A dataset containing occurrences of five virtual species throughout South America
#'
#' @format A data.frame with 215 rows and 3 variables:
#' \describe{
#'   \item{sp}{virtual species names}
#'   \item{x}{longitude of species occurrences}
#'   \item{y}{latitude of species occurrences}
#'   ...
#' }
#' @examples
#'
#' require(raster)
#' require(dplyr)
#' data("predictors")
#' data("occ")
#' plot(predictors[[1]])
#'
#' par(mfrow=c(2,3))
#' for(i in unique(occ[,1])){
#' plot(predictors[[1]], col='gray', legend=FALSE, main=i)
#' p <- occ %>% dplyr::filter(species==i)
#' points(p[,2:3], col='red', pch=19)
#' }
#' par(mfrow=c(1,1))
"occ"
