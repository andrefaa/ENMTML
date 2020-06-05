#' Virtual species occurrences.
#'
#' A dataset containing occurrences of five virtual species throughout South America.
#'
#' @format A data.frame with 215 rows and 3 variables:
#' \describe{
#'   \item{sp}{virtual species names}
#'   \item{x}{longitude of species occurrences}
#'   \item{y}{latitude of species occurrences}
#' }
#' @examples
#' require(ENMTML)
#' require(raster)
#' data("env")
#' data("occ")
#'
#' plot(env[[1]], col='gray')
#' points(points(occ[,2:3],
#'              col=rainbow(5)[as.factor(occ$species)],
#'              pch=19))
"occ"

#' Raster with 5 bioclimatic variables for current climate conditions
#'
#' A RasterBrick with Bioclimatic variables (Bio1, Bio3, Bio4, Bio12 and Bio15) from a portion of South America (source Worldclim: https://worldclim.org/version2).
#' @examples
#' require(ENMTML)
#' require(raster)
#' data("env")
#' env
#' plot(env)
"env"

#' Raster with 5 bioclimatic variables for 2080's climate conditions
#'
#' A list of RasterBrick with Bioclimatic variables (Bio1, Bio3, Bio4, Bio12 and Bio15) from a portion of South America (source GCM downscaled data portal: http://ccafs-climate.org/downscaling/). This dataset contains projections for 2080 based on MOHC HadGEM2-ES model for the representative concentration pathway of 4.5 and 8.5.
#'
#' @examples
#' require(ENMTML)
#' require(raster)
#' data("env_fut")
#' env_fut
#' plot(env_fut$`2080_4.5`)
#' plot(env_fut$`2080_8.5`)
"env_fut"

#' Polygon with terrestrial ecoregion
#'
#' This is a spatially simplified shapefile with terrestrial polygon ecoregion from a portion of South America from (Olson et al. 2001).
#'
#' @references
#'
#' Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell, G. V., Underwood, E. C., … others. (2001). Terrestrial ecoregions of the world: A new map of life on earth. BioScience, 51(11), 933–938.
#'
#' @examples
#' require(ENMTML)
#' require(raster)
#' data("ecoregions")
#' plot(ecoregions, col=rainbow(92)[as.factor(ecoregions@data$ECO_ID)])
"ecoregions"
