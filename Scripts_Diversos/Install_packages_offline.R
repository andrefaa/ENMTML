#########################################################################################
########################################## STEP 1 #######################################
#########################################################################################

#' Get package dependencies
#'
#' @param packs A string vector of package names
#'
#' @return A string vector with packs plus the names of any dependencies
getDependencies <- function(packs){
  dependencyNames <- unlist(
    tools::package_dependencies(packages = packs, db = available.packages(), 
                                which = c("Depends", "Imports"),
                                recursive = TRUE))
  packageNames <- union(packs, dependencyNames)
  packageNames
}
# Calculate dependencies
packages <- getDependencies(c("raster","sp","dismo","randomForest","kernlab","maxnet","RStoolbox","flexclust"))

#########################################################################################
########################################## STEP 2 #######################################
#########################################################################################

# Download the packages to the working directory.
# Package names and filenames are returned in a matrix.
setwd("D:/AulaPraticaII/packages/")
pkgInfo <- download.packages(pkgs = packages, destdir = getwd(), type = "win.binary")
# Save just the package file names (basename() strips off the full paths leaving just the filename)
write.csv(file = "pkgFilenames.csv", basename(pkgInfo[, 2]), row.names = FALSE)

#########################################################################################
########################################## STEP 3 #######################################
#########################################################################################

# Set working directory to the location of the package files
setwd("D:/AulaPraticaII/packages/")

# Read the package filenames and install
pkgFilenames <- read.csv("pkgFilenames.csv", stringsAsFactors = FALSE)[, 1]
install.packages(pkgFilenames, repos = NULL, type = "win.binary")
