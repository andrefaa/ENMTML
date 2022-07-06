# test ENMTML function

# Prepare data
require(ENMTML)
require(raster)

d_ex <- file.path(tempdir(), 'ENMTML_example')
d_ex
dir.create(d_ex)

# Virtual species occurrences
data("occ")
d_occ <- file.path(d_ex, 'occ.txt')
utils::write.table(occ, d_occ, sep = '\t', row.names = FALSE)
# Five bioclimatic variables for current conditions
data("env")
d_env <- file.path(d_ex, 'current_env_var')
dir.create(d_env)
raster::writeRaster(env, file.path(d_env, names(env)), bylayer=TRUE, format='GTiff', overwrite=TRUE)
# Five bioclimatic variables for future conditions
# (for more details see predictors_future help)
data("env_fut")
d_fut <- file.path(d_ex, 'future_env_var')
dir.create(d_fut)
d0 <- file.path(d_fut, names(env_fut))
sapply(d0, dir.create)

raster::writeRaster(env_fut$`2080_4.5`, file.path(d0[1],
                                                  names(env_fut$`2080_4.5`)), bylayer=TRUE, format='GTiff', overwrite=TRUE)
raster::writeRaster(env_fut$`2080_8.5`, file.path(d0[2],
                                                  names(env_fut$`2080_8.5`)), bylayer=TRUE, format='GTiff', overwrite=TRUE)

# Polygon of terrestrial ecoregions
data("ecoregions")
d_eco <- file.path(d_ex, 'ecoregions')
dir.create(d_eco)
d_eco <- file.path(d_eco, paste0('eco','.shp'))
shapefile(ecoregions, d_eco, overwrite=TRUE)

rm(list = c('d0', 'ecoregions', 'env', 'env_fut', 'occ'))


test_that('ENMTML: basic test works', {
  
  ENMTML(
    pred_dir = d_env,
    proj_dir = d_fut,
    result_dir = 'Test_01',
    occ_file = d_occ,
    sp = 'species',
    x = 'x',
    y = 'y',
    min_occ = 10,
    thin_occ = NULL,
    eval_occ = NULL,
    colin_var = c(method='PCA'),
    imp_var = FALSE,
    sp_accessible_area=c(method='BUFFER', type='2', width='300'),
    pseudoabs_method = c(method = 'RND'),
    pres_abs_ratio = 1,
    part=c(method= 'BLOCK'),
    save_part = FALSE,
    save_final = TRUE,
    algorithm = c("BIO","MAH","DOM","ENF","GLM","GAM","SVM","RDF","MXS","MXD","MLK","GAU"),
    thr = c(type='SORENSEN'),
    msdm = NULL,
    ensemble = NULL,
    extrapolation = FALSE,
    cores = 1
  )
    d_rslt <- file.path(dirname(d_env), 'Test_01')

  #Check if basic txt files were accurately created
  toMatch <- c("Evaluation_Table.txt", "Number_Unique_Occurrences.txt", "Thresholds_Algorithms.txt",
               "InfoModeling.txt","Occurrences_Cleaned.txt")
  expect_true(all(toMatch%in%list.files(d_rslt)))

  #Check if algorithms produced valid predictions
  d <- list.dirs(file.path(d_rslt, "Algorithm"), recursive = FALSE, full.names = TRUE)
  d <- sapply(d, function(x) length(list.files(x, pattern = '.tif')))
  expect_true(all(d==5) == TRUE)

  #Check if PCA folder was correctly created
  expect_true(any(grepl("PCA",list.files(d_env)) == TRUE))

  #Check if future predictions are correct
  f <- list.dirs(file.path(d_rslt, "Projection"), recursive = FALSE, full.names = TRUE)
  f_alg <- sapply(f, function(x) length(list.files(x)))
  expect_true(all(f_alg==12) == TRUE)
  fs <- sapply(f, function(x) list.files(x,full.names = T))
  fs <- sapply(fs, function(x) list.files(x,pattern = '.tif'))
  expect_true(all(sapply(fs, function(x) length(x==5)) == TRUE))

  #Check if masks and blocks were created for all species
  m <- length(list.files(file.path(d_rslt, "Extent_Masks")))
  expect_true(m==5)

  b <- length(list.files(file.path(d_rslt, "BLOCK"),pattern = ".tif$"))
  expect_true(b==5)
})


# test ENMTML function



test_that('ENMTML: test with future projection', {
  
  ENMTML(
    pred_dir = d_env,
    proj_dir = d_fut,
    result_dir = 'Test_01',
    occ_file = d_occ,
    sp = 'species',
    x = 'x',
    y = 'y',
    min_occ = 10,
    thin_occ = NULL,
    eval_occ = NULL,
    colin_var = c(method='PCA'),
    imp_var = FALSE,
    sp_accessible_area=c(method='BUFFER', type='2', width='300'),
    pseudoabs_method = c(method = 'RND'),
    pres_abs_ratio = 1,
    part=c(method= 'BLOCK'),
    save_part = FALSE,
    save_final = TRUE,
    algorithm = c("BIO","MAH","DOM","ENF","GLM","GAM","SVM","RDF","MXS","MXD","MLK","GAU"),
    thr = c(type='SORENSEN'),
    msdm = NULL,
    ensemble = NULL,
    extrapolation = FALSE,
    cores = 1
  )
  d_rslt <- file.path(dirname(d_env), 'Test_01')
  
  #Check if basic txt files were accurately created
  toMatch <- c("Evaluation_Table.txt", "Number_Unique_Occurrences.txt", "Thresholds_Algorithms.txt",
               "InfoModeling.txt","Occurrences_Cleaned.txt")
  expect_true(all(toMatch%in%list.files(d_rslt)))
  
  #Check if algorithms produced valid predictions
  d <- list.dirs(file.path(d_rslt, "Algorithm"), recursive = FALSE, full.names = TRUE)
  d <- sapply(d, function(x) length(list.files(x, pattern = '.tif')))
  expect_true(all(d==5) == TRUE)
  
  #Check if PCA folder was correctly created
  expect_true(any(grepl("PCA",list.files(d_env)) == TRUE))
  
  #Check if future predictions are correct
  f <- list.dirs(file.path(d_rslt, "Projection"), recursive = FALSE, full.names = TRUE)
  f_alg <- sapply(f, function(x) length(list.files(x)))
  expect_true(all(f_alg==12) == TRUE)
  fs <- sapply(f, function(x) list.files(x,full.names = T))
  fs <- sapply(fs, function(x) list.files(x,pattern = '.tif'))
  expect_true(all(sapply(fs, function(x) length(x==5)) == TRUE))
  
  #Check if masks and blocks were created for all species
  m <- length(list.files(file.path(d_rslt, "Extent_Masks")))
  expect_true(m==5)
  
  b <- length(list.files(file.path(d_rslt, "BLOCK"),pattern = ".tif$"))
  expect_true(b==5)
})

