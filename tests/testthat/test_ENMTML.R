# test ENMTML function

# Prepare data
require(ENMTML)
require(raster)
setwd('C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/65-ENM_TheMetaLand/ENMTML/tests/testthat')
d_ex <- file.path(getwd(), 'ENMTML_example')
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


test_that('ENMTML: basic test doesn\'t work', {
  ENMTML(
    pred_dir = d_env,
    proj_dir = NULL,
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
    sp_accessible_area = c(method='BUFFER', type='2', width='500'),
    pseudoabs_method = c(method = 'RND'),
    pres_abs_ratio = 1,
    part=c(method= 'KFOLD', folds='2'),
    save_part = FALSE,
    save_final = TRUE,
    algorithm = c("BIO", "GLM"),
                  # "MAH","DOM","ENF","GLM","GAM","SVM","BRT","RDF","MXS","MXD","MLK","GAU"),
    thr = c(type='MAX_TSS'),
    msdm = NULL,
    ensemble = NULL,
    extrapolation = FALSE,
    cores = 1
  )
  d_rslt <- file.path(dirname(d_env), 'Test_01')

  # Aqui temos que escreve o que esperamos que retorne o pacote ou seja
  # todos os exepct que achemos necessários deixei dois exemplos
  expect_true(any(grepl(
    "Evaluation_Table.txt", list.files(d_rslt)
  )) == TRUE)

  d <- list.dirs(file.path(d_rslt, "Algorithm"), recursive = FALSE, full.names = TRUE)
  d <- sapply(d, function(x) length(list.files(x, pattern = '.tif')))
  expect_true(all(d==5) == TRUE)
})

# Daqui para baixo podemos criar outros test


# Ao final de todo  escrevi esse código para eliminar as saidas dos testes. Caso contrário temos que estar apagando na mão o tempo todo
unlink(d_ex, recursive = TRUE)
