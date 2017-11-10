setwd ("C://Users//Andre//Desktop//Andre//Mestrado//Layers//IPCC//pH_RCP8")

library(raster)

r <- raster("ph_Oyr_MIROC-ESM_esmrcp85_r1i1p1_2006-2100.nc")
z <- raster("ph_Oyr_CMCC-CESM_rcp85_r1i1p1_2006-2100.nc", band="2", zvar="ph")
#band = variacao no tempo(1=2006, 2=2007,etc...)
extent(r) <- c(-180.0,180.0,-90.0,90.0)

s <- raster(ncol=4320, nrow=1680)
res(s) <- 0.08
#raster de resolucao 30 arc-sec

s <- resample(r, s, method='bilinear')
writeRaster(s, filename="rasterph2006.asc", format= "ascii", overwrite=T)

