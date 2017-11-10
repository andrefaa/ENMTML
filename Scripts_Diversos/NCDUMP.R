setwd ("C://Users//Andre//Desktop//Andre//Mestrado//Layers//IPCC//pH_RCP8")

library(ncdf)
nc <- open.ncdf("teste_o2")
myarray <- get.var.ncdf(nc, myvar)


writeRaster(o2_Oyr_CMCC-CESM_rcp85_r1i1p1_2006-2100, "o2_ascii")

library(ncdf)
library(maps)
library(fields)
library(scales)

nc <- open.ncdf("teste_o2")
lat <- get.var.ncdf(nc,"lat")
long <- get.var.ncdf(nc,"lon")
time <- get.var.ncdf(nc,"time")
var <- get.var.ncdf(nc, "o2")
close.ncdf(nc)

i=100
nom_var="Hgt"
var.slice <- var[,,i]

nlevel=64   #number of colors
image.plot(seq(from=-20, to=30, by=2),seq(from=28, to=60, by=2),var.slice, 
           xlab = "Longitude", ylab = "Latitude", legend.shrink = 0.9, 
           legend.width = 1.2, legend.mar = NULL, main =nom_var,
           graphics.reset = FALSE, horizontal = FALSE, bigplot = NULL
           , smallplot = NULL, legend.only = FALSE, col = tim.colors(nlevel),
           lab.breaks=NULL, axis.args=NULL)


map(add=TRUE, fill=TRUE, col = alpha("grey", 0.5))
abline(v=seq(-21, 31, 2), lty=2)
abline(h=seq(27, 61, 2), lty=2)