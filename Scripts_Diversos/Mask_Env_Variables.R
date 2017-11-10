#Funcao para cortar variaveis ambientais (.tif/.asc) a partir de uma mascara (.asc)

#Parametros
#dir_in: diretorio local dos arquivos asc
#dir_out:diretorio para salvar os arquivos cortados
#format: Formato das variaveis ambientais (.tif/.asc)

Mask_Env_Variables <- function (dir_in,dir_out,mask){

  setwd(dir_in)
  
  library(raster)
  library(SDMTools)
  
  #Le os arquivos asc e junta em um bloco
  
  bloco<-brick(stack(list.files(pattern= "\\.asc$")))

  #Carregar a mascara para corte
  
  if(mask=="raster"){
    print("Selecione a mascara(.asc/.tif/.bil):")
    masc<-raster(file.choose())
  }
  if(mask=="shp"){
    print("Selecione a mascara(.shp):")
    masc <- readOGR(file.choose())
  }
  e <- extent(masc)
  
  #Corta o RasterBrick pela extensao desejada
  
  corte<-mask(crop(bloco,e),masc)
  
  #Salva os novos ascs(sobreescrevendo os antigos)
  setwd(dir_out)
  
  writeRaster(corte,names(corte),bylayer=T,suffix='names',format="GTiff",overwrite=T)
}

Mask_Env_Variables("C:\\OneDrive\\Trabajos\\03 - MataAtlanticaDanira\\MSDM",
                   "C:\\OneDrive\\Trabajos\\03 - MataAtlanticaDanira\\MSDM\\Crop",
                   "shp")
