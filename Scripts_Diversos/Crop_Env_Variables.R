Crop_Env_Variables <- function (DirI,DirO,method=""){

#Funcao para cortar variaveis ambientais (tif)

  ###################################################
  ################### IMPORTANTE ####################
  ###################################################
  
  #sobreescreve os arquivos originais,caso DirI=DirO!

  #Parametros
    #DirI: diretorio input (arquivos em .tif)
    #DirO:diretorio para salvar os arquivos cortados
    #method: Metodo usado para cortar os arquivos em input
      #raster:usa um arquivo raster (.tif) como mascara
      #shape: usa um arquivo shape como mascara
      #extent: corta os layers a partir de uma caixa com limites a serem definidos pelo usuario
  
  setwd(DirI)
  
  library(raster)
  library(SDMTools)
  library(rgdal)
  
  #Read Env as a Raster Brick
  
  bloco<-brick(stack(list.files(pattern= ".tif")))
  
  #Selecionar o metodo de corte  
  
  if(method=="raster"){
    #Load Mask File
    print("Choose raster mask file:")
    mask<-raster(file.choose())
    e <- extent(mask)
  }
  
  if(method=="extent"){
    ans <- "N"
    
    while(ans=="N"){
      x11()
      plot(bloco[[1]])
      grid()
      #Selct Longitude Limits
        cat("Select Xmin:")
        x.min <- as.integer(readLines(n = 1))
        cat("Select Xmax:")
        x.max <- as.integer(readLines(n = 1))
        while(x.min>x.max){
          warning("Xmin is larger than Xmax! Selection invalid! Please select valid Xmin & Xmax limits")
          cat("Select Xmin:")
          x.min <- as.integer(readLines(n = 1))
          cat("Select Xmax:")
          x.max <- as.integer(readLines(n = 1))
        }
      #Selct Latitude Limits
        cat("Select Ymin:")
        y.min <- as.integer(readLines(n = 1))
        cat("Select Ymax:")
        y.max <- as.integer(readLines(n = 1))
        while(x.min>x.max){
          warning("Ymin is larger than Ymax! Selection invalid! Please select valid Ymin & Ymax limits")
          cat("Select Ymin:")
          y.min <- as.integer(readLines(n = 1))
          cat("Select Ymax:")
          y.max <- as.integer(readLines(n = 1))
        }
      e<-extent(x.min,x.max,y.min,y.max)
      plot(crop(bloco[[1]],e))
      cat("Are those specifications OK?")
      ans <- as.character(readLines(n = 1))
    }
  }
  
  if(method=="shape"){
    print("Choose shape mask file:")
    mask<-readOGR(file.choose())
    e <- extent(mask)
  }
  
  #Crop RasterBrick by the desired extent
  corte<- crop(bloco,e)
  
  #Extract by mask (if mask is provided)
  if(method=="shape" || method=="raster"){
    corte <- mask(corte,mask)
  }
  
  #Save Output
  setwd(DirO)
  
  if(DirO==DirI){
    warning("DirO=DirI, files will be overwritten! Are you sure? Y/N")
    res <- as.integer(readLines(n = 1))
    if(res=="N"){
      cat("Select a new output directory ")
      DirO <- as.character(readLines(n = 1))
    }
  }

  #Fix Problems caused by dots in filenames
  nomes <- gsub("\\.","_",names(corte))
  
  #Save Output
  writeRaster(corte,nomes,format="GTiff",bylayer=T,overwrite=F)
}

Crop_Env_Variables(DirI="C:\\OneDrive\\Trabajos\\05 - Poliquetas\\Art1-Lacunas\\BenthicLayers_MinDepth",
                   DirO="C:\\OneDrive\\Trabajos\\05 - Poliquetas\\Art1-Lacunas\\BenthicLayers_MinDepth\\Atlantic",
                   method="extent")


# Select Xmin:
#   -60
# Select Xmax:
#   -30
# Select Ymin:
#   -40
# Select Ymax:
#   5