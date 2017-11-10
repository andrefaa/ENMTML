MatrizComposicao <- function(dir.in,dir.out){
  
  library(raster)
  
  lista <- list.files(dir.in,patter="\\.asc$")
  ascs <- stack(paste(dir.in,lista,sep="/"))
  txt <- as.data.frame(rasterToPoints(ascs))
  write.table(txt,paste(dir.out,"Matriz_Composicao.txt",sep="/"),row.names=F,sep="\t")
}

MatrizComposicao("C:\\Users\\decoa\\Desktop\\OdonatasBrasil\\Result\\AMZ\\AMZ_BIN_CROP","C:\\Users\\decoa\\Desktop\\OdonatasBrasil\\Result\\AMZ\\AMZ_BIN_CROP")
