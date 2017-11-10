### By Caroline Nobrega

######## SDM no R!!! ##############
#install.packages ('rJava', 'maptools', 'raster')
library(maptools)
library(raster)
library(graphics)
#library(rJava)
library (dismo)
#memory.limit(size=8012) 

### Importar dados de ocorrencia ###

print ('Selecionar tabela dos pontos de ocorrencia - a mesma usada no Maxent (.csv).')
pts <- read.table(file=choose.files(), head=T, sep = ",")
# Validar pontos -> plotar pontos em um mapa 
data(wrld_simpl)
plot(wrld_simpl , xlim=c(-80,-35), ylim=c(-50,10), axes=TRUE, col='light yellow')
box() # restore the box around the map
points(pts$long, pts$lat, col='red', pch=20, cex=0.75) # plot points

### Importar variáveis ambientais (.asc) ###

print ('Selecionar diretório com as variáveis ambientais. Atenção: Nesse diretório não deve haver outros arquivos criados pelo ArcGis!')
files <- list.files(path=choose.dir(), pattern='.asc', full.names=TRUE )
mask <- raster(files[[1]]); 
ext<-extent(mask); predictors <- stack(files);

### Preparar os dados ###
# Selecionar pontos espacialmente únicos #
cell <- cellFromXY(mask, pts[,2:3]) # get the cell number for each point
dup <- duplicated(cbind(pts[,1],cell))
pts1 <- pts[!dup, ]# select the records that are not duplicated
#Criar pseudo-ausencias
pseudo0 <- randomPoints(predictors, n=10000, ext=ext, extf = 1.25)
pseudo<-data.frame(sp=rep(0,nrow(pseudo0)), long=pseudo0[,1], lat=pseudo0[,2])
names(pseudo)<-names(pts1)
# Extrair valores dos rasters para os pontos de ocorrencia #
pseudovals <- cbind(pseudo[,1],extract(predictors, pseudo[,2:3])); colnames(pseudovals)[1]<-'sp'
  # Pairs plot of the values of the climate data
  panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1));   r = (cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex )}
  #Plotar matriz com graficos e valores de correlacao das variaveis para os pontos de pseudo
#  pairs(pseudovals[,2:9], upper.panel=panel.cor)


### Model fitting ###

# Selecionar o diretório onde os modelos devem ser salvos
dir <- choose.dir()

# GLM
for (i in 1:length(unique(pts1[,1]))) {
#family = binomial(link = "logit"); family = gaussian(link = "identity"); family = poisson(link = "log")
  a<-subset(pts1, pts1[,1]==unique(pts1[,1])[i], c(2,3))
  presvals <- data.frame(sp=rep(1,nrow(a)),extract(predictors, a)) 
  sdmdata <- as.data.frame(rbind(presvals, pseudovals)); colnames(sdmdata)[1]<-'pb'
  mod <- glm(formula=pb~., family=binomial, data=sdmdata)
  adeq <- predict(predictors, mod, ext=ext, progress='text')
  writeRaster(adeq, paste(dir,"\\", unique(pts1[,1])[i],"_GLM",  ".asc",sep='') , format="ascii", overwrite=TRUE, NAflag=-9999)    
  }


