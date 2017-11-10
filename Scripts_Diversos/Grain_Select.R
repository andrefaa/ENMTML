Grain_Select <- function(RecordsData,
                         Variables){
  
  #1.Load occurrence data----
  RecordsData <- read.table(file.choose(),sep="\t",header=T)
  OccXY <- RecordsData[,2:3]
  OccXY2 <- xyFromCell(Env,cellFromXY(Env,OccXY))
  
  #2.Load Environmental Variables
  DirE <- "C:\\OneDrive\\Trabajos\\08 - Grain Size\\50km"
  Env <- brick(stack(file.path(DirE,list.files(DirE,pattern = ".tif"))))
  e <- extent(-120,-20,-60,30)
  Env <- crop(Env,e)

  #3.Extract EnvData & Cellsize
  OccE <- extract(Env,RecordsDataXY)
  OccE <- split(OccE, seq(nrow(OccE)))
  OccC <- cellFromXY(Env,RecordsDataXY)
  OccC <- split(OccC, seq(length(OccC)))
  
  Neig <- lapply(OccC, function (x) adjacent(Env, x, directions=8, pairs=FALSE, target=NULL, sorted=FALSE, 
         include=TRUE, id=FALSE))
  Neig <- lapply(Neig, function(x) xyFromCell(Env,x))
  
  #
  
  
  #Calculate Distances
  Dists <- list()
  for(x in 1:nrow(RecordsDataXY)){
    Dists[[x]] <- sqrt(((OccXY2[x,1]-Neig[[x]][,1])^2)+((OccXY2[x,2]-Neig[[x]][,2])^2))
  }
  Dists <- unlist(Dists)
  Dists <- ldply(Dists,data.frame,.id=NULL)

focal(Env,w=matrix(1,nrow=3,ncol=3),)
