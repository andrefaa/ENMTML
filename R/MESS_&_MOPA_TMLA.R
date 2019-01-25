MESS_and_MOP <- function(Fut,
                        spOcc,
                        VarCol,
                        ProjDir,
                        PAMethod=""){
  #Function to calculate MESS and MOPA metrics
  #Parameters:
    #envF: Raster stack in which the model will be extrapolated (other region/time period)
    #SpOccE: Species data with environmental information, created by ENM_TMLA
    #VarCol: Name of columns with the predictors information in species data
    #FoldF: Directory of envF variables (same structure as ENM_TMLA)
    #Background: Is the algorithm background-based 
  spN <- unique(spOcc$sp)
  #Initialisation
  for(i in 1:length(ProjDir)){
    if (!file.exists(file.path(ProjDir[i],"MESS"))){
      dir.create(file.path(ProjDir[i],"MESS"))
    }
    for(j in spN){
      spOccS <- spOcc[spOcc$sp==j,]
      if(PAMethod=="PresAbse"){
        MESS <- mess(Fut[[i]],spOccS[spOccS$PresAbse==1,VarCol])
        writeRaster(MESS,file.path(ProjDir[i],"MESS",paste0(j,"_MESS_Presence.tif")),format="GTiff",NAflag=-9999)
        MESS <- mess(Fut[[i]],spOccS[VarCol])
        writeRaster(MESS,file.path(ProjDir[i],"MESS",paste0(j,"_MESS_PresAbse.tif")),format="GTiff",NAflag=-9999)
      }else if(PAMethod=="Background"){
        writeRaster(MESS,file.path(ProjDir[i],"MESS",paste0(j,"_MESS_Background.tif")),format="GTiff",NAflag=-9999)
      }
    }
  }
}
# tes <- rasterToPoints(Fut[[1]])
# MOP_NB(spOccS, tes, 6:10, 3:7, "25", 0.1, 0.1, "x", "y", MxMESS="N")
##		m1 - Reference matrix (M), of variables in columns and 
## pixels in the rows, and Latitude and Longitude.
##		m2 - Extent matrix (G), of variables in columns and pixels 
## in the rows.
##      c1 - Environmental variables for matrix m1 (a vector describing 
## the positions in M of the variables to be used)
##      c2 - Environmental variables for matrix m2 (a vector describing 
## the positions in G of the variables to be used)
##      decil - A number which signifies what proportion of closest   
## points in M is to be compared. 
##      p1 - Subsampling percentage for M (proportion to be randomly
## sampled if M is too large)
##      p2 - Subsampling percentage for G (proportion to be randomly 
## sampled if G is too large)
##      Xcol - column containing the X coordinate (longitude)
##      Ycol - column containing the Y coordinate (latitude)
##      MxMESS - A boolean variable for generating Maxent MESS, 