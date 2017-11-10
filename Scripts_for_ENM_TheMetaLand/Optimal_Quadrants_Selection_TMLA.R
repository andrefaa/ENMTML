### This is the R code to conduct ENMs evaluation by the Odds-and-Evens Framework
# In the moment, the code is configured to create the splits only by its latitudinal axis
# Read carefully the Initialization and function parameters and the messages that appear on the console
# Developers: Paulo De Marco Junior & Andre Andrade

Optimal_Quadrants_Selection_TMLA <- function(env,occ,band){

  #Parameters
    #env: Predictors
    #occ: Occurrence List
    #band: Longitudinal(1) or Latitudinal(2) bands
  
  #Development
  res<-NULL

  #Separate data by groups
  occ <- lapply(occ, function(x) cbind(x,rep(0,nrow(x))))
  colnames <- c("x","y","Seg")
  occ <- lapply(occ, setNames, colnames)
  occ <- lapply(occ, function(x) round(x, digits=5))
  
  for(x in 1:length(occ)){
    opt <- NULL
    print(names(occ)[x])
    occ.s <- occ[[x]]

    for (quad in seq(2,20,2)){
      
      occ.st <- occ.s
      print(paste(quad,"Quadrants",sep=" "))
      
      #Bands----
      
      for (segm in 1:quad){
        axfin <- min(occ.st[,band])+((max(occ.st[,band])-min(occ.st[,band]))/quad)*segm
        axin <- min(occ.st[,band])+((max(occ.st[,band])-min(occ.st[,band]))/quad)*(segm-1)
        
        occ.st[occ.st[,band]>=axin & occ.st[,band]<=axfin, "Seg"] <- segm
      }
    
      occ.st <- cbind(occ.st,ifelse((occ.st$Seg/2)%%1,1,2))
      colnames(occ.st) <- c("x","y","Seg","Quad")
  
      #Moran's I----
      Moran<-Moran_for_Quadrants_Pair_TMLA(occ.st,env[[1]],quad)
      Moran <- data.frame(cbind(quad, Moran))
        
      #MESS----
      occ_e <- cbind(occ.st,extract(env, occ.st[,1:2]))
      occ_e <- split(occ_e, f=occ_e$Quad)
    
      mess <- MESS(occ_e[[1]][,-c(1:4)], occ_e[[2]][,-c(1:4)])
      mess <- mean(mess$TOTAL, na.rm = TRUE)
      
      #SD of number of records per band----
      #Sd <- sd(table(occ.st[occ.st$Quad == 1, ])) /mean(table(occ.st[occ.st$Quad == 1,]))
        
      res.t<-data.frame(cbind(Moran,mess))
      colnames(res.t) <- c("Quad","Moran","MESS")
      opt <- rbind(opt,res.t)
    }
    names(opt) <- names(res.t)
    
    # SELLECTION OF THE BEST NUMBER OF BANDS----
    Opt <- opt
    while (nrow(Opt) > 1) {
      # I MORAN
      if (nrow(Opt) == 1) break
      Opt <- Opt[which(Opt$Moran <= summary(Opt$Moran)[2]),]
      if (nrow(Opt) == 1) break
      # MESS
      Opt <- opt[which(Opt$MESS >= summary(Opt$MESS)[5]),]
      if (nrow(Opt) == 1) break
      # SD
      # Opt <- Opt[which(Opt$Sd <= summary(Opt$)[2]),]
    }

    #Create the ASC file to be used as mask
    
      msk<-env[[1]]
      msk[!is.na(msk[,])] <- 0

      quad <- Opt$Quad
        
      for (segm in 1:Opt$Quad){
        axfin <- min(occ.s[,band])+((max(occ.s[,band])-min(occ.s[,band]))/quad)*segm
        axin <- min(occ.s[,band])+((max(occ.s[,band])-min(occ.s[,band]))/quad)*(segm-1)
        
        occ.s[occ.s[,band]>=axin & occ.s[,band]<=axfin, "Seg"] <- segm
        
        if(band==1){
          if (((segm/2)-round(segm/2))!=0){
            y1<-ncol(msk)-floor((axfin-msk@extent[1])/xres(msk))
            y0<-ncol(msk)-floor((axin-msk@extent[1])/xres(msk))
            
            for (y in y1:y0){
              msk[,y] <- 1;
            }
            msk[env[[1]]%in%NA]<-NA
            
          }else{
            y1<-ncol(msk)-floor((axfin-msk@extent[1])/xres(msk))
            y0<-ncol(msk)-floor((axin-msk@extent[1])/xres(msk))
            
            for (y in y1:y0){
              msk[,y] <- 2;
            }
            msk[env[[1]]%in%NA] <- NA
          }
        }
        
        if(band==2){
          if (((segm/2)-round(segm/2))!=0){
            y1<-nrow(msk)-floor((axfin-msk@extent[3])/yres(msk))
            y0<-nrow(msk)-floor((axin-msk@extent[3])/yres(msk))
          
            for (y in y1:y0){
              msk[y,] <- 1;
            }
            msk[env[[1]]%in%NA]<-NA

          }else{
            y1<-nrow(msk)-floor((axfin-msk@extent[3])/yres(msk))
            y0<-nrow(msk)-floor((axin-msk@extent[3])/yres(msk))
          
            for (y in y1:y0){
              msk[y,] <- 2;
            }
            msk[env[[1]]%in%NA]<-NA
          }
        }
      }
      occ.s <- cbind(occ.s,ifelse((occ.s$Seg/2)%%1,1,2))
      occ.s <- cbind(rep(names(occ)[x],nrow(occ.s)),occ.s)
      colnames(occ.s) <- c("sp","x","y","Seg","Quad")
        
      msk[msk%in%0] <-  NA      #ASC with the ODD-EVEN quadrants
      writeRaster(msk,paste(getwd(),paste(names(occ)[x],".tif",sep=""),sep="\\"),format="GTiff",overwrite=T)
      
      res<-rbind(res,occ.s)
  }
    write.table(res,paste(getwd(),"OccBands.txt",sep="\\"),sep="\t",row.names=F)
    return(occ_t)
}