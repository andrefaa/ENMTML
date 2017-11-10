#Written by Andre Andrade

BandsPartition_TMLA <- function(Env,
                      Occ,
                      Band,
                      FromTo = FromTo,
                      pseudoabsencesMethod = PabM,
                      PrAbRatio = PabR,
                      DirSave = dirO){

  #Parameters
    #Env: Predictors
    #Occ: Occurrence List
    #Band: Longitudinal(1) or Latitudinal(2) bands
  
  #Development
  res<-NULL

  #Separate data by groups
  Occ <- lapply(Occ, function(x) cbind(x,rep(0,nrow(x))))
  colnames <- c("x","y","Seg")
  Occ <- lapply(Occ, setNames, colnames)
  Occ <- lapply(Occ, function(x) round(x, digits=5))
  
  for(x in FromTo:length(Occ)){
    opt <- NULL
    print(names(Occ)[x])
    Occ.s <- Occ[[x]]

    for (quad in seq(2,20,2)){
      
      Occ.st <- Occ.s
      print(paste(quad,"Quadrants",sep=" "))
      
      #Bands----
      
      for (segm in 1:quad){
        axfin <- min(Occ.st[,Band])+((max(Occ.st[,Band])-min(Occ.st[,Band]))/quad)*segm
        axin <- min(Occ.st[,Band])+((max(Occ.st[,Band])-min(Occ.st[,Band]))/quad)*(segm-1)
        
        Occ.st[Occ.st[,Band]>=axin & Occ.st[,Band]<=axfin, "Seg"] <- segm
      }
    
      Occ.st <- cbind(Occ.st,ifelse((Occ.st$Seg/2)%%1,1,2))
      colnames(Occ.st) <- c("x","y","Seg","Partition")
  
      #Moran's I----
      Moran<-Moran_for_Quadrants_Pair_TMLA(Occ.st,Env[[1]],quad)
      Moran <- data.frame(cbind(quad, Moran))
        
      #MESS----
      Occ_e <- cbind(Occ.st,extract(Env, Occ.st[,1:2]))
      Occ_e <- split(Occ_e, f=Occ_e$Partition)
    
      mess <- MESS(Occ_e[[1]][,-c(1:4)], Occ_e[[2]][,-c(1:4)])
      mess <- mean(mess$TOTAL, na.rm = TRUE)
      
      #SD of number of records per Band----
      #Sd <- sd(table(Occ.st[Occ.st$Quad == 1, ])) /mean(table(Occ.st[Occ.st$Quad == 1,]))
        
      res.t<-data.frame(cbind(Moran,mess))
      colnames(res.t) <- c("Partition","Moran","MESS")
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

    #Create Bands Mask
    
      msk<-Env[[1]]
      msk[!is.na(msk[,])] <- 0

      quad <- Opt$Partition
        
      for (segm in 1:quad){
        axfin <- min(Occ.s[,Band])+((max(Occ.s[,Band])-min(Occ.s[,Band]))/quad)*segm
        axin <- min(Occ.s[,Band])+((max(Occ.s[,Band])-min(Occ.s[,Band]))/quad)*(segm-1)
        
        Occ.s[Occ.s[,Band]>=axin & Occ.s[,Band]<=axfin, "Seg"] <- segm
        
        if(Band==1){
          if (((segm/2)-round(segm/2))!=0){
            y1<-ncol(msk)-floor((axfin-msk@extent[1])/xres(msk))
            y0<-ncol(msk)-floor((axin-msk@extent[1])/xres(msk))
            
            for (y in y1:y0){
              msk[,y] <- 1;
            }
            msk[Env[[1]]%in%NA]<-NA
            
          }else{
            y1<-ncol(msk)-floor((axfin-msk@extent[1])/xres(msk))
            y0<-ncol(msk)-floor((axin-msk@extent[1])/xres(msk))
            
            for (y in y1:y0){
              msk[,y] <- 2;
            }
            msk[Env[[1]]%in%NA] <- NA
          }
        }
        
        if(Band==2){
          if (((segm/2)-round(segm/2))!=0){
            y1<-nrow(msk)-floor((axfin-msk@extent[3])/yres(msk))
            y0<-nrow(msk)-floor((axin-msk@extent[3])/yres(msk))
          
            for (y in y1:y0){
              msk[y,] <- 1;
            }
            msk[Env[[1]]%in%NA]<-NA

          }else{
            y1<-nrow(msk)-floor((axfin-msk@extent[3])/yres(msk))
            y0<-nrow(msk)-floor((axin-msk@extent[3])/yres(msk))
          
            for (y in y1:y0){
              msk[y,] <- 2;
            }
            msk[Env[[1]]%in%NA]<-NA
          }
        }
      }
      Occ.s <- cbind(Occ.s,ifelse((Occ.s$Seg/2)%%1,1,2))
      Occ.s <- cbind(rep(names(Occ)[x],nrow(Occ.s)),Occ.s)
      colnames(Occ.s) <- c("sp","x","y","Seg","Partition")
      Occ.s <- Occ.s[,c("sp","x","y","Partition")]
      Occ.s <- cbind(Occ.s, rep(1,nrow(Occ.s)))
      colnames(Occ.s) <- c("sp","x","y","Partition","PresAbse")
      
      msk[msk%in%0] <-  NA      #ASC with the ODD-EVEN quadrants
      writeRaster(msk,paste(DirSave,paste(names(Occ)[x],".tif",sep=""),sep="\\"),
                  format="GTiff",NAflag = -9999,overwrite=T)
      
    # Pseudoabsences allocation-----
      
      # Random-----
      if(pseudoabsencesMethod=="rnd"){
        # Clip the mask raster to generate rando pseudoabsences
        pseudo.mask <- msk
        pseudo.mask2 <- list()

        for(i in 1:2){
          mask3 <- pseudo.mask
          mask3[!mask3[]==i] <- NA 
          pseudo.mask2[[i]] <- mask3
        }
        
        pseudo.mask <- brick(pseudo.mask2)
        rm(pseudo.mask2)

        # Random allocation of pseudoabsences 
        absences <- list()
        for (i in 1:2) {
          absences.0 <- randomPoints(pseudo.mask[[i]], (1 / PrAbRatio)*sum(Occ.s[,"Partition"]==i),
                                     ext = extent(pseudo.mask[[i]]),
                                     prob = FALSE)
          colnames(absences.0) <- c(x, y)
          absences[[i]] <- as.data.frame(absences.0)
          rm(absences.0)
        }
        for (i in 1:length(absences)){
          absences[[i]] <- cbind(absences[[i]], rep(i,nrow(absences[[i]])))  
        }
        absences <- lapply(absences, function(x) cbind(x, rep(0,nrow(x))))
        absences <- lapply(absences, function(y) cbind(rep(names(Occ)[x],nrow(y)),y))
        absences <- ldply(absences, data.frame)
        colnames(absences) <- colnames(Occ.s)
      }
      
      if(pseudoabsencesMethod=="constrain"){

        Model <- bioclim(Env, Occ.s[,c("x","y")])
        pseudo.mask <- dismo::predict(Model, Env, ext=extent(msk))
        names(pseudo.mask) <- "Group"
        pseudo.mask <- round(pseudo.mask, 5)
        pseudo.mask <-(pseudo.mask-minValue(pseudo.mask))/
          (maxValue(pseudo.mask)-minValue(pseudo.mask))
        pseudo.mask <-(1-pseudo.mask)>=0.99 #environmental constrain
        pseudo.mask[which(pseudo.mask[,]==FALSE)] <- NA
        
        # Split the raster of environmental layer with grids
        pseudo.mask2 <- list()

        for(i in 1:2){
          mask3 <- pseudo.mask
          mask3[!msk[]==i] <- NA 
          mask3[is.na(msk[])] <- NA 
          pseudo.mask2[[i]] <- mask3
        }
        pseudo.mask <- brick(pseudo.mask2)
        rm(pseudo.mask2)
        
        absences <- list()
        for (i in 1:2) {
          absences.0 <- randomPoints(pseudo.mask[[i]], (1 / PrAbRatio)*sum(Occ.s[,"Partition"]==i),
                                     ext = extent(pseudo.mask[[i]]),
                                     prob = FALSE)
          colnames(absences.0) <- c(x, y)
          absences[[i]] <- as.data.frame(absences.0)
          rm(absences.0)
        }
        for (i in 1:length(absences)){
          absences[[i]] <- cbind(absences[[i]], rep(i,nrow(absences[[i]])))  
        }
        absences <- lapply(absences, function(x) cbind(x, rep(0,nrow(x))))
        absences <- lapply(absences, function(y) cbind(rep(names(Occ)[x],nrow(y)),y))
        absences <- ldply(absences, data.frame)
        colnames(absences) <- colnames(Occ.s)
      }
      Occ.s <- rbind(Occ.s,absences)
      res<-rbind(res,Occ.s)
  }
    occ <- res[res[,"PresAbse"]==1,]
    write.table(occ,paste(DirSave,"OccBands.txt",sep="\\"),sep="\t",row.names=F)
    rm(occ)
    return(res)
}
