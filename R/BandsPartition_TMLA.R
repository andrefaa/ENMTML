#Written by Andre Andrade

BandsPartition_TMLA <- function(evnVariables,
                                RecordsData,
                                N,
                                pseudoabsencesMethod = PabM,
                                PrAbRatio = PabR,
                                DirSave = DirB,
                                DirM=DirM,
                                MRst=MRst,
                                type=TipoMoran){

  #Parameters
    #evnVariables: Predictors
    #RecordsData: Occurrence List
    #N: Longitudinal(1) or Latitudinal(2) bands
  
  #Development
  res<-NULL
  resOpt <- list()

  #Separate data by groups
  RecordsData <- lapply(RecordsData, function(x) cbind(x,rep(0,nrow(x))))
  colnames <- c("x","y","Seg")
  RecordsData <- lapply(RecordsData, setNames, colnames)
  RecordsData <- lapply(RecordsData, function(x) round(x, digits=5))
  
  grid <- seq(2,20,2)
  
  #Start species loop----
  for(x in 1:length(RecordsData)){
    opt <- NULL
    print(names(RecordsData)[x])
    RecordsData.s <- RecordsData[[x]]

    for (quad in seq(2,20,2)){
      
      RecordsData.st <- RecordsData.s
      print(paste(quad,"Quadrants",sep=" "))
      
      #Bands----
      
      for (segm in 1:quad){
        axfin <- min(RecordsData.st[,N])+((max(RecordsData.st[,N])-min(RecordsData.st[,N]))/quad)*segm
        axin <- min(RecordsData.st[,N])+((max(RecordsData.st[,N])-min(RecordsData.st[,N]))/quad)*(segm-1)
        
        RecordsData.st[RecordsData.st[,N]>=axin & RecordsData.st[,N]<=axfin, "Seg"] <- segm
      }
    
      RecordsData.st <- cbind(RecordsData.st,ifelse((RecordsData.st$Seg/2)%%1,1,2))
      colnames(RecordsData.st) <- c("x","y","Seg","Partition")
  
      #Moran's I----
      Moran<-Moran_for_Quadrants_Pair_TMLA(occ=RecordsData.st,pc1=evnVariables[[1]],
                                           quad=quad,type=type)
      Moran <- data.frame(cbind(quad, Moran))
        
      #MESS----
      RecordsData_e <- cbind(RecordsData.st,extract(evnVariables, RecordsData.st[,1:2]))
      RecordsData_e <- split(RecordsData_e, f=RecordsData_e$Partition)
    
      mess <- MESS(RecordsData_e[[1]][,-c(1:4)], RecordsData_e[[2]][,-c(1:4)])
      mess <- mean(mess$TOTAL, na.rm = TRUE)
      
      #SD of number of records per Band----
      pa <- RecordsData.st$Partition # Vector wiht presences and absences
      Sd <- sd(table(pa)) /
          mean(table(pa))
        
      #Combine calculations in a data frame
      res.t<-data.frame(cbind(Moran,mess,Sd))
      colnames(res.t) <- c("Partition","Moran","MESS","SD")
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
       Opt <- Opt[which(Opt$Sd <= summary(Opt$Sd)[2]),]
    }
    resOpt[[x]] <- cbind(names(RecordsData)[x],Opt)

    #Create Bands Mask
    
      msk<-evnVariables[[1]]
      msk[!is.na(msk[,])] <- 0

      quad <- Opt$Partition
        
      for (segm in 1:quad){
        axfin <- min(RecordsData.s[,N])+((max(RecordsData.s[,N])-min(RecordsData.s[,N]))/quad)*segm
        axin <- min(RecordsData.s[,N])+((max(RecordsData.s[,N])-min(RecordsData.s[,N]))/quad)*(segm-1)
        
        RecordsData.s[RecordsData.s[,N]>=axin & RecordsData.s[,N]<=axfin, "Seg"] <- segm
        
        if(N==1){
          if (((segm/2)-round(segm/2))!=0){
            y1<-ncol(msk)-floor((axfin-msk@extent[1])/xres(msk))
            y0<-ncol(msk)-floor((axin-msk@extent[1])/xres(msk))
            
            for (y in y1:y0){
              msk[,y] <- 1;
            }
            msk[evnVariables[[1]]%in%NA]<-NA
            
          }else{
            y1<-ncol(msk)-floor((axfin-msk@extent[1])/xres(msk))
            y0<-ncol(msk)-floor((axin-msk@extent[1])/xres(msk))
            
            for (y in y1:y0){
              msk[,y] <- 2;
            }
            msk[evnVariables[[1]]%in%NA] <- NA
          }
        }
        
        if(N==2){
          if (((segm/2)-round(segm/2))!=0){
            y1<-nrow(msk)-floor((axfin-msk@extent[3])/yres(msk))
            y0<-nrow(msk)-floor((axin-msk@extent[3])/yres(msk))
          
            for (y in y1:y0){
              msk[y,] <- 1;
            }
            msk[evnVariables[[1]]%in%NA]<-NA

          }else{
            y1<-nrow(msk)-floor((axfin-msk@extent[3])/yres(msk))
            y0<-nrow(msk)-floor((axin-msk@extent[3])/yres(msk))
          
            for (y in y1:y0){
              msk[y,] <- 2;
            }
            msk[evnVariables[[1]]%in%NA]<-NA
          }
        }
      }
      RecordsData.s <- cbind(RecordsData.s,ifelse((RecordsData.s$Seg/2)%%1,1,2))
      RecordsData.s <- cbind(rep(names(RecordsData)[x],nrow(RecordsData.s)),RecordsData.s)
      colnames(RecordsData.s) <- c("sp","x","y","Seg","Partition")
      RecordsData.s <- RecordsData.s[,c("sp","x","y","Partition")]
      RecordsData.s <- cbind(RecordsData.s, rep(1,nrow(RecordsData.s)))
      colnames(RecordsData.s) <- c("sp","x","y","Partition","PresAbse")
      
      msk[msk%in%0] <-  NA      #ASC with the ODD-EVEN quadrants
      writeRaster(msk,paste(DirSave,paste(names(RecordsData)[x],".tif",sep=""),sep="\\"),
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
          set.seed(x)
          if(MRst=="Y"){
            SpMask <- raster(file.path(DirM,paste0(names(RecordsData)[x],".tif")))
            pseudo.mask[[i]] <- pseudo.mask[[i]]*SpMask
            if(sum(is.na(SpMask[])==F)<(PrAbRatio*nrow(RecordsData[[x]]))){
              warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences") 
              stop("Please try again with another restriction type or without restricting the extent")
            }
          }
          absences.0 <- randomPoints(pseudo.mask[[i]], (1 / PrAbRatio)*sum(RecordsData.s[,"Partition"]==i),
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
        absences <- lapply(absences, function(y) cbind(rep(names(RecordsData)[x],nrow(y)),y))
        absences <- ldply(absences, data.frame)
        colnames(absences) <- colnames(RecordsData.s)
      }
      
      if(pseudoabsencesMethod=="const"){

        Model <- bioclim(evnVariables, RecordsData.s[,c("x","y")])
        pseudo.mask <- dismo::predict(Model, evnVariables, ext=extent(msk))
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
          set.seed(x)
          if(MRst=="Y"){
            SpMask <- raster(file.path(DirM,paste0(names(RecordsData)[x],".tif")))
            pseudo.mask[[i]] <- pseudo.mask[[i]]*SpMask
            if(sum(is.na(SpMask[])==F)<(PrAbRatio*nrow(RecordsData[[x]]))){
              warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences") 
              stop("Please try again with another restriction type or without restricting the extent")
            }
          }
          absences.0 <- randomPoints(pseudo.mask[[i]], (1 / PrAbRatio)*sum(RecordsData.s[,"Partition"]==i),
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
        absences <- lapply(absences, function(y) cbind(rep(names(RecordsData)[x],nrow(y)),y))
        absences <- ldply(absences, data.frame)
        colnames(absences) <- colnames(RecordsData.s)
      }
      RecordsData.s <- rbind(RecordsData.s,absences)
      res<-rbind(res,RecordsData.s)
  }
    resOpt <- ldply(resOpt,data.frame)
    write.table(resOpt,paste(DirSave,"Band_Moran_MESS.txt",sep="\\"),sep="\t",row.names=F)
    write.table(res,paste(DirSave,"OccBands.txt",sep="\\"),sep="\t",row.names=F)
    return(res)
}
