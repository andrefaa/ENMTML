### This is the R code to conduct ENMs evaluation by the Odds-and-Evens Framework
# In the moment, the code is configured to create the splits only by its latitudinal axis
# Read carefully the Initialization and function parameters and the messages that appear on the console
# Developers: Paulo De Marco Junior & Andre Andrade

Optimal_Quadrants_Selection_TML <- function(env,occ){
  
  #Function Parameters:
  #env:environmental variables stack
  #occ: occurrence data

  #Important Observations:
  #Occurence File MUST have the columns named Species (containing species name),Long (containing longitude),Lat (containing latitude)
  
  #Initialization
  
  library(SDMTools)
  library(raster)
  library(xlsx)
  library(biomod2)
  library(ade4)
  library(adehabitat)
  library(sp)
  library(gam)
  library(MASS)
  library(mvtnorm)
  library(gbm)
  library(dismo)
  library(ape)

  #Development
  res<-NULL
  
  dir.out<-"Optimal"
  if (file.exists(dir.out)){
  } else {
      dir.create(file.path(getwd(),dir.out))
  }
  
  ssp<-unique(occ[,1])
  
  for (a in 1:length(ssp)){
    tryCatch({
      
    Moran.Quad<-NULL
    Quadrant<-NULL
    Schoener<-NULL
    
    print(as.character(ssp[a]))
    occ.ssp<-occ[occ[,1]==ssp[a],]
    n.quad<-rep(2,nrow(occ.ssp))
    occ.ssp<-cbind(occ.ssp,n.quad)
    
  for (quad in seq(2,20,2)){
    print(paste(quad,"Quadrants",sep=" "))
    
    #Test different ASC file to be used as mask
    mask<-env[[1]]
    mask[is.na(mask)==F]<-0
    maskE<-mask
    
    for (segm in 1:quad){
      axfin<-min(occ.ssp[,3])+((max(occ.ssp[,3])-min(occ.ssp[,3]))/quad)*segm
      axin<-min(occ.ssp[,3])+((max(occ.ssp[,3])-min(occ.ssp[,3]))/quad)*(segm-1)
      
      if (((segm/2)-round(segm/2))!=0){
        y1<-nrow(mask)-floor((axfin-mask@extent[3])/yres(mask))
        y0<-nrow(mask)-floor((axin-mask@extent[3])/yres(mask))
        
        for (y in y1:y0){
              mask[y,]=1;
        }
        mask[env[[1]]%in%NA]<-NA
        
        if (segm==1){
          occ.ssp[occ.ssp[,3]<=axfin & occ.ssp[,3]>=axin,4]<-1
          occ.ssp[occ.ssp[,3]<=axfin & occ.ssp[,3]>=axin,5]<-segm
        }else{
          occ.ssp[occ.ssp[,3]<=axfin & occ.ssp[,3]>=axin,4]<-1
          occ.ssp[occ.ssp[,3]<=axfin & occ.ssp[,3]>=axin,5]<-segm
        }
      }else{
        y1<-nrow(mask)-floor((axfin-mask@extent[3])/yres(mask))
        y0<-nrow(mask)-floor((axin-mask@extent[3])/yres(mask))
        
        for (y in y1:y0){
              maskE[y,]=1;
        }
        maskE[env[[1]]%in%NA]<-NA
      }
      occ.ssp[occ.ssp[,3]<=axfin & occ.ssp[,3]>=axin,5]<-segm
    }
    
    mask[mask%in%0]<- NA      #ASC with the even quadrants
    maskE[maskE%in%0]<- NA    #ASC with the odd quadrants

    #Occurrence Data for Moran's I
    occ.ssp<-na.omit(occ.ssp[,c(1:3,5)])
    
    #Calculate Moran's I
    Moran<-Moran_for_Quadrants_Pair(occ.ssp,env[[1]])
    Moran.Quad<-c(Moran.Quad,Moran)
    Quadrant<-c(Quadrant,quad)

    #Calculate Schoener's D
    
    pcs1<-env
    pcs1[maskE%in%NA]<- NA

    pcs2<-env
    pcs2[mask%in%NA]<- NA

    clim1<-rasterToPoints(pcs1)
    colnames(clim1)<-c("long","lat","pc1","pc2")
    clim2<-rasterToPoints(pcs2)
    colnames(clim2)<-c("long","lat","pc1","pc2")
    clim12<-rbind(clim1,clim2)
    
    occ.sp1<-na.exclude(cbind(occ.ssp[,2:3],extract(pcs1,occ.ssp[,2:3])))
    colnames(occ.sp1)<-colnames(clim1)
    occ.sp2<-na.exclude(cbind(occ.ssp[,2:3],extract(pcs2,occ.ssp[,2:3])))
    colnames(occ.sp2)<-colnames(clim1)
    
    if (nrow(occ.sp1)<5 || nrow(occ.sp2)<5){
      D<-NA
      
    }else{
    
    # selection of the type of analysis.
    # If PROJ =F, the models are calibrated on both ranges.
    # If PROJ =T, the models are calibrated on species 1 range only and projected to range 2. 
    # Analyses where both ranges are needed (ex: LDA) are not done
    PROJ = T
    
    # selection of variables to include in the analyses
    colnames(clim12)
    Xvar<-c(3:4)
    nvar<-length(Xvar)
    
    #number of interation for the tests of equivalency and similarity
    iterations<-100
    
    #resolution of the gridding of the climate space
    R=1000
    
    # if PROJ = T
    row.w.occ.PROJT<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
    row.w.env.PROJT<-c(rep(1, nrow(clim1)),rep(0, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
    
    # global dataset for the analysis and rows for each sub dataset
    data.env.occ<-rbind(clim1,clim2,occ.sp1,occ.sp2)[Xvar]
    row.clim1<-1:nrow(clim1)
    row.clim2<-(nrow(clim1)+1):(nrow(clim1)+nrow(clim2))
    row.clim12<-1:(nrow(clim1)+nrow(clim2))
    row.sp1<-(nrow(clim1)+nrow(clim2)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1))
    row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+nrow(occ.sp2))
    
    #################################################################################################
    #################################### PCA-ENV ####################################################
    #################################################################################################
    
    # measures niche overlap along the two first axes of a PCA calibrated on all the pixels of the study areas
    
    if(PROJ == T){	#fit of the analyse using occurences from range 1		
      pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env.PROJT, center = T, scale = T, scannf = F, nf = 2)
    }
    # predict the scores on the axes
    scores.clim12<- pca.cal$li[row.clim12,]
    scores.clim1<- pca.cal$li[row.clim1,]
    scores.clim2<- pca.cal$li[row.clim2,]
    scores.sp1<- pca.cal$li[row.sp1,]
    scores.sp2<- pca.cal$li[row.sp2,]
    
    # calculation of occurence density and test of niche equivalency and similarity 
    z1<- grid.clim(scores.clim12,scores.clim1,scores.sp1,R)
    z2<- grid.clim(scores.clim12,scores.clim2,scores.sp2,R)
    D<-round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),3)
    }
    Schoener<-c(Schoener,D)
  }
    tst<-as.data.frame(cbind(Quadrant,Moran.Quad,Schoener))
    sp_col<-rep(as.character(ssp[a]),nrow(tst))
    sp_res<-cbind(sp_col,tst)
    
    rat<-(1-sqrt((sp_res[,3]-0)^2+(sp_res[,4]-1)^2))*((abs(sp_res[,3]))<0.11)
    sp_res.t<-cbind(sp_res,rat)
    max<-which(sp_res.t[,5]==max(na.omit(sp_res.t[1:nrow(sp_res.t),5])))
    opt<-as.numeric(sp_res.t[max,2])

    #Create the ASC file to be used as mask
      mask<-env[[1]]
      mask[is.na(mask)==F]<-0
      maskE<-mask
      
      for (segm.f in 1:opt){
        axfin<-min(occ.ssp[,3])+((max(occ.ssp[,3])-min(occ.ssp[,3]))/opt)*segm.f
        axin<-min(occ.ssp[,3])+((max(occ.ssp[,3])-min(occ.ssp[,3]))/opt)*(segm.f-1)
        
        if (((segm.f/2)-round(segm.f/2))!=0){
          y1<-nrow(mask)-floor((axfin-mask@extent[3])/yres(mask))
          y0<-nrow(mask)-floor((axin-mask@extent[3])/yres(mask))
          
          for (y in y1:y0){
            mask[y,]=1;
          }
          mask[env[[1]]%in%NA]<-NA
          
          if (segm.f==1){
            occ.ssp[occ.ssp[,3]<=axfin & occ.ssp[,3]>=axin,4]<-1
            occ.ssp[occ.ssp[,3]<=axfin & occ.ssp[,3]>=axin,5]<-segm.f
          }else{
            occ.ssp[occ.ssp[,3]<=axfin & occ.ssp[,3]>=axin,4]<-1
            occ.ssp[occ.ssp[,3]<=axfin & occ.ssp[,3]>=axin,5]<-segm.f
          }
        }else{
          y1<-nrow(mask)-floor((axfin-mask@extent[3])/yres(mask))
          y0<-nrow(mask)-floor((axin-mask@extent[3])/yres(mask))
          
          for (y in y1:y0){
            maskE[y,]=1;
          }
          maskE[env[[1]]%in%NA]<-NA
        }
        occ.ssp[occ.ssp[,3]<=axfin & occ.ssp[,3]>=axin,5]<-segm.f
      }
      
      mask[mask%in%0]<- NA      #ASC with the even quadrants
      maskE[maskE%in%0]<- NA    #ASC with the odd quadrants
      
      writeRaster(mask,paste(getwd(),dir.out,paste(ssp[a],"ODD.asc",sep=""),sep="\\"),format="ascii",overwrite=T)
      writeRaster(maskE,paste(getwd(),dir.out,paste(ssp[a],"EVEN.asc",sep=""),sep="\\"),format="ascii",overwrite=T)
    
    res<-rbind(res,sp_res)
    print(res)
    }, error=function(e){cat("ERROR : Deu ruim aqui!")})
  }
  write.table(res,paste(getwd(),dir.out,"Moran&Schoener_Quadrantes.txt",sep="\\"),sep="\t",row.names=F)
  dir.opt<-file.path(getwd(),dir.out)
  return(dir.opt)
}
