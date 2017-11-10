### This is the R code to conduct ENMs evaluation by the Odds-and-Evens Framework
# In the moment, the code is configured to create the splits only by its latitudinal axis
# Read carefully the Initialization and function parameters and the messages that appear on the console
# Developers: Paulo De Marco Junior & André Andrade

Odds.and.evens.Framework <- function(dir_input,dir_output,quad){
  #Function Parameters:
    #dir_input: The directory containing an ASC file that will be used as a mask(the whole desired extension) and the xls file with all occurrences
    #dir_output: The directory where will be saved the ASC files to be used as masks in the modelling proceddure
    #quad: Number of quadrants the area will be divided
  
  #Important Observations:
    #Occurence File MUST have the columns named Species (containing species name),Long (containing longitude),Lat (containing latitude)
  
  #Initialization
  
  require(SDMTools)
  require(raster)
  require(xlsx)
  
  setwd(dir_input)
  
  #Read the ASC file to be used as mask
  setwd(dir_input)
  mask<-read.asc(list.files(dir_input,pattern=".asc"))
  mask<-raster(mask)
  mask[is.na(mask)==F]<-0
  maskE<-mask

  #Import the xls files with species occurrences
  occ<-read.xlsx(list.files(dir_input,pattern=".xls"),sheetIndex = 1)
  occ<-occ[,which(names(occ) %in% c("Species","Long","Lat"))]
  occ.xy<-occ[,-1]
  
  #Identify unique species
  ssp<-unique(occ$Species)
  
  #Remove duplicate presences
  occ.var<-extract(mask,occ.xy,cellnumber=T)
  occ.var<-cbind(occ,occ.var)
  occ.var<-na.omit(occ.var)
    
  occ.unique<-NULL

  for (a in 1:length(ssp)){
      
    occ.ssp<-occ.var[occ.var[,1]==ssp[a],]
    duplicate<-which(duplicated(occ.ssp[,'cells'])==T)

    if (length(duplicate)!=0){
      occ.unique<-rbind(occ.unique,occ.ssp[-duplicate,])
    }else{
      occ.unique<-rbind(occ.unique,occ.ssp)
    }    
  }
  
  occ.unique<-occ.unique[,-4]
  
  #Development
  
  for (a in 1:length(ssp)){
    occ.ssp<-occ.unique[occ.unique[,1]==ssp[a],]
    n.quad<-rep(2,nrow(occ.ssp))
    occ.ssp<-cbind(occ.ssp,n.quad)
    
    #Re-read the ASC file to be used as mask
    setwd(dir_input)
    mask<-read.asc(list.files(dir_input,pattern=".asc"))
    mask<-raster(mask)
    mask[is.na(mask)==F]<-0
    maskE<-mask

    for (segm in 1:quad){
      axfin<-min(occ.ssp[,3])+((max(occ.ssp[,3])-min(occ.ssp[,3]))/quad)*segm
      axin<-min(occ.ssp[,3])+((max(occ.ssp[,3])-min(occ.ssp[,3]))/quad)*(segm-1)
      
      if (((segm/2)-round(segm/2))!=0){
        y1<-nrow(mask)-floor((axfin-mask@extent[3])/yres(mask))
        y0<-nrow(mask)-floor((axin-mask@extent[3])/yres(mask))
        
        for (y in y1:y0){
          print(y)
          for (cc in 1:ncol(mask)){
            if (is.na(mask[y,cc])==F){
              mask[y,cc]=1;
            }
          }
        }
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
          print(y)
          for (cc in 1:ncol(maskE)){
            if (is.na(maskE[y,cc])==F){
              maskE[y,cc]=1;
            }
          }
        }
      }
      occ.ssp[occ.ssp[,3]<=axfin & occ.ssp[,3]>=axin,5]<-segm
    }
    
    mask[mask%in%0]<- NA
    maskE[maskE%in%0]<- NA
    
    setwd(dir_output)
    writeRaster(mask,paste(ssp[a],"_evens.asc",sep=""),format="ascii")
    writeRaster(maskE,paste(ssp[a],"_odds.asc",sep=""),format="ascii")
    
    #Occurrence data of the even group
    occ.even<-occ.ssp[occ.ssp[,4]==1,1:3];
    write.table(occ.even,paste(ssp[a],"evens.txt",sep=""),row.names = F, sep="\t")
    
    #Occurrence data of the odds group
    occ.odds<-occ.ssp[occ.ssp[,4]==0,1:3];
    write.table(occ.odds,paste(ssp[a],"odds.txt",sep=""),row.names = F, sep="\t")
  }
}

Odds.and.evens.Framework("D:\\Users\\Andre\\Google Drive\\Mestrado\\Capitulo 1(Transferabilidade)\\Ocorrencias_20\\Conversao Script Matlab_R","D:\\Users\\Andre\\Google Drive\\Mestrado\\Capitulo 1(Transferabilidade)\\Ocorrencias_20\\Conversao Script Matlab_R\\TesteR",2)
