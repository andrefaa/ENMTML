Moran.Aleat <- function(dir_input,dir_output,m){
  #Parametros:
  
  #dir_input: input directory with xls files
  #dir_output: output directory
  
  #Occurrence data must have columns with the following names: Taxa,Long,Lat
  #Occurrence data must be a txt file
  
  library(SDMTools)
  library(raster)
  library(xlsx)
  
  setwd(dir_input)
  
  #Read the environmental asc file from where env_data will be extracted
  print("Select ASC file to be used to extrac environmental data from:")
  asc.base<-read.asc("C:\\Users\\Andre\\Google Drive\\Mestrado\\Capitulo 1(Transferabilidade)\\Ocorrencias_20\\Input\\calibrate\\PCA\\pca1.asc")
  asc.base<-raster(asc.base)
  
  #Create and divide xls files in odds and evens
  evens<-list.files(pattern="_par.xls")
  odds<-list.files(pattern="_impar.xls")
  
  Proportional.Moran<-NULL
  
  for (a in 1:length(odds)){
    
    #Import occurrence data
    print(substr(odds[a],1,nchar(odds[a])-11))
    
    #First Dataset
    occ1<-read.xlsx(evens[a],sheetIndex=1,header=T)
    
    subset<-rep(1,nrow(occ1))
    occ1<-cbind(occ1,subset)
    
    occ1.xy<-occ1[,which(names(occ1) %in% c("Long","Lat"))]
    occs.env<-extract(asc.base,occ1.xy)
    occ1<-cbind(occ1,occs.env)
    
    
    #Second Dataset
    occ2<-read.xlsx(odds[a],sheetIndex=1,header=T)
    
    subset<-rep(2,nrow(occ2))
    occ2<-cbind(occ2,subset)
    
    occ2.xy<-occ2[,which(names(occ2) %in% c("Long","Lat"))]
    occs.env<-extract(asc.base,occ2.xy)
    occ2<-cbind(occ2,occs.env)    
    
    #Calculate Z values for each observation
    occs<-rbind(occ1,occ2)
    occs.z<-occs$occs.env-mean(occs$occs.env)
    occs<-cbind(occs,occs.z)
    
    #Split matrixes
    occ1<-occs[occs[,4]==1,]
    occ2<-occs[occs[,4]==2,]

    #Calculate distance to nearest point in another quadrant
    Moran.temp<-NULL
    W<-0
    Div<-0
    Imax.temp<-0
        
    for (x in 1:nrow(occ2)){
      
      #Identify the focal test-occurrence
      
      Binary<-rep(0,nrow(occ2))
      occ2.focal<-cbind(occ2,Binary)
      occ2.focal[x,7]<-1
      
      dist.temp<-NULL
      quad.alvo<-NULL
      dist.quad<-NULL
        for (y in 1:nrow(occ1)){
          dist.temp[y]<-sqrt((occ2[x,2]-occ1[y,2])^2+(occ2[x,3]-occ1[y,3])^2)
          quad.alvo[y]<-occs[y,4]
        }
    dist.to<-cbind(occ1,dist.temp)    

    names(dist.to)<-c("Species","Lat","Long","Subset","Occ.Env","Occ.Z","Dist.to.Test X")
    
    Binary<-rep(0,nrow(dist.to))
    dist.to<-cbind(dist.to,Binary)
    
    #Create a Sub-matrix with only x% of the nearest points
    
    dist.to<-dist.to[order(dist.to[,7]),]
    dist.to[1:ceiling(m*nrow(dist.to)),8]<-1
   
      #Calculate Moran's I among subsets
      
      Moran.I<-0
      
      for (e in 1:nrow(occ2.focal)){
        for(b in 1:nrow(dist.to)){
          Moran.I<-Moran.I+((occ2.focal[e,6]*dist.to[b,6])*((occ2.focal[e,7]==1)&&(dist.to[b,8]==1)))
          W<-W+((occ2.focal[e,7]==1)&&(dist.to[b,8]==1))
        }
      }
      Moran.temp[x]<-Moran.I

      #Calculate I(max)
      
      Imax<-0
      
      for (c in 1:nrow(occ2.focal)){
        Imax.temp.sp<-0
        for (d in 1:nrow(dist.to)){
          Imax.temp.sp<-Imax.temp.sp+((dist.to[d,5]-mean(occs[,5]))*((occ2.focal[c,7]==1)&&(dist.to[d,8]==1)))
        }
        Imax<-Imax+(Imax.temp.sp^2)
        Div<-Div+(((occ2.focal[c,5]-mean(occs[,5]))^2)*(occ2.focal[c,7]==1))
      }
    Imax.temp[x]<-Imax

    
    } #Fecha o for, calculando o Moran I e Imax por especie
    
    W<-2*W
    Moran.I<-sum(Moran.temp)
    
    Moran.I<-(nrow(occs)*Moran.I)/(W*sum((occs$occs.z)^2))
    
    Imax<-sum(Imax.temp)
    
    Imax<-(nrow(occs)/W)*sqrt(Imax/Div)
    Imax<-abs(Imax)
    
    #Calculate Moran's I proportional to Imax
    Proportional.Moran[a]<-Moran.I/Imax
    print(Proportional.Moran[a])
    Sys.sleep(1)
  }
    setwd(dir_output)
    
    ssp<-substr(odds,1,nchar(odds)-11)
    Tabela.final<-cbind(ssp,Proportional.Moran)
    write.table(Tabela.final,paste("I_Moran_Entre_Quadrantes",m*100,".txt",sep=""),row.names=FALSE,sep="\t")
  }
  
  Moran.Aleat("C:\\Users\\Andre\\Google Drive\\Mestrado\\Capitulo 1(Transferabilidade)\\Ocorrencias_20\\Aleatorio\\Replicas\\Split Ocorrencias\\Rep10","C:\\Users\\Andre\\Google Drive\\Mestrado\\Capitulo 1(Transferabilidade)\\Ocorrencias_20\\Aleatorio\\Replicas\\Split Ocorrencias\\Rep10",0.3)
  