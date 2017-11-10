Spatial.Correlation <- function(dir_input,dir_output,m){
  #Parametros:
  
    #dir_input: input directory with xls files
    #dir_output: output directory
  
    #Occurrence data must have columns with the following names: Taxa,Long,Lat
    #Occurrence data must be a txt file
  
  library(SDMTools)
  library(raster)
  library(xlsx)
  
  setwd(dir_input)
  diretorio<-getwd()
  
  #Read the environmental asc file from where env_data will be extracted
    print("Select ASC file to be used to extrac environmental data from:")
    asc.base<-read.asc("D:\\Users\\Andre\\Google Drive\\Mestrado\\Capitulo 1(Transferabilidade)\\Ocorrencias_20\\Input\\calibrate\\PCA\\pca1.asc")
    asc.base<-raster(asc.base)
  
  #Create and divide xls files in odds and evens
    completo<-list.files(pattern=".xls")
    
  Proportional.Moran<-NULL
  
for (a in 1:length(completo)){
  
#Import occurrence data
print(substr(completo[a],1,nchar(completo[a])-4))
  
  #Full dataset
    occs<-read.xlsx(completo[a],sheetIndex=1,header=F)
    occs.xy<-occs[,c(3,2)]
    occs.env<-extract(asc.base,occs.xy)
    occs<-cbind(occs,occs.env)

  
  #Calculate Z values for each observation
    occs.z<-occs$occs.env-mean(occs$occs.env)
    occs<-cbind(occs,occs.z)

  #Remove Errors
  occs<-occs[occs[,4]!=0,]

  #Calculate distance to nearest point in another quadrant
  near.dist<-NULL
  
    for (x in 1:nrow(occs)){
      for (z in 1:max(occs[,4])){
      dist.temp<-NULL
      quad.alvo<-NULL
      dist.quad<-NULL
        for (y in 1:nrow(occs)){
          dist.temp[y]<-sqrt((occs[x,2]-occs[y,2])^2+(occs[x,3]-occs[y,3])^2)*((occs[y,4]==(occs[x,4]-z))||(occs[y,4]==(occs[x,4]+z)))
          quad.alvo[y]<-occs[y,4]
        }
        if ((any(dist.temp)!=0)){
          dist.quad<-cbind(dist.temp,quad.alvo)
          dist.quad<-as.matrix(dist.quad)
          dist.quad<-dist.quad[dist.quad[,1]!=0,]
          dist.quad<-matrix(dist.quad,ncol=2)
          near.dist<-rbind(near.dist,dist.quad[dist.quad[,1]==min(dist.quad[,1])])
          break
        }
      }
    }
  occs<-cbind(occs,near.dist)
  names(occs)<-c("Species","Lat","Long","Subset","Occ.Env","Occ.Z","Dist","Quad.Alvo")

  Binary<-rep(0,nrow(occs))
  occs<-cbind(occs,Binary)
  
  #Create a Sub-matrix with only x% of the nearest points
  occ.quad<-tabulate(occs[,4])
  subset.matrix<-NULL
  for (z in 1:max(occs[,4])){
    for (k in 1:max(occs[,4])){
    occs.quad<-occs[occs[,4]==z,]
    occs.adj1<-occs.quad[occs.quad[,8]==occs.quad[,4]+k,]
    occs.adj2<-occs.quad[occs.quad[,8]==occs.quad[,4]-k,]
    
    if (nrow(occs.adj1)!=0 || nrow(occs.adj2)!=0){
      occs.adj1<-occs.adj1[order(occs.adj1[,7]),]
      occs.adj2<-occs.adj2[order(occs.adj2[,7]),]
      occs.adj1[1:ceiling(m*occ.quad[z]),9]<-1
      occs.adj2[1:ceiling(m*occ.quad[z]),9]<-1
      subsets<-rbind(occs.adj1,occs.adj2)
      subsets<-na.omit(subsets)
      subset.matrix<-rbind(subset.matrix,subsets)
      break
    }
    }
  }
  occs<-subset.matrix
  
  #Calculate Moran's I among subsets
  
  W<-0
  Moran.I<-0
  
  for (x in 1:nrow(occs)){
    for(y in 1:nrow(occs)){
        Moran.I<-Moran.I+((occs[x,6]*occs[y,6])*((occs[y,4]==occs[x,8])&&(occs[y,8]==occs[x,4])&&(occs[x,9]==1)))
        W<-W+((occs[y,4]==occs[x,8])&&(occs[y,8]==occs[x,4]))
    }
  }
  W<-2*W
  Moran.I<-(nrow(occs)*Moran.I)/(W*sum((occs$Occ.Z)^2))
  
  #Calculate I(max)
  
  Imax<-0
  Div<-0
   
  for (x in 1:nrow(occs)){
    Imax.temp<-0
    for (y in 1:nrow(occs)){
        Imax.temp<-Imax.temp+((occs[y,5]-mean(occs[,5]))*((occs[y,4]==occs[x,8])&&(occs[y,8]==occs[x,4])&&(occs[x,9]==1)))
    }
    Imax<-Imax+(Imax.temp^2)
    Div<-Div+(occs[x,5]-mean(occs[,5]))^2
  }
  Imax<-(nrow(occs)/W)*sqrt(Imax/Div)
  Imax<-abs(Imax)
  
  #Calculate Moran's I proportional to Imax
  Proportional.Moran[a]<-Moran.I/Imax
  print(Proportional.Moran[a])
  Sys.sleep(1)
}

setwd(dir_output)

ssp<-substr(completo,1,nchar(completo)-4)
Tabela.final<-cbind(ssp,Proportional.Moran)
write.table(Tabela.final,paste("I_Moran_Entre_Quadrantes_",substr(diretorio, nchar(diretorio)-3+1, nchar(diretorio)),"_1classe",m*100,".txt",sep=""),row.names=FALSE,sep="\t")
}

Spatial.Correlation("D:\\Users\\Andre\\Google Drive\\Mestrado\\Capitulo 1(Transferabilidade)\\Ocorrencias_20\\I Moran x D Schoener\\Input2Q","D:\\Users\\Andre\\Google Drive\\Mestrado\\Capitulo 1(Transferabilidade)\\Ocorrencias_20\\I Moran x D Schoener",0.2)
