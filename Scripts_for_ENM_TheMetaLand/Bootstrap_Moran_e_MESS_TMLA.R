### This is the R code to conduct ENMs evaluation by the Odds-and-Evens Framework
# In the moment, the code is configured to create the splits only by its latitudinal axis
# Read carefully the Initialization and function parameters and the messages that appear on the console
# Developers: Paulo De Marco Junior & Andre Andrade

Bootstrap_Moran_e_MESS_TMLA <- function(Env,
                                        RecordsData,
                                        DirO,
                                        repl){
  
  #Function Parameters:
  #env:environmental variables stack
  #occ: occurrence data

  #Important Observations:
  #Occurence File MUST have the columns named Species (containing species name),Long (containing longitude),Lat (containing latitude)
  
  #Development
  occTR <- RecordsData[RecordsData$PresAbse==1,]
  
  # Lists for validation-----
  assign(paste("Bootstrap_Moran_MESS", sep=""),list())
  VALNAME <- paste0("Bootstrap_Moran_MESS",repl,".txt")
  
  #Calculate Moran's I
  Moran<-Moran_for_Bootstrap_TMLA(occTR=occTR,pc1=Env[[1]])

  #MESS----
  Occ_e <- cbind(occTR,extract(Env, occTR[,2:3]))
  Occ_e <- split(Occ_e,f=Occ_e$sp)
  Occ_e <- lapply(Occ_e, function(x) split(x,f=x$Partition))
    
  mess <- lapply(Occ_e,function(x) MESS(x[[1]][,-c(1:5)], x[[2]][,-c(1:5)]))
  mess <- ldply(lapply(mess, function(x) mean(x$TOTAL, na.rm = TRUE)),.id=NULL)
  Bootstrap_Moran_MESS[[1]]<-data.frame(Moran,MESS=mess)
  
  # Save .txt with models Moran&MESS---- 
  Obj <- ls(pattern = 'Bootstrap_Moran_MESS')
  res <- list()
  for(i in 1:length(Obj)){
    res[[i]] <- ldply(get(Obj[i]))}
  res <- ldply(res)
  colnames(res)[3] <- "MESS"
  if(is.null(repl)==F){
    res <- data.frame(res,Replicate=repl)
  }
  write.table(res,file.path(DirO,VALNAME),sep="\t",row.names=F,col.names=T)
}