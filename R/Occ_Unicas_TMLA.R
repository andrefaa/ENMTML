Occ_Unicas_TMLA<-function(env,
                          occ.xy,
                          DirO){
  spN <- names(occ.xy)
  occ.v <- lapply(occ.xy,function(x) extract(env,x,cellnumber=T))
  occ.v <- mapply(cbind, occ.xy, occ.v, SIMPLIFY=F)
  occ.v <- lapply(occ.v, function(x) na.omit(x))
  occ.v <- lapply(occ.v, function(x) x[!duplicated(x$cells),])
  uni <- data.frame(Species=spN,UniqueOcc=sapply(occ.v,function(x) nrow(x)))
  write.table(uni,file.path(DirO,"N_Occ_Unicas.txt"),sep="\t",row.names=F)
  return(occ.v)
}