Occ_Unicas_TMLA<-function(env,occ.xy){
  occ.v <- lapply(occ.xy,function(x) extract(env,x,cellnumber=T))
  occ.v <- mapply(cbind, occ.xy, occ.v, SIMPLIFY=F)
  occ.v <- lapply(occ.v, function(x) na.omit(x))
  occ.v <- lapply(occ.v, function(x) x[!duplicated(x$cells),])
  return(occ.v)
}