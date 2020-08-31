Occ_Unicas_TMLA <- function(env,
                            occ.xy,
                            DirO) {
  spN <- names(occ.xy)
  occ.v <-
    lapply(occ.xy, function(x)
      data.frame(raster::extract(env, x, cellnumber = T))[1])
  occ.v <- mapply(cbind, occ.xy, occ.v, SIMPLIFY = F)
  occ.v <- lapply(occ.v, function(x)
    stats::na.omit(x))
  occ.v <- lapply(occ.v, function(x)
     x[!duplicated(x$cells), ])
  uni <-
    data.frame(Species = spN,
               UniqueOcc = sapply(occ.v, function(x)
                 nrow(x)))
  utils::write.table(
    uni,
    file.path(DirO, "Number_Unique_Occurrences.txt"),
    sep = "\t",
    row.names = F
  )
  return(occ.v)
}
