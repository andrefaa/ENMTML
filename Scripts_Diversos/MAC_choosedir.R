choose.dir <- function() {
  system("osascript -e 'tell app \"R\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
         intern = FALSE, ignore.stderr = TRUE)
  p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
  return(ifelse(length(p), p, NA))
}