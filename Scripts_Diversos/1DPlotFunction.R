#Open R and change directory to the folder with data using the file, Change dir... option

setwd ("C://Users//Andre//Desktop//Andre//Mestrado//Teste pseudo-transferabilidade Mhispida//ASC & Suitability Plot//M1")

## This function requires RSAGA package. 
## Parameter passed....
##      ASCIIFile - Output file from modeling algorithms
##      Pmin - Minimum value of precipitation for the species
##      Pmax - Maximum value of precipitation for the species
##      Tmin - Minimum value of temperature for the species
##      Tmax - Maximum value of temperature for the species
##      PfileName - Plot file name for precipitation 
##      TfileName - Plot file name for temperature 

OneDPlot <- function(ASCIIFile, Pmin, Pmax, Tmin, Tmax, PfileName, TfileName)
{
   library(RSAGA)

   #Read asc file in as a matrix, excluding the header
   dat<-read.ascii.grid("Mussismilia_hispida", return.header = FALSE)
   #check the dimensions to make sure they are correct 15001 2001(varies with the map dimension)
   print(dim(dat))
   
   #Remove NA cells from a matrix(remove entire row or column if there is a NA)
   #library(functional)
   #datnasal<-dat[apply(dat, 2, Compose(is.finite, all)),]
   #print(dim(datnatemp))
   
#Average the suitability values in each column for temp and for each row for precipitation(salinity)
##Maxent multiply by 100 to scale between 0 - 100 in plot
#OMDK and OMOM   DO NOT MULTIPLY here
##DKDK multiply by 10

   AveT<-colMeans(dat)*100
   AveSal<-rowMeans(dat)*100


#Import the ASCS for salinity and temperature for plot comparison(ADICIONADO POR MIM, SUBSTITUTO AOS VETORES)

  sal<-read.ascii.grid("salinity", return.header = FALSE)
  temp<-read.ascii.grid("sstmax", return.header = FALSE)

  print(dim(sal))
  print(dim(temp))

#create a vector from -1000 to 1000 for the range of values for temp and 15000 to 0 for precip

   #temp<-seq(from = 0, to = 89, by = 1)
   #Sal<-seq(from = 50, to = 0, by = -1)

# Now create the plots and save them in file name supplied. 

   slimit<-c("35.31","35.76")
   windows()
   plot(sal,dat, main="Precipitation vs Suitability", xlab="Precipitation (in mm)", ylab="Average Suitability",pch=19,col = "blue",abline(v=slimit,col = "red"),xlim=c(35,37),ylim=c(0,1))
   savePlot(filename = "Mhispida_1_prec_suitability", type = c("tiff"), device = dev.cur(), restoreConsole = TRUE)

   tlimit<-c("25.92","29.35") 
   windows()
   plot(temp,dat, main="Temperature vs Suitability", xlab="Temperature", ylab="Average Suitability",pch=19,col = "blue",abline(v=tlimit,col = "red"), xlim=c(25,32), ylim=c(0,1))
   savePlot(filename = "Mhispida_1_temp_suitability", type = c("tiff"), device = dev.cur(), restoreConsole = TRUE)

}

