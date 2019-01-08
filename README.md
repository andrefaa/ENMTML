   # ENM_TheMetaLand

Welcome! This is the R Script from TheMetaLand Lab to create ENMs  
Overall, there is a main script (ENM_TheMetaLand) and a group of auxiliary functions  
Please follow the "installation" instructions

## Installation
**1.** install.packages("devtools")  
**2.** library(devtools)  
**3.** install_github("andrefaa/ENM_TheMetaLand")  
**4.** library(ENMTheMetaLand)  

**FOR NEW USERS**  
**5.** install_github('GRaF', 'goldingn')  
**6.** library(GRaF)

## Run
ENMs_TheMetaLand(Dir="",Sp="",x="",y="",NMin=,PCA="",Proj="",SetEval="",SpeciesM="",PabR=,PabM="",
                  Part="",SavePart="N",SaveFinal="Y",Alg=c(),Thr="",MSDM="",ENS=c())

**See possible input options below**

## What I can do with ENM_The_MetaLand?  
There are a couple of pre and post-processing available in the function, here is a list of what is currently available:  
**1.** PCA on environmental variables  
**2.** Project to other time/spatial locations (PCA included!)  
**3.** Automatically restrict the accessible area (M) before model fitting   
**4.** Specify an user defined dataset for evaluation    
**5.** Control Presence/Pseudo-absence Ratio    
**6.** Different pseudo-absence allocation methods    
**7.** Different data-partition methods for model evaluation (random or geographically structured)  
**8.** Nine different algorithms    
**9.** Create presence-absence maps (from dismo Thresholds)   
**10.** Incorporate spatial restrictions (M-SDM)     
**11.** Create Ensemble from the different algorithms  


## How to run ENM_TheMetaLand?  
The function has several input arguments, specify all of them as your desires.  

## Input Parameters:  
* **Dir**: Folder with predictors (4 file formats are supported: ASC/BILL/TIFF/TXT)  
* **Sp:** Name of the column with information about species names  
* **x:** Name of the column with information about longitude  
* **y:** Name of the column with information about latitude  
* **NMin:** Minimum number of unique occurrences (species with less than this number will be excluded)  
* **PCA:** Do you wish to perform a PC on your predictors?(Y/N) !Predictors will automatically be used for the modelling process!  
* **Proj:** Project the model onto another region or time period? (Y/N)  
* **SetEval:** Use an pre-determined set of occurrences for validation? (Y/N)
* **SpeciesM:** Restrict the acessible area M? (Species-specific) (Y/N)
* **PabR:** Presence-Absence Ratio  
* **PabM:** Pseudo-absence Selection Method  
  + **rnd:** Random  
  + **const:** Constrained by a Bioclim Model  
* **Part:** Data partition method  
  + **boot:** Random bootstrap partition (e.g. 70 training - 30% test)  
  + **cross:** Random partition in k-fold  
  + **band:** Geographic partition structured as bands (latitudinal or longitudinal)  
  + **check:** Geographic partition structured as a checkerboard
* **SavePart:** Save .tif files of the partitions? (Y/N)
* **SaveFinal:** Save .tif files of the final model (fitted with all the data)[Default="Y"]? (Y/N)
* **Alg:** List of available algorithms  
  + **BIO:** Bioclim  
  + **MXS:** Maxent Simple[only linear and quadratic features] (MaxNet)  
  + **MXD:** Maxent Default[all features] (MaxNet)  
  + **SVM:** Support Vector Machine  
  + **GLM:** Generalized Linear Model  
  + **GAM:** Generalizes Additive Model  
  + **RDF:** Random Forest  
  + **MLK:** Maximum Likelihood  
  + **GAU:** Gaussian   
* **Thr:** Threshold used for presence-absence maps (from package dismo)  
  + **no_omission:** The highest threshold at which there is no omission  
  + **spec_sens:** Threshold at which the sum of the sensitivity and specificity is highest
  + **kappa:** the threshold at which kappa is highest ("max kappa")
  + **prevalence:** modeled prevalence is closest to observed prevalence
  + **equal_sens_spec:** equal sensitivity and specificity
  + **sensitivty:** fixed (specified) sensitivity 
  + **Any number between 0-1**
* **MSDM:** Include Spatial Restrictions  
  + **N:** Do not include  
  + **XY:** Create two layers (Latitude and Longitude of each cell) [added as a predictor]  
  + **MIN:** Create a layer with information of the distance from each cell to the closest occurrence [added as a predictor]  
  + **CML:** Create a layer with information of the summed distance from each cell to ALL occurrences [added as a predictor]  
  + **KER:** Create a layer with a Gaussian-Kernel on the occurrence data [added as a predictor]  
  + **POST:** Posterior M-SDM Methods [NOT added as a predictor]  
    - **OBR:** Occurrence based restriciton, uses the distance between points to exclude far suitable patches (Mendes et al, in prep)
    - **LR:** Lower Quantile, select the nearest 25% patches (Mendes et al, in prep)  
    - **PRES:** Select only the patches with confirmed occurrence data (Mendes et al, in prep)  
    - **MCP:** Excludes suitable cells outside the Minimum Convex Polygon of the occurrence data (Kremen et al, 2008)  
    - **MCP-B:** Creates a Buffer around the MCP (distance defined by user; Kremen et al, 2008)  
* **ENS:** Ensemble of the different algorithms  
  + **N:** No Ensemble  
  + **Mean:** Simple average of the different models  
  + **Sup:** Average of the models with TSS value above the average of TSS values for all models  
  + **PCA:** Performs a PCA and returns the first axis  
  + **PCA_Sup:** Performs a PCA only with models with TSS value above the average of TSS values for all models  
  + **PCA_Thr:** Performs a PCA, but cells with suitability values under the threshold are set to 0  
      
## Where do I inform my Occurrence Data?  
Not everything will be input as arguments at the beggining!  
There will be some specific parameters which the user will need to "awnser" in the command during the execution  
Occurrence data will be selected at a given point during the process, the user will be asked to select the occurrence file [TXT FORMAT!!!]  

## What else do I need to awnser?  
Some specific questions will be asked during the process (e.g. Select the percentage of data for training(0-1) will be asked if you chose the "boot" partition)  
Pay attention and make your choices accordingly. There are specific valid awnsers, the procedure will NOT move ahead if the answer deviates from what is ecpected.  

## Where are my results?  
One level above the folder of your environmental variables will be created a **Result** folder, inside you will find a folder for each algorithm and a folder for Ensemble (ENS) and projections (FUT), if you chose to do so.  
There are also four (4) txt files:     
 **N_Unique_OCC:** Number of unique occurrences by species     
 **Info_Modelling:** Information of the modelling parameters       
 **Occ_filter:** Filtered occurrence with data used in the modelling (without the excluded species)        
 **Validation_Partition:** Information of the evaluation of partial models (e.g. while projecting the model onto the 30% left for test)       
 **Thresholds_Complete:** Information about the thresholds used to create the presence-absence maps (Presence-absence maps are created from the Threshold of complete models)    

## Last but not least  
There are no defaults! We believe every ENM experiment should be carefully planned and that every decision matters! There are some recommended parameters and literature on which those were based, but they were not included as a default option for the modelling routine.  
The function is still ongoing! So expect frequent updates and always check if you have the stable version!
Everytime a new version is updated it will accompanied by a log with the version number and what are the main changes/what was included

## And finally...  
Enjoy!Test the function and give us feedback! We want to make the world of ENM easier and more accessible to everyone!
