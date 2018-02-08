# ENM_TheMetaLand

Welcome! This is the R Script from TheMetaLand Lab to create ENMs  
Overall, there is a main script (ENM_TheMetaLand) and a group of auxiliary functions (Scripts_for_ENM_TheMetaLand)  
Please follow the "installation" instructions

## Installation  
1.Copy the folder "Scripts_for_ENM_TheMetaLand" to "C:" in your computer  
2.DONE!

## What I can do with ENM_The_MetaLand?  
There are a couple of pre and post-processing available in the function, here is a list of what is currently available:  
**1.** PCA on environmental variables  
**2.** Project to other time/spatial locations (PCA included!)  
**3.** Control Presence/Pseudo-absence Ratio  
**4.** Different pseudo-absence allocation methods  
**5.** Different data-partition methods for model evaluation (random or geographically structured)  
**6.** Nine different algorithms  
**7.** Create presence-absence maps  
**8.** Incorporate effects from spatial distribution (M-SDM)  
**9.** Create Ensemble from the different algorithms  


## How to run ENM_TheMetaLand?  
The function has several input arguments, specify all of them as your desires.  

## Input Parameters:  
* **Dir**: Folder with predictors (4 file formats are supported: ASC/BILL/TIFF/TXT)  
* **Sp:** Name of the column with information about species names  
* **x:** Name of the column with information about longitude  
* **y:** Name of the column with information about latitude  
* **NMin:** Minimum number of unique occurrences (species with less than this number will be excluded)  
* **PCA:** Do you wish to perform a PC on your predictors?(S/N) !Predictors will automatically be used for the modelling process!  
* **Proj:** Project the model onto another region or time period? (S/N)  
* **PabR:** Presence-Absence Ratio  
* **PabM:** Pseudo-absence Selection Method  
  + **rnd:** Random  
  + **const:** Constrained by a Bioclim Model  
  + **zoo:** Constrained by ZooRegions(Requires a ASC/TIFF file with information about the Regions, -Recommended Holt et al 2013-)  
* **Part:** Data partition method  
  + **boot:** Random bootstrap partition (e.g. 70 training - 30% test)  
  + **cross:** Random partition in k-fold  
  + **band:** Geographic partition structured as bands (latitudinal or longitudinal) (See Velazco et al 2017)  
  + **check:** Geographic partition structured as a checkerboard (See Velazco et al 2017)  
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
  + **ANN:** Artificial Neural Network  
* **Thr:** Available threshold for model evaluation and creation of presence-absence maps  
  + **LPT:** The highest threshold at which there is no omission  
  + **MAX:** Threshold at which the sum of the sensitivity and specificity is highest  
* **MSDM:** Include Spatial Restrictions  
  + **N:** Do not include  
  + **LatLong:** Create two layers (Latitude and Longitude of each cell) [added as a predictor]  
  + **Min:** Create a layer with information of the distance from each cell to the closest occurrence [added as a predictor]  
  + **Cum:** Create a layer with information of the summed distance from each cell to ALL occurrences [added as a predictor]  
  + **Kern:** Create a layer with a Gaussian-Kernel on the occurrence data [added as a predictor]  
  + **Land:** Spatial restriction based on adequability patches [NOT added as a predictor] (See Mendes et al 2017)  
* **ENS:** Ensemble of the different algorithms  
  + **N:** No Ensemble  
  + **Mean:** Simple average of the different models  
  + **Sup:** Average of the models with TSS value above the average of TSS values for all models  
  + **PCA:** Performs a PCA and returns the first axis  
  + **PCA_Sup:** Performs a PCA only with models with TSS value above the average of TSS values for all models  
  + **PCA_Thr:** Performs a PCA, but cells with suitability values under the threshold are set to 0  
      
## Where do I input my Occurrence Data?  
Not everything is initially input at the function!  
There will be some specific parameters in which the user will need to "awnser" in the command during the execution  
Occurrence data will be selected at a given point during the process, in which the user will be asked to select the occurrence file [TXT FORMAT!!!]  

## What else I will need to awnser?  
Some specific questions will be asked during the process (e.g. Select the percentage of data for training(0-1) will be asked if you chose the "boot" partition)  
Pay attention and make your choices accordingly. There are specific valid awnsers, the procedure will NOT move ahead if the awnser deviates from what is ecpected.  

## Where are my results?  
One level above the folder of your environmental variables will be created a Result folder, inside you will find a folder for each algorithm and a folder for Ensemble (ENS) and projections (FUT), if you chose to do so.  
There are also four (4) txt files:  
 **Info_Modelling:** Information of the modelling parameters  
 **Occ_filter:** Filtered occurrence with data used in the modelling (without the excluded species)  
 **Validation_Partial:** Information of the evaluation of partial models (e.g. while projecting the model onto the 30% left for test)  
 **Validation_Complete:** Information of the evaluation of models created with full data (Presence-absence maps are created from the Threshold of complete models)

## Last but not least  
There are no defaults! We believe every ENM experiment should be carefully planned and that every decision matters! There are some recommended parameters and literature on which those were based, but they were not included as a default option for the modelling routine.  
The function is still ongoing! So expect frequent updates and always check if you have the stable version!
Everytime a new version is updated it will accompanied by a log with the version number and what are the main changes/what was included

## And finally...  
Enjoy!Test the function and give us feedback! We want to make the world of ENM easier and more accessible to everyone!
