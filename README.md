   # ENMTML

Welcome! This is the R Script from TheMetaLand Lab to create ENMs  
Overall, there is a main script (ENMTML) and a group of auxiliary functions  
Please follow the "installation" instructions

## Installation
```ruby
install.packages("devtools")  
library(devtools)  
install_github("andrefaa/ENMTML")  
library(ENMTML)  
```

**FOR NEW USERS**  
```ruby
install_github('goldingn/GRaF')  
library(GRaF)
```

## Run
ENMTML(pred_dir, proj_dir = NULL, occ_file, sp, x, y, min_occ = 10,
  thin_occ = NULL, eval_occ = NULL, colin_var = NULL,
  imp_var = NULL, sp_accessible_area = NULL, pseudoabs_method,
  pres_abs_ratio = 1, part, save_part = FALSE, save_final = TRUE,
  algorithm, thr, msdm = NULL, ensemble = NULL,
  extrapolation = FALSE, cores = 1)

**See possible input options below**

## What I can do with ENM_The_MetaLand?  
There are a couple of pre and post-processing available in the function, here is a list of what is currently available:  
**1.** Control colinearity on environmental variables  
**2.** Project to other time/spatial locations (PCA included)  
**3.** Automatically restrict the accessible area (M) before model fitting   
**4.** Specify an user defined dataset for evaluation    
**5.** Control Presence/Pseudo-absence Ratio    
**6.** Different pseudo-absence allocation methods    
**7.** Different data-partition methods for model evaluation (random or geographically structured)  
**8.** Twelve different algorithms    
**9.** Create presence-absence maps   
**10.** Incorporate spatial restrictions (M-SDM)     
**11.** Create Ensemble from the different algorithms  
**12.** Create Stacked Species Distribution Modelling (S-SDM)


## How to run ENM_TheMetaLand?  
The function has several input arguments, specify all of them according to your desires.  

## Input Parameters:  
* **pred_dir**: Folder with predictors (4 file formats are supported: ASC/BILL/TIFF/TXT)  
* **occ_dir**: Path to occurrence file (TXT only!)  
* **sp:** Name of the column with information about species names  
* **x:** Name of the column with information about longitude  
* **y:** Name of the column with information about latitude  
* **min_occ:** Minimum number of unique occurrences (species with less than this number will be excluded)  
* **thin_occ:** Perform a spatial filtering (Thinning) on the presences? (Y/N)
* **colin_var:** Processes to reduce variable collinearity?  
  + **N**: Use original variables  
  + **Pearson:** Select variables by Pearson Correlation (Threshold specified by user)  
  + **VIF:** Variance Inflation Factor (Chatterjee and Hadi 2006)  
  + **PCA:** Perform a PCA on predictors and use PCs as environmental variables  
* **imp_var:** Calculate vairable importance and model response curve (Y/N)  
* **transfer:** Project the model onto another region or time period? (Y/N)  
* **eval_occ:** Specify a set of occurrences for validation? (Y/N)  
* **sp_accessible_area:** Create species-specific accessible areas (Y/N)  
* **pres_abs_ratio:** Presence-Absence Ratio  
* **pseudoabs_method:** Pseudo-absence Selection Method  
  + **RND:** Random allocation throughout area used to fit models  
  + **ENV_CONST:** Pseudo-absences are environmentally constrained to region with lower suitability values predicted by a Bioclim model  
  + **GEO_CONST:** Pseudo-absences are allocated far from occurrences, constrained by a geographical buffer  
  + **GEO_ENV_CONST:** Pseudo-absences are cosntrained both environmentally (Bioclim Model) and geographically (buffer)  
  + **GEO_ENV_KM_CONST:** Pseudo-absences constrained on a three-level proccedure  
* **part:** Data partition method for model evaluation  
  + **BOOT:** Random bootstrap partition (e.g. 70 training - 30% test)  
  + **KFOLD:** Random partition in k-fold cross validation  
  + **BAND:** Geographic partition structured as bands (latitudinal or longitudinal)  
  + **BLOCK:** Geographic partition structured as a checkerboard  
* **save_part:** Save .tif files of the partitions?[Default="N"] (Y/N)
* **save_final:** Save .tif files of the final model (fitted with all the data)[Default="Y"]? (Y/N)
* **algorithm:** List of available algorithms  
  + **BIO:** Bioclim  
  + **MAH:** Mahalanobis  
  + **DOM:** Domain  
  + **ENF:** ENFA  
  + **MXS:** Maxent Simple[only linear and quadratic features] (MaxNet)  
  + **MXD:** Maxent Default[all features] (MaxNet)  
  + **SVM:** Support Vector Machine  
  + **GLM:** Generalized Linear Model  
  + **GAM:** Generalizes Additive Model 
  + **BRT:** Boosted Regression Tree
  + **RDF:** Random Forest  
  + **MLK:** Maximum Likelihood  
  + **GAU:** Gaussian   
* **thr:** Threshold used for presence-absence maps (from package dismo)  
  + **LPT:** The highest threshold at which there is no omission  
  + **MAX_TSS:** Threshold at which the sum of the sensitivity and specificity is highest
  + **MAX_KAPPA:** the threshold at which kappa is highest ("max kappa")
  + **JACCARD:** the threshold at which Jaccard is highest  
  + **SORENSEN:** the threshold at which Sorensen is highest  
  + **SENSITIVITY:** fixed (specified) sensitivity 
* **msdm:** Include Spatial Restrictions  
  + **N:** Do not include  
  + **XY:** Create two layers (Latitude and Longitude of each cell) [added as a predictor]  
  + **MIN:** Create a layer with information of the distance from each cell to the closest occurrence [added as a predictor]  
  + **CML:** Create a layer with information of the summed distance from each cell to ALL occurrences [added as a predictor]  
  + **KER:** Create a layer with a Gaussian-Kernel on the occurrence data [added as a predictor]  
  + **POST:** Posterior M-SDM Methods (If chosen, prefered method will be asked later) [NOT added as a predictor]  
    + **OBR:** Occurrence based restriciton, uses the distance between points to exclude far suitable patches (Mendes et al, in prep)  
    + **LR:** Lower Quantile, select the nearest 25% patches (Mendes et al, in prep)  
    + **PRES:** Select only the patches with confirmed occurrence data (Mendes et al, in prep)  
    + **MCP:** Excludes suitable cells outside the Minimum Convex Polygon of the occurrence data (Kremen et al, 2008)  
    + **MCP-B:** Creates a Buffer around the MCP (distance defined by user; Kremen et al, 2008)  
* **ensemble:** Ensemble of the different algorithms  
  + **N:** No Ensemble  
  + **MEAN:** Simple average of the different models  
  + **W_MEAN:** Weighted Average  
  + **SUP:** Average of the best models (TSS over the average)  
  + **PCA:** Performs a PCA and returns the first axis  
  + **PCA_SUP:** PCA of the best models (TSS over the average)  
  + **PCA_THR:** PCA only with cells above the threshold  
* **s_sdm:** Stacked Species Distribution Model (Y/N)
      
## Is that all?  
Not everything will be input as arguments at the beggining!  

Some specific questions will be asked during the process (e.g. Select the percentage of data for training(0-1) will be asked if you chose the "boot" partition)  

## Where are my results?  
One level above the folder of your environmental variables will be created a **Result** folder, inside you will find a folder for algorithm results, ensemble results (if chosen), projections and maps of areas of extrapolation.  
There are also four (4) txt files:     
 **N_Unique_OCC:** Number of unique occurrences of each species    
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
