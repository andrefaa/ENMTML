# __ENMTML__ <img src="https://github.com/andrefaa/ENMTML/blob/master/ENMTML1.png?raw=true" align="right" width="170" />
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)


## An R package for an integrated construction of Ecological Niche Models



Please follow the installation instructions


### Installation
```ruby
# install.packages("devtools")  
library(devtools)  
install_github("andrefaa/ENMTML")  
library(ENMTML)  
```

**FOR NEW USERS**  
```ruby
if(!require(devtools)){
    install.packages("devtools")
}
if(!require(GRaF)){
    devtools::install_github("goldingn/GRaF")
}
if(!require(ellipsenm)){
    devtools::install_github("marlonecobos/ellipsenm")
}
devtools::install_github("andrefaa/ENMTML")
library(ENMTML)
```


### What I can do with __ENMTML__ package?  
There are a couple of pre and post-processing available in the package, here is a list of what is currently available:  
**1.** Control colinearity on environmental variables  
**2.** Project to other time/spatial locations (PCA included)  
**3.** Automatically restrict the accessible area (M) before model fitting   
**4.** Specify an user defined dataset for evaluation    
**5.** Control Presence/Pseudo-absence Ratio    
**6.** Different pseudo-absence allocation methods    
**7.** Different data-partition methods for model evaluation (random or geographically structured)  
**8.** Twelve different algorithms    
**9.** Create presence-absence maps   
**10.** Incorporate spatial restrictions (MSDM)     
**11.** Create Ensemble from the different algorithms  

For more information see 
https://www.sciencedirect.com/science/article/pii/S1364815219310424?via%3Dihub
### Please cite this pacakge as:
**Andrade AFA, Velazco SJE, De Marco Jr P. 2020. ENMTML: An R package for a straightforward construction of complex ecological niche models. Environmental Modelling & Software:104615.(https://doi.org/10.1016/j.envsoft.2019.104615)**


### How to run __ENMTML__?  
The package have in main function *ENMTML* with several arguments, specify all of them according to your modeling needs.  

```ruby
ENMTML(pred_dir, proj_dir = NULL, occ_file, sp, x, y, min_occ = 10,
  thin_occ = NULL, eval_occ = NULL, colin_var = NULL,
  imp_var = FALSE, sp_accessible_area = NULL, pseudoabs_method,
  pres_abs_ratio = 1, part, save_part = FALSE, save_final = TRUE,
  algorithm, thr, msdm = NULL, ensemble = NULL,
  extrapolation = FALSE, cores = 1)
```

**See possible input options below**

### ENMTML function arguments:  
* **pred_dir**: character. Directory path with predictors (file formats supported are ASC, BILL, TIFF or TXT)
* **proj_dir**: character. Directory path containing folders with predictors for different regions or time periods used to project models (file formats supported are ASC, BILL, TIFF, or TXT).
* **occ_file**: character. Directory path with the tab-delimited TXT file, which will contain at least three columns with information about species names, and the latitude and longitude of species occurrences.
* **sp**: character. Name of the column with information about species names.
* **x**: character. Name of the column with information about longitude.
* **y**: character. Name of the column with information about latitude.
* **min_occ**: integer. Minimum number of unique occurrences (species with less than this number will be excluded).
* **thin_occ**: character. Perform spatial filtering (Thinning, based on spThin package) on the presences. For this augment it is necessary provide a vector in which its elements need to have the names 'method' or 'method' and 'distance' (more information below). Three thinning methods are available (default `NULL`):
  + **MORAN** Distance defined by Moran Variogram. Usage </br> <b style='color:red'>`thin_occ=c(method='MORAN')`</b>.
  + **CELLSIZE** Distance defined by 2x cellsize (Haversine Transformation). Usage </br> 
  <b style='color:red'>`thin_occ=c(method='CELLSIZE')`</b>.
  + **USER-DEFINED** User defined distance. For this option it is neede provide a vector with two values. Usage </br>
  <b style='color:red'>`thin_occ=c(method='USER-DEFINED', ditance='300')`</b>. The second numeric value refers to the distance in km that will be used for thinning. So distance=300 means that all records within a radius of 300 km will be deleted.
* **eval_occ:** character. Directory path with tab-delimited TXT file with species names, latitude and longitude, these three columns must have the same columns names than the databased used in the `occ_file` argument. This external occurrence database will be used to external models validation (i.e., it will no be use to model fitting). (default `NULL`).
* **colin_var:** character. Method to reduce variable collinearity:
  + **PCA:** Perform a Principal Component Analysis on predictors and use Principal Components as environmental variables. Usage </br> <b style='color:red'>`colin_var=c(method='PCA')`</b>.
  + **VIF:** Variance Inflation Factor. Usage </br><b style='color:red'>`colin_var=c(method='VIF')`</b>.
  + **PEARSON**: Select variables by Pearson correlation, a threshold of maximum correlation must be specified by user. Usage </br><b style='color:red'>`colin_var=c(method='PEARSON', threshold='0.7')`</b>.
* **imp_var:** logical Perform variable importance and curves response for selected algorithms? (default `FALSE`)
* **sp_accessible_area:** character. Restrict for each species the accessible area, i.e., the area used to model fitting. It is necessary to provide a vector for this argument. Three methods were implemented
  + **BUFFER** area used to model fitting deliminted by a buffer with a width size equal to the maximum distance among pair of occurrences for each species. Usage </br><b style='color:red'>`sp_accessible_area=c(method='BUFFER', type='1')`</b>.
  + **BUFFER** area used to model fitting deliminted by a buffer with a width size difinited by the user in km. Note this width size of buffer will be used for all species. Usage </br><b style='color:red'>`sp_accessible_area=c(method='BUFFER', type='2', width='300')`</b>.
  + **MASK** this method consists in delimit the area used to model fitting based on the polygon where a species occurrences fall. For instance, it is possible delimit the calibration area based on ecoregion shapefile. For this option it is necessary inform the path to the file that will be used as mask. Next file format can be loaded '.bil', '.asc', '.tif', '.shp', and '.txt'. Usage </br><b style='color:red'>`sp_accessible_area=c(method='MASK', filepath='C:/Users/mycomputer/ecoregion/olson.shp')`</b>..
* **pseudoabs_method:** character. Pseudo-absence allocation method. It is necessary to provide a vector for this argument. Only one method can be chosen. The next methods are implemented:
  + **RND:** Random allocation of pseudo-absences throughout the area used for model fitting. Usage </br><b style='color:red'>`pseudoabs_method=c(method='RND')`</b>.
  + **ENV_CONST:** Pseudo-absences are environmentally constrained to a region with lower suitability values predicted by a Bioclim model. Usage pseudoabs_method=c(method='ENV_CONST'). Usage </br><b style='color:red'>`pseudoabs_method=c(method='ENV_CONST')`</b>.
  + **GEO_CONST:** Pseudo-absences are allocated far from occurrences based on a geographical buffer. For this method it is necessary provie a second value wich express the buffer width in km. Usage </br><b style='color:red'>`pseudoabs_method=c(method='GEO_CONST', width='50')`</b>.
  + **GEO_ENV_CONST:** Pseudo-absences are constrained environmentally (based on Bioclim model) but distributed geographically far from occurrences based on a geographical buffer. For this method it is necessary provide a second value which express the buffer width in km. Usage </br><b style='color:red'> `pseudoabs_method=c(method='GEO_ENV_CONST', width='50')`</b>
  + **GEO_ENV_KM_CONST:** Pseudo-absences are constrained on a three-level procedure; it is similar to the GEO_ENV_CONST with an additional step which distributes the pseudo-absences in the environmental space using k-means cluster analysis. For this method it is necessary provide a second value which express the buffer width in km. Usage </br><b style='color:red'>`pseudoabs_method=c(method='GEO_ENV_KM_CONST', width='50')`</b>.
* **pres_abs_ratio:** numeric. Presence-Absence ratio (values between 0 and 1).
* **part:** character. Partition method for model's validation. Only one method can be chosen. It is necessary to provide a vector for this argument. The next methods are implemented:
  + **BOOT:** Random bootstrap partition. Usage </br><b style='color:red'>`part=c(method='BOOT', replicates='2',  proportion='0.7')`</b>. `replicate` refers to the number of replicates, it assumes a value >=1. `proportion` refers to the proportion of occurrences used for model fitting, it assumes a value >0 and <=1. In this example `proportion='0.7'` mean that 70% of data will be used for model training, while 30% for model testing.
  + **KFOLD:** Random partition in k-fold cross-validation. Usage </br><b style='color:red'>`part=c(method= 'KFOLD', folds='5')`</b>. `folds` refers to the number of folds for data partitioning, it assumes value >=1.
  + **BANDS:** Geographic partition structured as bands arranged in a latitudinal way (type 1) or longitudinal way (type 2). Usage </br><b style='color:red'>`part=c(method= 'BANDS', type='1')`</b>. `type` refers to the bands disposition.
  + **BLOCK:** Geographic partition structured as a checkerboard (a.k.a. block cross-validation). Usage </br><b style='color:red'>`part=c(method= 'BLOCK')`</b>.
* **save_part:** logical. If `TRUE`, function will save .tif files of partial models, i.e. model created by each occurrence partitions. (default `FALSE`).
* **save_final:** logical. If `TRUE`, function will save .tif files of the final model, i.e. fitted with all occurrences data. (default `TRUE`)
* **algorithm:** character. Algorithm to construct ecological niche models (it is possible to use more than one method):
  + **BIO:** Bioclim  
  + **MAH:** Mahalanobis  
  + **DOM:** Domain  
  + **ENF:** Ecological-Niche Factor Analysis  
  + **MXS:** Maxent simple (only linear and quadratic features, based on MaxNet package)
  + **MXD:** Maxent default features (all features, based on MaxNet package)
  + **SVM:** Support Vector Machine  
  + **GLM:** Generalized Linear Model  
  + **GAM:** Generalized Additive Model 
  + **BRT:** Boosted Regression Tree
  + **RDF:** Random Forest  
  + **MLK:** Maximum Likelihood  
  + **GAU:** Gaussian Process</br> 
  Usage <b style='color:red'> algorithm=c('BIO', 'SVM', 'GLM', 'GAM', 'GAU').</b>
* **thr:** character. Threshold used for presence-absence predictions. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument:
  + **LPT:** The highest threshold at which there is no omission. Usage </br><b style='color:red'>`thr=c(type='LPT')`</b>.
  + **MAX_TSS:** Threshold at which the sum of the sensitivity and specificity is the highest. Usage </br><b style='color:red'> `thr=c(type='MAX_TSS')`</b>.
  + **MAX_KAPPA:** The threshold at which kappa is the highest ("max kappa"). Usage </br><b style='color:red'>`thr=c(type='MAX_KAPPA')`</b>.
  + **SENSITIVITY:** A threshold value specified by user. Usage </br><b style='color:red'>`thr=c(type='SENSITIVITY', sens='0.6')`</b>. 'sens' refers to models will be binarized using this suitability value. Note that this method assumes 'sens' value for all algorithm and species.
  + **JACCARD:** The threshold at which Jaccard is the highest. Usage </br><b style='color:red'>`thr=c(type='JACCARD')`</b>.
  + **SORENSEN:** The threshold at which Sorensen is highest. Usage </br><b style='color:red'>`thr=c(type='SORENSEN')`</b>.</br>
  
  In the case of use more than one threshold type it is necessary concatenate the names of threshold types, e.g., <b style='color:red'>`thr=c(type=c('LPT', 'MAX_TSS', 'JACCARD'))`</b>. When `SENSITIVITY` threshold is used in combination with other it is necessary specify the desired sensitivity value, e.g., <b style='color:red'>`thr=c(type=c('LPT', 'MAX_TSS', 'SENSITIVITY'), sens='0.8')`</b>.

* **msdm:** character. Include spatial restrictions to model projection. These methods restrict ecological niche models in order to have less potential prediction and turn models closer to species distribution models. They are classified in 'a Priori' and 'a Posteriori' methods. The first one encompasses method that include geographical layers as predictor of models' fitting, whereas a Posteriori constrain models based on occurrence and suitability patterns. This argument is filled only with a method, in the case of use MCP-B method msdm is filled in a different way se below (default `NULL`):</b>
  
  _a Priori_ methods (layer created area added as a predictor at moment of model fitting):</b>
  + **XY:** Create two layers latitude and longitude layer. Usage </br><b style='color:red'>`msdm=c(method='XY')`</b>.
  + **MIN:** Create a layer with information of the distance from each cell to the closest occurrence. Usage </br><b style='color:red'>`msdm=c(method='MIN')`</b>.
  + **CML:** Create a layer with information of the summed distance from each cell to all occurrences. Usage</br><b style='color:red'>`msdm=c(method='CML')`</b>.
  + **KER:** Create a layer with a Gaussian-Kernel on the occurrence data. Usage </br><b style='color:red'>`msdm=c(method='KER')`.</b>
  
  _a Posteriori_ methods:</b>
  + **OBR:** Occurrence based restriction, uses the distance between points to exclude far suitable patches (Mendes et al., in prep). Usage </br><b style='color:red'>`msdm=c(method='OBR')`</b>.
  + **LR:** Lower Quantile, select the nearest 25\% patches (Mendes et al., in prep). Usage </br><b style='color:red'>`msdm=c(method='LR')`</b>.
  + **PRES:** Select only the patches with confirmed occurrence data (Mendes et al, in prep). Usage </br><b style='color:red'>`msdm=c(method='PRES')`</b>.
  + **MCP:** Excludes suitable cells outside the Minimum Convex Polygon (MCP) built based on ocurrences data. Usage </br><b style='color:red'>`msdm=c(method='MCP')`</b>.
  + **MCP-B:** Creates a buffer (with a width size defined by user in km) around the MCP. Usage </br><b style='color:red'>`msdm=c(method='MCP-B', width=100)`</b>. In this case `width=100` means that a buffer with 100km of width will be created around the MCP.</b>

* **ensemble:** character. Method used to ensemble different algorithms. It is possible to use more than one method. A vector must be provided for this argument. For SUP, W_MEAN or PCA_SUP method it is necessary provide an evaluation metric to ensemble arguments (i.e., AUC, Kappa, TSS, Jaccard, Sorensen or Fpb) see below. (default NULL):
  + **MEAN:** Simple average of the different models. Usage </br><b style='color:red'>`ensemble=c(method='MEAN')`.</b>
  + **W_MEAN:** Weighted average of models based on their performance. An evaluation metric must be provided. Usage </br><b style='color:red'>`ensemble=c(method='W_MEAN', metric='TSS')`</b>.
  + **SUP:** Average of the best models (e.g., TSS over the average). An evaluation metric must be provided. Usage </br><b style='color:red'>`ensemble=c(method='SUP', metric='TSS')`</b>.
  + **PCA:** Performs a Principal Component Analysis (PCA) and returns the first axis. Usage </br><b style='color:red'>`ensemble=c(method='PCA')`</b>.
  + **PCA_SUP:** PCA of the best models (e.g., TSS over the average). An evaluation metric must be provided. Usage </br><b style='color:red'>`ensemble=c(method='PCA_SUP', metric='Fpb')`</b>.
  + **PCA_THR:** PCA performed only with those cells with suitability values above the selected threshold. Usage </br><b style='color:red'>`ensemble=c(method='PCA_THR')`</b>.</br>
  
  In the case of use more than one ensemble method it is necessary concatenate the names of ensemble methods within the argument, e.g., <b style='color:red'>`ensemble=c(type=c('MEAN', 'PCA'))`</b>, <b style='color:red'>`ensemble=c(method=c('MEAN, 'W_MEAN', 'PCA_SUP'), metric='Fpb')`</b>.

* **extrapolation**	logical. If TRUE the function will calculate extrapolation based on Mobility-Oriented Parity analysis (MOP) for current conditions. If the argument proj_dir is used, the extrapolation layers for other regions or time periods will also be calculated (default `FALSE`).

* **cores** numeric. Define the number of CPU cores to run modeling procedures in parallel (default `1`).

</br>
      

### Where are my results?  
One level above the folder of your environmental variables will be created a **Result** folder, inside you will find a folder for algorithm results, ensemble results (if chosen), projections and maps of areas of extrapolation.  
There are also some .txt files (note that some txt will be created under ceratin modeling settings, e.g., Thresholds_Ensemble.txt):     
 **Evaluation_Table.txt** Contains the results for model evaluation    
 ** InfoModeling.txt** Information of the modeling parameters       
 **Number_Unique_Occurrences.txt** Number of unique occurrences of each species    
 **Occurrences_Cleaned.txt** Datasets produced after occurrences were filtered by selecting a single occurrence per grid-cell       
**Occurrences_Filtered.txt** Datasets produced after occurrences were filtered by selecting corrected sampling spatial bias (thinned occurrences)       
 **Validation_Partition.txt** Information of the evaluation of partial models (e.g. while projecting the model onto the 30% left for test)       
**Thresholds_Algorithm.txt** Information about the thresholds used to create the presence-absence maps for each algorithm (Presence-absence maps are created from the Threshold of complete models)    
**Thresholds_Ensemble.txt** Information about the thresholds used to create the presence-absence maps for ensembled models    

### Last but not least  
There are no defaults! We believe every ecological niche model should be carefully planned and that every decision matters! There are some recommended parameters and literature on which those were based, but they were not included as a default option for the modelling routine.  
Everytime a new version is updated it will accompanied by a log with the version number and what are the main changes/what was included

> Test the package and give us feedback [here](https://github.com/andrefaa/ENMTML/issues) or send an e-mail to andrefaandrade@gmail.com or sjevelazco@gmail.com! 

> Soon we will provide ENMTML package vignette
