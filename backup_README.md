***********************
# Survey gaps
***********************
### A. Palermino, G. Coro, P. Bove, E.N. Armelloni, G. Scarcella
## Description
This repo contains the workflow developed to fill trawl survey haul gaps based on the paper: "Coro G, Bove P, Armelloni EN, Masnadi F, Scanu M and Scarcella G (2022) Filling Gaps in Trawl Surveys at Sea through Spatiotemporal and Environmental Modelling. Front. Mar. Sci. 9:919339. doi: 10.3389/fmars.2022.919339". 

The tool take into account the spatial mobility of the species throught BIMAC tool written(how far can a population go in one year and which is the mobility of the species?), the historical trends through SSA tool (what biomass do we expect in the missed hauls given the observation history?), the environmental aspect through MaxEnt tool (how inter-annual bio-geo-chemical change affect the species distribution in the survey year?). In the present repo DIVA tool described in he aforementioned paper has been replaced by BIMAC.

# Installation
* Download folder "R", "data" and "java" and store them in a unique folder.
* In "data" replace files with your survey data following "Data preparation".
* File "gebco_30sec_8.asc" contain world depth data needed for BIMAC computation. Do not remove, move or change this file.
* In the subfolder "Envirnmental inputs" replace environmental variables with your variable of interest for MaxEnt tool

# Data preparation
1) Build a csv file for each species on which you want to fill the missing hauls following this format:
      * "station","lat,"lon","year","species","biomindex"
      * e.g.: "1",44.55919167,12.33923333,"2005","MELIKER",35.23
      * name the csv = Genus_species.csv
2) Build a csv file with a list of the missing hauls that have to be calculated following this format:
      * "station","year","lat,"lon"
      * e.g.: "1","2005",44.55919167,12.33923333
      * name the csv = missing_hauls.csv
3) Build a csv file with the swept area by haul and strata for all the hauls: 
      * "station","Stratum","year,"SweptArea"
      * e.g.: "1","STR1_17",2019,0.0124
            "1","STR1_17",2020,0.0127
            "2","STR1_17",2019,0.0186
            "2","STR1_17",2020,0.0181
      * name the csv = HaulData.csv
4) Build a csv file with the strata weight in the area
      * "Stratum","Area","StratumWeight"
      * e.g.: "STR1_17",11361,0.268023969
      * name the csv = StrataWeight.csv
5) Make a folder named "data" and move files "Genus_species.csv"(for each interested species), "missing_hauls.csv", "HaulData.csv" and "StrataWeight.csv" inside this folder:
6) Download environmental variables for the interested years (monthly or daily etc...) in asc format (ALL THE VARIABLES NEED TO HAVE EQUAL RESOLUTION AND EXTENT)
      * e.g.: "CHL_summer_2019.asc"
7) Make a folder named "Environmental_inputs". Make subfolders named "MaxEnt_year" inside the folder "Environmental_inputs"
      * e.g.: "Environmental_inputs/MaxEnt_2019"
8) Move the downloaded environmental variables to the correct folder
      * e.g.: "Environmental_inputs/MaxEnt_2019/CHL_summer_2019.asc"
9) Move the folder "Environmental_inputs" to the folder "data"

## Required data and folder after data preparation
* Folder "data" containing: "missing_hauls.csv", "Genus_species.csv", "HaulData.csv", "StrataWeight.csv",
* Subfolder "Environmental_inputs" containing the asc files of environmental variables to be used for MaxEnt run each year, file "gebco_30sec_8.asc": world depth data for BIMAC computations (DO NOT EDIT)
* Folder "R" containing: Code "BIMAC_no_advection.R" for BIMAC computations (DO NOT EDIT), code "Workflow_Surveygaps.R", code "Workflow_Surveygaps_Solemon.R" with species and year to be filled for Solemon survey
* Folder "java" containing: Code "max_ent_cyb.jar" for MaxEnt computations (DO NOT EDIT), code "ssa.jar" for SSA computations (DO NOT EDIT), folder "cfg" containing: "operators.xlm" for ssa computations (DO NOT EDIT)

# Executing code
Run "Workflow_Surveygaps.R" in RStudio

## Install and load required libraries
library(raster)
library (R2jags)
library (coda)
library(plyr)
library(dplyr)
library(digest)
library(sqldf)

## Set inputs in "Workflow_Surveygaps.R"
* Line 3: feature selection = FALSE/TRUE:
  * if FALSE MaxEnt runs on all environmental variables. 
  * if TRUE MaxEnt select only the variables with a certain level of importance from the first run performingg a second run only on these variables
```
feature_selection=T
```
* line 4: generate a vector with the name of the species to compute.
  * The terms of the species have to be the same as the name of the files in the folder "data":
  * e.g.: folder= "data/Solea_solea.csv"
```
species<-c("Solea_solea","Sepia_officinalis", "Melicertus_keraturus", "Squilla_mantis", "Pecten_jacobeus") 
```
* line 5: generate vector with years.
  * The years have to be the same as the Environmental_inputs subfolders
  * e.g.: folder= "Environmental_inputs/MaxEnt_2019/CHL_summer_2019.asc"
```
years<-c(2019,2020,2021) 
```
# Predetermined settings
* Code "BIMAC_no_advection.R": smooth=F, resolution=0.1, SD compuation= 0.1 (DO NOT EDIT for small basins like the Adriatic Sea), alternative set for global computations: smooth = T and SD = 0.5

* Code "Workflow_Surveygaps.R": 
      * line 99: if haul biomass < 5% of the highest haul biomass level registered in that species and year the haul will be rejected for MaxEnt computation (EDITABLE)
   
      * line 146: if environmental variable x importance < 5% of the highest importance level among environmental variables the environmental variable is 
        rejected for second run maxent refine variables (EDITABLE)
   
      * lines 66-68 and 261-262: if SSA <0 -> SSA = spatial result; SSA = 0 -> SSA = 0; SSA = NaN -> SSA = 0 (DO NOT EDIT)
   
      * line 324: alpha = 1 (weight assigned to the spatial component for HBIE computation) (EDITABLE)
   
      * line 325: beta = 1 (weight assigned to the temporal component for HBIE computation) (EDITABLE)
   
      * line 326: penalty = 0.4 (penalty assign to the HBIE biomass results in case of hauls with ecological values < percent omission rate (EDITABLE)


