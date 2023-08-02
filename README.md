# Survey_gaps
This repo contains the workflow developed to easely fill the survey haul gaps

#DATA PREPARATION
1) Build a csv file for each species on which you want to fill the missing hauls following this format:
      "station","lat,"lon","year","species","biomindex"
e.g.: "1",44.55919167,12.33923333,"2005","MELIKER",35.23
      name the csv = Genus_species.csv
	  
2) Build a csv file with a list of the missing hauls that have to be calculated following this format:
      "station","year","lat,"lon"
e.g.: "1","2005",44.55919167,12.33923333
      name the csv = missing_hauls.csv
	  
3) Build a csv file with the swept area by haul and strata for all the hauls 
      "station","Stratum","year,"SweptArea"
e.g.: "1","STR1_17",2019,0.0124
      "1","STR1_17",2020,0.0127
      "2","STR1_17",2019,0.0186
      "2","STR1_17",2020,0.0181
      name the csv = HaulData.csv

4) Build a csv file with the strata weight in the area
      "Stratum","Area","StratumWeight"
e.g.: "STR1_17",11361,0.268023969
      name the csv = StrataWeight.csv
	  
5) Make a folder named "data" and move files "Genus_species.csv"(for each interested species), "missing_hauls.csv", "HaulData.csv" and "StrataWeight.csv" inside this folder:

6) Download environmental variables for the interested years (monthy or daily etc...) in asc format (ALL THE VARIABLES NEED TO HAVE EQUAL RESOLUTION AND EXTENT)
e.g.: "CHL_summer_2019.asc"

7) Make a folder named "Environmental_inputs". Make subfolders named "MaxEnt_year" inside the folder "Environmental_inputs"
e.g.: "Environmental_inputs/MaxEnt_2019"

8) Move the downloaded environmental variables in the correct folder
e.g.: "Environmental_inputs/MaxEnt_2019/CHL_summer_2019.asc"

9) Move the folder "Environmental_inputs" in the folder "data"

#DIRECTORY CONTENT AFTER DATA PREPARATION
-Folder "data" containing: "missing_hauls.csv", "Genus_species.csv", "HaulData.csv", "StrataWeight.csv", folder "Environmental_inputs"containing the asc files of environmental variables to be used for MaxEnt run each year, file "gebco_30sec_8.asc": world depth data for BIMAC computations (DO NOT EDIT)
-Folder "R" containing: Code "BIMAC_no_advection.R" for BIMAC computations (DO NOT EDIT), code "Workflow_Surveygaps.R", code "Workflow_Surveygaps_Solemon.R" with species and year to be filled for Solemon survey
-Folder "java" containing: Code "max_ent_cyb.jar" for MaxEnt computations (DO NOT EDIT), code "ssa.jar" for SSA computations (DO NOT EDIT), folder "cfg" containing: "operators.xlm" for ssa computations (DO NOT EDIT)

#REQUIRED LIBRARIES:
library(raster)
library (R2jags)
library (coda)
library(plyr)
library(dplyr)
library(digest)
library(sqldf)

#SET INPUTS IN "Workflow_Surveygaps.R":
-line 3: Feature selection = FALSE/TRUE: if FALSE MaxEnt run on all environmental variables. 
                                         if TRUE MaxEnt select only the variables with a certain level of importance from the first run performig a second run only on these variables
-line 4: generate vector with the name of the species to compute. The names of the species have to be the same of the name of the files in the folder "data"
         e.g.: folder= "data/Solea_solea.csv"
		       vector= species<-c("Solea_solea","Sepia_officinalis", "Melicertus_keraturus", "Squilla_mantis", "Pecten_jacobeus") 
-line 5: generate vector with years. The years have to be the same as the Environmental_inputs subfolders
         e.g.: folder= "Environmental_inputs/MaxEnt_2019/CHL_summer_2019.asc"
		       vector= years<-c(2019,2020,2021) 

#PREDETERMINED SETTINGS:
-Code "BIMAC_no_advection.R": smooth=F, resolution=0.1, SD compuation= 0.1 (DO NOT EDIT for small basins like Adriatic Sea), alternative set for global computations: smooth = T and SD = 0.5
-Code "Workflow_Surveygaps.R": 
   -line 99: if haul biomass < 5% of the highest haul biomass level registered in that species and year the haul will be rejected for MaxEnt computation (EDITABLE)
   -line 146: if environmental variable x importance < 5% of the highest importance level among environmental variables the environmental variable will be rejected for second run maxent refine variables (EDITABLE)
   -lines 66-68 and 261-262: if SSA <0 -> SSA = spatial result; SSA = 0 -> SSA = 0; SSA = NaN -> SSA = 0 (DO NOT EDIT)
   -line 324: alpha = 1 (weight assign to the spatial component for HBIE computation) (EDITABLE)
   -line 325: beta = 1 (weight assign to the temporal component for HBIE computation) (EDITABLE)
   -line 326: penalty = 0.4 (penalty assign to the HBIE biomass results in case of hauls with ecological values < percent omission rate (EDITABLE)

#RUN code "Workflow_Surveygaps.R"


ECC....

