Survey gaps
================
A. Palermino, G. Coro, P. Bove, E.N. Armelloni, G. Scarcella,
2023

## Overview

This repo contains the workflow developed to fill trawl survey haul gaps
based on the paper: “Coro G, Bove P, Armelloni EN, Masnadi F, Scanu M
and Scarcella G (2022) Filling Gaps in Trawl Surveys at Sea through
Spatiotemporal and Environmental Modelling. Front. Mar. Sci. 9:919339.
doi: 10.3389/fmars.2022.919339”.

The tool atake into account the spatial mobility of the species through
BIMAC tool (how far can a population go in one year and which is the
mobility of the species?), the historical trends through SSA tool (what
biomass do we expect in the missed hauls given the observation
history?), the environmental aspect through MaxEnt tool (how
inter-annual bio-geo-chemical change affect the species distribution in
the survey year?). Finally through the HBIE tool,the results of the
aforementioned tools are weighted and gather to obtain the final biomass
index of missing hauls. In the present repo DIVA tool described in Coro
et al. (2022) has been replaced by BIMAC.

# Installation

- Download folder “R”, “data” and “java” and store them in a unique
  folder.

## Data folders content (provided at installation)

- Solea_solea.csv; HaulData.csv; StrataWeight.csv: example of input data
- gebco_30sec_8.asc: world depth data needed for BIMAC computation. Do
  not remove, move or change this file.
- Subfolder “Envirnmental inputs”: contains environmental variables of
  interest for MaxEnt tool. Data provided refers to the Adriatic Sea
  from 2019 to 2022

## Data preparation guideline

To apply “Survey_gap” to your own case study you need to provide
following input files:

### Genus_species.csv

This file contains the available information regarding your target
species abundance and biomass for your study area. The file needs to
have the following structure:

    ##   station      lat      lon year species biomindex abunindex
    ## 1       1 44.55919 12.33923 2005 SOLEVUL 307.58783 6807.9443
    ## 2       2 44.58757 12.52665 2005 SOLEVUL  62.66961  716.2242
    ## 3       7 44.28683 12.38482 2005 SOLEVUL  55.14626 1246.4564
    ## 4       8 44.33746 12.50063 2005 SOLEVUL  94.46218 1208.7578
    ## 5      13 44.11446 12.54723 2005 SOLEVUL  54.91112 1402.5491
    ## 6      14 44.12538 12.70341 2005 SOLEVUL  84.19532  996.2042

where:

- station: identifier of the survey station. Needs to be consistent over
  the years.

- lat: Latitude of the station in decimal degree

- lon: Longitude of the station in decimal degree

- year: year of the survey during which the station was carried out

- species: name of the species in the preferred format. Needs to be
  consistent over the years

- biomindex: biomass index of the station for the interested species

- abunindex: abundance index of the station for the interested species

### missing_hauls.csv

This file contains a list of missing hauls that have to be computed with
relative information. The file needs to have the following structure:

    ##   station year      lat      lon
    ## 1       4 2022 45.09887 12.38338
    ## 2       5 2022 45.15225 12.50874
    ## 3       6 2022 45.25585 12.87548
    ## 4      10 2020 45.09946 13.25025
    ## 5      10 2022 45.09946 13.25025
    ## 6      14 2022 44.73999 13.08566

where:

- station: identifier of the survey station. Needs to be consistent over
  the years.

- year: year of the survey during which the station was carried out.

- lat: Latitude of the station in decimal degree

- lon: Longitude of the station in decimal degree

### HaulData.csv

This file contains the information on the surveyed area. In case of
missing haul use the swept area of the previous year:

    ##   Station Stratum year  SweptArea
    ## 1       1 STR1_17 2019 0.03599334
    ## 2       1 STR1_17 2020 0.02437610
    ## 3       1 STR1_17 2021 0.02349296
    ## 4       1 STR1_17 2022 0.02925132
    ## 5      10 STR2_17 2019 0.03656774
    ## 6      10 STR2_17 2020 0.03656774

where:

- Station: identifier of the survey station. Needs to be consistent over
  the years and includes all the station (missing and not missing
  station)

- Stratum: identifier of the stratum used for the computation of total
  biomass index of the area.

- year: year of the survey during which the station was carried out.
  Need to include the interested year to be filled

- SweptArea: Swept area of the station

### StrataWeight.csv

This file contains information on the surveyed area:

    ##    Stratum  Area StratumWeight
    ## 1  STR1_17 11361   0.268023969
    ## 2  STR2_17  8410   0.198405209
    ## 3  STR3_17 22466   0.530008493
    ## 4 STR1_SLO   151   0.003562329

where:

- Stratum: identifier of the stratum used for the computation of total
  biomass index of the area.

- Area: Total area covered by each stratum

- StratumWeight: Weight of the stratum on the total surveyd area

### Environmental_inputs/variable.asc

Download the raster layer file in asc format for the interested
environmental variables:

    ##            1         2         3         4         5         6         7
    ## 1         NA        NA        NA        NA        NA        NA        NA
    ## 2         NA        NA        NA        NA        NA 0.3571303 0.3275829
    ## 3         NA        NA        NA 0.4100096 0.3727627 0.3131888 0.2605357
    ## 4         NA 0.5278887 0.4521504 0.3670816 0.2999293 0.2466989 0.2129933
    ## 5  0.6281053 0.5507973 0.3928404 0.2796727 0.2383853 0.2186371 0.2056990
    ## 6  0.7750720 0.6094998 0.3925024 0.2695370 0.2414952 0.2330638 0.2263396
    ## 7         NA 0.7572239 0.5369453 0.3741204 0.3118363 0.2832875 0.2697312
    ## 8         NA        NA 0.7217675 0.6173192 0.4812652 0.3740655 0.3183495
    ## 9         NA        NA 1.1843615 0.8171226 0.6092055 0.4413580 0.3293392
    ## 10 1.9023403 1.5170109 0.9227595 0.7211673 0.5773040 0.4180491 0.2881048
    ##            8         9        10        11        12         13        14
    ## 1         NA        NA 0.3778360 0.3642117        NA 0.29353559 0.2814722
    ## 2  0.2734011 0.2493080 0.2782385 0.2499041 0.2034301 0.18954992 0.2050904
    ## 3  0.2213961 0.1925389 0.1782477 0.1599235 0.1415365 0.13839088        NA
    ## 4  0.1907711 0.1774113 0.1694480 0.1567162 0.1424684 0.13248047        NA
    ## 5  0.1954958 0.1904794 0.1831382 0.1694365 0.1530065 0.13775611        NA
    ## 6  0.2204139 0.2112434 0.1944607 0.1744659 0.1557369 0.13731483        NA
    ## 7  0.2518665 0.2263726 0.1929361 0.1627345 0.1440033 0.12913613 0.1172576
    ## 8  0.2695169 0.2184501 0.1724046 0.1402283 0.1272029 0.11738839 0.1105228
    ## 9  0.2474165 0.1844187 0.1420508 0.1192505 0.1108688 0.10523679 0.1021336
    ## 10 0.2017244 0.1489750 0.1246744 0.1086271 0.0999753 0.09423658 0.0924265
    ##            15         16         17         18         19         20
    ## 1  0.26235393         NA         NA         NA         NA         NA
    ## 2  0.21940364         NA         NA         NA         NA         NA
    ## 3          NA         NA         NA         NA         NA         NA
    ## 4          NA         NA         NA         NA         NA         NA
    ## 5          NA         NA         NA         NA         NA         NA
    ## 6          NA         NA         NA         NA         NA 0.07716765
    ## 7          NA         NA         NA         NA         NA 0.09092680
    ## 8  0.10219493         NA         NA         NA         NA 0.08775543
    ## 9  0.09902937 0.09241491         NA 0.08147497 0.08392087 0.08477256
    ## 10 0.09231476 0.09056097 0.08469907 0.07934242 0.08213522 0.08132850

    ## class      : RasterLayer 
    ## dimensions : 39, 62, 2418  (nrow, ncol, ncell)
    ## resolution : 0.1, 0.1  (x, y)
    ## extent     : 12.3, 18.5, 41.9, 45.8  (xmin, xmax, ymin, ymax)
    ## crs        : +proj=longlat +datum=WGS84 +no_defs 
    ## source     : CHL_summer_2019.asc 
    ## names      : CHL_summer_2019

where:

- file variable.asc is a matrix containing the values of environmental
  variable for each cell on the map. 0 values have to be set = -9999 in
  asc format in order to obtain NA.

- dimensions: identifier of matrix size. nrow, ncol and ncell need to be
  consistent across all the environmental variables

- resolution: identifier the resolution of the map determining the
  number of cells

- extent: extent of the map in longitude and latitude. Need to be
  consistent across all the environmental variables

## Required data and folder after data preparation

- Folder “data” containing: “missing_hauls.csv”, “Genus_species.csv”,
  “HaulData.csv”, “StrataWeight.csv”,
- Subfolder “Environmental_inputs” containing the asc files of
  environmental variables to be used for MaxEnt run each year, file
  “gebco_30sec_8.asc”: world depth data for BIMAC computations (DO NOT
  EDIT)
- Folder “R” containing: Code “BIMAC_no_advection.R” for BIMAC
  computations (DO NOT EDIT), code “Workflow_Surveygaps.R”, code
  “Workflow_Surveygaps_Solemon.R” with species and year to be filled for
  Solemon survey
- Folder “java” containing: Code “max_ent_cyb.jar” for MaxEnt
  computations (DO NOT EDIT), code “ssa.jar” for SSA computations (DO
  NOT EDIT), folder “cfg” containing: “operators.xlm” for ssa
  computations (DO NOT EDIT)

# Executing code

Open “Workflow_Surveygaps.R” in RStudio and set working directory to
“source file location”. In “Workflow_Surveygaps_Solemon.R”, data are
already selected for the unusable years (2005 and 2006) and hauls (“ms”
and “bis”)

## Install and load required libraries

library(raster) library (R2jags) library (coda) library(plyr)
library(dplyr) library(digest) library(sqldf)

## Set inputs in “Workflow_Surveygaps.R”

- Line 3: feature selection = FALSE/TRUE:

      # feature_selection=T

  - if FALSE all the environmental variables stored in the subfolder are
    used during MaxEnt computation.
  - if TRUE only the environmental variables with a certain level of
    importance (set at the 95% interval of the highest importance level
    among environmental variables) after the first run are included in
    MaxEnt computation

- line 4: generate a vector with the name of the species to compute.

<!-- -->

    ## [1] "Solea_solea"

    -   Terms of the species have to be consistent  with the file names in folder "data":

- line 5: generate vector with years.

<!-- -->

    ## [1] 2019 2020 2021 2022

    -   Years have to be consistent with Environmental_inputs subfolders

# Processing

Run the code from “Source”

README WORK IN PROGRESS

## SSA, BIMAC and MaxEnt computation

In this step several functions are applied to compute temporal, spatial
and ecological tools. The results are stored for HBIE computation

- Inputs are stored in folders named “toolname_input” by species and
  year

- Output are stored in folders named “toolname_output” by species and
  year

## HBIE computation

In this step the temporal (SSA), spatial (BIMAC) and ecological (MaxEnt)
results are pulled together and weighted to give the estimation of mean,
low and high biomass index of missing hauls.

- Inputs are stored in folder named “hbie_input” by species and year

- Output are stored in folder named “hbie_output” by species and year.
  The output csv file contain the final biomass index results per
  station in the column “biomass_est_mean”.

## Biomass computation

In this step the total biomass index of the surveyed area is computed
from the whole hauls including the filled missing hauls.

- Inputs are stored in folder named “biom_index_input” by species and
  year

- Output are stored in folder named “biom_index_output” by species and
  year.

README WORK IN PROGRESS
