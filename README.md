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

The tool take into account the spatial mobility of the species through
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

To apply “Survey_gap” to your own case study you need to provide the
following input files:

### Genus_species.csv

This file contains the available information regarding your target
species index (e.eg. abundance and biomass) for your study area. The
file needs to have the following structure:

    ##   station      lat      lon year species biomindex abunindex
    ## 1       1 44.55919 12.33923 2005 SOLEVUL 307.58783 6807.9443
    ## 2       2 44.58757 12.52665 2005 SOLEVUL  62.66961  716.2242
    ## 3       7 44.28683 12.38482 2005 SOLEVUL  55.14626 1246.4564
    ## 4       8 44.33746 12.50063 2005 SOLEVUL  94.46218 1208.7578
    ## 5      13 44.11446 12.54723 2005 SOLEVUL  54.91112 1402.5491
    ## 6      14 44.12538 12.70341 2005 SOLEVUL  84.19532  996.2042

where:

- **station**: identifier of the survey station. Needs to be consistent
  over the years, but admits gaps.
- **lat**: Latitude of the station in decimal degree
- **lon**: Longitude of the station in decimal degree
- **year**: year of the survey during which the station was carried out
- **species**: name of the species in the preferred format. Needs to be
  consistent over the years
- **biomindex**: biomass index of the station for the interested species
- **abunindex**: abundance index of the station for the interested
  species

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

- **station**: identifier of the survey station. Needs to be consistent
  over the years, but admits gaps.
- **year**: year of the survey during which the station was not carried
  out.
- **lat**: Latitude of the station in decimal degree
- **lon**: Longitude of the station in decimal degree

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

- **Station**: identifier of the survey station. Needs to be consistent
  over the years and includes all the stations (missing and not missing
  stations)
- **Stratum**: identifier of the stratum used for the computation of
  total biomass index of the area.
- **year**: year of the survey during which the station was carried out
  or not. Need to include the interested year to be filled
- **SweptArea**: Swept area of the station

### StrataWeight.csv

This file contains information on the surveyed area:

    ##    Stratum  Area StratumWeight
    ## 1  STR1_17 11361   0.268023969
    ## 2  STR2_17  8410   0.198405209
    ## 3  STR3_17 22466   0.530008493
    ## 4 STR1_SLO   151   0.003562329

where:

- **Stratum**: identifier of the stratum used for the computation of
  total biomass index of the area.
- **Area**: Total area covered by each stratum  
- **StratumWeight**: Weight of the stratum on the total surveyd area

### Environmental_inputs/variable.asc

Download the raster layer file in asc format for the interested
environmental variables:

    ##  [1]        NA        NA        NA        NA        NA        NA        NA
    ##  [8]        NA        NA        NA        NA        NA        NA 0.4100096
    ## [15] 0.3727627        NA 0.5278887 0.4521504 0.3670816 0.2999293 0.6281053
    ## [22] 0.5507973 0.3928404 0.2796727 0.2383853

    ## class      : RasterLayer 
    ## dimensions : 39, 62, 2418  (nrow, ncol, ncell)
    ## resolution : 0.1, 0.1  (x, y)
    ## extent     : 12.3, 18.5, 41.9, 45.8  (xmin, xmax, ymin, ymax)
    ## crs        : +proj=longlat +datum=WGS84 +no_defs 
    ## source     : CHL_summer_2019.asc 
    ## names      : CHL_summer_2019

where:

- **file “variable.asc”**: a matrix containing the values of
  environmental variable for each cell on the map. 0 values have to be
  set = -9999 in asc format in order to obtain NA.
- **dimensions**: identifier of matrix size. nrow, ncol and ncell need
  to be consistent across all the environmental variables
- **resolution**: identifier the resolution of the map determining the
  number of cells
- **extent**: extent of the map in longitude and latitude. Need to be
  consistent across all the environmental variables

## Required data and folder after data preparation

- Folder *“data”* containing: *“missing_hauls.csv”*,
  *“Genus_species.csv”*, *“HaulData.csv”*, *“StrataWeight.csv”*,
- Subfolder *“Environmental_inputs”* containing the asc files of
  environmental variables to be used for MaxEnt run each year, file
  *“gebco_30sec_8.asc”*: world depth data for BIMAC computations (DO NOT
  EDIT)
- Folder *“R”* containing: Code *“BIMAC_no_advection.R”* for BIMAC
  computations (DO NOT EDIT), code *“Workflow_Surveygaps.R”*, code
  *“Workflow_Surveygaps_Solemon.R”* with species and year to be filled
  for Solemon survey
- Folder *“java”* containing: Code *“max_ent_cyb.jar”* for MaxEnt
  computations (DO NOT EDIT), code *“ssa.jar”* for SSA computations (DO
  NOT EDIT), folder *“cfg”* containing: *“operators.xlm”* for ssa
  computations (DO NOT EDIT)

# Executing code

Open *“Workflow_Surveygaps.R”* in RStudio and set working directory to
“source file location”. In *“Workflow_Surveygaps_Solemon.R”*, data are
already selected for the unusable years (2005 and 2006) and hauls (“ms”
and “bis”)

## Set inputs in “Workflow_Surveygaps.R”

- Line 5: feature selection = FALSE/TRUE:

      # feature_selection=T

  - if FALSE all the environmental variables stored in the subfolder are
    used during MaxEnt computation.
  - if TRUE only the environmental variables with a certain level of
    importance (set at the 95% interval of the highest importance level
    among environmental variables) after the first run are included in
    MaxEnt computation

- line 6: generate a vector with the name of the species to compute.

<!-- -->

    ## [1] "Solea_solea"

    -   Terms of the species have to be consistent  with the file names in folder "data":

- line 7: generate vector with years.

<!-- -->

    ## [1] 2019 2020 2021 2022

    -   Years have to be consistent with Environmental_inputs subfolders

- line 8: select index to compute (e.g. abunindex or biomindex)

<!-- -->

    ## [1] "biomindex"

# Processing

Run the code from “Source”

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

## Index computation

In this step the total index of the surveyed area is computed through
the stratified mean method described in Coro et al. (2022) from the
whole hauls including the filled missing hauls.

- Inputs are stored in folder named “biom_index_input” by species and
  year
- Output are stored in folder named “biom_index_output” by species and
  year.

# Outputs and results

At the end of the run all inputs and outputs are stored in the following
folder tree structure per species

<figure>
<img
src="C:/Users/a.palermino/OneDrive%20-%20CNR/github/Survey_gaps/example_outputs/folder_structure.png"
alt="Folder tree structure" />
<figcaption aria-hidden="true">Folder tree structure</figcaption>
</figure>

The main results are reported in:

**hbie_output**

    ##   station      lon      lat ecological  spatial  temporal spatialerror
    ## 1      30 13.95125 43.94125 0.05007590 40.00000  0.000000    15.737645
    ## 2      69 14.24375 42.67375 0.06256380 66.83372  0.000000    16.638121
    ## 3      70 14.24375 42.57625 0.08174660 77.70821 15.008414    21.537151
    ## 4      71 14.63375 42.28375 0.04707840 52.24412  1.212261    12.507442
    ## 5      72 15.02375 42.28375 0.00784261 47.52965  0.000000     9.081131
    ## 6      73 14.92625 42.08875 0.04970930 48.73597 23.956892    19.635370
    ##   index_est_mean index_est_low index_est_high
    ## 1       20.00000      0.000000       43.60647
    ## 2       33.41686      8.459682       58.37404
    ## 3       46.35831     14.052583       78.66404
    ## 4       26.72819      7.967027       45.48935
    ## 5        9.50593      4.057252       14.95461
    ## 6       36.34643      6.893373       65.79948

where:

- **ecological** depicts the results of MaxEnt
- **spatial** and **spetail_error** depict the results of BIMAC
- **temporal** depicts the results of SSA
- **index_est_mean**, **index_est_low**, **index_est_high** depict the
  results of hbie

**biomindex_output**

    ## # A tibble: 4 × 2
    ##   method          biomindex
    ##   <chr>               <dbl>
    ## 1 biomindex_known      74.6
    ## 2 biomindex_mean       69.5
    ## 3 biomindex_low        65.2
    ## 4 biomindex_high       73.9

where:

- **biomass_known** depicts the biomass computation from the known hauls
  (excluding the missing hauls)
- **biomass_mean**, **biomass_mean**, **biomass_mean** depict the
  biomass computations from the estimated and known hauls (including the
  missing hauls)

At the end of each run the biomindex_outptut table is shown along with a
plot showing the index estimated per missing station (from hbie_output)
and Area ( from biomindex_output). Both are stored in biomindex_output
folder

<figure>
<img
src="C:/Users/a.palermino/OneDrive%20-%20CNR/github/Survey_gaps/example_outputs/Solea_solea_final_plot_2021.png"
alt="Example of final plot" />
<figcaption aria-hidden="true">Example of final plot</figcaption>
</figure>
