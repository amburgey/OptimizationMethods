# Optimization of Monitoring Methods

The goal of this project was to identify optimal monitoring scenarios for brown treesnakes while explicitly considering tradeoffs between the accuracy and precision of abundance estimates and the cost of conducting each monitoring scenario. This objective-driven optimization was integral to managers making cost-effective decisions regarding the suppression of this invasive species. We 1) analyzed real monitoring data collected in a population of brown treesnakes, 2) used estimates of key parameters describing the performance of these methods to seed monitoring simulations, 3) generated a list of monitoring scenarios of interest that varied the application of monitoring method by spatial and temporal deployment in different brown treesnake populations, and 4) assessed how each scenario performed via comparison by several metrics. We used spatial capture-recapture (SCR) analysis for analyzing real data and simulating and analyzing alternate scenarios.

Monitoring methods of interest were visual surveys (VIS; involving a surveyor hand-capturing snakes) or live trapping surveys (TRAP; involving a live mouse lure inside a snake trap). Snakes could be detected multiple times an evening with visual surveys but remained in traps overnight when using live trapping.

The alternative monitoring scenarios we assessed varied by:
Spatial application (half or the full number of transects)
Temporal application (14 or 30 days)
Monitoring method (VIS, TRAP, or VISTRAP - a combination of both)
Snake density (Normal [120] or reduced [60] snakes in the target area)
Snake size distribution (larger proportion of smaller- or larger-sized snakes)

## Table of Contents

Folder structure is divided into files pertaining to the analysis of real data ("Real Data Analysis") and the simulation and analysis of alternate scenarios ("Simulations").

### Within Real Data Analysis:
* *Figures* - supplementary figures that show information for both visual and trapping survey
* *Scripts* - scripts for creating supplementary figures and study area specification for both visual and trapping surveys
* *Trapping* - all files specific to trapping surveys
  + DataPrep - scripts to prepare individual project datasets for SCR analysis
  + Models - JAGS model script for trapping analysis
  + Results - need to be downloaded from ScienceBase, details in folder
  + Scripts - includes all scripts to pull data and check for errors (Select&PrepTrapData.R) and individual scripts for each file to format matrices for SCR analyses (SCRScriptPrep_projectname); the main file to focus on to run the analysis is **SCRScriptAnalysis_ALLTRAPCP.R**
* *Visual surveys* - all files specific to visual surveys
  + DataPrep - scripts to prepare individual project datasets for SCR analysis
  + Models - JAGS model script for visual survey analysis
  + Results - need to be downloaded from ScienceBase, details in folder
  + Scripts - includes all scripts to pull data and check for errors (Select&PrepVisualData.R) and individual scripts for each file to format matrices for SCR analyses (SCRScriptPrep_projectname); the main file to focus on to run the analysis is **SCRScriptAnalysis_ALLVISCP.R**

### Within Simulations:
* *Figures* - figures showing simulation transect designs and Pareto optimality for full, half, and no startup cost scenarios
* *Models* - JAGS model script for visual, trapping, and visual + trapping combined survey analyses
* *Scripts* - includes all scripts for simulating, checking, and plotting
  + BaseSimulationFramework.R - main script for specifying a monitoring scenario to simulate and analyze
  + CheckRhats.R - script to check for failed (Rhat > 1.1) simulations after analysis
  + FunctionsForSimulation_ClosedAndOneWayBarrier.R - functions used in BaseSimulationFramework
  + IdentifyingOptimParetoFrontier.R - script for highlighting Pareto efficient scenarios and plotting 4 panel figures for manuscript and supplementary material
  + PlottingSimulationStudyAreas.R - script for making figure of simulation transect designs
  + PostSimulationCalculations.R - script for calculating performance metrics (e.g., root mean square error) of all scenarios
* *simDat* - need to download cost model data file from ScienceBase, details in folder

## Required Packages and Versions Used
All code run in R Version 4.1.0 and on a PC

plyr_1.8.6
lubridate_1.7.10
dplyr_1.0.7
stringr_1.4.0
coda_0.19-4
parallel_4.1.0
runjags_2.2.0-3
tidyverse_1.3.1
tibble_3.1.2
reshape2_1.4.4
sf_1.0-1
raster_3.4-13
spatialEco_1.3-7
maptools_1.1-1
rgeos_0.5-5
sp_1.4-5
jagsUI_1.5.2
secr_4.4.5
ggplot2_3.3.5
HDInterval_0.2.2
ggpubr_0.4.0
RColorBrewer_1.1-2
forcats_0.5.1
purrr_0.3.4
scales_1.1.1
readr_2.0.0
patchwork_1.1.1

## Details of Article
Details of this work can be found in the published journal article on this topic:

Amburgey SA, Yackel Adams AA, Gardner B, Siers S, Converse SJ (in prep) Optimizing monitoring scenarios for an invasive predator. XXXX.

## How to Use This Repository
To run analyses of real data:
1. Go to Real Data Analysis folder
2. Go to Trapping folder
3. Go to Scripts folder
4. Run SCRScriptAnalysisALLTRAPCP.R
5. Repeat for Visual survey -> Scripts -> run SCRScriptAnalysisALLVISCP.R

To make figures of results of real data analyses:
1. Go to Real Data Analysis folder
2. Go to Scripts folder
3. Select and run plotting script

To simulate data and analyze simulated datasets:
1. Go to Simulations folder
2. Go to Scripts folder
3. Run BaseSimulationFramework.R with scenario specified (carefully read instructions in script)

After simulating and analyzing datasets:
1. Go to Simulations folder
2. Go to Scripts folder
3. Run CheckRhats.R to see which simulations need to be re-run or continue to fail
4. Run PostSimulationCalculations.R to calculate performance metrics

To plot results of simulations:
1. Go to Simulations folder
2. Go to Scripts folder
3. Run IdentifyingOptimParetoFrontier.R to make figures of Pareto optimality
4. Run PlottingSimulationStudyAreas.R to make figure of study areas
