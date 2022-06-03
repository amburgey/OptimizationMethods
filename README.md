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


## Details of Article
Details of this work can be found in the published journal article on this topic:

Amburgey SA, Yackel Adams AA, Gardner B, Siers S, Converse SJ (in prep) Optimizing monitoring scenarios for an invasive predator. XXXX.

## How to Use This Repository
