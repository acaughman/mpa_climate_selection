# Climate change reduces long-term population benefits from no-take marine protected areas through selective pressures on species movement

This repository contains code used in the paper: Caughman, A.M., Gaines, S.D., Bradley, D. (2024). Climate change reduces long-term population benefits from no-take marine protected areas through selective pressures on species movement. Global Change Biology. in Press.

The associated zenodo repository is located at https://doi.org/10.5281/zenodo.10778829

For any questions, comments, or concerns, please contact Alicia Caughman [acaughman@bren.ucsb.edu](acaughman@bren.ucsb.edu).

# Instructions

The order of running scripts should be as follows: 

1. The first script is the simulation model `01_final_model.R`. It will produce a CSV file with the model data.
2. The next script to be run after all simulation models is called `02_data_merge.Rmd`. This script merges all simulation data into one CSV for additional analysis and plotting.
3. The final script to be run is `03_analysis_figs.Rmd`. It contains all additional analysis and figure generation for the paper.
4. The fourth script in the folder `04_supplemental.Rmd` contains all code for analysis and figure generation in the supplemental.

the Data folder contains the RData (.rda) outputs from simulations used in the main text of Caughman et al. `model_list.xlsx` contains a list of parameters used in each model.

# Repository Structure

## Overview

```
data
  |__ model_list.xlsx
  |__ model1.rda
  .
  .
  .
  |__ model50.rda
  scripts
  |__ 01_final_model.R
  |__ 02_data_merge.Rmd
  |__ 03_analysis_figs.Rmd
  |__ 04_supplemental.Rmd
figs
  |__ supp
    |__ supplemental figures
  |__ fig1.pdf
  |__ fig2.pdf
  |__ fig3.pdf
  |__ fig4.pdf
```

# R Version

All code was run using R version 4.1.3

## Required Libraries

+ Data Ingestion, Cleaning, Harmonization, and Organization
  - `tidyverse` (version 1.3.2)
  - `here` (version 1.0.1)
+ Simulation
  - `pracma` (version 2.4.2)
+ Data Visualization
  - `patchwork` (version 1.1.2)
