# KC_NaChBac_Model

This repository accompanys "Ectopic sodium channel expression paradoxically decreases excitability of Drosophila Kenyon cells." (doi: https://doi.org/10.1113/JP288790)

It contains all code necessary for the development and simulations of the Kenyon cell model.

Scripts are written in matlab for data extraction from electrophysiology recorded data, simulations and skeletonization in python, and proofreading of skeleton model in R. The general outline and order of files is listed below. The packages and appropriate versions for each package is located in the packages file. General flow of modeling pipeline:

1. Generation of gamma-Kenyon cell skeleton model from full adult fly brain EM dataset (skeletonization.py) 
2. Manually proofread skeleton models using Neutube (https://neutracing.com/) / Meshlab (https://www.meshlab.net/) (Proofread.R)
3. Obtain ephys data --------> Process in matlab/extract data for comparisons with simulations (abfload.m to read data, and extracting_data.m to generate csv files)
4. Fitting of passive properties/fitting of NaChBac time constants/current injection simulations/plotting of simultion data (Final_KC_sim_file.py)


To run the scripts through the Windows Subsystem for Linux (WSL), install Ubuntu (https://ubuntu.com/desktop/wsl)
In the ubuntu terminal retrieve the files from Github. 
Change directory to the folder, and run the python files from the terminal. 
