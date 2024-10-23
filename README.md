# Subgroup_Analysis_Using_Exchangeability

This repository contains all of the necessary code to reproduce the results from "Subgroup Analysis Using Exchangeability Probabilities" as well as all simulated data and figures. 

- Test_Functions.R: This script contains all of the functions for calculating the weights for MEM methods as well as a wrapper function to run all three methods (MEM, MEMr, and a standard test) for any of the scenarios described in the paper. The main function is called RunModels(). Further explanation of this script can be found in its header. 

- Run_Simulations.R: This script runs the simulations in parallel using the snowfall package. The simulation function calls the RunModels() function from the Test_Functions.R script. The output of this script, simulation_results_raw.csv, is a data set containing one row for every randomly generated data set for every simulation scenario. The scenarios with unequal sample sizes are also found in this script. Further explanation of this script can be found in its header. 

- Figures.R: This script generates and saves all of the heatmaps based on the results in simulation_results_raw.csv. Further explanation of this script can be found in its header. 