# Subgroup_Analysis_Using_Exchangeability

This repository contains all of the necessary code to reproduce the results from "Subgroup Analysis Using Exchangeability Probabilities" as well as all simulated data and figures. 

- Test_Functions.R: This script contains all of the functions for calculating the weights for MEM methods as well as a wrapper function to run all three methods (MEM, MEMr, and a standard test) for any of the scenarios described in the paper. The main function is called RunModels(). Further explanation of this script can be found in its header. 

- Simulation_Function.R: This script contains the function, RunSim(), for parallelizing the simulations using the snowfall package. The simulation function calls the RunModels() function from the Test_Functions.R script. Further explanation of this script can be found in its header. 

- Run_Simulations.R: This script calls RunModels() for all simulation scenarios and outputs the results. The outputs of this script, simulation_results_raw.csv and simulation_results_unequaln_raw.csv, are data sets containing one row for every randomly generated data set for every simulation scenario. They are not available on GitHub because it is too large. The summary data sets are available on GitHub and sufficient to produce all of the heatmaps. Further explanation of this script can be found in its header. 

- Run_Additional_Simulations.R: This script is similar to Run_Simulations.R and runs the binary outcome settings with higher baseline probabilities (and thus higher variance).

- Figures.R: This script generates and saves all of the heatmaps based on the results in simulation_results_summary.csv and simulation_results_unequaln_summary.csv. Further explanation of this script can be found in its header. 

- RecommendedUse.R: This script creates the recommended use tables based on the summary data. 

- CaseStudies.R: This script contains the code used to run simulations for the case studies to get specific type 1 error rates and powers. It also runs the models using the RunModels() function to get p-values and probabilities of exchangeability for each case study. The functions in this script can be used to run any one/two-arm binary/continuous case study.  Further explanation of this script can be found in its header. 