#################################################################################
# TITLE: Run_Simulations.R
#
# PURPOSE: Run simulations in parallel using RunSim() func from 
#          Simulation_Function.R and create a master data output file. 
#
#
# OUTPUT: simulation_results_raw.csv - contains information for every simulation
#                                      run including means, probabilities of 
#                                      exchangeability, pvalues, number of arms,
#                                      outcome type, and all other relevant 
#                                      information
#
#         simulation_results_summary.csv - summary of power for each test for 
#                                          the main 0.2 MEM cutoff and 0.05 pval
#
#         simulation_results_summary_08cutoff.csv  - summary of power for each 
#                                                    test for a 0.8 MEM cutoff 
#                                                    and 0.05 pval
#         simulation_results_summary_calibratedcutoff.csv  - summary of power 
#                                                            for each test using
#                                                            MEM cutoffs calib-
#                                                            rated to 5% t1e
#         simulation_results_unequaln_raw.csv - same as above for unequal n
#
#         simulation_results_unequaln_summary.csv - same as above for unequal n
#
#         simulation_results_unequaln_summary_08cutoff.csv  - same as above for 
#                                                             unequal n
#
#         simulation_results_unequaln_summary_calibratedcutoff.csv - same as 
#                                                                    above for 
#                                                                    unequal n
#
# SECTIONS: Section 1
#              run simulations with various parameter settings 
#              1.1 parfunc() - function to create parameter list of lists
#              1.2 One-Arm Binary Outcome
#              1.3 One-Arm Continuous Outcome
#              1.4 Two-Arm Binary Outcome
#              1.5 Two-Arm Continuous Outcome
#           Section 2
#              aggregate data
#
#           Section 3
#              repeat section 1 for unequal sample sizes in groups
#           Section 4
#              repeat section 2 for unequal sample sizes in groups
#
# NOTES: Sections 1-4 took about 9.5 hours to run on my DELL laptop - 
#           Processor - Intel(R) Core(TM) i7-10750H CPU @ 2.60GHz, 2592 Mhz, 
#                       6 Core(s), 12 Logical Processor(s)
#           Installed RAM - 16.0 GB (15.8 GB usable)
#
# AUTHOR: Shannon Thomas
# DATE CREATED: OCT 22, 2024
#################################################################################



##########################################################
#### SECTION 0: LOAD DEPENDENCIES & SOURCE MODEL FUNC ####
##########################################################

source("Simulation_Function.R")


##########################################################
############### SECTION 1:  RUN SIMULATION ############### 
##########################################################


### SECTION 1.1: Function to create parameter list of lists 

parfunc <- function(means, sds = c(NA,NA), trt_effect = c(NA,NA)){
  pars=list(
    list(n=c(100, 100), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(150, 150), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(200, 200), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(250, 250), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(300, 300), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(350, 350), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(400, 400), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(450, 450), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(500, 500), means=means, sds = sds, trt_effect = trt_effect)
  )
  
  return(pars)
}


### SECTION 1.2: One-Arm Binary Outcome

#scenario 1: equal group probabilities (0.1 and 0.1)
system.time(OABin_11 <- RunSim(pars = parfunc(means = c(0.1,0.1)),
                               num_arms = 1, outcome_type = "binary"))

#scenario 2: unequal group probabilities (0.1 and 0.2)
system.time(OABin_12 <- RunSim(pars = parfunc(means = c(0.1,0.2)),
                               num_arms = 1, outcome_type = "binary"))

#scenario 3: unequal group probabilities (0.1 and 0.3)
system.time(OABin_13 <- RunSim(pars = parfunc(means = c(0.1,0.3)),
                               num_arms = 1, outcome_type = "binary"))

#scenario 4: unequal group probabilities (0.1 and 0.4)
system.time(OABin_14 <- RunSim(pars = parfunc(means = c(0.1,0.4)),
                               num_arms = 1, outcome_type = "binary"))

#scenario 5: unequal group probabilities (0.1 and 0.5)
system.time(OABin_15 <- RunSim(pars = parfunc(means = c(0.1,0.5)),
                               num_arms = 1, outcome_type = "binary"))


### SECTION 1.3: One-Arm Continuous Outcome

#scenario 1: equal group means (1 and 1) and sd = 1
system.time(OACont_var1_11 <- RunSim(pars = parfunc(means = c(1,1), sds = c(1,1)),
                                num_arms = 1, outcome_type = "continuous"))

#scenario 2: unequal group means (1 and 2) and sd = 1
system.time(OACont_var1_12 <- RunSim(pars = parfunc(means = c(1,2), sds = c(1,1)),
                                num_arms = 1, outcome_type = "continuous"))

#scenario 3: unequal group means (1 and 3) and sd = 1
system.time(OACont_var1_13 <- RunSim(pars = parfunc(means = c(1,3), sds = c(1,1)),
                                num_arms = 1, outcome_type = "continuous"))

#scenario 4: unequal group means (1 and 4) and sd = 1
system.time(OACont_var1_14 <- RunSim(pars = parfunc(means = c(1,4), sds = c(1,1)),
                                num_arms = 1, outcome_type = "continuous"))

#scenario 5: unequal group means (1 and 5) and sd = 1
system.time(OACont_var1_15 <- RunSim(pars = parfunc(means = c(1,5), sds = c(1,1)),
                                num_arms = 1, outcome_type = "continuous"))

#scenario 6: equal group means (1 and 1) and sd = sqrt(10)
system.time(OACont_var10_11 <- RunSim(pars = parfunc(means = c(1,1), c(sqrt(10),sqrt(10))),
                                num_arms = 1, outcome_type = "continuous"))

#scenario 7: unequal group means (1 and 2) and sd = sqrt(10)
system.time(OACont_var10_12 <- RunSim(pars = parfunc(means = c(1,2), sds = c(sqrt(10),sqrt(10))),
                                num_arms = 1, outcome_type = "continuous"))

#scenario 8: unequal group means (1 and 3) and sd = sqrt(10)
system.time(OACont_var10_13 <- RunSim(pars = parfunc(means = c(1,3), sds = c(sqrt(10),sqrt(10))),
                                num_arms = 1, outcome_type = "continuous"))

#scenario 9: unequal group means (1 and 4) and sd = sqrt(10)
system.time(OACont_var10_14 <- RunSim(pars = parfunc(means = c(1,4), sds = c(sqrt(10),sqrt(10))),
                                num_arms = 1, outcome_type = "continuous"))

#scenario 10: unequal group means (1 and 5) and sd = sqrt(10)
system.time(OACont_var10_15 <- RunSim(pars = parfunc(means = c(1,5), sds = c(sqrt(10),sqrt(10))),
                                num_arms = 1, outcome_type = "continuous"))

#scenario 11: equal group means (1 and 1) and sd = 10
system.time(OACont_var100_11 <- RunSim(pars = parfunc(means = c(1,1), c(10,10)),
                                      num_arms = 1, outcome_type = "continuous"))

#scenario 12: unequal group means (1 and 2) and sd = 10
system.time(OACont_var100_12 <- RunSim(pars = parfunc(means = c(1,2), sds = c(10,10)),
                                      num_arms = 1, outcome_type = "continuous"))

#scenario 13: unequal group means (1 and 3) and sd = 10
system.time(OACont_var100_13 <- RunSim(pars = parfunc(means = c(1,3), sds = c(10,10)),
                                      num_arms = 1, outcome_type = "continuous"))

#scenario 14: unequal group means (1 and 4) and sd = 10
system.time(OACont_var100_14 <- RunSim(pars = parfunc(means = c(1,4), sds = c(10,10)),
                                      num_arms = 1, outcome_type = "continuous"))

#scenario 15: unequal group means (1 and 5) and sd = 10
system.time(OACont_var100_15 <- RunSim(pars = parfunc(means = c(1,5), sds = c(10,10)),
                                      num_arms = 1, outcome_type = "continuous"))


### SECTION 1.4: Two-Arm Binary Outcome

#scenario 1: equal group probabilities (0.1 and 0.1) and treatment effects (0.1 and 0.1)
system.time(TABin_11 <- RunSim(pars = parfunc(means = c(0.1,0.1), trt_effect = c(0.1,0.1)),
                               num_arms = 2, outcome_type = "binary"))

#scenario 2: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.2)
system.time(TABin_12 <- RunSim(pars = parfunc(means = c(0.1,0.1), trt_effect = c(0.1,0.2)),
                               num_arms = 2, outcome_type = "binary"))

#scenario 3: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.3)
system.time(TABin_13 <- RunSim(pars = parfunc(means = c(0.1,0.1), trt_effect = c(0.1,0.3)),
                               num_arms = 2, outcome_type = "binary"))

#scenario 4: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.4)
system.time(TABin_14 <- RunSim(pars = parfunc(means = c(0.1,0.1), trt_effect = c(0.1,0.4)),
                               num_arms = 2, outcome_type = "binary"))

#scenario 5: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.5)
system.time(TABin_15 <- RunSim(pars = parfunc(means = c(0.1,0.1), trt_effect = c(0.1,0.5)),
                               num_arms = 2, outcome_type = "binary"))


### SECTION 1.5: Two-Arm Continuous Outcome

#scenario 1: equal group means (1 and 1) and treatment effects (1,1) and sd = 1
system.time(TACont_var1_11 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,1), sds = c(1,1)),
                                     num_arms = 2, outcome_type = "continuous"))

#scenario 2: unequal group means (1 and 1) and unequal treatment effects (1 and 2) and sd = 1
system.time(TACont_var1_12 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,2), sds = c(1,1)),
                                     num_arms = 2, outcome_type = "continuous"))

#scenario 3: unequal group means (1 and 1) and unequal treatment effects (1 and 3) and sd = 1
system.time(TACont_var1_13 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,3), sds = c(1,1)),
                                     num_arms = 2, outcome_type = "continuous"))

#scenario 4: unequal group means (1 and 1) and unequal treatment effects (1 and 4) and sd = 1
system.time(TACont_var1_14 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,4), sds = c(1,1)),
                                     num_arms = 2, outcome_type = "continuous"))

#scenario 5: unequal group means (1 and 1) and unequal treatment effects (1 and 5) and sd = 1
system.time(TACont_var1_15 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,5), sds = c(1,1)),
                                     num_arms = 2, outcome_type = "continuous"))

#scenario 6: equal group means (1 and 1) and treatment effects (1 and 1) and sd = sqrt(10)
system.time(TACont_var10_11 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,1), c(sqrt(10),sqrt(10))),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 7: unequal group means (1 and 1) and unequal treatment effects (1 and 2) and sd = sqrt(10)
system.time(TACont_var10_12 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,2), sds = c(sqrt(10),sqrt(10))),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 8: unequal group means (1 and 1) and unequal treatment effects (1 and 3) and sd = sqrt(10)
system.time(TACont_var10_13 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,3), sds = c(sqrt(10),sqrt(10))),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 9: unequal group means (1 and 1) and unequal treatment effects (1 and 4) and sd = sqrt(10)
system.time(TACont_var10_14 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,4), sds = c(sqrt(10),sqrt(10))),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 10: unequal group means (1 and 1) and unequal treatment effects (1 and 5) and sd = sqrt(10)
system.time(TACont_var10_15 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,5), sds = c(sqrt(10),sqrt(10))),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 11: equal group means (1 and 1) and treatment effects (1 and 1) and sd = 10
system.time(TACont_var100_11 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,1), c(10,10)),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 12: unequal group means (1 and 1) and unequal treatment effects (1 and 2) and sd = 10
system.time(TACont_var100_12 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,2), sds = c(10,10)),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 13: unequal group means (1 and 1) and unequal treatment effects (1 and 3) and sd = 10
system.time(TACont_var100_13 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,3), sds = c(10,10)),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 14: unequal group means (1 and 1) and unequal treatment effects (1 and 4) and sd = 10
system.time(TACont_var100_14 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,4), sds = c(10,10)),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 15: unequal group means (1 and 1) and unequal treatment effects (1 and 5) and sd = 10
system.time(TACont_var100_15 <- RunSim(pars = parfunc(means = c(1,1), trt_effect = c(1,5), sds = c(10,10)),
                                      num_arms = 2, outcome_type = "continuous"))





##########################################################
############### SECTION 2:  AGGREGATE DATA ############### 
##########################################################

allresults <- rbind(OABin_11, OABin_12, OABin_13, OABin_14, OABin_15,
                    OACont_var1_11, OACont_var1_12, OACont_var1_13, OACont_var1_14, OACont_var1_15,
                    OACont_var10_11, OACont_var10_12, OACont_var10_13, OACont_var10_14, OACont_var10_15,
                    OACont_var100_11, OACont_var100_12, OACont_var100_13, OACont_var100_14, OACont_var100_15,
                    TABin_11, TABin_12, TABin_13, TABin_14, TABin_15,
                    TACont_var1_11, TACont_var1_12, TACont_var1_13, TACont_var1_14, TACont_var1_15,
                    TACont_var10_11, TACont_var10_12, TACont_var10_13, TACont_var10_14, TACont_var10_15,
                    TACont_var100_11, TACont_var100_12, TACont_var100_13, TACont_var100_14, TACont_var100_15)


write.csv(allresults, "simulation_results_raw.csv")
#allresults <- read.csv("simulation_results_raw.csv")

#create summary data set using the following cutoffs
#PRIMARY CUTOFF = 0.2
MEM_cutoff <- 0.2
p_cutoff <- 0.05 
nsim <- 10000
allresults_summary <- allresults %>% 
                      group_by(num_arms, outcome_type, n1, n2, mean1, mean2, sd1, sd2, trteff1, trteff2) %>%
                      summarize(MEM_power = sum(MEMpexch < MEM_cutoff)/nsim,
                                MEMr_power = sum(MEMrpexch < MEM_cutoff)/nsim,
                                pval_power = sum(pval < p_cutoff)/nsim,
                                MEM_contapprox_power = sum(MEMapproxpexch < MEM_cutoff)/nsim,
                                pval_chisq_power = sum(pval_chisq < p_cutoff)/nsim)

allresults_summary$MEMcutoff <- MEM_cutoff
allresults_summary$pvalcutoff <- p_cutoff

write.csv(allresults_summary, "simulation_results_summary.csv")

#ALTERNATE CUTOFF = 0.8
MEM_cutoff <- 0.8
p_cutoff <- 0.05 
nsim <- 10000
allresults_summary_08 <- allresults %>% 
  group_by(num_arms, outcome_type, n1, n2, mean1, mean2, sd1, sd2, trteff1, trteff2) %>%
  summarize(MEM_power = sum(MEMpexch < MEM_cutoff)/nsim,
            MEMr_power = sum(MEMrpexch < MEM_cutoff)/nsim,
            pval_power = sum(pval < p_cutoff)/nsim,
            MEM_contapprox_power = sum(MEMapproxpexch < MEM_cutoff)/nsim,
            pval_chisq_power = sum(pval_chisq < p_cutoff)/nsim)

allresults_summary_08$MEMcutoff <- MEM_cutoff
allresults_summary_08$pvalcutoff <- p_cutoff

write.csv(allresults_summary_08, "simulation_results_summary_08cutoff.csv")

#T1E RATE CALIBRATED CUTOFF
MEM_cutoff <- "calibrated"
p_cutoff <- 0.05 
nsim <- 10000

allresults$effsize <- ifelse(allresults$num_arms == 1, 
                             allresults$mean2 - allresults$mean1, 
                             allresults$trteff2 - allresults$trteff1)
cutoffs <- allresults[allresults$effsize == 0,] %>% group_by(n1, n2, outcome_type, num_arms, sd1, sd2) %>%
  summarize(q5_MEM = quantile(MEMpexch, 0.05, na.rm = TRUE), 
            q5_MEMcontapprox = quantile(MEMapproxpexch, 0.05, na.rm = TRUE),
            q5_MEMr = quantile(MEMrpexch, 0.05)) %>% 
  as.data.frame()
allresults_cutoffs <- merge(allresults, cutoffs)

allresults_summary_calibrated <- allresults_cutoffs %>% 
  group_by(num_arms, outcome_type, n1, n2, mean1, mean2, sd1, sd2, trteff1, trteff2) %>%
  summarize(MEM_power = sum(MEMpexch < q5_MEM)/nsim,
            MEMr_power = sum(MEMrpexch < q5_MEMr)/nsim,
            pval_power = sum(pval < p_cutoff)/nsim,
            MEM_contapprox_power = sum(MEMapproxpexch < q5_MEMcontapprox)/nsim,
            pval_chisq_power = sum(pval_chisq < p_cutoff)/nsim)

allresults_summary_calibrated$MEMcutoff <- MEM_cutoff
allresults_summary_calibrated$pvalcutoff <- p_cutoff

write.csv(allresults_summary_calibrated, "simulation_results_summary_calibratedcutoff.csv")


##########################################################
########## SECTION 3: RUN SIMULATION (n1 != n2) ########## 
##########################################################


### SECTION 3.1: Function to create parameter list of lists 

parfunc_unequaln <- function(means, sds = c(NA,NA), trt_effect = c(NA,NA), n2mult){
  pars=list(
    #n = sample sizes, m0 = baseline mean for each source, 
    #p = probability of treatment being assigned, 
    #beta_t = treatment main effect 
    list(n=c(200*(1-n2mult), 200*n2mult), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(300*(1-n2mult), 300*n2mult), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(400*(1-n2mult), 400*n2mult), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(500*(1-n2mult), 500*n2mult), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(600*(1-n2mult), 600*n2mult), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(700*(1-n2mult), 700*n2mult), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(800*(1-n2mult), 800*n2mult), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(900*(1-n2mult), 900*n2mult), means=means, sds = sds, trt_effect = trt_effect),
    list(n=c(1000*(1-n2mult), 1000*n2mult), means=means, sds = sds, trt_effect = trt_effect)
  )
  
  return(pars)
}



### SECTION 3.2: One-Arm Binary Outcome

#scenario 1: equal group probabilities (0.1 and 0.1) and unequal sample sizes (n2 = 0.1*N)
system.time(OABin_11_n201 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.1), n2mult = 0.1),
                               num_arms = 1, outcome_type = "binary"))

#scenario 2: unequal group probabilities (0.1 and 0.2) and unequal sample sizes (n2 = 0.1*N)
system.time(OABin_12_n201 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.2), n2mult = 0.1),
                               num_arms = 1, outcome_type = "binary"))

#scenario 3: unequal group probabilities (0.1 and 0.3) and unequal sample sizes (n2 = 0.1*N)
system.time(OABin_13_n201 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.3), n2mult = 0.1),
                               num_arms = 1, outcome_type = "binary"))

#scenario 4: unequal group probabilities (0.1 and 0.4) and unequal sample sizes (n2 = 0.1*N)
system.time(OABin_14_n201 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.4), n2mult = 0.1),
                               num_arms = 1, outcome_type = "binary"))

#scenario 5: unequal group probabilities (0.1 and 0.5) and unequal sample sizes (n2 = 0.1*N)
system.time(OABin_15_n201 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.5), n2mult = 0.1),
                               num_arms = 1, outcome_type = "binary"))

#scenario 6: equal group probabilities (0.1 and 0.1) and unequal sample sizes (n2 = 0.25*N)
system.time(OABin_11_n2025 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.1), n2mult = 0.25),
                                    num_arms = 1, outcome_type = "binary"))

#scenario 7: unequal group probabilities (0.1 and 0.2) and unequal sample sizes (n2 = 0.25*N)
system.time(OABin_12_n2025 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.2), n2mult = 0.25),
                                    num_arms = 1, outcome_type = "binary"))

#scenario 8: unequal group probabilities (0.1 and 0.3) and unequal sample sizes (n2 = 0.25*N)
system.time(OABin_13_n2025 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.3), n2mult = 0.25),
                                    num_arms = 1, outcome_type = "binary"))

#scenario 9: unequal group probabilities (0.1 and 0.4) and unequal sample sizes (n2 = 0.25*N)
system.time(OABin_14_n2025 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.4), n2mult = 0.25),
                                    num_arms = 1, outcome_type = "binary"))

#scenario 10: unequal group probabilities (0.1 and 0.5) and unequal sample sizes (n2 = 0.25*N)
system.time(OABin_15_n2025 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.5), n2mult = 0.25),
                                    num_arms = 1, outcome_type = "binary"))



### SECTION 3.3: One-Arm Continuous Outcome

#scenario 1: equal group means (1 and 1) and sd = sqrt(10) and unequal sample sizes (n2 = 0.1*N)
system.time(OACont_var10_11_n201 <- RunSim(pars = parfunc_unequaln(means = c(1,1), c(sqrt(10),sqrt(10)), n2mult = 0.1),
                                      num_arms = 1, outcome_type = "continuous"))

#scenario 2: unequal group means (1 and 2) and sd = sqrt(10) and unequal sample sizes (n2 = 0.1*N)
system.time(OACont_var10_12_n201 <- RunSim(pars = parfunc_unequaln(means = c(1,2), sds = c(sqrt(10),sqrt(10)), n2mult = 0.1),
                                      num_arms = 1, outcome_type = "continuous"))

#scenario 3: unequal group means (1 and 3) and sd = sqrt(10) and unequal sample sizes (n2 = 0.1*N)
system.time(OACont_var10_13_n201 <- RunSim(pars = parfunc_unequaln(means = c(1,3), sds = c(sqrt(10),sqrt(10)), n2mult = 0.1),
                                      num_arms = 1, outcome_type = "continuous"))

#scenario 4: unequal group means (1 and 4) and sd = sqrt(10) and unequal sample sizes (n2 = 0.1*N)
system.time(OACont_var10_14_n201 <- RunSim(pars = parfunc_unequaln(means = c(1,4), sds = c(sqrt(10),sqrt(10)), n2mult = 0.1),
                                      num_arms = 1, outcome_type = "continuous"))

#scenario 5: unequal group means (1 and 5) and sd = sqrt(10) and unequal sample sizes (n2 = 0.1*N)
system.time(OACont_var10_15_n201 <- RunSim(pars = parfunc_unequaln(means = c(1,5), sds = c(sqrt(10),sqrt(10)), n2mult = 0.1),
                                      num_arms = 1, outcome_type = "continuous"))

#scenario 6: equal group means (1 and 1) and sd = sqrt(10) and unequal sample sizes (n2 = 0.25*N)
system.time(OACont_var10_11_n2025 <- RunSim(pars = parfunc_unequaln(means = c(1,1), c(sqrt(10),sqrt(10)), n2mult = 0.25),
                                           num_arms = 1, outcome_type = "continuous"))

#scenario 7: unequal group means (1 and 2) and sd = sqrt(10) and unequal sample sizes (n2 = 0.25*N)
system.time(OACont_var10_12_n2025 <- RunSim(pars = parfunc_unequaln(means = c(1,2), sds = c(sqrt(10),sqrt(10)), n2mult = 0.25),
                                           num_arms = 1, outcome_type = "continuous"))

#scenario 8: unequal group means (1 and 3) and sd = sqrt(10) and unequal sample sizes (n2 = 0.25*N)
system.time(OACont_var10_13_n2025 <- RunSim(pars = parfunc_unequaln(means = c(1,3), sds = c(sqrt(10),sqrt(10)), n2mult = 0.25),
                                           num_arms = 1, outcome_type = "continuous"))

#scenario 9: unequal group means (1 and 4) and sd = sqrt(10) and unequal sample sizes (n2 = 0.25*N)
system.time(OACont_var10_14_n2025 <- RunSim(pars = parfunc_unequaln(means = c(1,4), sds = c(sqrt(10),sqrt(10)), n2mult = 0.25),
                                           num_arms = 1, outcome_type = "continuous"))

#scenario 10: unequal group means (1 and 5) and sd = sqrt(10) and unequal sample sizes (n2 = 0.25*N)
system.time(OACont_var10_15_n2025 <- RunSim(pars = parfunc_unequaln(means = c(1,5), sds = c(sqrt(10),sqrt(10)), n2mult = 0.25),
                                           num_arms = 1, outcome_type = "continuous"))



### SECTION 3.4: Two-Arm Binary Outcome

#scenario 1: equal group probabilities (0.1 and 0.1) and treatment effects (0.1 and 0.1) and unequal sample sizes (n2 = 0.1*N)
system.time(TABin_11_n201 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.1), trt_effect = c(0.1,0.1), n2mult = 0.1),
                               num_arms = 2, outcome_type = "binary"))

#scenario 2: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.2) and unequal sample sizes (n2 = 0.1*N)
system.time(TABin_12_n201 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.1), trt_effect = c(0.1,0.2), n2mult = 0.1),
                               num_arms = 2, outcome_type = "binary"))

#scenario 3: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.3) and unequal sample sizes (n2 = 0.1*N)
system.time(TABin_13_n201 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.1), trt_effect = c(0.1,0.3), n2mult = 0.1),
                               num_arms = 2, outcome_type = "binary"))

#scenario 4: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.4) and unequal sample sizes (n2 = 0.1*N)
system.time(TABin_14_n201 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.1), trt_effect = c(0.1,0.4), n2mult = 0.1),
                               num_arms = 2, outcome_type = "binary"))

#scenario 5: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.5) and unequal sample sizes (n2 = 0.1*N)
system.time(TABin_15_n201 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.1), trt_effect = c(0.1,0.5), n2mult = 0.1),
                               num_arms = 2, outcome_type = "binary"))

#scenario 6: equal group probabilities (0.1 and 0.1) and treatment effects (0.1 and 0.1) and unequal sample sizes (n2 = 0.25*N)
system.time(TABin_11_n2025 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.1), trt_effect = c(0.1,0.1), n2mult = 0.25),
                                    num_arms = 2, outcome_type = "binary"))

#scenario 7: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.2) and unequal sample sizes (n2 = 0.25*N)
system.time(TABin_12_n2025 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.1), trt_effect = c(0.1,0.2), n2mult = 0.25),
                                    num_arms = 2, outcome_type = "binary"))

#scenario 8: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.3) and unequal sample sizes (n2 = 0.25*N)
system.time(TABin_13_n2025 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.1), trt_effect = c(0.1,0.3), n2mult = 0.25),
                                    num_arms = 2, outcome_type = "binary"))

#scenario 9: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.4) and unequal sample sizes (n2 = 0.25*N)
system.time(TABin_14_n2025 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.1), trt_effect = c(0.1,0.4), n2mult = 0.25),
                                    num_arms = 2, outcome_type = "binary"))

#scenario 10: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.5) and unequal sample sizes (n2 = 0.25*N)
system.time(TABin_15_n2025 <- RunSim(pars = parfunc_unequaln(means = c(0.1,0.1), trt_effect = c(0.1,0.5), n2mult = 0.25),
                                    num_arms = 2, outcome_type = "binary"))


### SECTION 3.5: Twp-Arm Continuous Outcome

#scenario 1: equal group means (1 and 1) and treatment effects (1 and 1) and sd = sqrt(10) and unequal sample sizes (n2 = 0.1*N)
system.time(TACont_var10_11_n201 <- RunSim(pars = parfunc_unequaln(means = c(1,1), trt_effect = c(1,1), c(sqrt(10),sqrt(10)), n2mult = 0.1),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 2: unequal group means (1 and 1) and unequal treatment effects (1 and 2) and sd = sqrt(10) and unequal sample sizes (n2 = 0.1*N) 
system.time(TACont_var10_12_n201 <- RunSim(pars = parfunc_unequaln(means = c(1,1), trt_effect = c(1,2), sds = c(sqrt(10),sqrt(10)), n2mult = 0.1),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 3: unequal group means (1 and 1) and unequal treatment effects (1 and 3) and sd = sqrt(10) and unequal sample sizes (n2 = 0.1*N)
system.time(TACont_var10_13_n201 <- RunSim(pars = parfunc_unequaln(means = c(1,1), trt_effect = c(1,3), sds = c(sqrt(10),sqrt(10)), n2mult = 0.1),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 4: unequal group means (1 and 1) and unequal treatment effects (1 and 4) and sd = sqrt(10) and unequal sample sizes (n2 = 0.1*N)
system.time(TACont_var10_14_n201 <- RunSim(pars = parfunc_unequaln(means = c(1,1), trt_effect = c(1,4), sds = c(sqrt(10),sqrt(10)), n2mult = 0.1),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 5: unequal group means (1 and 1) and unequal treatment effects (1 and 5) and sd = sqrt(10) and unequal sample sizes (n2 = 0.1*N)
system.time(TACont_var10_15_n201 <- RunSim(pars = parfunc_unequaln(means = c(1,1), trt_effect = c(1,5), sds = c(sqrt(10),sqrt(10)), n2mult = 0.1),
                                      num_arms = 2, outcome_type = "continuous"))

#scenario 6: equal group means (1 and 1) and treatment effects (1 and 1) and sd = sqrt(10) and unequal sample sizes (n2 = 0.25*N)
system.time(TACont_var10_11_n2025 <- RunSim(pars = parfunc_unequaln(means = c(1,1), trt_effect = c(1,1), c(sqrt(10),sqrt(10)), n2mult = 0.25),
                                           num_arms = 2, outcome_type = "continuous"))

#scenario 7: unequal group means (1 and 1) and unequal treatment effects (1 and 2) and sd = sqrt(10) and unequal sample sizes (n2 = 0.25*N) 
system.time(TACont_var10_12_n2025 <- RunSim(pars = parfunc_unequaln(means = c(1,1), trt_effect = c(1,2), sds = c(sqrt(10),sqrt(10)), n2mult = 0.25),
                                           num_arms = 2, outcome_type = "continuous"))

#scenario 8: unequal group means (1 and 1) and unequal treatment effects (1 and 3) and sd = sqrt(10) and unequal sample sizes (n2 = 0.25*N)
system.time(TACont_var10_13_n2025 <- RunSim(pars = parfunc_unequaln(means = c(1,1), trt_effect = c(1,3), sds = c(sqrt(10),sqrt(10)), n2mult = 0.25),
                                           num_arms = 2, outcome_type = "continuous"))

#scenario 9: unequal group means (1 and 1) and unequal treatment effects (1 and 4) and sd = sqrt(10) and unequal sample sizes (n2 = 0.25*N)
system.time(TACont_var10_14_n2025 <- RunSim(pars = parfunc_unequaln(means = c(1,1), trt_effect = c(1,4), sds = c(sqrt(10),sqrt(10)), n2mult = 0.25),
                                           num_arms = 2, outcome_type = "continuous"))

#scenario 10: unequal group means (1 and 1) and unequal treatment effects (1 and 5) and sd = sqrt(10) and unequal sample sizes (n2 = 0.25*N)
system.time(TACont_var10_15_n2025 <- RunSim(pars = parfunc_unequaln(means = c(1,1), trt_effect = c(1,5), sds = c(sqrt(10),sqrt(10)), n2mult = 0.25),
                                           num_arms = 2, outcome_type = "continuous"))





##########################################################
########## SECTION 4:  AGGREGATE UNEQUAL N DATA ########## 
##########################################################

allresults_unequaln <- rbind(OABin_11_n201, OABin_12_n201, OABin_13_n201, OABin_14_n201, OABin_15_n201,
                    OABin_11_n2025, OABin_12_n2025, OABin_13_n2025, OABin_14_n2025, OABin_15_n2025,
                    OACont_var10_11_n201, OACont_var10_12_n201, OACont_var10_13_n201, OACont_var10_14_n201, OACont_var10_15_n201,
                    OACont_var10_11_n2025, OACont_var10_12_n2025, OACont_var10_13_n2025, OACont_var10_14_n2025, OACont_var10_15_n2025,
                    TABin_11_n201, TABin_12_n201, TABin_13_n201, TABin_14_n201, TABin_15_n201,
                    TABin_11_n2025, TABin_12_n2025, TABin_13_n2025, TABin_14_n2025, TABin_15_n2025,
                    TACont_var10_11_n201, TACont_var10_12_n201, TACont_var10_13_n201, TACont_var10_14_n201, TACont_var10_15_n201,
                    TACont_var10_11_n2025, TACont_var10_12_n2025, TACont_var10_13_n2025, TACont_var10_14_n2025, TACont_var10_15_n2025)


write.csv(allresults_unequaln, "simulation_results_unequaln_raw.csv")
#allresults_unequaln <- read.csv("simulation_results_unequaln_raw.csv")

#create summary data set using the following cutoffs
#PRIMARY CUTOFF = 0.2
MEM_cutoff <- 0.2
p_cutoff <- 0.05 
nsim <- 10000
allresults_unequaln_summary <- allresults_unequaln %>% 
  group_by(num_arms, outcome_type, n1, n2, mean1, mean2, sd1, sd2, trteff1, trteff2) %>%
  summarize(MEM_power = sum(MEMpexch < MEM_cutoff, na.rm = TRUE)/nsim,
            MEMr_power = sum(MEMrpexch < MEM_cutoff)/nsim,
            pval_power = sum(pval < p_cutoff)/nsim,
            MEM_contapprox_power = sum(MEMapproxpexch < MEM_cutoff)/nsim,
            pval_chisq_power = sum(pval_chisq < p_cutoff)/nsim)

allresults_unequaln_summary$MEMcutoff <- MEM_cutoff
allresults_unequaln_summary$pvalcutoff <- p_cutoff

write.csv(allresults_unequaln_summary, "simulation_results_unequaln_summary.csv")

#ALTERNATE CUTOFF = 0.8
MEM_cutoff <- 0.8
p_cutoff <- 0.05 
nsim <- 10000
allresults_unequaln_summary_08 <- allresults_unequaln %>% 
  group_by(num_arms, outcome_type, n1, n2, mean1, mean2, sd1, sd2, trteff1, trteff2) %>%
  summarize(MEM_power = sum(MEMpexch < MEM_cutoff, na.rm = TRUE)/nsim,
            MEMr_power = sum(MEMrpexch < MEM_cutoff)/nsim,
            pval_power = sum(pval < p_cutoff)/nsim,
            MEM_contapprox_power = sum(MEMapproxpexch < MEM_cutoff)/nsim,
            pval_chisq_power = sum(pval_chisq < p_cutoff)/nsim)

allresults_unequaln_summary_08$MEMcutoff <- MEM_cutoff
allresults_unequaln_summary_08$pvalcutoff <- p_cutoff

write.csv(allresults_unequaln_summary_08, "simulation_results_unequaln_summary_08cutoff.csv")

#T1E RATE CALIBRATED CUTOFF
MEM_cutoff <- "calibrated"
p_cutoff <- 0.05 
nsim <- 10000

allresults_unequaln$effsize <- ifelse(allresults_unequaln$num_arms == 1, 
                                      allresults_unequaln$mean2 - allresults_unequaln$mean1, 
                                      allresults_unequaln$trteff2 - allresults_unequaln$trteff1)
cutoffs <- allresults_unequaln[allresults_unequaln$effsize == 0,] %>% group_by(n1, n2, outcome_type, num_arms, sd1, sd2) %>%
  summarize(q5_MEM = quantile(MEMpexch, 0.05, na.rm = TRUE), 
            q5_MEMcontapprox = quantile(MEMapproxpexch, 0.05, na.rm = TRUE), 
            q5_MEMr = quantile(MEMrpexch, 0.05)) %>% 
  as.data.frame()
allresults_unequaln_cutoffs <- merge(allresults_unequaln, cutoffs)

allresults_unequaln_summary_calibrated <- allresults_unequaln_cutoffs %>% 
  group_by(num_arms, outcome_type, n1, n2, mean1, mean2, sd1, sd2, trteff1, trteff2) %>%
  summarize(MEM_power = sum(MEMpexch < q5_MEM, na.rm = TRUE)/nsim,
            MEMr_power = sum(MEMrpexch < q5_MEMr)/nsim,
            pval_power = sum(pval < p_cutoff)/nsim,
            MEM_contapprox_power = sum(MEMapproxpexch < q5_MEMcontapprox)/nsim,
            pval_chisq_power = sum(pval_chisq < p_cutoff)/nsim)

allresults_unequaln_summary_calibrated$MEMcutoff <- MEM_cutoff
allresults_unequaln_summary_calibrated$pvalcutoff <- p_cutoff

write.csv(allresults_unequaln_summary_calibrated, "simulation_results_unequaln_summary_calibratedcutoff.csv")


