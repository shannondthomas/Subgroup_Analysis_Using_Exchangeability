#################################################################################
# TITLE: Run_Simulations.R
#
# PURPOSE: Run simulations in parallel using the snowfall package and create
#          a master data output file. 
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
#         simulation_results_unequaln_raw.csv - same as above for unequal n
#
#         simulation_results_unequaln_summary.csv - same as above for unequal n
#
#         simulation_results_unequaln_summary_08cutoff.csv  - same as above for 
#                                                             unequal n
#
# SECTIONS: Section 1
#              RunSim() - function to run the simulation in parallel 
#                         given a number of simulations per setting and 
#                         a list with parameter settings containing 
#                         n, means, sds (if continuous),& trt_effect (if two-arm),
#                         num_arms, outcome_type, and marginal
#           Section 2
#              run simulations with various parameter settings 
#              2.1 parfunc() - function to create parameter list of lists
#              2.2 One-Arm Binary Outcome
#              2.3 One-Arm Continuous Outcome
#              2.4 Two-Arm Binary Outcome
#              2.5 Two-Arm Continuous Outcome
#           Section 3
#              aggregate data
#
#           Section 4
#              repeat section 2 for unequal sample sizes in groups
#           Section 5
#              repeat section 3 for unequal sample sizes in groups
#
# NOTES: Sections 1-5 took 16 hours to run on my DELL laptop - 
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

source("Test_Functions.R")
#R version 4.3.2 was used for this project.
library(tidyverse) #version 2.0.0
library(snowfall)  #version 1.84-6.3




##########################################################
############# SECTION 1: SIMULATION FUNCTION ############# 
##########################################################

RunSim <- function(nsim = 10000, pars, num_arms, outcome_type, marginal = "BIC") {
  
  #summary matrix: first weight 1 (non-exchangeable) and second is weight 2 (exchangeable)
  summ <- data.frame(n1 = rep(NA, length(pars)*nsim),
                     n2 = rep(NA, length(pars)*nsim),
                     mean1 = rep(NA, length(pars)*nsim),
                     mean2 = rep(NA, length(pars)*nsim),                     
                     sd1 = rep(NA, length(pars)*nsim),
                     sd2 = rep(NA, length(pars)*nsim),
                     trteff1 = rep(NA, length(pars)*nsim),
                     trteff2 = rep(NA, length(pars)*nsim),
                     MEMpexch = rep(NA, length(pars)*nsim), 
                     MEMpnexch = rep(NA, length(pars)*nsim),
                     MEMrpexch = rep(NA, length(pars)*nsim), 
                     MEMrpnexch = rep(NA, length(pars)*nsim), 
                     pval = rep(NA, length(pars)*nsim))
  
  
  pb= txtProgressBar(min = 0, max = length(pars), style = 3, char=":)")
  for(i in 1:length(pars)) {
    ## Define the parameters
    par <- pars[[i]] 
    
    #### Do the simulation with these settings
    sfInit(parallel=T, cpus=12, type='SOCK')
    sfLibrary(mvtnorm)
    sfLibrary(gsDesign)
    sfLibrary(car)
    sfLibrary(nlme)
    sfLibrary(HDInterval)
    sfLibrary(lmtest)
    sfLibrary(matrixStats)
    sfClusterSetupRNG(seed=12345)
    
    sfExportAll()
    sim.results <- t(sfLapply(rep(1,nsim),
                              function(x) RunModels(n=par$n, means=par$means, sds=par$sds,
                                                    trt_effect = par$trt_effect,
                                                    marginal = marginal, 
                                                    outcome_type = outcome_type,
                                                    num_arms = num_arms)))
    sfStop()
    
    result <- as.data.frame(unlist(sim.results))
    
    summ[(nsim*(i-1) + 1):(nsim*i),1] <- (result %>% filter(row_number() %% 13 == 1))
    summ[(nsim*(i-1) + 1):(nsim*i),2] <- (result %>% filter(row_number() %% 13 == 2))
    summ[(nsim*(i-1) + 1):(nsim*i),3] <- (result %>% filter(row_number() %% 13 == 3))
    summ[(nsim*(i-1) + 1):(nsim*i),4] <- (result %>% filter(row_number() %% 13 == 4))
    summ[(nsim*(i-1) + 1):(nsim*i),5] <- (result %>% filter(row_number() %% 13 == 5))
    summ[(nsim*(i-1) + 1):(nsim*i),6] <- (result %>% filter(row_number() %% 13 == 6))    
    summ[(nsim*(i-1) + 1):(nsim*i),7] <- (result %>% filter(row_number() %% 13 == 7))
    summ[(nsim*(i-1) + 1):(nsim*i),8] <- (result %>% filter(row_number() %% 13 == 8))
    summ[(nsim*(i-1) + 1):(nsim*i),9] <- (result %>% filter(row_number() %% 13 == 9))
    summ[(nsim*(i-1) + 1):(nsim*i),10] <- (result %>% filter(row_number() %% 13 == 10))
    summ[(nsim*(i-1) + 1):(nsim*i),11] <- (result %>% filter(row_number() %% 13 == 11))    
    summ[(nsim*(i-1) + 1):(nsim*i),12] <- (result %>% filter(row_number() %% 13 == 12))
    summ[(nsim*(i-1) + 1):(nsim*i),13] <- (result %>% filter(row_number() %% 13 == 0))
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  summ$outcome_type <- outcome_type
  summ$num_arms <- num_arms
  
  return(summ)
  
}



##########################################################
############### SECTION 2:  RUN SIMULATION ############### 
##########################################################


### SECTION 2.1: Function to create parameter list of lists 

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


### SECTION 2.2: One-Arm Binary Outcome

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


### SECTION 2.3: One-Arm Continuous Outcome

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


### SECTION 2.4: Two-Arm Binary Outcome

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


### SECTION 2.5: Two-Arm Continuous Outcome

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
############### SECTION 3:  AGGREGATE DATA ############### 
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
                                pval_power = sum(pval < p_cutoff)/nsim)

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
            pval_power = sum(pval < p_cutoff)/nsim)

allresults_summary_08$MEMcutoff <- MEM_cutoff
allresults_summary_08$pvalcutoff <- p_cutoff

write.csv(allresults_summary_08, "simulation_results_summary_08cutoff.csv")




##########################################################
########## SECTION 4: RUN SIMULATION (n1 != n2) ########## 
##########################################################


### SECTION 2.1: Function to create parameter list of lists 

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



### SECTION 2.2: One-Arm Binary Outcome

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



### SECTION 2.2: One-Arm Continuous Outcome

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



### SECTION 2.2: Two-Arm Binary Outcome

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


### SECTION 2.2: Twp-Arm Continuous Outcome

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
########## SECTION 5:  AGGREGATE UNEQUAL N DATA ########## 
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
  summarize(MEM_power = sum(MEMpexch < MEM_cutoff)/nsim,
            MEMr_power = sum(MEMrpexch < MEM_cutoff)/nsim,
            pval_power = sum(pval < p_cutoff)/nsim)

allresults_unequaln_summary$MEMcutoff <- MEM_cutoff
allresults_unequaln_summary$pvalcutoff <- p_cutoff

write.csv(allresults_unequaln_summary, "simulation_results_unequaln_summary.csv")

#ALTERNATE CUTOFF = 0.8
MEM_cutoff <- 0.8
p_cutoff <- 0.05 
nsim <- 10000
allresults_unequaln_summary_08 <- allresults_unequaln %>% 
  group_by(num_arms, outcome_type, n1, n2, mean1, mean2, sd1, sd2, trteff1, trteff2) %>%
  summarize(MEM_power = sum(MEMpexch < MEM_cutoff)/nsim,
            MEMr_power = sum(MEMrpexch < MEM_cutoff)/nsim,
            pval_power = sum(pval < p_cutoff)/nsim)

allresults_unequaln_summary_08$MEMcutoff <- MEM_cutoff
allresults_unequaln_summary_08$pvalcutoff <- p_cutoff

write.csv(allresults_unequaln_summary_08, "simulation_results_unequaln_summary_08cutoff.csv")




