#################################################################################
# TITLE: Run_Additional_Simulations.R
#
# PURPOSE: Run simulations in parallel using RunSim() func from 
#          Simulation_Function.R and create a master data output file
#          for additional binary scenarios with higher baseline prevalence
#          and thus higher variance. 
#
# OUTPUT: higherbaselineprobs_simulation_results_raw.csv
#         higherbaselineprobs_simulation_results_summary.csv
#         higherbaselineprobs_simulation_results_summary_08cutoff.csv
#         higherbaselineprobs_simulation_results_summary_calibratedcutoff.csv
#         
#
# SECTIONS: Section 1
#              run simulations with various parameter settings 
#              1.1 parfunc() - function to create parameter list of lists
#              1.2 One-Arm Binary Outcome
#              1.3 Two-Arm Binary Outcome
#           Section 2
#              aggregate data
#
# NOTES: Sections 1-2 took about 1 hour to run on my DELL laptop - 
#           Processor - Intel(R) Core(TM) i7-10750H CPU @ 2.60GHz, 2592 Mhz, 
#                       6 Core(s), 12 Logical Processor(s)
#           Installed RAM - 16.0 GB (15.8 GB usable)
#
# AUTHOR: Shannon Thomas
# DATE CREATED: DEC 17, 2024
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
system.time(OABin_55 <- RunSim(pars = parfunc(means = c(0.5,0.5)),
                               num_arms = 1, outcome_type = "binary"))

#scenario 2: unequal group probabilities (0.1 and 0.2)
system.time(OABin_56 <- RunSim(pars = parfunc(means = c(0.5,0.6)),
                               num_arms = 1, outcome_type = "binary"))

#scenario 3: unequal group probabilities (0.1 and 0.3)
system.time(OABin_57 <- RunSim(pars = parfunc(means = c(0.5,0.7)),
                               num_arms = 1, outcome_type = "binary"))

#scenario 4: unequal group probabilities (0.1 and 0.4)
system.time(OABin_58 <- RunSim(pars = parfunc(means = c(0.5,0.8)),
                               num_arms = 1, outcome_type = "binary"))

#scenario 5: unequal group probabilities (0.1 and 0.5)
system.time(OABin_59 <- RunSim(pars = parfunc(means = c(0.5,0.9)),
                               num_arms = 1, outcome_type = "binary"))


### SECTION 1.3: Two-Arm Binary Outcome

#scenario 1: equal group probabilities (0.1 and 0.1) and treatment effects (0.1 and 0.1)
system.time(TABin_11 <- RunSim(pars = parfunc(means = c(0.5,0.5), trt_effect = c(0.1,0.1)),
                               num_arms = 2, outcome_type = "binary"))

#scenario 2: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.2)
system.time(TABin_12 <- RunSim(pars = parfunc(means = c(0.5,0.5), trt_effect = c(0.1,0.2)),
                               num_arms = 2, outcome_type = "binary"))

#scenario 3: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.3)
system.time(TABin_13 <- RunSim(pars = parfunc(means = c(0.5,0.5), trt_effect = c(0.1,0.3)),
                               num_arms = 2, outcome_type = "binary"))

#scenario 4: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.4)
system.time(TABin_14 <- RunSim(pars = parfunc(means = c(0.5,0.5), trt_effect = c(0.1,0.4)),
                               num_arms = 2, outcome_type = "binary"))

#scenario 5: equal group probabilities (0.1 and 0.1) and unequal treatment effects (0.1 and 0.5)
system.time(TABin_15 <- RunSim(pars = parfunc(means = c(0.5,0.5), trt_effect = c(0.1,0.5)),
                               num_arms = 2, outcome_type = "binary"))




##########################################################
############### SECTION 2:  AGGREGATE DATA ############### 
##########################################################

allresults <- rbind(OABin_55, OABin_56, OABin_57, OABin_58, OABin_59,
                    TABin_11, TABin_12, TABin_13, TABin_14, TABin_15)


write.csv(allresults, "higherbaselineprobs_simulation_results_raw.csv")

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

write.csv(allresults_summary, "higherbaselineprobs_simulation_results_summary.csv")

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

write.csv(allresults_summary_08, "higherbaselineprobs_simulation_results_summary_08cutoff.csv")

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

write.csv(allresults_summary_calibrated, "higherbaselineprobs_simulation_results_summary_calibratedcutoff.csv")






