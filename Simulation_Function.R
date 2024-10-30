#################################################################################
# TITLE: Simulation_Function.R
#
# PURPOSE: Function to run simulations in parallel using the snowfall package.
#
# OUTPUT: data frame with simulation results
#
# SECTIONS: Section 1
#              RunSim() - function to run the simulation in parallel 
#                         given a number of simulations per setting and 
#                         a list with parameter settings containing 
#                         n, means, sds (if continuous),& trt_effect (if two-arm),
#                         num_arms, outcome_type, and marginal
#           
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
