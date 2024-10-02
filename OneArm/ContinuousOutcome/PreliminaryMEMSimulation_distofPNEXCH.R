################################################################################
# TITLE: PreliminaryMEMSimulation_distofPNEXCH.R
#
# PURPOSE: Preliminary MEM Simulation to Determine Cut-Offs for Exchangeability 
#       Probabilities Comparable to Two-Sample T-Test Power Values (look at distribution of P(Not Exchangeable))
#
# AUTHOR: Shannon Murphy
# DATE CREATED: NOV 11, 2023
################################################################################

# This simulation will compare the exchangeability probabilities given by MEMs to 
# the power of a two-sample t.test. 

library(ggplot2)
library(MESS) #power_t_test


#load in MEM function from Joe's GitHub
source("C:/Users/mushanno/Desktop/Work/Dissertation1/CodeFromAlex/Kaizer_MEM_Weights.R")

#test MEM function
calc.weights_MEM(xvec = c(1,1),svec = c(10,10),nvec = c(200,200),prior = 'pi_e')

## Define function to simulate data for all parameter combinations
sim_samp <- function(params, nsamp, MEMprior, power){

  tt_pval <- data.frame(matrix(nrow = dim(params), ncol = nsamp))
  MEM_pnexch <- data.frame(matrix(nrow = dim(params), ncol = nsamp))
  tt_pval_null <- data.frame(matrix(nrow = dim(params), ncol = nsamp))
  MEM_pnexch_null <- data.frame(matrix(nrow = dim(params), ncol = nsamp))
  
  for (i in 1:dim(params)[[1]]){
    
    for (j in 1:nsamp){
      dat1 <- rnorm((1-params[i,2])*params[i,1], mean = params[i,4], sd = sqrt(params[i,3])) #new data
      dat2 <- rnorm(params[i,2]*params[i,1], mean = params[i,5], sd = sqrt(params[i,3]))     #historic data
      datnull <- rnorm(params[i,2]*params[i,1], mean = params[i,4], sd = sqrt(params[i,3]))  #null case of historic data
      
      tt_pval[i,j] <- t.test(dat1, dat2, var.equal = TRUE)$p.value
      tt_pval_null[i,j] <- t.test(dat1, dat2, var.equal = TRUE)$p.value
      
      
      #store probability that sources are NOT exchangeable
      MEM_pnexch[i,j] <- calc.weights_MEM(xvec = c(mean(dat1), mean(dat2)), 
                                        nvec = c((1-params[i,2])*params[i,1], params[i,2]*params[i,1]),
                                        svec = c(sd(dat1), sd(dat2)), 
                                        prior = MEMprior)[1]
      MEM_pnexch_null[i,j] <- calc.weights_MEM(xvec = c(mean(dat1), mean(datnull)), 
                                        nvec = c((1-params[i,2])*params[i,1], params[i,2]*params[i,1]),
                                        svec = c(sd(dat1), sd(datnull)), 
                                        prior = MEMprior)[1]
    }
    
    params$power[i] <- sum(tt_pval[i,] < 0.05)/nsamp          #power for t-test
    params$tttypeIerror <- sum(tt_pval_null[i,] < 0.05)/nsamp #type I error for t-test
    
    params$pexch_cutoff[i] <- round(quantile((MEM_pnexch[i,]), 1-power, na.rm = TRUE),5) #get cut off value for given power
    
    params$pexch_typeIerror[i] <- round(sum(MEM_pnexch_null[i,] >= params$pexch_cutoff[i])/nsamp,5) #type I error for MEM
  }
  
  params[,(dim(params)[[2]]+1):(dim(params)[[2]]+nsamp+1)] <- MEM_pnexch
  
  return(params)
  
}



# First, the subgroups will be equal sizes and have equal variance (n1 = n2, 
# v1 = v2), priors will be equal and fixed, and the variances and difference 
# in effect size will vary.

##Set Simulation Parameters for equal sample sizes
N = seq(200,1000,200) ## total sample sizes
n2mult = 0.5 ##percentage of total sample size that group 2 will be
v1 = v2 = c(1,10,100) ## variance
p_exch <- 'pi_e' ## prior for MEM
effsize1 <- 1 ## Effect Size for group 1
effsize2 <- seq(1,5,0.5) ## Effect Sizes for group 2
paramgrid <- expand.grid(N=N, n2=n2mult, v1=v1, es1=effsize1, es2=effsize2, power = NA, tttypeIerror = NA, pexch_cutoff = NA, pexch_typeIerror = NA)

ns <- 10000 ## Set number of samples per parameter combination

## run simulation
set.seed(1000)
pnexch_distdata <- sim_samp(paramgrid, ns, p_exch, 0.8)
write.csv(pnexch_distdata, "C:\\Users\\mushanno\\Desktop\\Work\\Dissertation1\\OneArm\\ContinuousOutcome\\resultsdata\\pnexch_dist.csv", row.names = FALSE)


##Set Simulation Parameters for unequal sample size
n2mult = c(0.1, 0.2, 0.5) ##percentage of total sample size that group 2 will be
v1 = v2 = 100 ## constant variance for easier visualization
paramgrid2 <- expand.grid(N=N, n2=n2mult, v1=v1, es1=effsize1, es2=effsize2, power = NA, pexch = NA)

## run simulation
set.seed(1001)
pnexch_distdata_unequaln <- sim_samp(paramgrid2,ns,p_exch, power = 0.8)
write.csv(pnexch_distdata_unequaln, "C:\\Users\\mushanno\\Desktop\\Work\\Dissertation1\\OneArm\\ContinuousOutcome\\resultsdata\\pnexch_dist_unequaln.csv", row.names = FALSE)

