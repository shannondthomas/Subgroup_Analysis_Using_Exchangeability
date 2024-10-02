################################################################################
# TITLE: PowerSimulation_OA_Binary_MEM.R
#
# PURPOSE: Power simulation for LRT and MEM for one-arm binary outcome
#          with multiple exchangeability cutoffs possible. 
#
#
# AUTHOR: Shannon Murphy
# DATE CREATED: FEB 29, 2024
################################################################################


source("C:/Users/mushanno/Desktop/Work/Dissertation1/CodeFromAlex/GitHub_Functions_AdaptivePlatformDesign.R")


# calc.MEM.betabin <- function(xvec, nvec, avec, bvec, prior, constraint=1){
  ###Function to calculate the MEM model weights for binomial data with beta(alpha,beta) prior
  #xvec: vector with counts of those with event (primary source first, supplemental afterwards)
  #nvec: vector of total number of subjects
  #avec: vector of alpha parameters for beta priors for each source
  #bvec: vector of beta parameters for beta priors for each source
  #prior: prior to use on MEM source inclusion probability: equal (pi_e), pool (naive pooling), opt.source (pi_EB), opt.source_constrain (pi_EBc)
  #constraint: value to limit maximum value of empirical Bayes prior to for pi_EBc
  
# test function
# calc.MEM.betabin(xvec = c(8,1), nvec = c(10,10), avec = c(1,1), bvec = c(1,1), prior = 'equal',
#                  constraint = 1)$q
#$q[1] is prob not exchangeable, $q[2] is prob exchangeable


################################################################################
# TITLE: PreliminaryMEMSimulation3.R
#
# PURPOSE: Preliminary MEM Simulation to Determine Cut-Offs for Exchangeability 
#       Probabilities Comparable to Two-Sample T-Test Power Values
#
# AUTHOR: Shannon Murphy
# DATE CREATED: NOV 27, 2023
################################################################################

# This simulation will compare the exchangeability probabilities given by MEMs to 
# the power of a two-sample t.test. 

library(ggplot2)

## Define function to simulate data for all parameter combinations
sim_samp_bin <- function(params, nsamp, alpha, beta, MEMprior, MEMcutoffs){

  pb <- txtProgressBar(min = 0, max = dim(params)[[1]],style=3)
  
  for (i in 1:dim(params)[[1]]){
    
    pt_pval <- rep(NA,nsamp)
    MEM_pnexch <- rep(NA, nsamp)
    
    for (j in 1:nsamp){
      ss <- c((1-params[i,2])*params[i,1], 
              params[i,2]*params[i,1])
      
      dat1 <- rbinom(n=1, size = ss[1], prob = params[i,4])
      dat2 <- rbinom(n=1, size = ss[2], prob = params[i,5])
      
      pt_pval[j] <- prop.test(c(dat1, dat2), ss)$p.value
      #store probability that sources are NOT exchangeable
      MEM_pnexch[j] <- calc.MEM.betabin(xvec = c(dat1, dat2), 
                                        nvec = c(ss[1], ss[2]),
                                        avec = alpha, bvec = beta, prior = MEMprior,
                                        constraint = 1)$q[1]
      
    }
    
    params$power[i] <- sum(pt_pval < 0.05)/nsamp
    params$pnexch[i] <- round(mean(MEM_pnexch),5)
    
    for (k in 1:length(MEMcutoffs)){
      params[i,8+k] <- round(sum(MEM_pnexch >= MEMcutoffs[k])/nsamp,5)
      colnames(params)[8+k] <- paste('pnexch_cutoff',substr(as.character(MEMcutoffs[k]),start = 3,stop = 5),'_power',sep='')
    }
    
    setTxtProgressBar(pb, i)
    
  }
  
  close(pb)
  
  return(params)
  
}



# First, the subgroups will be equal sizes and have equal variance (n1 = n2, 
# v1 = v2), priors will be equal and fixed, and the difference 
# in effect size will vary.

##Set Simulation Parameters
N = seq(200,1000,200) ## total sample sizes
n2mult = 0.5 ##percentage of total sample size that group 2 will be
p_exch <- 'equal' ## prior for MEM
effsize1 <- 0.1 ## Effect Size for group 1
effsize2 <- seq(0.1,0.9,0.1) ## Effect Sizes for group 2
paramgrid <- expand.grid(N=N, n2=n2mult, v1=NA, es1=effsize1, es2=effsize2, power = NA, pnexch = NA, pnexch_cutoff = NA)

ns <- 10000 ## Set number of samples per parameter combination

## run simulation
set.seed(1000)
cutoffs <- c(0.2,0.7976455,0.8,0.9) 
OA_binary_MEM_result <- sim_samp_bin(paramgrid, ns, c(1,1), c(1,1), p_exch, cutoffs )

write.csv(OA_binary_MEM_result, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/OA_binary_MEM_result.csv")

# ## Plot results
# ## simulated sample plots
# ggplot(OA_binary_MEM_result, aes(es2-es1, N, fill=power)) +
#   geom_tile() + geom_text(aes(label=power)) +
#   labs(title = "Heatmap of prop.test Power for Varying Sample Sizes and Differences in Effect Sizes",
#        fill = 'Type I Error/Power') +
#   xlab("Difference in Effect Sizes") +
#   ylab("Total Sample Size (N)")
# 
# for(i in 1:length(cutoffs)){
#   fillvar <- colnames(OA_binary_MEM_result)[8+i]
#   
#   print(
#     ggplot(OA_binary_MEM_result, aes(es2-es1, N, fill=.data[[fillvar]])) + 
#       geom_tile() + geom_text(aes(label=.data[[fillvar]])) + 
#       labs(title = paste("Heatmap of P(Not Exchangeable) >=", cutoffs[i], 
#                          "for Varying Sample Sizes and Differences in Effect Sizes",
#                          sep = ' '),
#            fill = "Type I Error/Power") +
#       xlab("Difference in Effect Sizes") +
#       ylab("Total Sample Size (N)")
#   )
# }




####################################################################################
# Now, the subgroups will be unequal sizes (n1 != n2), variance and priors will be #
# equal and fixed, and the difference in effect size will vary.                    #
####################################################################################

##Set Simulation Parameters
n2mult = c(0.1, 0.25) ##percentage of total sample size that group 2 will be
v1 = v2 = 10 ## constant variance for easier visualization
paramgrid <- expand.grid(N=N, n2=n2mult, v1=NA, es1=effsize1, es2=effsize2, power = NA, pnexch = NA, pnexch_cutoff = NA)

## run simulation
set.seed(1001)
OA_binary_MEM_diffn_result <- sim_samp_bin(paramgrid, ns, c(1,1), c(1,1), p_exch, cutoffs )

#create n2mult character variable for plotting
OA_binary_MEM_diffn_result$n2mult_c <- NA
OA_binary_MEM_diffn_result$n2mult_c[OA_binary_MEM_diffn_result$n2 == 0.1] <- 'n2 = 0.1*N'
OA_binary_MEM_diffn_result$n2mult_c[OA_binary_MEM_diffn_result$n2 == 0.25] <- 'n2 = 0.25*N'

write.csv(OA_binary_MEM_diffn_result, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/OA_binary_MEM_diffn_result.csv")

# ## Plot results
# ## simulated sample results
# ggplot(OA_binary_MEM_diffn_result, aes(es2-es1, N, fill=power)) +
#   geom_tile() + geom_text(aes(label=power)) + facet_grid(n2mult_c~.) +
#   labs(title = "Heatmap of prop.test Power for Varying Sample Sizes per Group and Differences in Effect Sizes",
#        fill = 'Type I Error/Power') +
#   xlab("Difference in Effect Sizes") +
#   ylab("Total Sample Size (N)")
# 
# for(i in 1:length(cutoffs)){
#   fillvar <- colnames(OA_binary_MEM_diffn_result)[8+i]
#   
#   print(
#     ggplot(OA_binary_MEM_diffn_result, aes(es2-es1, N, fill=.data[[fillvar]])) + 
#       geom_tile() + geom_text(aes(label=.data[[fillvar]])) + facet_grid(n2mult_c~.) +
#       labs(title = paste("Heatmap of P(Not Exchangeable) >=", cutoffs[i], 
#                          "for Varying Sample Sizes per Group and Differences in Effect Sizes",
#                          sep = ' '),
#            fill = "Type I Error/Power") +
#       xlab("Difference in Effect Sizes") +
#       ylab("Total Sample Size (N)")
#   )
# }
