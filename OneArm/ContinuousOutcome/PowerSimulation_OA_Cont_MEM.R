################################################################################
# TITLE: PowerSimulation_OA_Cont_MEM.R
#
# PURPOSE: Power simulation for T-test and MEM for one-arm continuous outcome
#          with multiple exchangeability cutoffs possible. 
#
# NOTE: The simulation also outputs information on calculated cutoffs for given
#       type I error and power levels as an artifact of previous ideation.
#
# AUTHOR: Shannon Murphy
# DATE CREATED: FEB 29, 2024
################################################################################

 

library(ggplot2)
library(MESS) #power_t_test


#MEM function downloaded from Joe Koopmeiner's GitHub
source("C:/Users/mushanno/Desktop/Work/Dissertation1/CodeFromAlex/Kaizer_MEM_Weights.R")

# #test MEM function
# calc.weights_MEM(xvec = c(1,1),svec = c(10,10),nvec = c(200,200),prior = 'pi_e')


####################################################################
## Define function to simulate data for all parameter combinations #
####################################################################

sim_samp <- function(params, nsamp, MEMprior, MEMcutoffs, typeierror, typeiierror){
  
  pb= txtProgressBar(min = 0, max = dim(params)[[1]], style = 3, char=":)")
  
  for (i in 1:dim(params)[[1]]){
    
    tt_pval <- rep(NA,nsamp)
    MEM_pnexch <- rep(NA, nsamp)
    
    for (j in 1:nsamp){
      dat1 <- rnorm((1-params[i,2])*params[i,1], mean = params[i,4], sd = sqrt(params[i,3]))
      dat2 <- rnorm(params[i,2]*params[i,1], mean = params[i,5], sd = sqrt(params[i,3]))
      
      tt_pval[j] <- t.test(dat1, dat2, var.equal = TRUE)$p.value
      #store probability that sources are NOT exchangeable
      MEM_pnexch[j] <- calc.weights_MEM(xvec = c(mean(dat1), mean(dat2)), 
                                               nvec = c((1-params[i,2])*params[i,1], params[i,2]*params[i,1]),
                                               svec = c(sd(dat1), sd(dat2)), 
                                               prior = MEMprior)[1]

    }
    
    params$power[i] <- sum(tt_pval < 0.05)/nsamp
    params$pnexch[i] <- round(mean(MEM_pnexch),5)
    if (params[i,4] == params[i,5]){
      params$pnexch_power[i] <- round(quantile(MEM_pnexch, 1-typeierror),5) #get cut off value for given Type I error (for null)
    }
    else {
      params$pnexch_power[i] <- round(quantile(MEM_pnexch, typeiierror),5) #get cut off value for given Type II error (for alternative)
      
    }
    
    for (k in 1:length(MEMcutoffs)){
      params[i,8+k] <- round(sum(MEM_pnexch >= MEMcutoffs[k])/nsamp,5)
      colnames(params)[8+k] <- paste('pnexch_cutoff',substr(as.character(MEMcutoffs[k]),start = 3,stop = 5),'_power',sep='')
    }

    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  return(params)
  
}




##############################################################################
# First, the subgroups will be equal sizes and have equal variance (n1 = n2, # 
# v1 = v2), priors will be equal and fixed, and the variances and difference # 
# in effect size will vary.                                                  #
##############################################################################

##Set Simulation Parameters
N = seq(200,1000,200) ## total sample sizes
n2mult = 0.5 ##percentage of total sample size that group 2 will be
v1 = v2 = c(1,10,100) ## variance
p_exch <- 'pi_e' ## prior for MEM
effsize1 <- 1 ## Effect Size for group 1
effsize2 <- seq(1,5,0.5) ## Effect Sizes for group 2
paramgrid <- expand.grid(N=N, n2=n2mult, v1=v1, es1=effsize1, es2=effsize2, power = NA, pnexch = NA, pnexch_power = NA)

ns <- 10000 ## Set number of samples per parameter combination

## run simulation
set.seed(1000)
cutoffs <- c(0.2,0.7976455,0.8,0.9) 
OA_cont_MEM_result <- sim_samp(paramgrid, ns, p_exch, cutoffs, 0.05, 0.2)

#create var character variable for plotting
OA_cont_MEM_result$var_c <- NA
OA_cont_MEM_result$var_c[OA_cont_MEM_result$v1 == 1] <- 'Var = 1'
OA_cont_MEM_result$var_c[OA_cont_MEM_result$v1 == 10] <- 'Var = 10'
OA_cont_MEM_result$var_c[OA_cont_MEM_result$v1 == 100] <- 'Var = 100'

write.csv(OA_cont_MEM_result, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/OA_cont_MEM_result.csv")
    

######expand var = 1 case
  N = seq(200,1000,200) ## total sample sizes
  n2mult = 0.5 ##percentage of total sample size that group 2 will be
  v1 = v2 = c(1) ## variance
  p_exch <- 'pi_e' ## prior for MEM
  effsize1 <- 1 ## Effect Size for group 1
  effsize2 <- seq(1,2,0.1) ## Effect Sizes for group 2
  paramgrid <- expand.grid(N=N, n2=n2mult, v1=v1, es1=effsize1, es2=effsize2, power = NA, pnexch = NA, pnexch_power = NA)
  
  ns <- 10000 ## Set number of samples per parameter combination
  
  ## run simulation
  set.seed(1010)
  cutoffs <- c(0.2,0.7976455,0.8,0.9) 
  OA_cont_MEM_result_expandvar1 <- sim_samp(paramgrid, ns, p_exch, cutoffs, 0.05, 0.2)
  write.csv(OA_cont_MEM_result_expandvar1, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/OA_cont_MEM_result_expandvar1.csv")

## Plot results
## simulated sample plots
ggplot(OA_cont_MEM_result, aes(es2-es1, N, fill=power)) +
  geom_tile() + geom_text(aes(label=power)) + facet_grid(var_c~.) +
  labs(title = "Heatmap of T-test Power for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = 'Type I Error/Power') +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)")

for(i in 1:length(cutoffs)){
  fillvar <- colnames(OA_cont_MEM_result)[8+i]
  
  print(
    ggplot(OA_cont_MEM_result, aes(es2-es1, N, fill=.data[[fillvar]])) + 
      geom_tile() + geom_text(aes(label=.data[[fillvar]])) + facet_grid(var_c~.) + 
      labs(title = paste("Heatmap of P(Not Exchangeable) >=", cutoffs[i], 
                         "for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
                         sep = ' '),
           fill = "Type I Error/Power") +
      xlab("Difference in Effect Sizes") +
      ylab("Total Sample Size (N)")
  )
}

# ggplot(OA_cont_MEM_result, aes(es2-es1, N, fill=pnexch_power)) + 
#   geom_tile() + geom_text(aes(label=pnexch_power)) + facet_grid(var_c~.) + 
#   labs(title = "Heatmap of P(Not Exchangeable) Cutoff Needed for 5% Type I Error and 20% Type II Error for Varying Sample Sizes, Variances, and Differences in Effect Sizes")





####################################################################################
# Now, the subgroups will be unequal sizes (n1 != n2), variance and priors will be #
# equal and fixed, and the difference in effect size will vary.                    #
####################################################################################

##Set Simulation Parameters
n2mult = c(0.1, 0.25) ##percentage of total sample size that group 2 will be
v1 = v2 = 10 ## constant variance for easier visualization
paramgrid <- expand.grid(N=N, n2=n2mult, v1=v1, es1=effsize1, es2=effsize2, power = NA, pexch = NA)

## run simulation
set.seed(1001)
OA_cont_MEM_diffn_result <- sim_samp(paramgrid,ns,p_exch, c(0.2,0.7976455,0.8,0.9), 0.05, 0.2)
write.csv(OA_cont_MEM_diffn_result, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/OA_cont_MEM_diffn_result.csv")

#create n2mult character variable for plotting
OA_cont_MEM_diffn_result$n2mult_c <- NA
OA_cont_MEM_diffn_result$n2mult_c[OA_cont_MEM_diffn_result$n2 == 0.1] <- 'n2 = 0.1*N'
OA_cont_MEM_diffn_result$n2mult_c[OA_cont_MEM_diffn_result$n2 == 0.25] <- 'n2 = 0.25*N'

## Plot results
## simulated sample results
ggplot(OA_cont_MEM_diffn_result, aes(es2-es1, N, fill=power)) +
  geom_tile() + geom_text(aes(label=power)) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of T-test Power for Varying Sample Sizes per Group and Differences in Effect Sizes",
       fill = 'Type I Error/Power') +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)")

for(i in 1:length(cutoffs)){
  fillvar <- colnames(OA_cont_MEM_diffn_result)[8+i]
  
  print(
    ggplot(OA_cont_MEM_diffn_result, aes(es2-es1, N, fill=.data[[fillvar]])) + 
      geom_tile() + geom_text(aes(label=.data[[fillvar]])) + facet_grid(n2mult_c~.) + 
      labs(title = paste("Heatmap of P(Not Exchangeable) >=", cutoffs[i], 
                         "for Varying Sample Sizes per Group and Differences in Effect Sizes",
                         sep = ' '),
           fill = "Type I Error/Power") +
      xlab("Difference in Effect Sizes") +
      ylab("Total Sample Size (N)")
  )
}

# ggplot(OA_cont_MEM_diffn_result, aes(es2-es1, N, fill=pnexch_power)) +
#   geom_tile() + geom_text(aes(label=pnexch_power)) + facet_grid(n2mult_c~.) +
#   labs(title = "Heatmap of Mean P(Not Exchangeable) Cutoff Needed for 5% Type I Error and 20% Type II Error for Varying Sample Sizes and Differences in Effect Sizes")
