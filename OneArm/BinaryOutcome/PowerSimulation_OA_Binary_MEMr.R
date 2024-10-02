library(tidyverse)
library(snowfall)

#interaction lm only functional for two sources/groups
source("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/AlesandLMSimFunc_OneArm_BinaryOutcome.R")

# #function to input beta values and output control/treatment probabilities
# invlogit <- function(mu0, beta_t){
#   pc <- exp(mu0)/(1+exp(mu0)) #control group event probability
#   pt <- exp(mu0 + beta_t)/(1+exp(mu0+beta_t)) #treatment group event probability
#   
#   OR <- exp(mu0 + beta_t)/exp(mu0)
#   
#   return(data.frame(prob_cont = pc, prob_treat = pt, oddsratio = OR))
# }
# invlogit(mu0 = c(5,5), beta_t = c(1,1))
# invlogit(mu0 = c(5,5), beta_t = c(1,2))
# invlogit(mu0 = c(5,5), beta_t = c(1,3))
# invlogit(mu0 = c(-5,-5), beta_t = c(1,1))
# 
# #function to input probabilities and output beta coefficients and odds ratio
# logitfunc <- function(pc, pt){
#   mu0 <- -1*log(-1*(1-(1/pc)))
#   beta_t <- -1*mu0 -log(-1*(1-(1/pt)))
#   
#   OR <- exp(mu0 + beta_t)/exp(mu0)
#   
#   return(data.frame(mu0 = mu0, beta_t = beta_t, oddsratio = OR))
# }
# logitfunc(pc = 0.1, pt = 0.2)
# logitfunc(pc = 0.1, pt = 0.3)
# logitfunc(pc = 0.1, pt = 0.4)
# logitfunc(pc = 0.1, pt = 0.5) #use pt = 0.1-0.5 for reasonable odds ratios
# logitfunc(pc = 0.1, pt = 0.6)
# logitfunc(pc = 0.1, pt = 0.7)
# logitfunc(pc = 0.1, pt = 0.8)
# logitfunc(pc = 0.1, pt = 0.9)

#function to run simulation in parallel and output results nicely
runsim_MEM1_onearm_binary <- function(nsim, pars) {
  
  summ <- data.frame(totsampsize = rep(NA, length(pars)*nsim),
                     pc = rep(NA,length(pars)*nsim),
                     pt = rep(NA, length(pars)*nsim),
                     pexch = rep(NA, length(pars)*nsim), 
                     pnexch = rep(NA, length(pars)*nsim), 
                     pval = rep(NA, length(pars)*nsim))
  
  
  pb= txtProgressBar(min = 0, max = length(pars), style = 3, char=":)")
  for(i in 1:length(pars)) {
    ## Define the parameters
    par <- pars[[i]] 

    #### Do the simulation with these settings
    sfInit(parallel=T, cpus=12, type='SOCK')
    #sfLibrary(R2jags)
    sfLibrary(mvtnorm)
    sfLibrary(gsDesign)
    sfLibrary(car)
    sfLibrary(nlme)
    sfLibrary(HDInterval)
    sfLibrary(lmtest)
    sfClusterSetupRNG(seed=12345)
    
    sfExportAll()
    sim.results <- t(sfLapply(rep(1,nsim),
                              function(x) doSim_MEM1_OneArm_binary(n=par$n, probs=par$probs, 
                                                BICiterations=par$BICiterations, 
                                                marginal=par$marginal)))
    sfStop()
    
    result <- as.data.frame(unlist(sim.results))
    
    summ[(nsim*(i-1) + 1):(nsim*i),1] <- (result %>% filter(row_number() %% 6 == 1))
    summ[(nsim*(i-1) + 1):(nsim*i),2] <- (result %>% filter(row_number() %% 6 == 2))
    summ[(nsim*(i-1) + 1):(nsim*i),3] <- (result %>% filter(row_number() %% 6 == 3))
    summ[(nsim*(i-1) + 1):(nsim*i),4] <- (result %>% filter(row_number() %% 6 == 4))
    summ[(nsim*(i-1) + 1):(nsim*i),5] <- (result %>% filter(row_number() %% 6 == 5))
    summ[(nsim*(i-1) + 1):(nsim*i),6] <- (result %>% filter(row_number() %% 6 == 0))

    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(summ)
  
}

#############################################
##############EQUAL SAMPLE###################
#############################################

#function to format inputs 
parfunc_onearm_binary <- function(eventprobs){
  pars=list(
    #n = sample sizes, m0 = baseline mean for each source, 
    #p = probability of treatment being assigned, 
    #beta_t = treatment main effect 
    list(n=c(100, 100), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(150, 150), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(200, 200), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(250, 250), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(300, 300), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(350, 350), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(400, 400), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(450, 450), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(500, 500), probs=eventprobs, BICiterations=1000, 
         marginal="BIC")
  )
  
  return(pars)
}


#set number of samples per setting
nsim=1000



#scenario 1: equal source means (1 and 1)
system.time(equalsourceeff11 <- runsim_MEM1_onearm_binary(nsim, 
                                              parfunc_onearm_binary(eventprobs = c(0.1,0.1))))

equalsourceeff11$plt05 <- equalsourceeff11$pval < 0.05
equalsourceeff11$pnexchgt090 <- equalsourceeff11$pnexch > 0.90
equalsourceeff11$pnexchgt080 <- equalsourceeff11$pnexch > 0.80
equalsourceeff11$pnexchgt0798 <- equalsourceeff11$pnexch > 0.7976455
equalsourceeff11$pnexchgt020 <- equalsourceeff11$pnexch > 0.2


#scenario 2: similar source means (1 and 2)
system.time(diffsourceeff12 <- runsim_MEM1_onearm_binary(nsim, 
                                             parfunc_onearm_binary(eventprobs = c(0.1,0.2))))
diffsourceeff12$plt05 <- diffsourceeff12$pval < 0.05
diffsourceeff12$pnexchgt090 <- diffsourceeff12$pnexch > 0.90
diffsourceeff12$pnexchgt080 <- diffsourceeff12$pnexch > 0.80
diffsourceeff12$pnexchgt0798 <- diffsourceeff12$pnexch > 0.7976455
diffsourceeff12$pnexchgt020 <- diffsourceeff12$pnexch > 0.2

#scenario 3: different source means (1 and 3)
system.time(diffsourceeff13 <- runsim_MEM1_onearm_binary(nsim, 
                                             parfunc_onearm_binary(eventprobs = c(0.1,0.3))))
diffsourceeff13$plt05 <- diffsourceeff13$pval < 0.05
diffsourceeff13$pnexchgt090 <- diffsourceeff13$pnexch > 0.90
diffsourceeff13$pnexchgt080 <- diffsourceeff13$pnexch > 0.80
diffsourceeff13$pnexchgt0798 <- diffsourceeff13$pnexch > 0.7976455
diffsourceeff13$pnexchgt020 <- diffsourceeff13$pnexch > 0.2

#scenario 4: different source means (1 and 4)
system.time(diffsourceeff14 <- runsim_MEM1_onearm_binary(nsim, 
                                             parfunc_onearm_binary(eventprobs = c(0.1,0.4))))
diffsourceeff14$plt05 <- diffsourceeff14$pval < 0.05
diffsourceeff14$pnexchgt090 <- diffsourceeff14$pnexch > 0.90
diffsourceeff14$pnexchgt080 <- diffsourceeff14$pnexch > 0.80
diffsourceeff14$pnexchgt0798 <- diffsourceeff14$pnexch > 0.7976455
diffsourceeff14$pnexchgt020 <- diffsourceeff14$pnexch > 0.2

#scenario 5: different source means (1 and 5)
system.time(diffsourceeff15 <- runsim_MEM1_onearm_binary(nsim, 
                                             parfunc_onearm_binary(eventprobs = c(0.1,0.5))))
diffsourceeff15$plt05 <- diffsourceeff15$pval < 0.05
diffsourceeff15$pnexchgt090 <- diffsourceeff15$pnexch > 0.90
diffsourceeff15$pnexchgt080 <- diffsourceeff15$pnexch > 0.80
diffsourceeff15$pnexchgt0798 <- diffsourceeff15$pnexch > 0.7976455
diffsourceeff15$pnexchgt020 <- diffsourceeff15$pnexch > 0.2



#combine data sets
allresults <- rbind(equalsourceeff11,
                    diffsourceeff12, diffsourceeff13, diffsourceeff14, diffsourceeff15)

all_summ0 <- allresults %>% group_by(pc, pt, totsampsize)

all_summ <- all_summ0 %>% summarize(plt05 = mean(plt05), pnexchgt090 = mean(pnexchgt090),
                                    pnexchgt080 = mean(pnexchgt080),
                                    pnexchgt0798 = mean(pnexchgt0798), 
                                    pnexchgt020 = mean(pnexchgt020))

write.csv(all_summ, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/OA_binary_MEMr_result_summary.csv")
write.csv(allresults, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/OA_binary_MEMr_result_raw.csv")



#############################################
#############UNEQUAL SAMPLE##################
#############################################

#function to format inputs 
parfunc_onearm_binary_unequalsamp <- function(eventprobs, n2mult){
  pars=list(
    #n = sample sizes, m0 = baseline mean for each source, 
    #p = probability of treatment being assigned, 
    #beta_t = treatment main effect 
    list(n=c(200*(1-n2mult), 200*n2mult), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(300*(1-n2mult), 300*n2mult), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(400*(1-n2mult), 400*n2mult), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(500*(1-n2mult), 500*n2mult), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(600*(1-n2mult), 600*n2mult), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(700*(1-n2mult), 700*n2mult), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(800*(1-n2mult), 800*n2mult), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(900*(1-n2mult), 900*n2mult), probs=eventprobs, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(1000*(1-n2mult), 1000*n2mult), probs=eventprobs, BICiterations=1000, 
         marginal="BIC")
  )
  
  return(pars)
}


#set number of samples per setting
nsim=1000


#scenario 1: equal source means (1 and 1)
system.time(equalsourceeff11_unequalsamp01 <- runsim_MEM1_onearm_binary(nsim, 
                                                                        parfunc_onearm_binary_unequalsamp(eventprobs = c(0.1,0.1), 
                                                                                                          n2mult = 0.1)))

equalsourceeff11_unequalsamp01$n2mult <- 0.1
equalsourceeff11_unequalsamp01$plt05 <- equalsourceeff11_unequalsamp01$pval < 0.05
equalsourceeff11_unequalsamp01$pnexchgt090 <- equalsourceeff11_unequalsamp01$pnexch > 0.90
equalsourceeff11_unequalsamp01$pnexchgt080 <- equalsourceeff11_unequalsamp01$pnexch > 0.80
equalsourceeff11_unequalsamp01$pnexchgt0798 <- equalsourceeff11_unequalsamp01$pnexch > 0.7976455
equalsourceeff11_unequalsamp01$pnexchgt020 <- equalsourceeff11_unequalsamp01$pnexch > 0.2

system.time(equalsourceeff11_unequalsamp025 <- runsim_MEM1_onearm_binary(nsim, 
                                                                        parfunc_onearm_binary_unequalsamp(eventprobs = c(0.1,0.1), 
                                                                                                          n2mult = 0.25)))
equalsourceeff11_unequalsamp025$n2mult <- 0.25
equalsourceeff11_unequalsamp025$plt05 <- equalsourceeff11_unequalsamp025$pval < 0.05
equalsourceeff11_unequalsamp025$pnexchgt090 <- equalsourceeff11_unequalsamp025$pnexch > 0.90
equalsourceeff11_unequalsamp025$pnexchgt080 <- equalsourceeff11_unequalsamp025$pnexch > 0.80
equalsourceeff11_unequalsamp025$pnexchgt0798 <- equalsourceeff11_unequalsamp025$pnexch > 0.7976455
equalsourceeff11_unequalsamp025$pnexchgt020 <- equalsourceeff11_unequalsamp025$pnexch > 0.2


#scenario 2: similar source means (1 and 2)
system.time(diffsourceeff12_unequalsamp01  <- runsim_MEM1_onearm_binary(nsim, 
                                                         parfunc_onearm_binary_unequalsamp(eventprobs = c(0.1,0.2), 
                                                                                           n2mult = 0.1)))
diffsourceeff12_unequalsamp01$n2mult <- 0.1
diffsourceeff12_unequalsamp01$plt05 <- diffsourceeff12_unequalsamp01$pval < 0.05
diffsourceeff12_unequalsamp01$pnexchgt090 <- diffsourceeff12_unequalsamp01$pnexch > 0.90
diffsourceeff12_unequalsamp01$pnexchgt080 <- diffsourceeff12_unequalsamp01$pnexch > 0.80
diffsourceeff12_unequalsamp01$pnexchgt0798 <- diffsourceeff12_unequalsamp01$pnexch > 0.7976455
diffsourceeff12_unequalsamp01$pnexchgt020 <- diffsourceeff12_unequalsamp01$pnexch > 0.20

system.time(diffsourceeff12_unequalsamp025  <- runsim_MEM1_onearm_binary(nsim, 
                                                                        parfunc_onearm_binary_unequalsamp(eventprobs = c(0.1,0.2), 
                                                                                                          n2mult = 0.25)))
diffsourceeff12_unequalsamp025$n2mult <- 0.25
diffsourceeff12_unequalsamp025$plt05 <- diffsourceeff12_unequalsamp025$pval < 0.05
diffsourceeff12_unequalsamp025$pnexchgt090 <- diffsourceeff12_unequalsamp025$pnexch > 0.90
diffsourceeff12_unequalsamp025$pnexchgt080 <- diffsourceeff12_unequalsamp025$pnexch > 0.80
diffsourceeff12_unequalsamp025$pnexchgt0798 <- diffsourceeff12_unequalsamp025$pnexch > 0.7976455
diffsourceeff12_unequalsamp025$pnexchgt020 <- diffsourceeff12_unequalsamp025$pnexch > 0.20

#scenario 3: different source means (1 and 3)
system.time(diffsourceeff13_unequalsamp01  <- runsim_MEM1_onearm_binary(nsim, 
                                                         parfunc_onearm_binary_unequalsamp(eventprobs = c(0.1,0.3), 
                                                                                           n2mult = 0.1)))
diffsourceeff13_unequalsamp01$n2mult <- 0.1
diffsourceeff13_unequalsamp01$plt05 <- diffsourceeff13_unequalsamp01$pval < 0.05
diffsourceeff13_unequalsamp01$pnexchgt090 <- diffsourceeff13_unequalsamp01$pnexch > 0.90
diffsourceeff13_unequalsamp01$pnexchgt080 <- diffsourceeff13_unequalsamp01$pnexch > 0.80
diffsourceeff13_unequalsamp01$pnexchgt0798 <- diffsourceeff13_unequalsamp01$pnexch > 0.7976455
diffsourceeff13_unequalsamp01$pnexchgt020 <- diffsourceeff13_unequalsamp01$pnexch > 0.20

system.time(diffsourceeff13_unequalsamp025  <- runsim_MEM1_onearm_binary(nsim, 
                                                                        parfunc_onearm_binary_unequalsamp(eventprobs = c(0.1,0.3), 
                                                                                                          n2mult = 0.25)))
diffsourceeff13_unequalsamp025$n2mult <- 0.25
diffsourceeff13_unequalsamp025$plt05 <- diffsourceeff13_unequalsamp025$pval < 0.05
diffsourceeff13_unequalsamp025$pnexchgt090 <- diffsourceeff13_unequalsamp025$pnexch > 0.90
diffsourceeff13_unequalsamp025$pnexchgt080 <- diffsourceeff13_unequalsamp025$pnexch > 0.80
diffsourceeff13_unequalsamp025$pnexchgt0798 <- diffsourceeff13_unequalsamp025$pnexch > 0.7976455
diffsourceeff13_unequalsamp025$pnexchgt020 <- diffsourceeff13_unequalsamp025$pnexch > 0.20

#scenario 4: different source means (1 and 4)
system.time(diffsourceeff14_unequalsamp01  <- runsim_MEM1_onearm_binary(nsim, 
                                                         parfunc_onearm_binary_unequalsamp(eventprobs = c(0.1,0.4), 
                                                                                           n2mult = 0.1)))
diffsourceeff14_unequalsamp01$n2mult <- 0.1
diffsourceeff14_unequalsamp01$plt05 <- diffsourceeff14_unequalsamp01$pval < 0.05
diffsourceeff14_unequalsamp01$pnexchgt090 <- diffsourceeff14_unequalsamp01$pnexch > 0.90
diffsourceeff14_unequalsamp01$pnexchgt080 <- diffsourceeff14_unequalsamp01$pnexch > 0.80
diffsourceeff14_unequalsamp01$pnexchgt0798 <- diffsourceeff14_unequalsamp01$pnexch > 0.7976455
diffsourceeff14_unequalsamp01$pnexchgt020 <- diffsourceeff14_unequalsamp01$pnexch > 0.20

system.time(diffsourceeff14_unequalsamp025  <- runsim_MEM1_onearm_binary(nsim, 
                                                                        parfunc_onearm_binary_unequalsamp(eventprobs = c(0.1,0.4), 
                                                                                                          n2mult = 0.25)))
diffsourceeff14_unequalsamp025$n2mult <- 0.25
diffsourceeff14_unequalsamp025$plt05 <- diffsourceeff14_unequalsamp025$pval < 0.05
diffsourceeff14_unequalsamp025$pnexchgt090 <- diffsourceeff14_unequalsamp025$pnexch > 0.90
diffsourceeff14_unequalsamp025$pnexchgt080 <- diffsourceeff14_unequalsamp025$pnexch > 0.80
diffsourceeff14_unequalsamp025$pnexchgt0798 <- diffsourceeff14_unequalsamp025$pnexch > 0.7976455
diffsourceeff14_unequalsamp025$pnexchgt020 <- diffsourceeff14_unequalsamp025$pnexch > 0.20

#scenario 5: different source means (1 and 5)
system.time(diffsourceeff15_unequalsamp01  <- runsim_MEM1_onearm_binary(nsim, 
                                                         parfunc_onearm_binary_unequalsamp(eventprobs = c(0.1,0.5), 
                                                                                           n2mult = 0.1)))
diffsourceeff15_unequalsamp01$n2mult <- 0.1
diffsourceeff15_unequalsamp01$plt05 <- diffsourceeff15_unequalsamp01$pval < 0.05
diffsourceeff15_unequalsamp01$pnexchgt090 <- diffsourceeff15_unequalsamp01$pnexch > 0.90
diffsourceeff15_unequalsamp01$pnexchgt080 <- diffsourceeff15_unequalsamp01$pnexch > 0.80
diffsourceeff15_unequalsamp01$pnexchgt0798 <- diffsourceeff15_unequalsamp01$pnexch > 0.7976455
diffsourceeff15_unequalsamp01$pnexchgt020 <- diffsourceeff15_unequalsamp01$pnexch > 0.20

system.time(diffsourceeff15_unequalsamp025  <- runsim_MEM1_onearm_binary(nsim, 
                                                                        parfunc_onearm_binary_unequalsamp(eventprobs = c(0.1,0.5), 
                                                                                                          n2mult = 0.25)))
diffsourceeff15_unequalsamp025$n2mult <- 0.25
diffsourceeff15_unequalsamp025$plt05 <- diffsourceeff15_unequalsamp025$pval < 0.05
diffsourceeff15_unequalsamp025$pnexchgt090 <- diffsourceeff15_unequalsamp025$pnexch > 0.90
diffsourceeff15_unequalsamp025$pnexchgt080 <- diffsourceeff15_unequalsamp025$pnexch > 0.80
diffsourceeff15_unequalsamp025$pnexchgt0798 <- diffsourceeff15_unequalsamp025$pnexch > 0.7976455
diffsourceeff15_unequalsamp025$pnexchgt020 <- diffsourceeff15_unequalsamp025$pnexch > 0.20

#combine data sets
allresults_unequalsamp <- rbind(equalsourceeff11_unequalsamp01,equalsourceeff11_unequalsamp025,
                                diffsourceeff12_unequalsamp01, diffsourceeff12_unequalsamp025,
                                diffsourceeff13_unequalsamp01, diffsourceeff13_unequalsamp025, 
                                diffsourceeff14_unequalsamp01, diffsourceeff14_unequalsamp025, 
                                diffsourceeff15_unequalsamp01, diffsourceeff15_unequalsamp025)

all_summ0_unequalsamp <- allresults_unequalsamp %>% group_by(pc, pt, totsampsize, n2mult)

all_summ_unequalsamp <- all_summ0_unequalsamp %>% summarize(plt05 = mean(plt05), pnexchgt090 = mean(pnexchgt090), 
                                                            pnexchgt080 = mean(pnexchgt080),
                                                            pnexchgt0798 = mean(pnexchgt0798), 
                                                            pnexchgt020 = mean(pnexchgt020))

write.csv(all_summ_unequalsamp, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/OA_binary_MEMr_unequalsamp_result_summary.csv")
write.csv(allresults_unequalsamp, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/OA_binary_MEMr_unequalsamp_result_raw.csv")
