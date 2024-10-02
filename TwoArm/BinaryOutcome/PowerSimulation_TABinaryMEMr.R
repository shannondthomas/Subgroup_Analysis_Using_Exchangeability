library(tidyverse)
library(snowfall)

#interaction lm only functional for two sources/groups
source("C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/AlesandLMSimFunc_BinaryOutcome.R")

# #function to input probabilities and output beta coefficients and odds ratio
# logitfunc <- function(pc1,pc2,pt1,pt2){
#   b0 <- log(pc1/(1-pc1)) #log odds for control group in group 1
#   b1 <- -1*b0 + log(pc2/(1-pc2)) # log odds for control group in group 2
#   b2 <- -1*b0 + log(pt1/(1-pt1)) # log odds for treatment group in group 1
#   b3 <- -1*b0 + (-1)*b1 + (-1)*b2 + log(pt2/(1-pt2)) # interaction effect (change in log odds for group 2 vs group 1)
#     
#   OR <- exp(b3+b1)
#     
#   return(data.frame(b0=b0,b1=b1,b2=b2,b3=b3,OR=OR))
# }
# logitfunc(pc1=0.1,pc2=0.1,pt1=0.2,pt2=0.2)
# logitfunc(pc1=0.1,pc2=0.1,pt1=0.2,pt2=0.3)
# logitfunc(pc1=0.1,pc2=0.1,pt1=0.2,pt2=0.4)
# logitfunc(pc1=0.1,pc2=0.1,pt1=0.2,pt2=0.5) 
# logitfunc(pc1=0.1,pc2=0.1,pt1=0.2,pt2=0.6)
# logitfunc(pc1=0.1,pc2=0.1,pt1=0.2,pt2=0.7) #use pt2 = 0.1-0.7 for reasonable odds ratios
# logitfunc(pc1=0.1,pc2=0.1,pt1=0.2,pt2=0.8)
# logitfunc(pc1=0.1,pc2=0.1,pt1=0.2,pt2=0.9)


#function to run simulation in parallel and output results nicely
runsim_MEM1_binary <- function(nsim, pars) {
  
  #summary matrix: first weight 1 (non-exchangeable) and second is weight 2 (exchangeable)
  summ <- data.frame(totsampsize = rep(NA, length(pars)*nsim),
                     pc1 = rep(NA, length(pars)*nsim),
                     pc2 = rep(NA, length(pars)*nsim),                     
                     pt1 = rep(NA, length(pars)*nsim),
                     pt2 = rep(NA, length(pars)*nsim),
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
                              function(x) doSim_MEM1_binary(n=par$n, pc=par$pc, pt=par$pt,
                                                p = par$p,
                                                BICiterations=par$BICiterations, 
                                                marginal=par$marginal)))
    sfStop()
    
    result <- as.data.frame(unlist(sim.results))
    
    summ[(nsim*(i-1) + 1):(nsim*i),1] <- (result %>% filter(row_number() %% 8 == 1))
    summ[(nsim*(i-1) + 1):(nsim*i),2] <- (result %>% filter(row_number() %% 8 == 2))
    summ[(nsim*(i-1) + 1):(nsim*i),3] <- (result %>% filter(row_number() %% 8 == 3))
    summ[(nsim*(i-1) + 1):(nsim*i),4] <- (result %>% filter(row_number() %% 8 == 4))
    summ[(nsim*(i-1) + 1):(nsim*i),5] <- (result %>% filter(row_number() %% 8 == 5))
    summ[(nsim*(i-1) + 1):(nsim*i),6] <- (result %>% filter(row_number() %% 8 == 6))    
    summ[(nsim*(i-1) + 1):(nsim*i),7] <- (result %>% filter(row_number() %% 8 == 7))
    summ[(nsim*(i-1) + 1):(nsim*i),8] <- (result %>% filter(row_number() %% 8 == 0))
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(summ)
  
}

#############################################
##############EQUAL SAMPLE###################
#############################################

#function to format inputs 
parfunc_twoarm <- function(probc, probt, effectsizes){
  pars=list(
    #n = sample sizes, m0 = baseline mean for each source, 
    #p = probability of treatment being assigned, 
    #beta_t = treatment main effect 
    list(n=c(100, 100), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(150, 150), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(200, 200), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(250, 250), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(300, 300), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(350, 350), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(400, 400), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(450, 450), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(500, 500), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC")
  )
  
  return(pars)
}


#set number of samples per setting
nsim=1000


#scenario 1: equal treatment effects (1 and 1)
system.time(equaltrteff11 <- runsim_MEM1_binary(nsim, 
                                              parfunc_twoarm(probc = c(0.1,0.1), probt = c(0.2,0.2),
                                                               effectsizes = c(1,1))))

equaltrteff11$plt05 <- equaltrteff11$pval < 0.05
equaltrteff11$pnexchgt090 <- equaltrteff11$pnexch > 0.9
equaltrteff11$pnexchgt080 <- equaltrteff11$pnexch > 0.8
equaltrteff11$pnexchgt0798 <- equaltrteff11$pnexch > 0.7976455
equaltrteff11$pnexchgt020 <- equaltrteff11$pnexch > 0.2

table(equaltrteff11$pnexchgt0798, equaltrteff11$plt05)



#scenario 2: similar treatment effects (1 and 2)
system.time(difftrteff12 <- runsim_MEM1_binary(nsim, 
                                             parfunc_twoarm(probc = c(0.1,0.1), probt = c(0.2,0.3),
                                                           effectsizes = c(1,2))))
difftrteff12$plt05 <- difftrteff12$pval < 0.05
difftrteff12$pnexchgt090 <- difftrteff12$pnexch > 0.9
difftrteff12$pnexchgt080 <- difftrteff12$pnexch > 0.8
difftrteff12$pnexchgt0798 <- difftrteff12$pnexch > 0.7976455
difftrteff12$pnexchgt020 <- difftrteff12$pnexch > 0.2

table(difftrteff12$pnexchgt0798, difftrteff12$plt05)

#scenario 3: different treatment effects (1 and 3)
system.time(difftrteff13 <- runsim_MEM1_binary(nsim, 
                                             parfunc_twoarm(probc = c(0.1,0.1), probt = c(0.2,0.4),
                                                           effectsizes = c(1,3))))
difftrteff13$plt05 <- difftrteff13$pval < 0.05
difftrteff13$pnexchgt090 <- difftrteff13$pnexch > 0.9
difftrteff13$pnexchgt080 <- difftrteff13$pnexch > 0.8
difftrteff13$pnexchgt0798 <- difftrteff13$pnexch > 0.7976455
difftrteff13$pnexchgt020 <- difftrteff13$pnexch > 0.2

table(difftrteff13$pnexchgt0798, difftrteff13$plt05)

#scenario 4: different treatment effects (1 and 4)
system.time(difftrteff14 <- runsim_MEM1_binary(nsim, 
                                             parfunc_twoarm(probc = c(0.1,0.1), probt = c(0.2,0.5),
                                                           effectsizes = c(1,4))))
difftrteff14$plt05 <- difftrteff14$pval < 0.05
difftrteff14$pnexchgt090 <- difftrteff14$pnexch > 0.9
difftrteff14$pnexchgt080 <- difftrteff14$pnexch > 0.8
difftrteff14$pnexchgt0798 <- difftrteff14$pnexch > 0.7976455
difftrteff14$pnexchgt020 <- difftrteff14$pnexch > 0.2

table(difftrteff14$pnexchgt0798, difftrteff14$plt05)

#scenario 5: different treatment effects (1 and 5)
system.time(difftrteff15 <- runsim_MEM1_binary(nsim, 
                                             parfunc_twoarm(probc = c(0.1,0.1), probt = c(0.2,0.6),
                                                           effectsizes = c(1,5))))
difftrteff15$plt05 <- difftrteff15$pval < 0.05
difftrteff15$pnexchgt090 <- difftrteff15$pnexch > 0.9
difftrteff15$pnexchgt080 <- difftrteff15$pnexch > 0.8
difftrteff15$pnexchgt0798 <- difftrteff15$pnexch > 0.7976455
difftrteff15$pnexchgt020 <- difftrteff15$pnexch > 0.2

table(difftrteff15$pnexchgt0798, difftrteff15$plt05)

#scenario 6: different treatment effects (1 and 6)
system.time(difftrteff16 <- runsim_MEM1_binary(nsim, 
                                             parfunc_twoarm(probc = c(0.1,0.1), probt = c(0.2,0.7),
                                                           effectsizes = c(1,6))))
difftrteff16$plt05 <- difftrteff16$pval < 0.05
difftrteff16$pnexchgt090 <- difftrteff16$pnexch > 0.9
difftrteff16$pnexchgt080 <- difftrteff16$pnexch > 0.8
difftrteff16$pnexchgt0798 <- difftrteff16$pnexch > 0.7976455
difftrteff16$pnexchgt020 <- difftrteff16$pnexch > 0.2

table(difftrteff16$pnexchgt0798, difftrteff16$plt05)



#combine data sets
allresults <- rbind(equaltrteff11,
                    difftrteff12, difftrteff13, difftrteff14, difftrteff15, difftrteff16)

all_summ0 <- allresults %>% group_by(pc1,pc2,pt1,pt2, totsampsize)

all_summ <- all_summ0 %>% summarize(plt05 = mean(plt05), pnexchgt090 = mean(pnexchgt090), 
                                    pnexchgt080 = mean(pnexchgt080), pnexchgt0798 = mean(pnexchgt0798), 
                                    pnexchgt020 = mean(pnexchgt020))

write.csv(all_summ, "C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/resultsdata/TA_binary_MEMr_result_summary.csv")
write.csv(allresults, "C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/resultsdata/TA_binary_MEMr_result_raw.csv")


#############################################
#############UNEQUAL SAMPLE##################
#############################################

#function to format inputs for unequal samples
parfunc_twoarm_unequalsamp <- function(probc, probt, effectsizes, n2mult){
  pars=list(
    #n = sample sizes, m0 = baseline mean for each source, 
    #p = probability of treatment being assigned, 
    #beta_t = treatment main effect 
    list(n=c(200*(1-n2mult), 200*n2mult), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(300*(1-n2mult), 300*n2mult), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(400*(1-n2mult), 400*n2mult), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(500*(1-n2mult), 500*n2mult), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(600*(1-n2mult), 600*n2mult), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(700*(1-n2mult), 700*n2mult), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(800*(1-n2mult), 800*n2mult), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(900*(1-n2mult), 900*n2mult), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(1000*(1-n2mult), 1000*n2mult), pc = probc, pt = probt,   p=c(0.5, 0.5),
         beta_t=effectsizes, BICiterations=1000, 
         marginal="BIC")
  )
}

#set number of samples per setting
nsim=1000


#scenario 1: equal treatment effects (1 and 1)
system.time(equaltrteff11_n2mult01 <- runsim_MEM1_binary(nsim, 
                                                         parfunc_twoarm_unequalsamp(probc = c(0.1,0.1), probt = c(0.2,0.2),
                                                               effectsizes = c(1,1),
                                                               n2mult = 0.1)))

equaltrteff11_n2mult01$plt05 <- equaltrteff11_n2mult01$pval < 0.05
equaltrteff11_n2mult01$pnexchgt090 <- equaltrteff11_n2mult01$pnexch > 0.9
equaltrteff11_n2mult01$pnexchgt080 <- equaltrteff11_n2mult01$pnexch > 0.8
equaltrteff11_n2mult01$pnexchgt0798 <- equaltrteff11_n2mult01$pnexch > 0.7976455
equaltrteff11_n2mult01$pnexchgt020 <- equaltrteff11_n2mult01$pnexch > 0.2

table(equaltrteff11_n2mult01$pnexchgt0798, equaltrteff11_n2mult01$plt05)

system.time(equaltrteff11_n2mult025 <- runsim_MEM1_binary(nsim, 
                                                          parfunc_twoarm_unequalsamp(probc = c(0.1,0.1), probt = c(0.2,0.2),
                                                                        effectsizes = c(1,1),
                                                                        n2mult = 0.25)))

equaltrteff11_n2mult025$plt05 <- equaltrteff11_n2mult025$pval < 0.05
equaltrteff11_n2mult025$pnexchgt090 <- equaltrteff11_n2mult025$pnexch > 0.9
equaltrteff11_n2mult025$pnexchgt080 <- equaltrteff11_n2mult025$pnexch > 0.8
equaltrteff11_n2mult025$pnexchgt0798 <- equaltrteff11_n2mult025$pnexch > 0.7976455
equaltrteff11_n2mult025$pnexchgt020 <- equaltrteff11_n2mult025$pnexch > 0.2

table(equaltrteff11_n2mult025$pnexchgt0798, equaltrteff11_n2mult025$plt05)

#scenario 2: similar treatment effects (1 and 2)
system.time(difftrteff12_n2mult01 <- runsim_MEM1_binary(nsim, 
                                                        parfunc_twoarm_unequalsamp(probc = c(0.1,0.1), probt = c(0.2,0.3),
                                                              effectsizes = c(1,2),
                                                              n2mult = 0.1)))
difftrteff12_n2mult01$plt05 <- difftrteff12_n2mult01$pval < 0.05
difftrteff12_n2mult01$pnexchgt090 <- difftrteff12_n2mult01$pnexch > 0.9
difftrteff12_n2mult01$pnexchgt080 <- difftrteff12_n2mult01$pnexch > 0.8
difftrteff12_n2mult01$pnexchgt0798 <- difftrteff12_n2mult01$pnexch > 0.7976455
difftrteff12_n2mult01$pnexchgt020 <- difftrteff12_n2mult01$pnexch > 0.2

table(difftrteff12_n2mult01$pnexchgt0798, difftrteff12_n2mult01$plt05)


system.time(difftrteff12_n2mult025 <- runsim_MEM1_binary(nsim, 
                                                        parfunc_twoarm_unequalsamp(probc = c(0.1,0.1), probt = c(0.2,0.3),
                                                                                   effectsizes = c(1,2),
                                                                                   n2mult = 0.25)))
difftrteff12_n2mult025$plt05 <- difftrteff12_n2mult025$pval < 0.05
difftrteff12_n2mult025$pnexchgt090 <- difftrteff12_n2mult025$pnexch > 0.9
difftrteff12_n2mult025$pnexchgt080 <- difftrteff12_n2mult025$pnexch > 0.8
difftrteff12_n2mult025$pnexchgt0798 <- difftrteff12_n2mult025$pnexch > 0.7976455
difftrteff12_n2mult025$pnexchgt020 <- difftrteff12_n2mult025$pnexch > 0.2

table(difftrteff12_n2mult025$pnexchgt0798, difftrteff12_n2mult025$plt05)

#scenario 3: different treatment effects (1 and 3)
system.time(difftrteff13_n2mult01 <- runsim_MEM1_binary(nsim, 
                                                        parfunc_twoarm_unequalsamp(probc = c(0.1,0.1), probt = c(0.2,0.4),
                                                              effectsizes = c(1,3),
                                                              n2mult = 0.1)))
difftrteff13_n2mult01$plt05 <- difftrteff13_n2mult01$pval < 0.05
difftrteff13_n2mult01$pnexchgt090 <- difftrteff13_n2mult01$pnexch > 0.9
difftrteff13_n2mult01$pnexchgt080 <- difftrteff13_n2mult01$pnexch > 0.8
difftrteff13_n2mult01$pnexchgt0798 <- difftrteff13_n2mult01$pnexch > 0.7976455
difftrteff13_n2mult01$pnexchgt020 <- difftrteff13_n2mult01$pnexch > 0.2

table(difftrteff13_n2mult01$pnexchgt0798, difftrteff13_n2mult01$plt05)


system.time(difftrteff13_n2mult025 <- runsim_MEM1_binary(nsim, 
                                                        parfunc_twoarm_unequalsamp(probc = c(0.1,0.1), probt = c(0.2,0.4),
                                                                                   effectsizes = c(1,3),
                                                                                   n2mult = 0.25)))
difftrteff13_n2mult025$plt05 <- difftrteff13_n2mult025$pval < 0.05
difftrteff13_n2mult025$pnexchgt090 <- difftrteff13_n2mult025$pnexch > 0.9
difftrteff13_n2mult025$pnexchgt080 <- difftrteff13_n2mult025$pnexch > 0.8
difftrteff13_n2mult025$pnexchgt0798 <- difftrteff13_n2mult025$pnexch > 0.7976455
difftrteff13_n2mult025$pnexchgt020 <- difftrteff13_n2mult025$pnexch > 0.2

table(difftrteff13_n2mult025$pnexchgt0798, difftrteff13_n2mult025$plt05)


#scenario 4: different treatment effects (1 and 4)
system.time(difftrteff14_n2mult01 <- runsim_MEM1_binary(nsim, 
                                                        parfunc_twoarm_unequalsamp(probc = c(0.1,0.1), probt = c(0.2,0.5),
                                                              effectsizes = c(1,4),
                                                              n2mult = 0.1)))
difftrteff14_n2mult01$plt05 <- difftrteff14_n2mult01$pval < 0.05
difftrteff14_n2mult01$pnexchgt090 <- difftrteff14_n2mult01$pnexch > 0.9
difftrteff14_n2mult01$pnexchgt080 <- difftrteff14_n2mult01$pnexch > 0.8
difftrteff14_n2mult01$pnexchgt0798 <- difftrteff14_n2mult01$pnexch > 0.7976455
difftrteff14_n2mult01$pnexchgt020 <- difftrteff14_n2mult01$pnexch > 0.2

table(difftrteff14_n2mult01$pnexchgt0798, difftrteff14_n2mult01$plt05)

system.time(difftrteff14_n2mult025 <- runsim_MEM1_binary(nsim, 
                                                        parfunc_twoarm_unequalsamp(probc = c(0.1,0.1), probt = c(0.2,0.5),
                                                                                   effectsizes = c(1,4),
                                                                                   n2mult = 0.25)))
difftrteff14_n2mult025$plt05 <- difftrteff14_n2mult025$pval < 0.05
difftrteff14_n2mult025$pnexchgt090 <- difftrteff14_n2mult025$pnexch > 0.9
difftrteff14_n2mult025$pnexchgt080 <- difftrteff14_n2mult025$pnexch > 0.8
difftrteff14_n2mult025$pnexchgt0798 <- difftrteff14_n2mult025$pnexch > 0.7976455
difftrteff14_n2mult025$pnexchgt020 <- difftrteff14_n2mult025$pnexch > 0.2

table(difftrteff14_n2mult025$pnexchgt0798, difftrteff14_n2mult025$plt05)

#scenario 5: different treatment effects (1 and 5)
system.time(difftrteff15_n2mult01 <- runsim_MEM1_binary(nsim, 
                                                        parfunc_twoarm_unequalsamp(probc = c(0.1,0.1), probt = c(0.2,0.6),
                                                              effectsizes = c(1,5),
                                                              n2mult = 0.1)))
difftrteff15_n2mult01$plt05 <- difftrteff15_n2mult01$pval < 0.05
difftrteff15_n2mult01$pnexchgt090 <- difftrteff15_n2mult01$pnexch > 0.9
difftrteff15_n2mult01$pnexchgt080 <- difftrteff15_n2mult01$pnexch > 0.8
difftrteff15_n2mult01$pnexchgt0798 <- difftrteff15_n2mult01$pnexch > 0.7976455
difftrteff15_n2mult01$pnexchgt020 <- difftrteff15_n2mult01$pnexch > 0.2

table(difftrteff15_n2mult01$pnexchgt0798, difftrteff15_n2mult01$plt05)

system.time(difftrteff15_n2mult025 <- runsim_MEM1_binary(nsim, 
                                                        parfunc_twoarm_unequalsamp(probc = c(0.1,0.1), probt = c(0.2,0.6),
                                                                                   effectsizes = c(1,5),
                                                                                   n2mult = 0.25)))
difftrteff15_n2mult025$plt05 <- difftrteff15_n2mult025$pval < 0.05
difftrteff15_n2mult025$pnexchgt090 <- difftrteff15_n2mult025$pnexch > 0.9
difftrteff15_n2mult025$pnexchgt080 <- difftrteff15_n2mult025$pnexch > 0.8
difftrteff15_n2mult025$pnexchgt0798 <- difftrteff15_n2mult025$pnexch > 0.7976455
difftrteff15_n2mult025$pnexchgt020 <- difftrteff15_n2mult025$pnexch > 0.2

table(difftrteff15_n2mult025$pnexchgt0798, difftrteff15_n2mult025$plt05)

#scenario 6: different treatment effects (1 and 6)
system.time(difftrteff16_n2mult01 <- runsim_MEM1_binary(nsim, 
                                                        parfunc_twoarm_unequalsamp(probc = c(0.1,0.1), probt = c(0.2,0.7),
                                                              effectsizes = c(1,6),
                                                              n2mult = 0.1)))
difftrteff16_n2mult01$plt05 <- difftrteff16_n2mult01$pval < 0.05
difftrteff16_n2mult01$pnexchgt090 <- difftrteff16_n2mult01$pnexch > 0.9
difftrteff16_n2mult01$pnexchgt080 <- difftrteff16_n2mult01$pnexch > 0.8
difftrteff16_n2mult01$pnexchgt0798 <- difftrteff16_n2mult01$pnexch > 0.7976455
difftrteff16_n2mult01$pnexchgt020 <- difftrteff16_n2mult01$pnexch > 0.2

table(difftrteff16_n2mult01$pnexchgt0798, difftrteff16_n2mult01$plt05)

system.time(difftrteff16_n2mult025 <- runsim_MEM1_binary(nsim, 
                                                        parfunc_twoarm_unequalsamp(probc = c(0.1,0.1), probt = c(0.2,0.7),
                                                                                   effectsizes = c(1,6),
                                                                                   n2mult = 0.25)))
difftrteff16_n2mult025$plt05 <- difftrteff16_n2mult025$pval < 0.05
difftrteff16_n2mult025$pnexchgt090 <- difftrteff16_n2mult025$pnexch > 0.9
difftrteff16_n2mult025$pnexchgt080 <- difftrteff16_n2mult025$pnexch > 0.8
difftrteff16_n2mult025$pnexchgt0798 <- difftrteff16_n2mult025$pnexch > 0.7976455
difftrteff16_n2mult025$pnexchgt020 <- difftrteff16_n2mult025$pnexch > 0.2

table(difftrteff16_n2mult025$pnexchgt0798, difftrteff16_n2mult025$plt05)

#combine data sets
equaltrteff11_n2mult01$n2mult <- rep(0.1, length(equaltrteff11_n2mult01$pnexch))

difftrteff12_n2mult01$n2mult <- rep(0.1, length(difftrteff12_n2mult01$pnexch))
difftrteff13_n2mult01$n2mult <- rep(0.1, length(difftrteff13_n2mult01$pnexch))
difftrteff14_n2mult01$n2mult <- rep(0.1, length(difftrteff14_n2mult01$pnexch))
difftrteff15_n2mult01$n2mult <- rep(0.1, length(difftrteff15_n2mult01$pnexch))
difftrteff16_n2mult01$n2mult <- rep(0.1, length(difftrteff16_n2mult01$pnexch))

equaltrteff11_n2mult025$n2mult <- rep(0.25, length(equaltrteff11_n2mult025$pnexch))

difftrteff12_n2mult025$n2mult <- rep(0.25, length(difftrteff12_n2mult025$pnexch))
difftrteff13_n2mult025$n2mult <- rep(0.25, length(difftrteff13_n2mult025$pnexch))
difftrteff14_n2mult025$n2mult <- rep(0.25, length(difftrteff14_n2mult025$pnexch))
difftrteff15_n2mult025$n2mult <- rep(0.25, length(difftrteff15_n2mult025$pnexch))
difftrteff16_n2mult025$n2mult <- rep(0.25, length(difftrteff16_n2mult025$pnexch))

allresults_unequalsamp <- rbind(equaltrteff11_n2mult01, equaltrteff11_n2mult025,
                    difftrteff12_n2mult01, difftrteff12_n2mult025,
                    difftrteff13_n2mult01, difftrteff13_n2mult025,
                    difftrteff14_n2mult01, difftrteff14_n2mult025,
                    difftrteff15_n2mult01, difftrteff15_n2mult025,
                    difftrteff16_n2mult01, difftrteff16_n2mult025)

all_summ0_unequalsamp <- allresults_unequalsamp %>% group_by(pc1,pc2,pt1,pt2, n2mult, totsampsize)

all_summ_unequalsamp <- all_summ0_unequalsamp %>% summarize(plt05 = mean(plt05), 
                                                            pnexchgt090 = mean(pnexchgt090), 
                                                            pnexchgt080 = mean(pnexchgt080), 
                                                            pnexchgt0798 = mean(pnexchgt0798), 
                                                            pnexchgt020 = mean(pnexchgt020))

write.csv(all_summ_unequalsamp, "C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/resultsdata/TA_binary_MEMr_result_unequalsamp_summary.csv")
write.csv(allresults_unequalsamp, "C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/resultsdata/TA_binary_MEMr_result_unequalsamp_raw.csv")
