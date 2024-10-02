library(tidyverse)
library(snowfall)

#interaction lm only functional for two sources/groups
source("C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/AlesandLMSimFunc_MEM1Only.R")

#function to run simulation in parallel and output results nicely
runsim_MEM1 <- function(nsim, pars) {
  
  #summary matrix: first weight 1 (non-exchangeable) and second is weight 2 (exchangeable)
  summ <- data.frame(totsampsize = rep(NA, length(pars)*nsim),
                     source1te = rep(NA, length(pars)*nsim),
                     source2te = rep(NA, length(pars)*nsim),
                     pnexch = rep(NA, length(pars)*nsim), 
                     pexch = rep(NA, length(pars)*nsim), 
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
    sfClusterSetupRNG(seed=12345)
    
    sfExportAll()
    sim.results <- t(sfLapply(rep(1,nsim),
                              function(x) doSim_MEM1(n=par$n, mu0=par$mu0, p=par$p,
                                                beta_t=par$beta_t, 
                                                sd=par$sd, 
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
parfunc_twoarm <- function(baselinemeans, effectsizes, sds){
  pars=list(
    #n = sample sizes, m0 = baseline mean for each source, 
    #p = probability of treatment being assigned, 
    #beta_t = treatment main effect 
    list(n=c(100, 100), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(150, 150), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(200, 200), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(250, 250), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(300, 300), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(350, 350), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(400, 400), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(450, 450), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(500, 500), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC")
  )
  
  return(pars)
}


#set number of samples per setting
nsim=1000


######SD = 1

#scenario 1: equal treatment effects (1 and 1), sd = 1
system.time(t1_equaltrteff11_var1 <- runsim_MEM1(nsim, 
                                              parfunc_twoarm(baselinemeans = c(5,5),
                                                               effectsizes = c(1,1),
                                                               sds = c(1,1))))

t1_equaltrteff11_var1$plt05 <- t1_equaltrteff11_var1$pval < 0.05
t1_equaltrteff11_var1$pnexchgt080 <- t1_equaltrteff11_var1$pnexch > 0.7976455
t1_equaltrteff11_var1$pnexchgt020 <- t1_equaltrteff11_var1$pnexch > 0.2
t1_equaltrteff11_var1$var <- 1

table(t1_equaltrteff11_var1$pnexchgt080, t1_equaltrteff11_var1$plt05)



#scenario 2: similar treatment effects (1 and 2), sd = 1
system.time(difftrteff12_var1 <- runsim_MEM1(nsim, 
                                             parfunc_twoarm(baselinemeans = c(5,5),
                                                           effectsizes = c(1,2),
                                                           sds = c(1,1))))
difftrteff12_var1$plt05 <- difftrteff12_var1$pval < 0.05
difftrteff12_var1$pnexchgt080 <- difftrteff12_var1$pnexch > 0.7976455
difftrteff12_var1$pnexchgt020 <- difftrteff12_var1$pnexch > 0.2
difftrteff12_var1$var <- 1

table(difftrteff12_var1$pnexchgt080, difftrteff12_var1$plt05)

#scenario 3: different treatment effects (1 and 3), sd = 1
system.time(difftrteff13_var1 <- runsim_MEM1(nsim, 
                                             parfunc_twoarm(baselinemeans = c(5,5),
                                                           effectsizes = c(1,3),
                                                           sds = c(1,1))))
difftrteff13_var1$plt05 <- difftrteff13_var1$pval < 0.05
difftrteff13_var1$pnexchgt080 <- difftrteff13_var1$pnexch > 0.7976455
difftrteff13_var1$pnexchgt020 <- difftrteff13_var1$pnexch > 0.2
difftrteff13_var1$var <- 1

table(difftrteff13_var1$pnexchgt080, difftrteff13_var1$plt05)

#scenario 4: different treatment effects (1 and 4), sd = 1
system.time(difftrteff14_var1 <- runsim_MEM1(nsim, 
                                             parfunc_twoarm(baselinemeans = c(5,5),
                                                           effectsizes = c(1,4),
                                                           sds = c(1,1))))
difftrteff14_var1$plt05 <- difftrteff14_var1$pval < 0.05
difftrteff14_var1$pnexchgt080 <- difftrteff14_var1$pnexch > 0.7976455
difftrteff14_var1$pnexchgt020 <- difftrteff14_var1$pnexch > 0.2
difftrteff14_var1$var <- 1

table(difftrteff14_var1$pnexchgt080, difftrteff14_var1$plt05)

#scenario 5: different treatment effects (1 and 5), sd = 1
system.time(difftrteff15_var1 <- runsim_MEM1(nsim, 
                                             parfunc_twoarm(baselinemeans = c(5,5),
                                                           effectsizes = c(1,5),
                                                           sds = c(1,1))))
difftrteff15_var1$plt05 <- difftrteff15_var1$pval < 0.05
difftrteff15_var1$pnexchgt080 <- difftrteff15_var1$pnexch > 0.7976455
difftrteff15_var1$pnexchgt020 <- difftrteff15_var1$pnexch > 0.2
difftrteff15_var1$var <- 1

table(difftrteff15_var1$pnexchgt080, difftrteff15_var1$plt05)

#scenario 6: different treatment effects (1 and 6), sd = 1
system.time(difftrteff16_var1 <- runsim_MEM1(nsim, 
                                             parfunc_twoarm(baselinemeans = c(5,5),
                                                           effectsizes = c(1,6),
                                                           sds = c(1,1))))
difftrteff16_var1$plt05 <- difftrteff16_var1$pval < 0.05
difftrteff16_var1$pnexchgt080 <- difftrteff16_var1$pnexch > 0.7976455
difftrteff16_var1$pnexchgt020 <- difftrteff16_var1$pnexch > 0.2
difftrteff16_var1$var <- 1

table(difftrteff16_var1$pnexchgt080, difftrteff16_var1$plt05)




######SD = sqrt(10)

#scenario 7: equal treatment effects (1 and 1), sd = sqrt(10)
system.time(equaltrteff11_var10 <- runsim_MEM1(nsim, 
                                               parfunc_twoarm(baselinemeans = c(5,5),
                                                             effectsizes = c(1,1),
                                                             sds = c(sqrt(10),sqrt(10)))))
equaltrteff11_var10$plt05 <- equaltrteff11_var10$pval < 0.05
equaltrteff11_var10$pnexchgt080 <- equaltrteff11_var10$pnexch > 0.7976455
equaltrteff11_var10$pnexchgt020 <- equaltrteff11_var10$pnexch > 0.2
equaltrteff11_var10$var <- 10

table(equaltrteff11_var10$pnexchgt080, equaltrteff11_var10$plt05)

#scenario 8: different treatment effects (1 and 2), sd = sqrt(10)
system.time(difftrteff12_var10 <- runsim_MEM1(nsim,
                                              parfunc_twoarm(baselinemeans = c(5,5),
                                                            effectsizes = c(1,2),
                                                            sds = c(sqrt(10),sqrt(10)))))
difftrteff12_var10$plt05 <- difftrteff12_var10$pval < 0.05
difftrteff12_var10$pnexchgt080 <- difftrteff12_var10$pnexch > 0.7976455
difftrteff12_var10$pnexchgt020 <- difftrteff12_var10$pnexch > 0.2
difftrteff12_var10$var <- 10

table(difftrteff12_var10$pnexchgt080, difftrteff12_var10$plt05)


#scenario 9: different treatment effects (1 and 3), sd = sqrt(10)
system.time(difftrteff13_var10 <- runsim_MEM1(nsim, 
                                              parfunc_twoarm(baselinemeans = c(5,5),
                                                            effectsizes = c(1,3),
                                                            sds = c(sqrt(10),sqrt(10)))))
difftrteff13_var10$plt05 <- difftrteff13_var10$pval < 0.05
difftrteff13_var10$pnexchgt080 <- difftrteff13_var10$pnexch > 0.7976455
difftrteff13_var10$pnexchgt020 <- difftrteff13_var10$pnexch > 0.2
difftrteff13_var10$var <- 10

table(difftrteff13_var10$pnexchgt080, difftrteff13_var10$plt05)


#scenario 10: different treatment effects (1 and 4), sd = sqrt(10)
system.time(difftrteff14_var10 <- runsim_MEM1(nsim, 
                                              parfunc_twoarm(baselinemeans = c(5,5),
                                                            effectsizes = c(1,4),
                                                            sds = c(sqrt(10),sqrt(10)))))
difftrteff14_var10$plt05 <- difftrteff14_var10$pval < 0.05
difftrteff14_var10$pnexchgt080 <- difftrteff14_var10$pnexch > 0.7976455
difftrteff14_var10$pnexchgt020 <- difftrteff14_var10$pnexch > 0.2
difftrteff14_var10$var <- 10

table(difftrteff14_var10$pnexchgt080, difftrteff14_var10$plt05)


#scenario 11: different treatment effects (1 and 5), sd = sqrt(10)
system.time(difftrteff15_var10 <- runsim_MEM1(nsim, 
                                              parfunc_twoarm(baselinemeans = c(5,5),
                                                            effectsizes = c(1,5),
                                                            sds = c(sqrt(10),sqrt(10)))))
difftrteff15_var10$plt05 <- difftrteff15_var10$pval < 0.05
difftrteff15_var10$pnexchgt080 <- difftrteff15_var10$pnexch > 0.7976455
difftrteff15_var10$pnexchgt020 <- difftrteff15_var10$pnexch > 0.2
difftrteff15_var10$var <- 10

table(difftrteff15_var10$pnexchgt080, difftrteff15_var10$plt05)


#scenario 12: different treatment effects (1 and 6), sd = sqrt(10)
system.time(difftrteff16_var10 <- runsim_MEM1(nsim, 
                                              parfunc_twoarm(baselinemeans = c(5,5),
                                                            effectsizes = c(1,6),
                                                            sds = c(sqrt(10),sqrt(10)))))
difftrteff16_var10$plt05 <- difftrteff16_var10$pval < 0.05
difftrteff16_var10$pnexchgt080 <- difftrteff16_var10$pnexch > 0.7976455
difftrteff16_var10$pnexchgt020 <- difftrteff16_var10$pnexch > 0.2
difftrteff16_var10$var <- 10

table(difftrteff16_var10$pnexchgt080, difftrteff16_var10$plt05)




######SD = 10

#scenario 13: equal treatment effects (1 and 1), sd = 10
system.time(equaltrteff11_var100 <- runsim_MEM1(nsim, 
                                                parfunc_twoarm(baselinemeans = c(5,5),
                                                              effectsizes = c(1,1),
                                                              sds = c(10,10))))
equaltrteff11_var100$plt05 <- equaltrteff11_var100$pval < 0.05
equaltrteff11_var100$pnexchgt080 <- equaltrteff11_var100$pnexch > 0.7976455
equaltrteff11_var100$pnexchgt020 <- equaltrteff11_var100$pnexch > 0.2
equaltrteff11_var100$var <- 100

table(equaltrteff11_var100$pnexchgt080, equaltrteff11_var100$plt05)

#scenario 14: different treatment effects (1 and 2), sd = 10
system.time(difftrteff12_var100 <- runsim_MEM1(nsim, 
                                               parfunc_twoarm(baselinemeans = c(5,5),
                                                             effectsizes = c(1,2),
                                                             sds = c(10,10))))
difftrteff12_var100$plt05 <- difftrteff12_var100$pval < 0.05
difftrteff12_var100$pnexchgt080 <- difftrteff12_var100$pnexch > 0.7976455
difftrteff12_var100$pnexchgt020 <- difftrteff12_var100$pnexch > 0.2
difftrteff12_var100$var <- 100

table(difftrteff12_var100$pnexchgt080, difftrteff12_var100$plt05)

#scenario 15: different treatment effects (1 and 3), sd = 10
system.time(difftrteff13_var100 <- runsim_MEM1(nsim, 
                                               parfunc_twoarm(baselinemeans = c(5,5),
                                                             effectsizes = c(1,3),
                                                             sds = c(10,10))))
difftrteff13_var100$plt05 <- difftrteff13_var100$pval < 0.05
difftrteff13_var100$pnexchgt080 <- difftrteff13_var100$pnexch > 0.7976455
difftrteff13_var100$pnexchgt020 <- difftrteff13_var100$pnexch > 0.2
difftrteff13_var100$var <- 100

table(difftrteff13_var100$pnexchgt080, difftrteff13_var100$plt05)


#scenario 16: different treatment effects (1 and 4), sd = 10
system.time(difftrteff14_var100 <- runsim_MEM1(nsim, 
                                               parfunc_twoarm(baselinemeans = c(5,5),
                                                             effectsizes = c(1,4),
                                                             sds = c(10,10))))
difftrteff14_var100$plt05 <- difftrteff14_var100$pval < 0.05
difftrteff14_var100$pnexchgt080 <- difftrteff14_var100$pnexch > 0.7976455
difftrteff14_var100$pnexchgt020 <- difftrteff14_var100$pnexch > 0.2
difftrteff14_var100$var <- 100

table(difftrteff14_var100$pnexchgt080, difftrteff14_var100$plt05)


#scenario 17: different treatment effects (1 and 5), sd = 10
system.time(difftrteff15_var100 <- runsim_MEM1(nsim, 
                                               parfunc_twoarm(baselinemeans = c(5,5),
                                                             effectsizes = c(1,5),
                                                             sds = c(10,10))))
difftrteff15_var100$plt05 <- difftrteff15_var100$pval < 0.05
difftrteff15_var100$pnexchgt080 <- difftrteff15_var100$pnexch > 0.7976455
difftrteff15_var100$pnexchgt020 <- difftrteff15_var100$pnexch > 0.2
difftrteff15_var100$var <- 100

table(difftrteff15_var100$pnexchgt080, difftrteff15_var100$plt05)


#scenario 18: different treatment effects (1 and 6), sd = 10
system.time(difftrteff16_var100 <- runsim_MEM1(nsim, 
                                               parfunc_twoarm(baselinemeans = c(5,5),
                                                             effectsizes = c(1,6),
                                                             sds = c(10,10))))
difftrteff16_var100$plt05 <- difftrteff16_var100$pval < 0.05
difftrteff16_var100$pnexchgt080 <- difftrteff16_var100$pnexch > 0.7976455
difftrteff16_var100$pnexchgt020 <- difftrteff16_var100$pnexch > 0.2
difftrteff16_var100$var <- 100

table(difftrteff16_var100$pnexchgt080, difftrteff16_var100$plt05)




#combine data sets
allresults <- rbind(t1_equaltrteff11_var1, equaltrteff11_var10, equaltrteff11_var100,
                    difftrteff12_var1, difftrteff13_var1, difftrteff14_var1, difftrteff15_var1, difftrteff16_var1,
                    difftrteff12_var10, difftrteff13_var10, difftrteff14_var10, difftrteff15_var10, difftrteff16_var10,
                    difftrteff12_var100, difftrteff13_var100, difftrteff14_var100, difftrteff15_var100, difftrteff16_var100)

all_summ0 <- allresults %>% group_by(source1te, source2te, var, totsampsize)

all_summ <- all_summ0 %>% summarize(plt05 = mean(plt05), pnexchgt080 = mean(pnexchgt080), pnexchgt020 = mean(pnexchgt020))

write.csv(all_summ, "C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/resultsdata/results.csv")
write.csv(allresults, "C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/resultsdata/results_raw.csv")




#############################################
#############UNEQUAL SAMPLE##################
#############################################

#function to format inputs 
parfunc_twoarm_unequalsamp <- function(baselinemeans, effectsizes, sds, n2mult){
  pars=list(
    #n = sample sizes, m0 = baseline mean for each source, 
    #p = probability of treatment being assigned, 
    #beta_t = treatment main effect 
    list(n=c(200*(1-n2mult), 200*n2mult), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(300*(1-n2mult), 300*n2mult), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(400*(1-n2mult), 400*n2mult), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(500*(1-n2mult), 500*n2mult), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(600*(1-n2mult), 600*n2mult), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(700*(1-n2mult), 700*n2mult), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(800*(1-n2mult), 800*n2mult), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(900*(1-n2mult), 900*n2mult), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(1000*(1-n2mult), 1000*n2mult), mu0=baselinemeans,   p=c(0.5, 0.5),
         beta_t=effectsizes, sd=sds, BICiterations=1000, 
         marginal="BIC")
  )
  
  return(pars)
}


#set number of samples per setting
nsim=1000

######SD = sqrt(10)

#scenario 1: equal treatment effects (1 and 1), sd = sqrt(10)
system.time(equaltrteff11_var10_n2mult01 <- runsim_MEM1(nsim, 
                                               parfunc_twoarm_unequalsamp(baselinemeans = c(5,5),
                                                              effectsizes = c(1,1),
                                                              sds = c(sqrt(10),sqrt(10)),
                                                              n2mult = 0.1)))

equaltrteff11_var10_n2mult01$plt05 <- equaltrteff11_var10_n2mult01$pval < 0.05
equaltrteff11_var10_n2mult01$pnexchgt080 <- equaltrteff11_var10_n2mult01$pnexch > 0.7976455
equaltrteff11_var10_n2mult01$pnexchgt020 <- equaltrteff11_var10_n2mult01$pnexch > 0.2
equaltrteff11_var10_n2mult01$var <- 10
equaltrteff11_var10_n2mult01$n2mult <- 0.1

table(equaltrteff11_var10_n2mult01$pnexchgt080, equaltrteff11_var10_n2mult01$plt05)

system.time(equaltrteff11_var10_n2mult025 <- runsim_MEM1(nsim, 
                                                        parfunc_twoarm_unequalsamp(baselinemeans = c(5,5),
                                                                                   effectsizes = c(1,1),
                                                                                   sds = c(sqrt(10),sqrt(10)),
                                                                                   n2mult = 0.25)))

equaltrteff11_var10_n2mult025$plt05 <- equaltrteff11_var10_n2mult025$pval < 0.05
equaltrteff11_var10_n2mult025$pnexchgt080 <- equaltrteff11_var10_n2mult025$pnexch > 0.7976455
equaltrteff11_var10_n2mult025$pnexchgt020 <- equaltrteff11_var10_n2mult025$pnexch > 0.2
equaltrteff11_var10_n2mult025$var <- 10
equaltrteff11_var10_n2mult025$n2mult <- 0.25

table(equaltrteff11_var10_n2mult025$pnexchgt080, equaltrteff11_var10_n2mult025$plt05)

#scenario 2: different treatment effects (1 and 2), sd = sqrt(10)
system.time(difftrteff12_var10_n2mult01 <- runsim_MEM1(nsim,
                                              parfunc_twoarm_unequalsamp(baselinemeans = c(5,5),
                                                             effectsizes = c(1,2),
                                                             sds = c(sqrt(10),sqrt(10)),
                                                             n2mult = 0.1)))

difftrteff12_var10_n2mult01$plt05 <- difftrteff12_var10_n2mult01$pval < 0.05
difftrteff12_var10_n2mult01$pnexchgt080 <- difftrteff12_var10_n2mult01$pnexch > 0.7976455
difftrteff12_var10_n2mult01$pnexchgt020 <- difftrteff12_var10_n2mult01$pnexch > 0.2
difftrteff12_var10_n2mult01$var <- 10
difftrteff12_var10_n2mult01$n2mult <- 0.1

table(difftrteff12_var10_n2mult01$pnexchgt080, difftrteff12_var10_n2mult01$plt05)

system.time(difftrteff12_var10_n2mult025 <- runsim_MEM1(nsim,
                                                       parfunc_twoarm_unequalsamp(baselinemeans = c(5,5),
                                                                                  effectsizes = c(1,2),
                                                                                  sds = c(sqrt(10),sqrt(10)),
                                                                                  n2mult = 0.25)))

difftrteff12_var10_n2mult025$plt05 <- difftrteff12_var10_n2mult025$pval < 0.05
difftrteff12_var10_n2mult025$pnexchgt080 <- difftrteff12_var10_n2mult025$pnexch > 0.7976455
difftrteff12_var10_n2mult025$pnexchgt020 <- difftrteff12_var10_n2mult025$pnexch > 0.2
difftrteff12_var10_n2mult025$var <- 10
difftrteff12_var10_n2mult025$n2mult <- 0.25

table(difftrteff12_var10_n2mult025$pnexchgt080, difftrteff12_var10_n2mult025$plt05)

#scenario 3: different treatment effects (1 and 3), sd = sqrt(10)
system.time(difftrteff13_var10_n2mult01 <- runsim_MEM1(nsim, 
                                              parfunc_twoarm_unequalsamp(baselinemeans = c(5,5),
                                                             effectsizes = c(1,3),
                                                             sds = c(sqrt(10),sqrt(10)),
                                                             n2mult = 0.1)))

difftrteff13_var10_n2mult01$plt05 <- difftrteff13_var10_n2mult01$pval < 0.05
difftrteff13_var10_n2mult01$pnexchgt080 <- difftrteff13_var10_n2mult01$pnexch > 0.7976455
difftrteff13_var10_n2mult01$pnexchgt020 <- difftrteff13_var10_n2mult01$pnexch > 0.2
difftrteff13_var10_n2mult01$var <- 10
difftrteff13_var10_n2mult01$n2mult <- 0.1

table(difftrteff13_var10_n2mult01$pnexchgt080, difftrteff13_var10_n2mult01$plt05)

system.time(difftrteff13_var10_n2mult025 <- runsim_MEM1(nsim, 
                                                       parfunc_twoarm_unequalsamp(baselinemeans = c(5,5),
                                                                                  effectsizes = c(1,3),
                                                                                  sds = c(sqrt(10),sqrt(10)),
                                                                                  n2mult = 0.25)))

difftrteff13_var10_n2mult025$plt05 <- difftrteff13_var10_n2mult025$pval < 0.05
difftrteff13_var10_n2mult025$pnexchgt080 <- difftrteff13_var10_n2mult025$pnexch > 0.7976455
difftrteff13_var10_n2mult025$pnexchgt020 <- difftrteff13_var10_n2mult025$pnexch > 0.2
difftrteff13_var10_n2mult025$var <- 10
difftrteff13_var10_n2mult025$n2mult <- 0.25

table(difftrteff13_var10_n2mult025$pnexchgt080, difftrteff13_var10_n2mult025$plt05)

#scenario 4: different treatment effects (1 and 4), sd = sqrt(10)
system.time(difftrteff14_var10_n2mult01 <- runsim_MEM1(nsim, 
                                              parfunc_twoarm_unequalsamp(baselinemeans = c(5,5),
                                                             effectsizes = c(1,4),
                                                             sds = c(sqrt(10),sqrt(10)),
                                                             n2mult = 0.1)))

difftrteff14_var10_n2mult01$plt05 <- difftrteff14_var10_n2mult01$pval < 0.05
difftrteff14_var10_n2mult01$pnexchgt080 <- difftrteff14_var10_n2mult01$pnexch > 0.7976455
difftrteff14_var10_n2mult01$pnexchgt020 <- difftrteff14_var10_n2mult01$pnexch > 0.2
difftrteff14_var10_n2mult01$var <- 10
difftrteff14_var10_n2mult01$n2mult <- 0.1

table(difftrteff14_var10_n2mult01$pnexchgt080, difftrteff14_var10_n2mult01$plt05)

system.time(difftrteff14_var10_n2mult025 <- runsim_MEM1(nsim, 
                                                       parfunc_twoarm_unequalsamp(baselinemeans = c(5,5),
                                                                                  effectsizes = c(1,4),
                                                                                  sds = c(sqrt(10),sqrt(10)),
                                                                                  n2mult = 0.25)))

difftrteff14_var10_n2mult025$plt05 <- difftrteff14_var10_n2mult025$pval < 0.05
difftrteff14_var10_n2mult025$pnexchgt080 <- difftrteff14_var10_n2mult025$pnexch > 0.7976455
difftrteff14_var10_n2mult025$pnexchgt020 <- difftrteff14_var10_n2mult025$pnexch > 0.2
difftrteff14_var10_n2mult025$var <- 10
difftrteff14_var10_n2mult025$n2mult <- 0.25

table(difftrteff14_var10_n2mult025$pnexchgt080, difftrteff14_var10_n2mult025$plt05)

#scenario 5: different treatment effects (1 and 5), sd = sqrt(10)
system.time(difftrteff15_var10_n2mult01 <- runsim_MEM1(nsim, 
                                              parfunc_twoarm_unequalsamp(baselinemeans = c(5,5),
                                                             effectsizes = c(1,5),
                                                             sds = c(sqrt(10),sqrt(10)),
                                                             n2mult = 0.1)))

difftrteff15_var10_n2mult01$plt05 <- difftrteff15_var10_n2mult01$pval < 0.05
difftrteff15_var10_n2mult01$pnexchgt080 <- difftrteff15_var10_n2mult01$pnexch > 0.7976455
difftrteff15_var10_n2mult01$pnexchgt020 <- difftrteff15_var10_n2mult01$pnexch > 0.2
difftrteff15_var10_n2mult01$var <- 10
difftrteff15_var10_n2mult01$n2mult <- 0.1

table(difftrteff15_var10_n2mult01$pnexchgt080, difftrteff15_var10_n2mult01$plt05)

system.time(difftrteff15_var10_n2mult025 <- runsim_MEM1(nsim, 
                                                       parfunc_twoarm_unequalsamp(baselinemeans = c(5,5),
                                                                                  effectsizes = c(1,5),
                                                                                  sds = c(sqrt(10),sqrt(10)),
                                                                                  n2mult = 0.25)))

difftrteff15_var10_n2mult025$plt05 <- difftrteff15_var10_n2mult025$pval < 0.05
difftrteff15_var10_n2mult025$pnexchgt080 <- difftrteff15_var10_n2mult025$pnexch > 0.7976455
difftrteff15_var10_n2mult025$pnexchgt020 <- difftrteff15_var10_n2mult025$pnexch > 0.2
difftrteff15_var10_n2mult025$var <- 10
difftrteff15_var10_n2mult025$n2mult <- 0.25

table(difftrteff15_var10_n2mult025$pnexchgt080, difftrteff15_var10_n2mult025$plt05)

#scenario 6: different treatment effects (1 and 6), sd = sqrt(10)
system.time(difftrteff16_var10_n2mult01 <- runsim_MEM1(nsim, 
                                              parfunc_twoarm_unequalsamp(baselinemeans = c(5,5),
                                                             effectsizes = c(1,6),
                                                             sds = c(sqrt(10),sqrt(10)),
                                                             n2mult = 0.1)))

difftrteff16_var10_n2mult01$plt05 <- difftrteff16_var10_n2mult01$pval < 0.05
difftrteff16_var10_n2mult01$pnexchgt080 <- difftrteff16_var10_n2mult01$pnexch > 0.7976455
difftrteff16_var10_n2mult01$pnexchgt020 <- difftrteff16_var10_n2mult01$pnexch > 0.2
difftrteff16_var10_n2mult01$var <- 10
difftrteff16_var10_n2mult01$n2mult <- 0.1

table(difftrteff16_var10_n2mult01$pnexchgt080, difftrteff16_var10_n2mult01$plt05)

system.time(difftrteff16_var10_n2mult025 <- runsim_MEM1(nsim, 
                                                       parfunc_twoarm_unequalsamp(baselinemeans = c(5,5),
                                                                                  effectsizes = c(1,6),
                                                                                  sds = c(sqrt(10),sqrt(10)),
                                                                                  n2mult = 0.25)))

difftrteff16_var10_n2mult025$plt05 <- difftrteff16_var10_n2mult025$pval < 0.05
difftrteff16_var10_n2mult025$pnexchgt080 <- difftrteff16_var10_n2mult025$pnexch > 0.7976455
difftrteff16_var10_n2mult025$pnexchgt020 <- difftrteff16_var10_n2mult025$pnexch > 0.2
difftrteff16_var10_n2mult025$var <- 10
difftrteff16_var10_n2mult025$n2mult <- 0.25

table(difftrteff16_var10_n2mult025$pnexchgt080, difftrteff16_var10_n2mult025$plt05)



#combine data sets
allresults_unequalsamp <- rbind(equaltrteff11_var10_n2mult01,equaltrteff11_var10_n2mult025,
                    difftrteff12_var10_n2mult01, difftrteff12_var10_n2mult025, 
                    difftrteff13_var10_n2mult01, difftrteff13_var10_n2mult025,
                    difftrteff14_var10_n2mult01, difftrteff14_var10_n2mult025,
                    difftrteff15_var10_n2mult01, difftrteff15_var10_n2mult025,
                    difftrteff16_var10_n2mult01, difftrteff16_var10_n2mult025)

all_summ0_unequalsamp <- allresults_unequalsamp %>% group_by(source1te, source2te, n2mult, var, totsampsize)

all_summ_unequalsamp <- all_summ0_unequalsamp %>% summarize(plt05 = mean(plt05), pnexchgt080 = mean(pnexchgt080), pnexchgt020 = mean(pnexchgt020))

write.csv(all_summ_unequalsamp, "C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/resultsdata/TA_cont_MEMr_result_unequalsamp_summary.csv")
write.csv(allresults_unequalsamp, "C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/resultsdata/TA_cont_MEMr_result_unequalsamp_raw.csv")
