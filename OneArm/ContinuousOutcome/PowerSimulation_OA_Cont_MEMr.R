library(tidyverse)
library(snowfall)

#interaction lm only functional for two sources/groups
source("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/AlesandLMSimFunc_OneArm.R")

#function to run simulation in parallel and output results nicely
runsim_OneArm <- function(nsim, pars) {
  
  #summary matrix: first weight 1 (non-exchangeable) and second is weight 2 (exchangeable)
  summ <- data.frame(totsampsize = rep(NA, length(pars)*nsim),
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
    sfClusterSetupRNG(seed=1234)
    
    sfExportAll()

    #set.seed(1023) #NOT REPRODUCIBLE
    sim.results <- t(sfLapply(rep(1,nsim),
                              function(x) doSim_MEM1_OneArm(n=par$n, mu0=par$mu0,
                                                            sd=par$sd, 
                                                            BICiterations=par$BICiterations, 
                                                            marginal=par$marginal)))
    sfStop()
    
    result <- as.data.frame(unlist(sim.results))
    
    # USE THIS ONCE PVAL FOR LM IS INCLUDED
    summ[(nsim*(i-1) + 1):(nsim*i),1] <- (result %>% filter(row_number() %% 4 == 1))
    summ[(nsim*(i-1) + 1):(nsim*i),2] <- (result %>% filter(row_number() %% 4 == 2))
    summ[(nsim*(i-1) + 1):(nsim*i),3] <- (result %>% filter(row_number() %% 4 == 3))
    summ[(nsim*(i-1) + 1):(nsim*i),4] <- (result %>% filter(row_number() %% 4 == 0))

    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(summ)
  
}



#############################################
##############EQUAL SAMPLE###################
#############################################

#function to format inputs 
parfunc <- function(effectsizes, sds){
  pars=list(
    #n = sample sizes, m0 = mean for each source, 
    list(n=c(100, 100), mu0=effectsizes,
         sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(150, 150), mu0=effectsizes,
         sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(200, 200), mu0=effectsizes,
         sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(250, 250), mu0=effectsizes,
         sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(300, 300), mu0=effectsizes,
         sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(350, 350), mu0=effectsizes,
         sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(400, 400), mu0=effectsizes,
         sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(450, 450), mu0=effectsizes,
         sd=sds, BICiterations=1000, 
         marginal="BIC"),
    list(n=c(500, 500), mu0=effectsizes,
         sd=sds, BICiterations=1000, 
         marginal="BIC")
  )
  
  return(pars)
}

nsim=1000

#scenario 1: equal source effects (1 and 1), sd = 1
system.time(equalsourceeff11_var1 <- runsim_OneArm(nsim, parfunc(effectsizes=c(1,1),sds=c(1,1))))
equalsourceeff11_var1$var <- 1
equalsourceeff11_var1$efdiff <- 1-1
equalsourceeff11_var1$plt05 <- equalsourceeff11_var1$pval < 0.05
equalsourceeff11_var1$pnexchlt080 <- equalsourceeff11_var1$pnexch > 0.7976455
equalsourceeff11_var1$pnexchlt020 <- equalsourceeff11_var1$pnexch > 0.20

table(equalsourceeff11_var1$pnexchlt080, equalsourceeff11_var1$plt05)
table(equalsourceeff11_var1$pnexchlt020, equalsourceeff11_var1$plt05)

#scenario 2: different source effects (1 and 2), sd = 1
system.time(diffsourceeff12_var1 <- runsim_OneArm(nsim, parfunc(c(1,2),c(1,1))))
diffsourceeff12_var1$var <- 1
diffsourceeff12_var1$efdiff <- 2-1
diffsourceeff12_var1$plt05 <- diffsourceeff12_var1$pval < 0.05
diffsourceeff12_var1$pnexchlt080 <- diffsourceeff12_var1$pnexch > 0.7976455
diffsourceeff12_var1$pnexchlt020 <- diffsourceeff12_var1$pnexch > 0.20

table(diffsourceeff12_var1$pnexchlt080, diffsourceeff12_var1$plt05)
table(diffsourceeff12_var1$pnexchlt020, diffsourceeff12_var1$plt05)

#scenario 3: different source effects (1 and 3), sd = 1
system.time(diffsourceeff13_var1 <- runsim_OneArm(nsim, parfunc(c(1,3),c(1,1))))
diffsourceeff13_var1$var <- 1
diffsourceeff13_var1$efdiff <- 3-1
diffsourceeff13_var1$plt05 <- diffsourceeff13_var1$pval < 0.05
diffsourceeff13_var1$pnexchlt080 <- diffsourceeff13_var1$pnexch > 0.7976455
diffsourceeff13_var1$pnexchlt020 <- diffsourceeff13_var1$pnexch > 0.20

table(diffsourceeff13_var1$pnexchlt080, diffsourceeff13_var1$plt05)
table(diffsourceeff13_var1$pnexchlt020, diffsourceeff13_var1$plt05)

#scenario 4: different source effects (1 and 4), sd = 1
system.time(diffsourceeff14_var1 <- runsim_OneArm(nsim, parfunc(c(1,4),c(1,1))))
diffsourceeff14_var1$var <- 1
diffsourceeff14_var1$efdiff <- 4-1
diffsourceeff14_var1$plt05 <- diffsourceeff14_var1$pval < 0.05
diffsourceeff14_var1$pnexchlt080 <- diffsourceeff14_var1$pnexch > 0.7976455
diffsourceeff14_var1$pnexchlt020 <- diffsourceeff14_var1$pnexch > 0.20

table(diffsourceeff14_var1$pnexchlt080, diffsourceeff14_var1$plt05)
table(diffsourceeff14_var1$pnexchlt020, diffsourceeff14_var1$plt05)

#scenario 5: different source effects (1 and 5), sd = 1
system.time(diffsourceeff15_var1 <- runsim_OneArm(nsim, parfunc(c(1,5),c(1,1))))
diffsourceeff15_var1$var <- 1
diffsourceeff15_var1$efdiff <- 5-1
diffsourceeff15_var1$plt05 <- diffsourceeff15_var1$pval < 0.05
diffsourceeff15_var1$pnexchlt080 <- diffsourceeff15_var1$pnexch > 0.7976455
diffsourceeff15_var1$pnexchlt020 <- diffsourceeff15_var1$pnexch > 0.20

table(diffsourceeff15_var1$pnexchlt080, diffsourceeff15_var1$plt05)
table(diffsourceeff15_var1$pnexchlt020, diffsourceeff15_var1$plt05)





#scenario 6: equal source effects (1 and 1), sd = sqrt(10)
system.time(equalsourceeff11_var10 <- runsim_OneArm(nsim, parfunc(c(1,1),c(sqrt(10),sqrt(10)))))
equalsourceeff11_var10$var <- 10
equalsourceeff11_var10$efdiff <- 1-1
equalsourceeff11_var10$plt05 <- equalsourceeff11_var10$pval < 0.05
equalsourceeff11_var10$pnexchlt080 <- equalsourceeff11_var10$pnexch > 0.7976455
equalsourceeff11_var10$pnexchlt020 <- equalsourceeff11_var10$pnexch > 0.20

table(equalsourceeff11_var10$pnexchlt080, equalsourceeff11_var10$plt05)
table(equalsourceeff11_var10$pnexchlt020, equalsourceeff11_var10$plt05)

#scenario 7: different source effects (1 and 2), sd = sqrt(10)
system.time(diffsourceeff12_var10 <- runsim_OneArm(nsim, parfunc(c(1,2),c(sqrt(10),sqrt(10)))))
diffsourceeff12_var10$var <- 10
diffsourceeff12_var10$efdiff <- 2-1
diffsourceeff12_var10$plt05 <- diffsourceeff12_var10$pval < 0.05
diffsourceeff12_var10$pnexchlt080 <- diffsourceeff12_var10$pnexch > 0.7976455
diffsourceeff12_var10$pnexchlt020 <- diffsourceeff12_var10$pnexch > 0.20

table(diffsourceeff12_var10$pnexchlt080, diffsourceeff12_var10$plt05)
table(diffsourceeff12_var10$pnexchlt020, diffsourceeff12_var10$plt05)

#scenario 8: different source effects (1 and 3), sd = sqrt(10)
system.time(diffsourceeff13_var10 <- runsim_OneArm(nsim, parfunc(c(1,3),c(sqrt(10),sqrt(10)))))
diffsourceeff13_var10$var <- 10
diffsourceeff13_var10$efdiff <- 3-1
diffsourceeff13_var10$plt05 <- diffsourceeff13_var10$pval < 0.05
diffsourceeff13_var10$pnexchlt080 <- diffsourceeff13_var10$pnexch > 0.7976455
diffsourceeff13_var10$pnexchlt020 <- diffsourceeff13_var10$pnexch > 0.20

table(diffsourceeff13_var10$pnexchlt080, diffsourceeff13_var10$plt05)
table(diffsourceeff13_var10$pnexchlt020, diffsourceeff13_var10$plt05)

#scenario 9: different source effects (1 and 4), sd = sqrt(10)
system.time(diffsourceeff14_var10 <- runsim_OneArm(nsim, parfunc(c(1,4),c(sqrt(10),sqrt(10)))))
diffsourceeff14_var10$var <- 10
diffsourceeff14_var10$efdiff <- 4-1
diffsourceeff14_var10$plt05 <- diffsourceeff14_var10$pval < 0.05
diffsourceeff14_var10$pnexchlt080 <- diffsourceeff14_var10$pnexch > 0.7976455
diffsourceeff14_var10$pnexchlt020 <- diffsourceeff14_var10$pnexch > 0.20

table(diffsourceeff14_var10$pnexchlt080, diffsourceeff14_var10$plt05)
table(diffsourceeff14_var10$pnexchlt020, diffsourceeff14_var10$plt05)

#scenario 10: different source effects (1 and 5), sd = sqrt(10)
system.time(diffsourceeff15_var10 <- runsim_OneArm(nsim, parfunc(c(1,5),c(sqrt(10),sqrt(10)))))
diffsourceeff15_var10$var <- 10
diffsourceeff15_var10$efdiff <- 5-1
diffsourceeff15_var10$plt05 <- diffsourceeff15_var10$pval < 0.05
diffsourceeff15_var10$pnexchlt080 <- diffsourceeff15_var10$pnexch > 0.7976455
diffsourceeff15_var10$pnexchlt020 <- diffsourceeff15_var10$pnexch > 0.20

table(diffsourceeff15_var10$pnexchlt080, diffsourceeff15_var10$plt05)
table(diffsourceeff15_var10$pnexchlt020, diffsourceeff15_var10$plt05)





#scenario 11: equal source effects (1 and 1), sd = 10
system.time(equalsourceeff11_var100 <- runsim_OneArm(nsim, parfunc(c(1,1),c(10,10))))
equalsourceeff11_var100$var <- 100
equalsourceeff11_var100$efdiff <- 1-1
equalsourceeff11_var100$plt05 <- equalsourceeff11_var100$pval < 0.05
equalsourceeff11_var100$pnexchlt080 <- equalsourceeff11_var100$pnexch > 0.7976455
equalsourceeff11_var100$pnexchlt020 <- equalsourceeff11_var100$pnexch > 0.20

table(equalsourceeff11_var100$pnexchlt080, equalsourceeff11_var100$plt05)
table(equalsourceeff11_var100$pnexchlt020, equalsourceeff11_var100$plt05)

#scenario 12: different source effects (1 and 2), sd = 10
system.time(diffsourceeff12_var100 <- runsim_OneArm(nsim, parfunc(c(1,2),c(10,10))))
diffsourceeff12_var100$var <- 100
diffsourceeff12_var100$efdiff <- 2-1
diffsourceeff12_var100$plt05 <- diffsourceeff12_var100$pval < 0.05
diffsourceeff12_var100$pnexchlt080 <- diffsourceeff12_var100$pnexch > 0.7976455
diffsourceeff12_var100$pnexchlt020 <- diffsourceeff12_var100$pnexch > 0.20

table(diffsourceeff12_var100$pnexchlt080, diffsourceeff12_var100$plt05)
table(diffsourceeff12_var100$pnexchlt020, diffsourceeff12_var100$plt05)

#scenario 13: different source effects (1 and 3), sd = 10
system.time(diffsourceeff13_var100 <- runsim_OneArm(nsim, parfunc(c(1,3),c(10,10))))
diffsourceeff13_var100$var <- 100
diffsourceeff13_var100$efdiff <- 3-1
diffsourceeff13_var100$plt05 <- diffsourceeff13_var100$pval < 0.05
diffsourceeff13_var100$pnexchlt080 <- diffsourceeff13_var100$pnexch > 0.7976455
diffsourceeff13_var100$pnexchlt020 <- diffsourceeff13_var100$pnexch > 0.20

table(diffsourceeff13_var100$pnexchlt080, diffsourceeff13_var100$plt05)
table(diffsourceeff13_var100$pnexchlt020, diffsourceeff13_var100$plt05)

#scenario 14: different source effects (1 and 4), sd = 10
system.time(diffsourceeff14_var100 <- runsim_OneArm(nsim, parfunc(c(1,4),c(10,10))))
diffsourceeff14_var100$var <- 100
diffsourceeff14_var100$efidff <- 4-1
diffsourceeff14_var100$plt05 <- diffsourceeff14_var100$pval < 0.05
diffsourceeff14_var100$pnexchlt080 <- diffsourceeff14_var100$pnexch > 0.7976455
diffsourceeff14_var100$pnexchlt020 <- diffsourceeff14_var100$pnexch > 0.20

table(diffsourceeff14_var100$pnexchlt080, diffsourceeff14_var100$plt05)
table(diffsourceeff14_var100$pnexchlt020, diffsourceeff14_var100$plt05)

#scenario 15: different source effects (1 and 5), sd = 10
system.time(diffsourceeff15_var100 <- runsim_OneArm(nsim, parfunc(c(1,5),c(10,10))))
diffsourceeff15_var100$var <- 100
diffsourceeff15_var100$efdiff <- 5-1
diffsourceeff15_var100$plt05 <- diffsourceeff15_var100$pval < 0.05
diffsourceeff15_var100$pnexchlt080 <- diffsourceeff15_var100$pnexch > 0.7976455
diffsourceeff15_var100$pnexchlt020 <- diffsourceeff15_var100$pnexch > 0.20

table(diffsourceeff15_var100$pnexchlt080, diffsourceeff15_var100$plt05)
table(diffsourceeff15_var100$pnexchlt020, diffsourceeff15_var100$plt05)







#combine data sets
allresults2 <- rbind(equalsourceeff11_var1,equalsourceeff11_var10,equalsourceeff11_var100,
                     diffsourceeff12_var1,diffsourceeff12_var10,diffsourceeff12_var100,
                     diffsourceeff13_var1,diffsourceeff13_var10,diffsourceeff13_var100,
                     diffsourceeff14_var1,diffsourceeff14_var10,diffsourceeff14_var100,
                     diffsourceeff15_var1,diffsourceeff15_var10,diffsourceeff15_var100)

all_summ02 <- allresults2 %>% group_by(efdiff, var, totsampsize)

all_summ2 <- all_summ02 %>% summarize(plt05 = mean(plt05), 
                                      pnexchlt080 = mean(pnexchlt080), 
                                      pnexchlt020 = mean(pnexchlt020))


#write.csv(allresults2, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/lmresults_raw.csv")
#write.csv(all_summ2, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/lmresults.csv")


#############################################
#############UNEQUAL SAMPLE##################
#############################################
#function to format inputs 
parfunc_unequalsamp <- function(effectsizes, n2mult){
  pars=list(
    #n = sample sizes, m0 = mean for each source, 
    list(n=c(200*(1-n2mult), 200*n2mult), mu0=effectsizes,
         sd=c(sqrt(10),sqrt(10)), BICiterations=1000, 
         marginal="BIC"),
    list(n=c(300*(1-n2mult), 300*n2mult), mu0=effectsizes,
         sd=c(sqrt(10),sqrt(10)), BICiterations=1000, 
         marginal="BIC"),
    list(n=c(400*(1-n2mult), 400*n2mult), mu0=effectsizes,
         sd=c(sqrt(10),sqrt(10)), BICiterations=1000, 
         marginal="BIC"),
    list(n=c(500*(1-n2mult), 500*n2mult), mu0=effectsizes,
         sd=c(sqrt(10),sqrt(10)), BICiterations=1000, 
         marginal="BIC"),
    list(n=c(600*(1-n2mult), 600*n2mult), mu0=effectsizes,
         sd=c(sqrt(10),sqrt(10)), BICiterations=1000, 
         marginal="BIC"),
    list(n=c(700*(1-n2mult), 700*n2mult), mu0=effectsizes,
         sd=c(sqrt(10),sqrt(10)), BICiterations=1000, 
         marginal="BIC"),
    list(n=c(800*(1-n2mult), 800*n2mult), mu0=effectsizes,
         sd=c(sqrt(10),sqrt(10)), BICiterations=1000, 
         marginal="BIC"),
    list(n=c(900*(1-n2mult), 900*n2mult), mu0=effectsizes,
         sd=c(sqrt(10),sqrt(10)), BICiterations=1000, 
         marginal="BIC"),
    list(n=c(1000*(1-n2mult), 1000*n2mult), mu0=effectsizes,
         sd=c(sqrt(10),sqrt(10)), BICiterations=1000, 
         marginal="BIC")
  )
  
  return(pars)
}

nsim=1000

#scenario 1: equal source effects (1 and 1), n2mult = 0.1
system.time(equalsourceeff11_n2mult01 <- runsim_OneArm(nsim, parfunc_unequalsamp(effectsizes=c(1,1),n2mult = 0.1)))
equalsourceeff11_n2mult01$var <- 10
equalsourceeff11_n2mult01$efdiff <- 1-1
equalsourceeff11_n2mult01$n2mult <- 0.1
equalsourceeff11_n2mult01$plt05 <- equalsourceeff11_n2mult01$pval < 0.05
equalsourceeff11_n2mult01$pnexchlt080 <- equalsourceeff11_n2mult01$pnexch > 0.7976455
equalsourceeff11_n2mult01$pnexchlt020 <- equalsourceeff11_n2mult01$pnexch > 0.20

table(equalsourceeff11_n2mult01$pnexchlt080, equalsourceeff11_n2mult01$plt05)
table(equalsourceeff11_n2mult01$pnexchlt020, equalsourceeff11_n2mult01$plt05)

#scenario 2: different source effects (1 and 2), n2mult = 0.1
system.time(diffsourceeff12_n2mult01 <- runsim_OneArm(nsim, parfunc_unequalsamp(c(1,2),n2mult = 0.1)))
diffsourceeff12_n2mult01$var <- 10
diffsourceeff12_n2mult01$efdiff <- 2-1
diffsourceeff12_n2mult01$n2mult <- 0.1
diffsourceeff12_n2mult01$plt05 <- diffsourceeff12_n2mult01$pval < 0.05
diffsourceeff12_n2mult01$pnexchlt080 <- diffsourceeff12_n2mult01$pnexch > 0.7976455
diffsourceeff12_n2mult01$pnexchlt020 <- diffsourceeff12_n2mult01$pnexch > 0.20

table(diffsourceeff12_n2mult01$pnexchlt080, diffsourceeff12_n2mult01$plt05)
table(diffsourceeff12_n2mult01$pnexchlt020, diffsourceeff12_n2mult01$plt05)

#scenario 3: different source effects (1 and 3), n2mult = 0.1
system.time(diffsourceeff13_n2mult01 <- runsim_OneArm(nsim, parfunc_unequalsamp(c(1,3),n2mult = 0.1)))
diffsourceeff13_n2mult01$var <- 10
diffsourceeff13_n2mult01$efdiff <- 3-1
diffsourceeff13_n2mult01$n2mult <- 0.1
diffsourceeff13_n2mult01$plt05 <- diffsourceeff13_n2mult01$pval < 0.05
diffsourceeff13_n2mult01$pnexchlt080 <- diffsourceeff13_n2mult01$pnexch > 0.7976455
diffsourceeff13_n2mult01$pnexchlt020 <- diffsourceeff13_n2mult01$pnexch > 0.20

table(diffsourceeff13_n2mult01$pnexchlt080, diffsourceeff13_n2mult01$plt05)
table(diffsourceeff13_n2mult01$pnexchlt020, diffsourceeff13_n2mult01$plt05)

#scenario 4: different source effects (1 and 4), n2mult = 0.1
system.time(diffsourceeff14_n2mult01 <- runsim_OneArm(nsim, parfunc_unequalsamp(c(1,4),n2mult = 0.1)))
diffsourceeff14_n2mult01$var <- 10
diffsourceeff14_n2mult01$efdiff <- 4-1
diffsourceeff14_n2mult01$n2mult <- 0.1
diffsourceeff14_n2mult01$plt05 <- diffsourceeff14_n2mult01$pval < 0.05
diffsourceeff14_n2mult01$pnexchlt080 <- diffsourceeff14_n2mult01$pnexch > 0.7976455
diffsourceeff14_n2mult01$pnexchlt020 <- diffsourceeff14_n2mult01$pnexch > 0.20

table(diffsourceeff14_n2mult01$pnexchlt080, diffsourceeff14_n2mult01$plt05)
table(diffsourceeff14_n2mult01$pnexchlt020, diffsourceeff14_n2mult01$plt05)

#scenario 5: different source effects (1 and 5), n2mult = 0.1
system.time(diffsourceeff15_n2mult01 <- runsim_OneArm(nsim, parfunc_unequalsamp(c(1,5),n2mult = 0.1)))
diffsourceeff15_n2mult01$var <- 10
diffsourceeff15_n2mult01$efdiff <- 5-1
diffsourceeff15_n2mult01$n2mult <- 0.1
diffsourceeff15_n2mult01$plt05 <- diffsourceeff15_n2mult01$pval < 0.05
diffsourceeff15_n2mult01$pnexchlt080 <- diffsourceeff15_n2mult01$pnexch > 0.7976455
diffsourceeff15_n2mult01$pnexchlt020 <- diffsourceeff15_n2mult01$pnexch > 0.20

table(diffsourceeff15_n2mult01$pnexchlt080, diffsourceeff15_n2mult01$plt05)
table(diffsourceeff15_n2mult01$pnexchlt020, diffsourceeff15_n2mult01$plt05)





#scenario 6: equal source effects (1 and 1), n2mult = 0.25
system.time(equalsourceeff11_n2mult025 <- runsim_OneArm(nsim, parfunc_unequalsamp(c(1,1),n2mult = 0.25)))
equalsourceeff11_n2mult025$var <- 10
equalsourceeff11_n2mult025$efdiff <- 1-1
equalsourceeff11_n2mult025$n2mult <- 0.25
equalsourceeff11_n2mult025$plt05 <- equalsourceeff11_n2mult025$pval < 0.05
equalsourceeff11_n2mult025$pnexchlt080 <- equalsourceeff11_n2mult025$pnexch > 0.7976455
equalsourceeff11_n2mult025$pnexchlt020 <- equalsourceeff11_n2mult025$pnexch > 0.20

table(equalsourceeff11_n2mult025$pnexchlt080, equalsourceeff11_n2mult025$plt05)
table(equalsourceeff11_n2mult025$pnexchlt020, equalsourceeff11_n2mult025$plt05)

#scenario 7: different source effects (1 and 2), n2mult = 0.25
system.time(diffsourceeff12_n2mult025 <- runsim_OneArm(nsim, parfunc_unequalsamp(c(1,2),n2mult = 0.25)))
diffsourceeff12_n2mult025$var <- 10
diffsourceeff12_n2mult025$efdiff <- 2-1
diffsourceeff12_n2mult025$n2mult <- 0.25
diffsourceeff12_n2mult025$plt05 <- diffsourceeff12_n2mult025$pval < 0.05
diffsourceeff12_n2mult025$pnexchlt080 <- diffsourceeff12_n2mult025$pnexch > 0.7976455
diffsourceeff12_n2mult025$pnexchlt020 <- diffsourceeff12_n2mult025$pnexch > 0.20

table(diffsourceeff12_n2mult025$pnexchlt080, diffsourceeff12_n2mult025$plt05)
table(diffsourceeff12_n2mult025$pnexchlt020, diffsourceeff12_n2mult025$plt05)

#scenario 8: different source effects (1 and 3), n2mult = 0.25
system.time(diffsourceeff13_n2mult025 <- runsim_OneArm(nsim, parfunc_unequalsamp(c(1,3),n2mult = 0.25)))
diffsourceeff13_n2mult025$var <- 10
diffsourceeff13_n2mult025$efdiff <- 3-1
diffsourceeff13_n2mult025$n2mult <- 0.25
diffsourceeff13_n2mult025$plt05 <- diffsourceeff13_n2mult025$pval < 0.05
diffsourceeff13_n2mult025$pnexchlt080 <- diffsourceeff13_n2mult025$pnexch > 0.7976455
diffsourceeff13_n2mult025$pnexchlt020 <- diffsourceeff13_n2mult025$pnexch > 0.20

table(diffsourceeff13_n2mult025$pnexchlt080, diffsourceeff13_n2mult025$plt05)
table(diffsourceeff13_n2mult025$pnexchlt020, diffsourceeff13_n2mult025$plt05)

#scenario 9: different source effects (1 and 4), n2mult = 0.25
system.time(diffsourceeff14_n2mult025 <- runsim_OneArm(nsim, parfunc_unequalsamp(c(1,4),n2mult = 0.25)))
diffsourceeff14_n2mult025$var <- 10
diffsourceeff14_n2mult025$efdiff <- 4-1
diffsourceeff14_n2mult025$n2mult <- 0.25
diffsourceeff14_n2mult025$plt05 <- diffsourceeff14_n2mult025$pval < 0.05
diffsourceeff14_n2mult025$pnexchlt080 <- diffsourceeff14_n2mult025$pnexch > 0.7976455
diffsourceeff14_n2mult025$pnexchlt020 <- diffsourceeff14_n2mult025$pnexch > 0.20

table(diffsourceeff14_n2mult025$pnexchlt080, diffsourceeff14_n2mult025$plt05)
table(diffsourceeff14_n2mult025$pnexchlt020, diffsourceeff14_n2mult025$plt05)

#scenario 10: different source effects (1 and 5), n2mult = 0.25
system.time(diffsourceeff15_n2mult025 <- runsim_OneArm(nsim, parfunc_unequalsamp(c(1,5),n2mult = 0.25)))
diffsourceeff15_n2mult025$var <- 10
diffsourceeff15_n2mult025$efdiff <- 5-1
diffsourceeff15_n2mult025$n2mult <- 0.25
diffsourceeff15_n2mult025$plt05 <- diffsourceeff15_n2mult025$pval < 0.05
diffsourceeff15_n2mult025$pnexchlt080 <- diffsourceeff15_n2mult025$pnexch > 0.7976455
diffsourceeff15_n2mult025$pnexchlt020 <- diffsourceeff15_n2mult025$pnexch > 0.20

table(diffsourceeff15_n2mult025$pnexchlt080, diffsourceeff15_n2mult025$plt05)
table(diffsourceeff15_n2mult025$pnexchlt020, diffsourceeff15_n2mult025$plt05)


#combine data sets
allresults2_unequalsamp <- rbind(equalsourceeff11_n2mult01,equalsourceeff11_n2mult025,
                     diffsourceeff12_n2mult01,diffsourceeff12_n2mult025,
                     diffsourceeff13_n2mult01,diffsourceeff13_n2mult025,
                     diffsourceeff14_n2mult01,diffsourceeff14_n2mult025,
                     diffsourceeff15_n2mult01,diffsourceeff15_n2mult025)

all_summ02_unequalsamp <- allresults2_unequalsamp %>% group_by(efdiff, var, totsampsize, n2mult)

all_summ2_unequalsamp <- all_summ02_unequalsamp %>% summarize(plt05 = mean(plt05), 
                                      pnexchlt080 = mean(pnexchlt080), 
                                      pnexchlt020 = mean(pnexchlt020))


write.csv(allresults2_unequalsamp, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/OA_cont_MEMr_result_raw.csv")
write.csv(all_summ2_unequalsamp, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/OA_cont_MEMr_result.csv")
