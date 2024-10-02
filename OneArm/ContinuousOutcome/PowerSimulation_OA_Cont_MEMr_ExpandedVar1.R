library(tidyverse)
library(snowfall)
library(ggplot2)
library(tidyverse)

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
system.time(equalsourceeff0_var1 <- runsim_OneArm(nsim, parfunc(effectsizes=c(1,1),sds=c(1,1))))
equalsourceeff0_var1$plt05 <- equalsourceeff0_var1$pval < 0.05
equalsourceeff0_var1$pnexchlt080 <- equalsourceeff0_var1$pnexch > 0.7976455
equalsourceeff0_var1$pnexchlt020 <- equalsourceeff0_var1$pnexch > 0.20

table(equalsourceeff0_var1$pnexchlt080, equalsourceeff0_var1$plt05)
table(equalsourceeff0_var1$pnexchlt020, equalsourceeff0_var1$plt05)

#scenario 2: different source effects (1 and 1.1), sd = 1
system.time(diffsourceeff01_var1 <- runsim_OneArm(nsim, parfunc(c(1,1.1),c(1,1))))
diffsourceeff01_var1$plt05 <- diffsourceeff01_var1$pval < 0.05
diffsourceeff01_var1$pnexchlt080 <- diffsourceeff01_var1$pnexch > 0.7976455
diffsourceeff01_var1$pnexchlt020 <- diffsourceeff01_var1$pnexch > 0.20

table(diffsourceeff01_var1$pnexchlt080, diffsourceeff01_var1$plt05)
table(diffsourceeff01_var1$pnexchlt020, diffsourceeff01_var1$plt05)

#scenario 3: different source effects (1 and 1.2), sd = 1
system.time(diffsourceeff02_var1 <- runsim_OneArm(nsim, parfunc(c(1,1.2),c(1,1))))
diffsourceeff02_var1$plt05 <- diffsourceeff02_var1$pval < 0.05
diffsourceeff02_var1$pnexchlt080 <- diffsourceeff02_var1$pnexch > 0.7976455
diffsourceeff02_var1$pnexchlt020 <- diffsourceeff02_var1$pnexch > 0.20

table(diffsourceeff02_var1$pnexchlt080, diffsourceeff02_var1$plt05)
table(diffsourceeff02_var1$pnexchlt020, diffsourceeff02_var1$plt05)

#scenario 4: different source effects (1 and 1.3), sd = 1
system.time(diffsourceeff03_var1 <- runsim_OneArm(nsim, parfunc(c(1,1.3),c(1,1))))
diffsourceeff03_var1$plt05 <- diffsourceeff03_var1$pval < 0.05
diffsourceeff03_var1$pnexchlt080 <- diffsourceeff03_var1$pnexch > 0.7976455
diffsourceeff03_var1$pnexchlt020 <- diffsourceeff03_var1$pnexch > 0.20

table(diffsourceeff03_var1$pnexchlt080, diffsourceeff03_var1$plt05)
table(diffsourceeff03_var1$pnexchlt020, diffsourceeff03_var1$plt05)

#scenario 5: different source effects (1 and 1.4), sd = 1
system.time(diffsourceeff04_var1 <- runsim_OneArm(nsim, parfunc(c(1,1.4),c(1,1))))
diffsourceeff04_var1$plt05 <- diffsourceeff04_var1$pval < 0.05
diffsourceeff04_var1$pnexchlt080 <- diffsourceeff04_var1$pnexch > 0.7976455
diffsourceeff04_var1$pnexchlt020 <- diffsourceeff04_var1$pnexch > 0.20

table(diffsourceeff04_var1$pnexchlt080, diffsourceeff04_var1$plt05)
table(diffsourceeff04_var1$pnexchlt020, diffsourceeff04_var1$plt05)

#scenario 6: different source effects (1 and 1.5), sd = 1
system.time(diffsourceeff05_var1 <- runsim_OneArm(nsim, parfunc(c(1,1.5),c(1,1))))
diffsourceeff05_var1$plt05 <- diffsourceeff05_var1$pval < 0.05
diffsourceeff05_var1$pnexchlt080 <- diffsourceeff05_var1$pnexch > 0.7976455
diffsourceeff05_var1$pnexchlt020 <- diffsourceeff05_var1$pnexch > 0.20

table(diffsourceeff05_var1$pnexchlt080, diffsourceeff05_var1$plt05)
table(diffsourceeff05_var1$pnexchlt020, diffsourceeff05_var1$plt05)

#scenario 7: different source effects (1 and 1.6), sd = 1
system.time(diffsourceeff06_var1 <- runsim_OneArm(nsim, parfunc(c(1,1.6),c(1,1))))
diffsourceeff06_var1$plt05 <- diffsourceeff06_var1$pval < 0.05
diffsourceeff06_var1$pnexchlt080 <- diffsourceeff06_var1$pnexch > 0.7976455
diffsourceeff06_var1$pnexchlt020 <- diffsourceeff06_var1$pnexch > 0.20

table(diffsourceeff06_var1$pnexchlt080, diffsourceeff06_var1$plt05)
table(diffsourceeff06_var1$pnexchlt020, diffsourceeff06_var1$plt05)

#scenario 8: different source effects (1 and 1.7), sd = 1
system.time(diffsourceeff07_var1 <- runsim_OneArm(nsim, parfunc(c(1,1.7),c(1,1))))
diffsourceeff07_var1$plt05 <- diffsourceeff07_var1$pval < 0.05
diffsourceeff07_var1$pnexchlt080 <- diffsourceeff07_var1$pnexch > 0.7976455
diffsourceeff07_var1$pnexchlt020 <- diffsourceeff07_var1$pnexch > 0.20

table(diffsourceeff07_var1$pnexchlt080, diffsourceeff07_var1$plt05)
table(diffsourceeff07_var1$pnexchlt020, diffsourceeff07_var1$plt05)

#scenario 9: different source effects (1 and 1.8), sd = 1
system.time(diffsourceeff08_var1 <- runsim_OneArm(nsim, parfunc(c(1,1.8),c(1,1))))
diffsourceeff08_var1$plt05 <- diffsourceeff08_var1$pval < 0.05
diffsourceeff08_var1$pnexchlt080 <- diffsourceeff08_var1$pnexch > 0.7976455
diffsourceeff08_var1$pnexchlt020 <- diffsourceeff08_var1$pnexch > 0.20

table(diffsourceeff08_var1$pnexchlt080, diffsourceeff08_var1$plt05)
table(diffsourceeff08_var1$pnexchlt020, diffsourceeff08_var1$plt05)

#scenario 10: different source effects (1 and 1.9), sd = 1
system.time(diffsourceeff09_var1 <- runsim_OneArm(nsim, parfunc(c(1,1.9),c(1,1))))
diffsourceeff09_var1$plt05 <- diffsourceeff09_var1$pval < 0.05
diffsourceeff09_var1$pnexchlt080 <- diffsourceeff09_var1$pnexch > 0.7976455
diffsourceeff09_var1$pnexchlt020 <- diffsourceeff09_var1$pnexch > 0.20

table(diffsourceeff09_var1$pnexchlt080, diffsourceeff09_var1$plt05)
table(diffsourceeff09_var1$pnexchlt020, diffsourceeff09_var1$plt05)

#scenario 11: different source effects (1 and 2), sd = 1
system.time(diffsourceeff10_var1 <- runsim_OneArm(nsim, parfunc(c(1,2),c(1,1))))
diffsourceeff10_var1$plt05 <- diffsourceeff10_var1$pval < 0.05
diffsourceeff10_var1$pnexchlt080 <- diffsourceeff10_var1$pnexch > 0.7976455
diffsourceeff10_var1$pnexchlt020 <- diffsourceeff10_var1$pnexch > 0.20

table(diffsourceeff10_var1$pnexchlt080, diffsourceeff10_var1$plt05)
table(diffsourceeff10_var1$pnexchlt020, diffsourceeff10_var1$plt05)







#combine data sets
equalsourceeff0_var1$var <- rep(1, length(equalsourceeff0_var1$pnexch))
diffsourceeff01_var1$var <- rep(1, length(diffsourceeff01_var1$pnexch))
diffsourceeff02_var1$var <- rep(1, length(diffsourceeff02_var1$pnexch))
diffsourceeff03_var1$var <- rep(1, length(diffsourceeff03_var1$pnexch))
diffsourceeff04_var1$var <- rep(1, length(diffsourceeff04_var1$pnexch))
diffsourceeff05_var1$var <- rep(1, length(diffsourceeff05_var1$pnexch))
diffsourceeff06_var1$var <- rep(1, length(diffsourceeff06_var1$pnexch))
diffsourceeff07_var1$var <- rep(1, length(diffsourceeff07_var1$pnexch))
diffsourceeff08_var1$var <- rep(1, length(diffsourceeff08_var1$pnexch))
diffsourceeff09_var1$var <- rep(1, length(diffsourceeff09_var1$pnexch))
diffsourceeff10_var1$var <- rep(1, length(diffsourceeff10_var1$pnexch))

equalsourceeff0_var1$efdiff <- rep(0, length(equalsourceeff0_var1$pnexch))
diffsourceeff01_var1$efdiff <- rep(0.1, length(diffsourceeff01_var1$pnexch))
diffsourceeff02_var1$efdiff <- rep(0.2, length(diffsourceeff02_var1$pnexch))
diffsourceeff03_var1$efdiff <- rep(0.3, length(diffsourceeff03_var1$pnexch))
diffsourceeff04_var1$efdiff <- rep(0.4, length(diffsourceeff04_var1$pnexch))
diffsourceeff05_var1$efdiff <- rep(0.5, length(diffsourceeff05_var1$pnexch))
diffsourceeff06_var1$efdiff <- rep(0.6, length(diffsourceeff06_var1$pnexch))
diffsourceeff07_var1$efdiff <- rep(0.7, length(diffsourceeff07_var1$pnexch))
diffsourceeff08_var1$efdiff <- rep(0.8, length(diffsourceeff08_var1$pnexch))
diffsourceeff09_var1$efdiff <- rep(0.9, length(diffsourceeff09_var1$pnexch))
diffsourceeff10_var1$efdiff <- rep(1, length(diffsourceeff10_var1$pnexch))

allresults2 <- rbind(equalsourceeff0_var1,
                     diffsourceeff01_var1,
                     diffsourceeff02_var1,
                     diffsourceeff03_var1,
                     diffsourceeff04_var1,
                     diffsourceeff05_var1,
                     diffsourceeff06_var1,
                     diffsourceeff07_var1,
                     diffsourceeff08_var1,
                     diffsourceeff09_var1,
                     diffsourceeff10_var1)

all_summ02 <- allresults2 %>% group_by(efdiff, var, totsampsize)

all_summ2 <- all_summ02 %>% summarize(plt05 = mean(plt05), 
                                      pnexchlt080 = mean(pnexchlt080), 
                                      pnexchlt020 = mean(pnexchlt020))


write.csv(allresults2, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/lmresults_expandvar1_raw.csv")
write.csv(all_summ2, "C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/lmresults_expandvar1.csv")


ggplot(all_summ2, aes(efdiff, totsampsize, fill=plt05)) + 
  geom_tile() + geom_text(aes(label=plt05)) + 
  labs(title = "Heatmap of p-value < 0.05 for Varying Sample Sizes and Differences in Effect Sizes")

ggplot(all_summ2, aes(efdiff, totsampsize, fill=pnexchlt080)) + 
  geom_tile() + geom_text(aes(label=pnexchlt080)) +
  labs(title = "Heatmap of P(Not Exchangeable) > 0.7976455 for Varying Sample Sizes and Differences in Effect Sizes")

ggplot(all_summ2, aes(efdiff, totsampsize, fill=pnexchlt020)) + 
  geom_tile() + geom_text(aes(label=pnexchlt020)) +
  labs(title = "Heatmap of P(Not Exchangeable) > 0.2 for Varying Sample Sizes and Differences in Effect Sizes")



