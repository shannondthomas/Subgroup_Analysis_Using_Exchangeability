#source('H:/Advising/Murphy_Shannon/Subgroup Detection/Kaizer_MEM_Weights.R')
source("C:/Users/mushanno/Desktop/Work/Dissertation1/CodeFromAlex/GitHub_Functions_AdaptivePlatformDesign.R")
library(pROC)
library(pracma)


# calc.MEM.betabin <- function(xvec, nvec, avec, bvec, prior, constraint=1){
###Function to calculate the MEM model weights for binomial data with beta(alpha,beta) prior
#xvec: vector with counts of those with event (primary source first, supplemental afterwards)
#nvec: vector of total number of subjects
#avec: vector of alpha parameters for beta priors for each source
#bvec: vector of beta parameters for beta priors for each source
#prior: prior to use on MEM source inclusion probability: equal (pi_e), pool (naive pooling), opt.source (pi_EB), opt.source_constrain (pi_EBc)
#constraint: value to limit maximum value of empirical Bayes prior to for pi_EBc


#calc.MEM.betabin(xvec = c(8,1), nvec = c(10,10), avec = c(1,1), bvec = c(1,1), prior = 'equal',
#                 constraint = 1)$q



getZ <- function(Z, p1, n1, n2){
  
  pb <- txtProgressBar(min = 0, max = length(Z),style=3)
  
  p2 <- rep(NA, length(p1))
  
  for (i in 1:length(p1)){
    zerofunc <- function(p2){
      Z[i] - (p1[i]-p2)/sqrt(((p1[i]*n1[i] + p2*n2[i])/(n1[i] + n2[i]))*(1-(p1[i]*n1[i] + p2*n2[i])/(n1[i] + n2[i]))*((1/n1[i]) + (1/n2[i])))
    }
    

    if (length(findzeros(zerofunc,a = 0, b = 1) != 0)){
      p2[i] <- findzeros(zerofunc, a = 0, b = 1)
    }
    
    setTxtProgressBar(pb, i)

  }
  
  close(pb)
  
  df <- data.frame(p2 = p2, x2 = p2*n2, n2 = n2, p1 = p1, x1 = p1*n1, n1 = n1, Z = Z)
  return(df)
}
# 
# p1 <- seq(0.1,0.9,0.001)
# n1 <- seq(50,200,50)
# n2 <- seq(50,200,50)
# Z <- seq(-4,4,0.1)
# 
# testvals <- expand.grid(p1=p1,n1=n1,n2=n2,Z=Z)
# 
# dftest <- getZ(Z = testvals$Z, p1 = testvals$p1, n1 = testvals$n1, n2 = testvals$n2)
# 
# #get exchangeability probabilities 
# pb <- txtProgressBar(min = 0, max = length(dftest),style=3)
# for(i in 1:nrow(dftest)){
#   dftest$p1e[i] <- calc.MEM.betabin(xvec = c(dftest$x1[i],dftest$x2[i]), 
#                                     nvec = c(dftest$n1[i],dftest$n2[i]), 
#                                     avec = c(1,1), bvec = c(1,1), prior = 'equal',
#                                     constraint = 1)$q[1]
# 
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# 
# write.csv(dftest, "C:/Users/mushanno/Desktop/Work/BinaryOutcome/binaryZ.csv", row.names = FALSE)

dftest <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/binaryZ.csv")
dftest_sub <- dftest[!is.na(dftest$p2),]

plot(x=abs(dftest_sub$Z), y=dftest_sub$p1e, ylim=c(0,1), xlab='|Z|', ylab='Pr(Not Exchangeable)', main='pi_e')

# Logistic regression to predict |Z|>=1
dftest_sub$z_gte1 <- abs(dftest_sub$Z)>=1
dftest_sub$perc_p1e <- dftest_sub$p1e*100

dftest_sub$col <- 'black'
dftest_sub$col[which(dftest_sub$z_gte1==T)] <- 'blue'
dftest_sub$pch <- 3
dftest_sub$pch[which(dftest_sub$z_gte1==T)] <- 3

mod3 <- glm(z_gte1 ~ perc_p1e, data=dftest_sub, family='binomial')
summary(mod3)
mod3pred <- predict(mod3, type='response')
mean( dftest_sub$z_gte1[which(abs(dftest_sub$z) <= 2)] == (mod3pred >0.5))


roc( dftest_sub$z_gte1 ~ mod3pred)
coordout <- coords(roc( dftest_sub$z_gte1 ~ mod3pred), best.method='youden', x='best')
coordout #threshold = 0.8082359 
table(zgte1=dftest_sub$z_gte1, pred=(mod3pred>coordout[[1]] ))
mean( dftest_sub$z_gte1 == (mod3pred>coordout[[1]] )) 

plot(x=dftest_sub$perc_p1e, y=mod3pred, col=dftest_sub$col, pch=dftest_sub$pch)


#visualize and summarize results to look for other potential cutoffs
boxplot(perc_p1e ~ z_gte1, data=dftest_sub[dftest_sub$Z < 2,], main = "|Z| < 2 Subset")
boxplot(perc_p1e ~ z_gte1, data=dftest_sub, main = "Full Data Set")

max(dftest_sub$perc_p1e[dftest_sub$Z == 0]) #19.6897
quantile(dftest_sub$perc_p1e[dftest_sub$Z == 0], 0.975) #18.86296
max(dftest_sub$perc_p1e[dftest_sub$z_gte1 == FALSE]) #26.73151
quantile(dftest_sub$perc_p1e[dftest_sub$z_gte1 == FALSE], 0.975) #22.55774

# Strategy 1: simulation calibration, potentially using p*C_5 + (1-p)*C_80
# Strategy 2: general calibration above using range of values/combos:
#- Logistic regression with just pi_e
#- Logistic regression with pi_e + n1, n2, s1, s2, ratio, etc.
#- These are affected by the choice of threshold (e.g., >0.5 versus >Youden's J)
# Strategy 3: general calibration but focused moreso on range of values for trial
#- Could focus on power calculation assumptions
#- Could assume range of prevalence in each subgroup
# Strategy 4: consider different models for Z<0.5; 0.5 to 1.5; etc.




