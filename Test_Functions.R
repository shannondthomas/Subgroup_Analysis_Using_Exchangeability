#################################################################################
# TITLE: Test_Functions.R
#
# PURPOSE: Store functions for MEM, MEMr, and standard test given 
#          a (simulated) data set. The main function, RunModels(), takes the
#          the following inputs and gives the outputs described below. Other 
#          functions needed to calculate pure MEM weights are also included
#          in this script. 
#
#
#
# INPUT: num_arms - number of treatment arms (1 or 2)   
#        outcome_type - "binary" or "continuous"
#                        NOTE IF data IS PROVIDED THIS MUST MATCH data$Y FORM, 
#                        THE FUNCTION WILL NOT CHECK IF THEY MATCH
#        n_sources - number of subgroup levels (must be 2)
#        marginal - "BIC" or "AIC" model selection criteria 
#                   (only BIC was used in this project)
#        data - (OPTIONAL) must have the following columns
#           Y: outcome variable (can be binary or continuous)
#           S2: indicator for subgroup level (source) 2
#           df: subgroup level (1 if level 1 and 2 if level 2)
#           IF TWO-ARM SETTING, THE FOLLOWING TWO COLUMNS ALSO MUST BE INCLUDED
#           trt: treatment indicator 
#                (1 for trt, 0 for placebo)
#           trt_S2: treatment AND level 2 indicator (trt*S2)  
#
#        *THE FOLLOWING PARAMETERS MUST BE USED IF data IS NOT PROVIDED.
#         ALL OF THE FOLLOWING ARE VECTORS OF LENGTH 2 (n_sources)
#        n - sample sizes for each subgroup level
#        means - means for each subgroup 
#                if outcome is binary, these values must be between 0 and 1
#        sds - standard deviations for each subgroup
#              only needed for continuous outcome
#        trt_effect - difference in means for treatment vs placebo for each
#                     subgroup level
#                     only needed for two arm setting
#
#
# OUTPUT: list with the following elements
#            $sampsize: total sample size of data
#            $weights_MEM: double with first element giving the 
#                          MEM probability of exchangeability and 
#                          second element giving the probability of 
#                          nonexchangeability
#            $weights_MEMr: double with first element giving the 
#                           MEMr probability of exchangeability and 
#                           second element giving the probability of 
#                           nonexchangeability
#            $pval: p-value from standard test
#                   -binary standard test is LRT
#                   -continuous standard test is a partial F-test
#
#
# FUNCTIONS: Section 1
#              calc.weights_MEM_cont() - calc. weights for MEM w/ continuous Y
#              BinaryESS2() - calculate ESS (used in cal.weights_MEM_bin)
#              calc.weights_MEM_bin() - calc. weights for MEM w/ binary Y
#            Section 2
#              SimulateData() - function used by RunModels() if no data is input
#            Section 3
#              RunModels() - run standard test and MEMs
#                            THIS IS THE MAIN FUNCTION CALLED IN OTHER SCRIPTS
#                            AND GIVES THE OUTPUTS DESCRIBED ABOVE.
#
# AUTHOR: Shannon Thomas
# DATE CREATED: OCT 22, 2024
#################################################################################


##########################################################
############## SECTION 0: LOAD DEPENDENCIES ############## 
##########################################################

#R version 4.3.2 was used for this project.
library(mvtnorm)      #version 1.2-4
library(invgamma)     #version 1.1
library(lmtest)       #version 0.9-40
library(matrixStats)  #version 1.3.0
library(xtable)       #version 1.8-4



##########################################################
########## SECTION 1: PURE MEM WEIGHT FUNCTIONS ########## 
##########################################################

#### SECTION 1.1: CONTINUOUS OUTCOME

calc.weights_MEM_cont = function(xvec,svec,nvec,prior){
  ###function to calculate model weights for MEM approach with a continuous outcome 
  #xvec: means for sources
  #svec: standard deviation for sources
  #nvec: sample size for sources
  #prior: prior to use for calculations
  
  if(length(xvec)==2){
    xbar = xvec[1] 
    xbar01 = xvec[2]
    
    s = svec[1]^2
    s01 = svec[2]^2
    
    n = nvec[1]
    n01 = nvec[2]
    
    v = s/n
    v01 = s01/n01
    
    ###Calculating the weights
    #Writing out the marginal models
    m1 = sqrt(2*pi)^2 / sqrt(1/(v*v01))
    
    m2 = sqrt(2*pi)^1 / sqrt((1/v + 1/v01)) * exp(-0.5 * ((xbar-xbar01)^2/(v + v01)) )
    
    if(prior=='pi_e'){
      w1 <- 0.5*m1; w2 <- 0.5*m2
      m.sum.prior <- w1 + w2
      
      q1 <- w1/m.sum.prior; q2 <- w2/m.sum.prior
    }else if(prior=='pi_n'){
      m1inv.noCurN = sqrt(1/(s*v01))/(2*pi)^(2/2) #only current study excluded
      m2inv.noCurN = sqrt(1/s + 1/v01)/(2*pi)^(1/2)
      
      prior1.1 = m2inv.noCurN
      prior1.0 = m1inv.noCurN
      
      m1pr <- prior1.0; m2pr <- prior1.1
      
      w1 = m1*m1pr; w2 = m2*m2pr
      m.sum.prior = w1 + w2
      
      #Final weights:
      q1 = w1 / m.sum.prior; q2 = w2 / m.sum.prior
    }else{print('Only pi_e and pi_n currently supported.')}
    
    return(c(q1,q2))
  }
  
}


#### SECTION 1.2: BINARY OUTCOME

BinaryESS2 <- function(V,M){
  ###Function to calculate ESS for binary data
  #M: posterior mean
  #V: posterior variance
  
  u <- (M*(1-M))/V - 1
  a <- M*u
  b <- (1-M)*u
  return(a+b)
  
}

calc.weights_MEM_bin <- function(xvec, nvec, avec = c(1,1), bvec = c(1,1), prior, constraint=1){
  ###Function to calculate the MEM model weights for binomial data with beta(alpha,beta) prior
  #xvec: vector with counts of those with event (primary source first, supplemental afterwards)
  #nvec: vector of total number of subjects
  #avec: vector of alpha parameters for beta priors for each source (default is (1,1) -> Uniform(0,1))
  #bvec: vector of beta parameters for beta priors for each source (default is (1,1) -> Uniform(0,1))
  #prior: prior to use on MEM source inclusion probability: equal (pi_e), pool (naive pooling), opt.source (pi_EB), opt.source_constrain (pi_EBc)
  #constraint: value to limit maximum value of empirical Bayes prior to for pi_EBc
  
  mod.mat <- as.matrix(expand.grid( rep(list(c(0,1)), length(xvec)-1) )) #create matrix of all models with potential combinations of supplemental sources
  mod.mat <- mod.mat[order(rowSums(mod.mat)),] #group by number of supplemental sources
  mod.mat <- cbind(1, mod.mat) #add column before 1st column for primary source indicator
  colnames(mod.mat) <- c('p',paste0('s',seq(1,length(xvec)-1)))
  
  prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) #calculate the product portion of integrated marginal likelihood
  p.vec <- apply( t(sapply(1:dim(mod.mat)[1], function(x) prod.vec^(1-mod.mat[x,]))), 1, prod) #calculate the product portion of the integrated marginal likelihood corresponding to each model by taking power w/indicator for each source in a model and then using apply over rows of resulting matrix
  
  ###Calculate the integrated marginal likelihood given data:
  marg.vec <- (beta(avec[1] + mod.mat%*%xvec, bvec[1] + mod.mat%*%(nvec-xvec)) / beta(avec[1],bvec[1]) ) * p.vec
  
  ###Calculate prior:
  if(prior=='equal'){
    prior1 <- rep(.5, length(xvec)-1)
    prior0 <- rep(.5, length(xvec)-1)
  }else if(prior=='pool'){
    prior1 <- rep(1, length(xvec)-1)
    prior0 <- rep(0, length(xvec)-1)
  }else if(prior=='opt.source'){
    #this prior identifies the optimal MEM and gives it a weight of 1
    
    s.mat_opt.source <- mod.mat[,-1]
    min.vec <- -log( marg.vec ) 
    if(length(xvec)>2){
      if( length(which(min.vec == min(min.vec))) > 1){
        prior1.sum <- colSums(s.mat_opt.source[ which(min.vec == min(min.vec)) ,])
        prior1.sum[which(prior1.sum != 0)] <- 1
        prior1 <- prior1.sum
      }else{
        prior1 <- s.mat_opt.source[ which(min.vec == min(min.vec)) ,]
      }
      prior0 <- 1 - prior1
    }else{
      prior1 <- s.mat_opt.source[ which(min.vec == min(min.vec)) ]
      prior0 <- 1 - prior1
    }
  }else if(prior=='opt.source_constrain'){
    #this prior identifies the optimal MEM and gives it a weight of 1
    
    s.mat_opt.source <- mod.mat[,-1]
    min.vec <- -log( marg.vec ) 
    if(length(xvec)>2){
      if( length(which(min.vec == min(min.vec))) > 1){
        prior1.sum <- colSums(s.mat_opt.source[ which(min.vec == min(min.vec)) ,])
        prior1.sum[which(prior1.sum != 0)] <- 1
        prior1 <- prior1.sum*constraint
      }else{
        prior1 <- s.mat_opt.source[ which(min.vec == min(min.vec)) ,]*constraint
      }
      prior0 <- 1 - prior1
    }else{
      prior1 <- s.mat_opt.source[ which(min.vec == min(min.vec)) ]*constraint
      prior0 <- 1 - prior1
    }
  }else{print('Prior not valid, please enter a valid prior option.')}
  
  ###Calculate model weights given priors
  
  if(length(xvec)==2){mps <- matrix( rbind(prior0,prior1)[ paste0('prior',(mod.mat[,2:length(xvec)])), 1], ncol=1)} #extract priors for sources for each MEM
  if(length(xvec)>2){mps <- sapply(1:(length(xvec)-1), function(x) rbind(prior0,prior1)[ paste0('prior',(mod.mat[,2:length(xvec)])[,x]), x])} #extract priors for sources for each MEM
  
  q.vec <- marg.vec * ( rowProds(mps)/sum(rowProds(mps)) ) / sum(marg.vec * ( rowProds(mps)/sum(rowProds(mps)) ) ) #model weights, note: need to include sum(rowProds(mps)) otherwise some cases result in NaN because values are so small
  
  ret <- list(q = q.vec, mod.mat = mod.mat, prior = prior )
  
  return(ret)
}





##########################################################
################ SECTION 2: SIMULATE DATA ################ 
##########################################################

SimulateData <- function(n, means, sds = NULL, trt_effect = NULL, num_arms, outcome_type){
  #n: vector of sample sizes for each subgroup level (length 2)
  #means: vector of means for each subgroup (should be between 0 and 1 for binary outcome)
  #sd: vector of sds for each subgroup (only needed for continuous outcome)
  #trt_effect: difference in means/probabilities for treatment vs placebo each subgroup 
  #            (only needed when num_arms = 2)
  #outcome_type: "binary" or "continuous"
  
  
  n_sources <- length(n)
  
  ###ONE-ARM
  
  if((num_arms == 1) && (outcome_type == "binary")){
    
    #simulate data:
    #primary/level 1 of grouping variable
    Y <- rbinom(n[1], 1, prob=means[1]) 
    primary <- data.frame(Y)
    primary$df <- 1
    primary$secondary <- 0
    
    #secondary
    secondary <- NULL
    for(i in 2:n_sources) {
      Y <- rbinom(n[i], 1, prob= means[i])
      df <- i
      secondary <- rbind(secondary, data.frame(Y,df))
    }
    secondary$secondary <- 1
    
    
    #data frame together
    data <- data.frame(rbind(primary, secondary))
    
    
    #dummies for source
    for(i in 2:n_sources) {
      x <- 1*(data$df==i)
      data <- cbind(data, x)
    }
    
    names(data) <- c("Y", "df", "secondary", paste("S", 2:n_sources, sep=""))
    
    return(data)
  }
  
  if((num_arms == 1) && (outcome_type == "binary")){
    
    #simulate data:
    #primary/level 1 of grouping variable
    Y <- rbinom(n[1], 1, prob=means[1]) 
    primary <- data.frame(Y)
    primary$df <- 1
    primary$secondary <- 0
    
    #secondary
    secondary <- NULL
    for(i in 2:n_sources) {
      Y <- rbinom(n[i], 1, prob= means[i])
      df <- i
      secondary <- rbind(secondary, data.frame(Y,df))
    }
    secondary$secondary <- 1
    
    
    #data frame together
    data <- data.frame(rbind(primary, secondary))
    
    
    #dummies for source
    for(i in 2:n_sources) {
      x <- 1*(data$df==i)
      data <- cbind(data, x)
    }
    
    names(data) <- c("Y", "df", "secondary", paste("S", 2:n_sources, sep=""))
    
    return(data)
  }
  
  
  if((num_arms == 1) && (outcome_type == "continuous")){
    
    #simulate data:
    #primary
    Y <- rnorm(n[1], mean= means[1], sd=sds[1])
    primary <- data.frame(Y)
    primary$df <- 1
    primary$secondary <- 0
    
    #secondary
    secondary <- NULL
    for(i in 2:n_sources) {
      Y <- rnorm(n[i], mean= means[i], sd=sds[i])
      df <- i
      secondary <- rbind(secondary, data.frame(Y,df))
    }
    secondary$secondary <- 1
    
    #data frame together
    data <- data.frame(rbind(primary, secondary))
    
    
    #dummies for source
    for(i in 2:n_sources) {
      x <- 1*(data$df==i)
      data <- cbind(data, x)
    }
    
    names(data) <- c("Y", "df", "secondary", paste("S", 2:n_sources, sep=""))
    
    return(data)
  }
  
  
  ###TWO-ARM
  
  if((num_arms == 2) && (outcome_type == "binary")){
    #simulate data:
    #primary
    trt <- rbinom(n[1], 1, prob = 0.5) #ALWAYS 1:1 TREATMENT:CONTROL SPLIT
    Y <- rbinom(n[1], 1, prob= (means[1] + trt_effect[1]*(trt))) 
    primary <- data.frame(Y, trt)
    primary$df <- 1
    primary$secondary <- 0
    
    #secondary
    secondary <- NULL
    for(i in 2:n_sources) {
      trt <- rbinom(n[i], 1, 0.5) #ALWAYS 1:1 TREATMENT:CONTROL SPLIT
      Y <- rbinom(n[i], 1, prob= (means[i] + trt_effect[i]*(trt)))
      df <- i
      secondary <- rbind(secondary, data.frame(Y,trt,df))
    }
    secondary$secondary <- 1
    
    #data frame together
    data <- data.frame(rbind(primary, secondary))
    
    
    #dummies for source
    for(i in 2:n_sources) {
      x <- 1*(data$df==i)
      data <- cbind(data, x)
    }
    
    
    #create interactions for source*trt which will decide whether I borrow or not on the trt effect
    for(i in 2:n_sources) {
      x <- 1*(data$df==i) * data$trt
      data <- cbind(data, x)
    }
    names(data) <- c("Y", "trt", "df", "secondary", paste("S", 2:n_sources, sep=""), paste("trt_S", 2:n_sources, sep=""))
    
    data$type <- 1*(data$df==1 & data$trt==0) + 2*(data$df==1 & data$trt==1) + 3*(data$df==2 & data$trt==0) + 4*(data$df==2 & data$trt==1) + 5*(data$df==3 & data$trt==0) + 6*(data$df==3 & data$trt==1)
    
    return(data)
  }
  
  
  
  if((num_arms == 2) && (outcome_type == "continuous")){
    #simulate data:
    #primary
    trt <- rbinom(n[1], 1, 0.5)
    Y <- rnorm(n[1], mean= means[1] + trt_effect[1]*trt, sd=sds[1])
    primary <- data.frame(Y, trt)
    primary$df <- 1
    primary$secondary <- 0
    
    #secondary
    secondary <- NULL
    for(i in 2:n_sources) {
      trt <- rbinom(n[i], 1, 0.5)
      Y <- rnorm(n[i], mean= means[i] + trt_effect[i]*trt, sd=sds[i])
      df <- i
      secondary <- rbind(secondary, data.frame(Y,trt,df))
    }
    secondary$secondary <- 1
    
    #data frame together
    data <- data.frame(rbind(primary, secondary))
    
    
    #dummies for source
    for(i in 2:n_sources) {
      x <- 1*(data$df==i)
      data <- cbind(data, x)
    }
    
    #create interactions for source*trt which will decide whether I borrow or not on the trt effect
    for(i in 2:n_sources) {
      x <- 1*(data$df==i) * data$trt
      data <- cbind(data, x)
    }
    names(data) <- c("Y", "trt", "df", "secondary", paste("S", 2:n_sources, sep=""), paste("trt_S", 2:n_sources, sep=""))
    
    data$type <- 1*(data$df==1 & data$trt==0) + 2*(data$df==1 & data$trt==1) + 3*(data$df==2 & data$trt==0) + 4*(data$df==2 & data$trt==1) + 5*(data$df==3 & data$trt==0) + 6*(data$df==3 & data$trt==1)
    
    return(data)
  }
}


# # TEST FUNCTION
# SimulateData(n = c(100,100), means = c(0.1,0.2), num_arms = 1, outcome_type = "binary")
# SimulateData(n = c(100,100), means = c(0.1,100), sds = c(1,1), num_arms = 1, outcome_type = "continuous")
# SimulateData(n = c(100,100), means = c(0.1,0.1), trt_effect = c(0.1,0.5), num_arms = 2, outcome_type = "binary")
# SimulateData(n = c(100,100), means = c(0.1,15), sds = c(1,1), trt_effect = c(1,5), num_arms = 2, outcome_type = "continuous")



###########################################################
################## SECTION 3: RUN MODELS ################## 
###########################################################

RunModels <- function(data = NULL, n = NULL, means = NULL, sds = NULL, trt_effect = NULL, num_arms, n_sources = 2, outcome_type, marginal) {
  
  if((outcome_type != "binary") && (outcome_type != "continuous")){
    return("Error: outcome_type must be exactly `binary` or `continuous`")
  }
  
  #if no data set is provided, simulate the data using the SimulateData() func above
  if(is.null(data)){
    data <- SimulateData(n = n, means = means, sds = sds, 
                         trt_effect = trt_effect, num_arms = num_arms, 
                         outcome_type = outcome_type)
  }
  
  n_tot <- nrow(data)
  data$secondary <- data$S2
  data$sourcefactor <- as.factor(data$df)
  
  if((num_arms == 1) && (outcome_type == "binary")){
    #Logistic Regression Model
    lmfit <- glm(Y ~ sourcefactor, family = binomial(link = "logit"), data = data) 
    
    mainfit <- glm(Y~1, family = binomial(link = "logit"), data = data)
    standardtestresult <- lmtest::lrtest(mainfit, lmfit)
    
    standardtestpval <- standardtestresult$`Pr(>Chisq)`[2] #get p-value from LRT for flexibility in # of sources
    
    #MEM weights
    #get summary stats of simulated data
    count_g1 <- sum(data$Y[data$secondary == 0])
    count_g2 <- sum(data$Y[data$secondary == 1])
    
    wts_MEM <- calc.weights_MEM_bin(xvec = c(count_g1, count_g2),
                                    nvec=c(sum(data$S2 == 0), sum(data$S2 == 1)), 
                                    avec = c(1,1), 
                                    bvec = c(1,1),
                                    prior = "equal",
                                    constraint = 1)$q
    
    #MEMr weights
    ls <- vector("list", n_sources-1)
    for(i in 1:(n_sources-1)) {
      x<-c(TRUE, FALSE)
      ls[[i]] <- x
    }
    
    regMat <- expand.grid(ls)
    sources <- paste("S", 2:n_sources, sep="")
    
    allModelsList <- apply(regMat, 1, function(x) as.formula(
      paste(c("Y ~ 1 ",sources[x]),
            collapse=" + ")) )
    
    allModelsResults <- lapply(allModelsList,
                               function(x) glm(x, family = binomial(link = "logit"), data=data))
    
  }
  
  if((num_arms == 1) && (outcome_type == "continuous")){
    #Linear Regression Model
    lmfit <- lm(Y ~ sourcefactor, data = data) 
    
    mainfit <- lm(Y~1, data = data)
    standardtestresult <- anova(mainfit, lmfit)
    
    standardtestpval <- standardtestresult$`Pr(>F)`[2] #get p-value from F-test for flexibility in # of sources
    
    #MEM weights
    #get summary stats of simulated data
    simmean_g1 <- mean(data$Y[data$secondary == 0])
    simmean_g2 <- mean(data$Y[data$secondary == 1])
    simsd_g1 <- sd(data$Y[data$secondary == 0])
    simsd_g2 <- sd(data$Y[data$secondary == 1])
    
    wts_MEM <- calc.weights_MEM_cont(xvec = c(simmean_g1, simmean_g2),
                                     svec = c(simsd_g1, simsd_g2),
                                     nvec=c(sum(data$S2 == 0), sum(data$S2 == 1)), 
                                     prior = "pi_e")
    
    #MEMr weights
    ls <- vector("list", n_sources-1)
    for(i in 1:(n_sources-1)) {
      x<-c(TRUE, FALSE)
      ls[[i]] <- x
    }
    
    regMat <- expand.grid(ls)
    sources <- paste("S", 2:n_sources, sep="")
    
    allModelsList <- apply(regMat, 1, function(x) as.formula(
      paste(c("Y ~ 1 ",sources[x]),
            collapse=" + ")) )
    
    allModelsResults <- lapply(allModelsList,
                               function(x) lm(x, data=data))
    
  }
  
  if((num_arms == 2) && (outcome_type == "binary")){
    data$trtfactor <- as.factor(data$trt)
    
    #Logistic Regression Interaction Model
    lmfit <- glm(Y ~ trtfactor*sourcefactor, family = binomial(link = "logit"), data = data) 
    
    mainfit <- glm(Y~sourcefactor + trtfactor, family = binomial(link = "logit"), data = data)
    standardtestresult <- lrtest(mainfit, lmfit)
    
    standardtestpval <- standardtestresult$`Pr(>Chisq)`[2] #get p-value from LRT for flexibility in # of sources
    
    #MEM weights (NA for two arm)
    wts_MEM <- c(NA,NA)
    
    #MEMr weights
    #first marginal MEM: only considers main effects
    ls <- vector("list", n_sources-1)
    for(i in 1:(n_sources-1)) {
      x<-c(TRUE, FALSE)
      ls[[i]] <- x
    }
    
    regMat <- expand.grid(ls)
    regressors <- paste("trt_S", 2:n_sources, sep="")
    
    allModelsList <- apply(regMat, 1, function(x) as.formula(
      paste(c("Y ~ 1 + trt",paste("S", 2:n_sources, sep=""), regressors[x]),
            collapse=" + ")) )
    
    allModelsResults <- lapply(allModelsList,
                               function(x) glm(x, family = binomial(link = "logit"), data=data))
    
  }
  
  if((num_arms == 2) && (outcome_type == "continuous")){
    data$trtfactor <- as.factor(data$trt)
    
    #Linear Regression Interaction Model
    lmfit <- lm(Y ~ trtfactor*sourcefactor, data = data) 
    
    mainfit <- lm(Y~sourcefactor + trtfactor, data = data)
    standardtestresult <- anova(mainfit, lmfit)
    
    standardtestpval <- standardtestresult$`Pr(>F)`[2] #get p-value from F-test for flexibility in # of sources
    
    #MEM weights (NA for two arm)
    wts_MEM <- c(NA,NA)
    
    #MEMr weights
    #first marginal MEM: only considers main effects
    ls <- vector("list", n_sources-1)
    for(i in 1:(n_sources-1)) {
      x<-c(TRUE, FALSE)
      ls[[i]] <- x
    }
    
    regMat <- expand.grid(ls)
    regressors <- paste("trt_S", 2:n_sources, sep="")
    
    allModelsList <- apply(regMat, 1, function(x) as.formula(
      paste(c("Y ~ 1 + trt",paste("S", 2:n_sources, sep=""), regressors[x]),
            collapse=" + ")) )
    
    allModelsResults <- lapply(allModelsList,
                               function(x) lm(x, data=data))
    
  }
  
  
  #finish MEMr weight calculations (same for all settings from here)
  bic <- NULL
  aic <- NULL
  betas <- NULL
  for(i in 1:length(allModelsResults)) {
    fit <- allModelsResults[[i]]
    X <- as.matrix(model.matrix(fit))
    
    if(marginal == "BIC" | marginal == "AIC") {
      
      bic <- c(bic, BIC(fit))
      aic <- c(aic, AIC(fit))
      
    }
    
  }
  
  if(marginal == "BIC") {
    wts_MEMr <- exp(0.5*(bic-max(bic)))/sum(exp(0.5*(bic-max(bic))))
  }else if(marginal == "AIC") {
    wts_MEMr <- exp(0.5*(aic-max(aic)))/sum(exp(0.5*(aic-max(aic))))
  }
  
  
  #return relevant output
  return(list(n_g1 = n[1],
              n_g2 = n[2], 
              mean_g1 = means[1],
              mean_g2 = means[2],
              sd_g1 = sds[1], 
              sd_g2 = sds[2], 
              trteff_g1 = trt_effect[1],
              trteff_g2 = trt_effect[2],
              weights_MEM = c(wts_MEM[2],wts_MEM[1]), #change MEM order because it is opposite of MEMr order
                                                      #makes first weight probability of exchangeability and second nonexchangeability
              weights_MEMr = wts_MEMr,                #first weight is probability of exchangeability and second is nonexchangeability
              pval = standardtestpval))
  
}



# #TEST FUNCTION
# 
# #OA CONTINUOUS
# #run with data set
# set.seed(10)
# y1 <- as.data.frame(rnorm(n = 100, mean = 10, sd = 1))
# colnames(y1) <- "Y"
# y2 <- as.data.frame(rnorm(n = 100, mean = 25, sd = 1))
# colnames(y2) <- "Y"
# 
# testdat <- data.frame(Y = rbind(y1,y2),
#                       S2 = rep(c(0,1), each = 100),
#                       df = rep(c(0,1), each = 100) + 1)
# 
# RunModels(data = testdat, num_arms = 1, outcome_type = "continuous", marginal = "BIC")
# 
# #run without data set
# set.seed(15)
# RunModels(n = c(100,100), means = c(10,25), sds = c(1,1), num_arms = 1, outcome_type = "continuous", marginal = "BIC")
# 
# 
# #OA BINARY
# #run with data set
# set.seed(11)
# y1 <- as.data.frame(rbinom(n = 100, size = 1, prob = 0.1))
# colnames(y1) <- "Y"
# y2 <- as.data.frame(rbinom(n = 100, size = 1, prob = 0.7))
# colnames(y2) <- "Y"
# 
# testdat <- data.frame(Y = rbind(y1,y2),
#                       S2 = rep(c(0,1), each = 100),
#                       df = rep(c(0,1), each = 100) + 1)
# 
# RunModels(data = testdat, num_arms = 1, outcome_type = "binary", marginal = "BIC")
# 
# #run without data set
# set.seed(16)
# RunModels(n = c(100,100), means = c(0.1,0.7), num_arms = 1, outcome_type = "binary", marginal = "BIC")
# 
# #TA CONTINUOUS
# #run with data set
# set.seed(20)
# y1 <- as.data.frame(rnorm(n = 100, mean = 10, sd = 1))
# colnames(y1) <- "Y"
# y2 <- as.data.frame(rnorm(n = 100, mean = 12, sd = 1))
# colnames(y2) <- "Y"
# y3 <- as.data.frame(rnorm(n = 100, mean = 10, sd = 1))
# colnames(y3) <- "Y"
# y4 <- as.data.frame(rnorm(n = 100, mean = 25, sd = 1))
# colnames(y4) <- "Y"
# 
# testdat <- data.frame(Y = rbind(y1,y2,y3,y4),
#                       S2 = rep(c(0,1), each = 200),
#                       trt = rep(c(0,1,0,1), each = 100),
#                       df = rep(c(0,1), each = 200) + 1)
# testdat$trt_S2 <- testdat$trt*testdat$S2
# 
# RunModels(data = testdat, num_arms = 2, outcome_type = "continuous", marginal = "BIC")
# 
# #run without data set
# set.seed(25)
# RunModels(n = c(200,200), means = c(10,10), sds = c(1,1), trt_effect = c(2,15), num_arms = 2, outcome_type = "continuous", marginal = "BIC")
# 
# 
# #TA BINARY
# #run with data set
# set.seed(21)
# y1 <- as.data.frame(rbinom(n = 100, size = 1, prob = 0.1))
# colnames(y1) <- "Y"
# y2 <- as.data.frame(rbinom(n = 100, size = 1, prob = 0.2))
# colnames(y2) <- "Y"
# y3 <- as.data.frame(rbinom(n = 100, size = 1, prob = 0.1))
# colnames(y3) <- "Y"
# y4 <- as.data.frame(rbinom(n = 100, size = 1, prob = 0.7))
# colnames(y4) <- "Y"
# 
# testdat <- data.frame(Y = rbind(y1,y2,y3,y4),
#                       S2 = rep(c(0,1), each = 200),
#                       trt = rep(c(0,1,0,1), each = 100),
#                       df = rep(c(0,1), each = 200) + 1)
# testdat$trt_S2 <- testdat$trt*testdat$S2
# 
# RunModels(data = testdat, num_arms = 2, outcome_type = "binary", marginal = "BIC")
# 
# #run without data set
# set.seed(26)
# RunModels(n = c(200,200), means = c(0.1,0.1), trt_effect = c(0.1,0.6), num_arms = 2, outcome_type = "binary", marginal = "BIC")

