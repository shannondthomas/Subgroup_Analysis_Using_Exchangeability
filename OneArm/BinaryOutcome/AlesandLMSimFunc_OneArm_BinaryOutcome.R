#manuscript simulation code

# libraries used
library(mvtnorm)
library(ggplot2)
library(invgamma)


#simulation function
doSim_MEM1_OneArm_binary <- function(n, probs, BICiterations=NULL, marginal=NULL) {
    n_tot=sum(n)
    n_sources = length(n)
    
    #simulate data:
    #primary
    Y <- rbinom(n[1], 1, prob= probs[1]) 
    primary <- data.frame(Y)
    primary$df <- 1
    primary$secondary <- 0
    
    #secondary
    secondary <- NULL
    for(i in 2:n_sources) {
      Y <- rbinom(n[i], 1, prob= probs[i])
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
  
  
  #Logistic Regression Interaction Model
  data$sourcefactor <- as.factor(data$df)
  lmfit <- glm(Y ~ sourcefactor, family = binomial(link = "logit"), data = data) 

  mainfit <- glm(Y~1, family = binomial(link = "logit"), data = data)
  lrtestresult <- lrtest(mainfit, lmfit)
  
  lrtestpval <- lrtestresult$`Pr(>Chisq)`[2] #get p-value from LRT for flexibility in # of sources
  
  
  #first marginal MEM: only considers main effects
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
  
  bic <- NULL
  aic <- NULL
  betas <- NULL
  for(i in 1:length(allModelsResults)) {
    fit <- allModelsResults[[i]]
    X <- as.matrix(model.matrix(fit))
    
    #################EDITED THIS FOR BINARY OUTCOME############################################
    if(marginal == "BIC" | marginal == "AIC") {
        # #BIC weights:
        # #iterated weighted least squares
        # cf <- coef(fit)
        # s <- rep(NA, n_sources)
        # for(k in 1:BICiterations) {
        #   old <- cf
        #   for(j in 1:n_sources){
        #     s[j] <- sum((data$Y[data$df==j] - X[data$df==j,] %*% cf)^2)/(length(data$Y[data$df==j]))
        #   }
        #   #W <- diag(rep(s, times=n))
        #   #cf <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% data$Y
        # 
        #   cf <- coef(lm(data$Y ~ X + 0, weights = sqrt(rep(s, times=n))))
        #   if (sum(abs(old - cf)) < 0.0000000000001)
        #   {
        #     break;
        #   }
        # }
        
        #BIC:
        bic <- c(bic, BIC(fit))
        aic <- c(aic, AIC(fit))
      
    }
    
  }
  
  if(marginal == "BIC") {
    wts <- exp(0.5*(bic-max(bic)))/sum(exp(0.5*(bic-max(bic))))
  }else if(marginal == "AIC") {
    wts <- exp(0.5*(aic-max(aic)))/sum(exp(0.5*(aic-max(aic))))
  }
  
  
  return(list(sampsize = n_tot,
              pc = probs[1],
              pt = probs[2],
              weights= wts, pval = lrtestpval))
  
}


