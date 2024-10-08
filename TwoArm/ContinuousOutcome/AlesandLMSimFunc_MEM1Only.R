#manuscript simulation code

# libraries used
library(mvtnorm)
library(ggplot2)
library(invgamma)


#simulation function

#n = sample sizes, m0 = baseline mean for each source, beta_z = covariate effects on outcome,
#p = probability of treatment being assigned, beta_t = treatment main effect, 
#beta_m = interaction term coefficient (if = 0, then no trt eff modifier), muz = mean for each covariate

doSim_MEM1 <- function(n, mu0, p, beta_t, sd, BICiterations=NULL, marginal=NULL) {
  n_tot=sum(n)
  n_sources = length(n)
  
  #simulate data:
  #primary
  trt <- rbinom(n[1], 1, p[1])
  Y <- rnorm(n[1], mean= mu0[1] + beta_t[1]*trt, sd=sd[1])
  primary <- data.frame(Y, trt)
  primary$df <- 1
  primary$secondary <- 0
  
  #secondary
  secondary <- NULL
  for(i in 2:n_sources) {
    trt <- rbinom(n[i], 1, p[i])
    Y <- rnorm(n[i], mean= mu0[i] + beta_t[i]*trt, sd=sd[i])
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
  
  
  
  #Linear Regression Interaction Model
  data$sourcefactor <- as.factor(data$df)
  data$trtfactor <- as.factor(data$trt)
  lmfit <- lm(Y ~ trtfactor*sourcefactor, data = data) 
  #pval <- summary(lmfit)$coefficients[4,4] #extract p-value for interaction term

  mainfit <- lm(Y~sourcefactor + trtfactor, data = data)
  ftestresult <- anova(mainfit, lmfit)
  
  ftestpval <- ftestresult$`Pr(>F)`[2] #get p-value from F-test for flexibility in # of sources
  
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
    
    bic <- NULL
    aic <- NULL
    betas <- NULL
    for(i in 1:length(allModelsResults)) {
      fit <- allModelsResults[[i]]
      X <- as.matrix(model.matrix(fit))

      if(marginal == "BIC" | marginal == "AIC") {
        #BIC weights:
        #iterated weighted least squares
        cf <- coef(fit)
        s <- rep(NA, n_sources)
        for(k in 1:BICiterations) {
          old <- cf
          for(j in 1:n_sources){
            s[j] <- sum((data$Y[data$df==j] - X[data$df==j,] %*% cf)^2)/(length(data$Y[data$df==j]))
          }
          #W <- diag(rep(s, times=n))
          #cf <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% data$Y
          
          cf <- coef(lm(data$Y ~ X + 0, weights = sqrt(rep(s, times=n))))
          if (sum(abs(old - cf)) < 0.0000000000001)
          {
            break;
          }
        }
        
        #BIC:
        bic <- c(bic, - log(nrow(data))*(ncol(X)+n_sources) + 2*dmvnorm(data$Y, mean= X %*% cf, sigma = diag(rep(s, times=n)), log=TRUE))
        aic <- c(aic, - 2*(ncol(X)+n_sources) + 2*dmvnorm(data$Y, mean= X %*% cf, sigma = diag(rep(s, times=n)), log=TRUE))
        
      }
      
    }
    
    if(marginal == "BIC") {
      wts <- exp(0.5*(bic-max(bic)))/sum(exp(0.5*(bic-max(bic))))
    }else if(marginal == "AIC") {
      wts <- exp(0.5*(aic-max(aic)))/sum(exp(0.5*(aic-max(aic))))
    }
    
  
  return(list(sampsize = n_tot, te = beta_t,
              weights= wts, pval = ftestpval))
  
}


