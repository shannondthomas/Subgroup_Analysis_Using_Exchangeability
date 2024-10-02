calc.weights_MEM = function(xvec,svec,nvec,prior){
  ###function to calculate model weights for MEM approach with "correct" calculations which don't assume conditional independence
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