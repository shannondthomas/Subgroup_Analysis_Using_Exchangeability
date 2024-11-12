#################################################################################
# TITLE: CaseStudies.R
#
# PURPOSE: Run all models based on summary statistics from case studies and 
#          evaluate model performance in those settings. 
#
# INPUT: summary statistics from clinical trials
#
# OUTPUT: p-values, probabilities of exchangeability, power, and t1e
#         for each case study
#
# SECTIONS: Section 0 - source functions and load packages
#           Section 1 - function to run simulations given summary stats and   
#                       summarize t1e rates and power
#           Section 2 - function to create binary data set given summary stats
#                       to match the exact counts rather than simulating 
#                       random counts (RunModels() would simulate random counts
#                       if no data was given)
#           Section 3 - run binary case study
#           Section 4 - function to create continuous data set given summary 
#                       stats (n for each group x trt), means, and sds
#                       (RunModels() would simulate random data if no data 
#                        was given)
#           Section 5 - run continuous case study
#
# AUTHOR: Shannon Thomas
# DATE CREATED: OCT 29, 2024
#################################################################################


##########################################################
#### SECTION 0: LOAD DEPENDENCIES & SOURCE MODEL FUNC ####
##########################################################

source("Test_Functions.R")
source("Simulation_Function.R")


##########################################################
####### SECTION 1:  SIMULATION RESULT SUMMARY FUNC ####### 
##########################################################

casestudy_sim_summary <- function(nsim = 10000, num_arms, outcome_type, n_sources = 2, n, events = c(NA,NA), means = c(NA,NA), sds = c(NA,NA), MEMcutoffs = c(0.2), pvalcutoffs = c(0.05)){
  
  if((num_arms == 2) & (outcome_type == "binary")){
    #run simulation to get performance in null case and given data case
    n_sim <- c(n[1]+n[2], n[3]+n[4])               #subgroup sample sizes
    means_sim <- c(events[2]/n[2], events[4]/n[4]) #placebo event rates
    trteff_sim <- c(events[1]/n[1] - events[2]/n[2], #trt effect is difference in effect size between trt and placebo groups
                    events[3]/n[3] - events[4]/n[4])
    trteff_sim_null <- rep((events[1] + events[3])/(n[1] + n[3]) - (events[2] + events[4])/(n[2] + n[4]), 2)      #assume null case is the pooled trt effect
    
    params <- list(
      list(n=n_sim, means=means_sim, sds = NA, trt_effect = trteff_sim),
      list(n=n_sim, means=means_sim, sds = NA, trt_effect = trteff_sim_null))
    
  }
  
  if((num_arms == 2) & (outcome_type == "continuous")){
    #run simulation to get performance in null case and given data case
    n_sim <- c(n[1]+n[2], n[3]+n[4])                    #subgroup sample sizes
    means_sim <- c(means[2], means[4])                  #placebo means
    sds_sim <- c(max(sds[1],sds[2]),max(sds[3],sds[4])) #define group sd as the maximum of the trt and placebo arms 
    trteff_sim <- c(means[1] - means[2], #trt effect is difference in effect size between trt and placebo groups
                    means[3] - means[4])
    trteff_sim_null <- rep((n[1]*means[1] + n[3]*means[3])/(n[1] + n[3]) - 
                             (n[2]*means[2] + n[4]*means[4])/(n[2] + n[4]), 2)  #assume null case is the pooled trt effect
    
    params <- list(
      list(n=n_sim, means=means_sim, sds = sds_sim, trt_effect = trteff_sim),
      list(n=n_sim, means=means_sim, sds = sds_sim, trt_effect = trteff_sim_null))
    
  }
  
  if((num_arms == 1) & (outcome_type == "binary")){
    #run simulation to get performance in null case and given data case
    n_sim <- n                                     #subgroup sample sizes
    means_sim <- events/n                          #event rates
    means_sim_null <- rep(sum(events)/sum(n),2)        #assume null event rate is pooled event rate
    params <- list(
      list(n=n_sim, means=means_sim, sds = c(NA,NA), trt_effect = c(NA,NA)),
      list(n=n_sim, means=means_sim_null, sds = c(NA,NA), trt_effect = c(NA,NA)))
  }
  
  if((num_arms == 1) & (outcome_type == "continuous")){
    #run simulation to get performance in null case and given data case
    means_sim_null <- rep((n[1]*means[1] + n[2]*means[2])/sum(n),2)        #assume null event rate is pooled mean
    params <- list(
      list(n=n_sim, means=means, sds = sds, trt_effect = c(NA,NA)),
      list(n=n_sim, means=means_sim_null, sds = sds, trt_effect = c(NA,NA)))
  }
  
  allsimresults <- RunSim(nsim = nsim, pars = params,num_arms = num_arms, outcome_type = outcome_type)
  
  #get summary info
  for(i in 1:length(MEMcutoffs)){
    cutoffi_summary <- allsimresults %>% 
      group_by(num_arms, outcome_type, n1, n2, mean1, mean2, sd1, sd2, trteff1, trteff2) %>%
      summarize(MEM_power = sum(MEMpexch < MEMcutoffs[i])/nsim,
                MEMr_power = sum(MEMrpexch < MEMcutoffs[i])/nsim,
                pval_power = sum(pval < pvalcutoffs[i])/nsim)
    
    cutoffi_summary$MEMcutoff <- MEMcutoffs[i]
    cutoffi_summary$pvalcutoff <- pvalcutoffs[i]
    
    if(i == 1){
      allsimresults_summary <- cutoffi_summary
    }
    else{allsimresults_summary <- rbind(allsimresults_summary, cutoffi_summary)}
  }

  if(num_arms == 2){
    MEMt1e <- allsimresults_summary$MEM_power[allsimresults_summary$trteff2 == allsimresults_summary$trteff1]
    MEMrt1e <- allsimresults_summary$MEMr_power[allsimresults_summary$trteff2 == allsimresults_summary$trteff1]
    standardtestt1e <- allsimresults_summary$pval_power[allsimresults_summary$trteff2 == allsimresults_summary$trteff1]
    
    MEMpower <- allsimresults_summary$MEM_power[allsimresults_summary$trteff2 != allsimresults_summary$trteff1]
    MEMrpower <- allsimresults_summary$MEMr_power[allsimresults_summary$trteff2 != allsimresults_summary$trteff1]
    standardtestpower <- allsimresults_summary$pval_power[allsimresults_summary$trteff2 != allsimresults_summary$trteff1]  
  }
  if(num_arms == 1){
    MEMt1e <- allsimresults_summary$MEM_power[allsimresults_summary$mean2 == allsimresults_summary$mean1]
    MEMrt1e <- allsimresults_summary$MEMr_power[allsimresults_summary$mean2 == allsimresults_summary$mean1]
    standardtestt1e <- allsimresults_summary$pval_power[allsimresults_summary$mean2 == allsimresults_summary$mean1]
    
    MEMpower <- allsimresults_summary$MEM_power[allsimresults_summary$mean2 != allsimresults_summary$mean1]
    MEMrpower <- allsimresults_summary$MEMr_power[allsimresults_summary$mean2 != allsimresults_summary$mean1]
    standardtestpower <- allsimresults_summary$pval_power[allsimresults_summary$mean2 != allsimresults_summary$mean1]  
  }
  
  return(data.frame(test = c(rep("MEM",length(MEMcutoffs)),
                             rep("MEMr",length(MEMcutoffs)),
                             rep("Standard",length(MEMcutoffs))), 
                    MEMcutoff = c(rep(MEMcutoffs,3)),
                    pvalcutoff = c(rep(pvalcutoffs,3)),
                    t1e = c(MEMt1e, MEMrt1e, standardtestt1e), 
                    power = c(MEMpower, MEMrpower, standardtestpower)))
}





##########################################################
########## SECTION 2:  BINARY DATA SET FUNCTION ########## 
##########################################################

binarydataset <- function(n, events, num_arms){
  #function to create a data set with exact counts based on given information rather than random counts
  
  #n = vector of sample sizes for each group (length = number of treatment arms x number of subgroups)
  #If two arm, n[1] = samp size of trt and subgroup level 1,
  #            n[2] = samp size of placebo and subgroup level 1,
  #            n[3] = samp size of trt and subgroup level 2,
  #            n[4] = samp size of placebo and subgroup level 2.
  #
  #If one arm, n[1] = samp size of subgroup level 1,
  #            n[2] = samp size of subgroup level 2
  
  #events = vector of event counts for each group (same length and ordering as n)
  #num_arms = number of treatment arms (1 or 2)
  
  n_tot=sum(n) #total sample size
  n_levels = length(n)/num_arms #number of subgroup levels
  
  #create data frame
  data <- data.frame(df = rep(NA, n_tot), trt = rep(NA, n_tot), Y = rep(NA, n_tot))
  
  if(num_arms == 2){
    data$df[1:(n[1] + n[2])] <- 1       #label level 1, first n[1] + n[2] 
    data$df[(n[1]+n[2]+1):n_tot] <- 2   #label level 2, last n[3] + n[4] 
    
    data$trt[1:n[1]] <- 1                                       #label treatment group, first n[1]
    data$trt[(n[1] + 1):(n[1] + n[2])] <- 0                     #label placebo group, next n[2]
    data$trt[(n[1] + n[2] + 1):(n[1] + n[2] + n[3])] <- 1       #label treatment group, next n[3]
    data$trt[(n[1] + n[2] + n[3] + 1):n_tot] <- 0               #label placebo group, last n[4]
    
    data$Y <- 0                                                 #create binary response variable
    data$Y[(data$df == 1) & (data$trt == 1)][1:events[1]] <- 1  #create events[1] events in the treatment x level 1 group
    data$Y[(data$df == 1) & (data$trt == 0)][1:events[2]] <- 1  #create events[2] events in the placebo x level 1 group
    data$Y[(data$df == 2) & (data$trt == 1)][1:events[3]] <- 1  #create events[3] events in the treatment x level 2 group
    data$Y[(data$df == 2) & (data$trt == 0)][1:events[4]] <- 1  #create events[4] events in the placebo x level 2 group
    
    data$S2 <- data$df -1                            #create source (level) 2 indicator variable 
    data$trt_S2 <- (data$df == 2) & (data$trt == 1)  #create  level 2 and treatment indicator variable
  }
  if(num_arms == 1){
    data$df[1:n[1]] <- 1          #label level 1, first n[1]
    data$df[(n[1]+1):n_tot] <- 2  #label level 2, last n[2] 
    
    data$Y <- 0                               #create binary response variable
    data$Y[(data$df == 1)][1:events[1]] <- 1  #create events[1] events in level 1 group
    data$Y[(data$df == 2)][1:events[2]] <- 1  #create events[2] events in level 2 group
    
    data$S2 <- data$df -1 #create source (level) 2 indicator variable 
  }
  
  return(data)
  
}


##########################################################
############# SECTION 3: BINARY CASE STUDIES ############# 
##########################################################


# #TESTING FOR DIFFERENCE IN PROPORTION PA+ AT WEEK 70
# #source = cycle-based OR culture-based
# #trt = ciprofloxacin
# n <- c(67,65,66,69) #n[1] is samp size for cycle based TIS-ciprofloxacin
#                     #n[2] is samp size for cycle based TIS-placebo
#                     #n[3] is samp size for culture based TIS-ciprofloxacin
#                     #n[4] is samp size for culture based TIS-placebo
# 
# events <- c(7,6,10,6)
# 
# #WEEK 46 VALUES 
# # n <- c(87,88,88,88) #n[1] is samp size for cycle based TIS-ciprofloxacin
# # #n[2] is samp size for cycle based TIS-placebo
# # #n[3] is samp size for culture based TIS-ciprofloxacin
# # #n[4] is samp size for culture based TIS-placebo
# # 
# # events <- c(4,11,12,6)
# 
# #test two-arm setting
# effsize <- (events[1]/n[1] - events[2]/n[2]) - (events[3]/n[3] - events[4]/n[4]) #estimated difference in effect size
# effsize
#   
# casestudy_sim_summary(num_arms = 2, outcome_type = "binary", n = n, events = events, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
# RunModels(num_arms = 2,
#           outcome_type = "binary",
#           marginal = "BIC",
#           data = binarydataset(n = n, events = events, num_arms = 2))
# 
# #test one-arm setting
# #pooled
# n1 <- c(n[1]+n[3], n[2] + n[4])
# events1 <- c(events[1]+events[3],events[2]+events[4])
# effsize1 <- (events[1] + events[3])/(n[1]+n[3]) - (events[2] + events[4])/(n[2]+n[4]) #estimated difference in effect size
# effsize1
# casestudy_sim_summary(num_arms = 1, outcome_type = "binary", n = n1, events = events1, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
# RunModels(num_arms = 1,
#           outcome_type = "binary",
#           marginal = "BIC",
#           data = binarydataset(n=n1, 
#                                events = events1, 
#                                num_arms = 1))
# 
# #treatment group only (NOT THIS ONE - WAY WORSE)
# n1 <- c(n[1],n[3])
# events1 <- c(events[1],events[3])
# effsize1 <- events[1]/n[1] - events[3]/n[3] #estimated difference in effect size
# effsize1
# casestudy_sim_summary(num_arms = 1, outcome_type = "binary", n = n1, events = events1, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
# RunModels(num_arms = 1,
#           outcome_type = "binary",
#           marginal = "BIC",
#           data = binarydataset(n=n1, 
#                                events = events1, 
#                                num_arms = 1))
# 

#TESTING FOR DIFFERENCE IN SAE RESPIRATORY EVENTS 
#(significant difference was found by trt/placebo groups)
#source/trt = ciprofloxacin

n <- c(152,152) #n[1] is samp size for TIS-ciprofloxacin
                #n[2] is samp size for TIS-placebo

events <- c(22,9)

#test one-arm setting
(events[1]/n[1] - events[2]/n[2])  #estimated difference in effect size
casestudy_sim_summary(num_arms = 1, outcome_type = "binary", n = n, events = events, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
RunModels(num_arms = 1,
          outcome_type = "binary",
          marginal = "BIC",
          data = binarydataset(n = n, events = events, num_arms = 1))




##########################################################
######## SECTION 4:  CONTINUOUS DATA SET FUNCTION ######## 
##########################################################

continuousdataset <- function(n, means, sds, num_arms){
  #function to create data set with summary statistics matching those given
  
  #n = vector of sample sizes for each group (length = number of treatment arms x number of subgroups)
  #If two arm, n[1] = samp size of trt and subgroup level 1,
  #            n[2] = samp size of placebo and subgroup level 1,
  #            n[3] = samp size of trt and subgroup level 2,
  #            n[4] = samp size of placebo and subgroup level 2.
  #
  #If one arm, n[1] = samp size of subgroup level 1,
  #            n[2] = samp size of subgroup level 2
  
  #means: vector of means for each subgroup x treatment group (same ordering as n above)
  #sd: vector of sds for each subgroup x treatment group (same ordering as n above)
  
  tolerance <- 0.1
  
  n_sources <- 2
  
  if(num_arms == 2){
    truetrteff <- (means[1] - means[2]) - (means[3] - means[4])
    simtrteff <- rep(1000000,2)
    
  }
  if(num_arms == 1){
    truetrteff <- means[1] - means[2]
    simtrteff <- 1000000
  }
  
  simmeans <- rep(1000000, num_arms*2)
  simsds <- rep(1000000, num_arms*2)
  
  while((sum(abs(simtrteff - truetrteff)) > tolerance) || 
        (sum(abs(simmeans - means)) > tolerance) || 
        (sum(abs(simsds - sds)) > tolerance)){
  
        ###ONE-ARM
        if((num_arms == 1)){
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
        }
        
        
        ###TWO-ARM
        if((num_arms == 2)){
          #simulate data:
          #primary
          trt <- rep(c(1, 0), times = c(n[1],n[2]))
          Y <- rnorm(n[1]+n[2], mean= means[1]*trt + means[2]*(1-trt), sd=sds[1]*trt + sds[2]*(1-trt))
          primary <- data.frame(Y, trt)
          primary$df <- 1
          primary$secondary <- 0
          
          #secondary
          secondary <- NULL
            trt <- rep(c(1, 0), times = c(n[3],n[4]))
            Y <- rnorm(n[3]+n[4], mean= means[3]*trt + means[4]*(1-trt), sd=sds[3]*trt + sds[4]*(1-trt))
            df <- 2
            secondary <- rbind(secondary, data.frame(Y,trt,df))

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
          
         }
    
        #recalculate summary stats 
        if(num_arms == 2){
          simmeans <- c(mean(data$Y[(data$trt == 1) & (data$secondary == 0)]),
                        mean(data$Y[(data$trt == 0) & (data$secondary == 0)]),
                        mean(data$Y[(data$trt == 1) & (data$secondary == 1)]),
                        mean(data$Y[(data$trt == 0) & (data$secondary == 1)]))
          
          simsds <- c(sd(data$Y[(data$trt == 1) & (data$secondary == 0)]),
                      sd(data$Y[(data$trt == 0) & (data$secondary == 0)]),
                      sd(data$Y[(data$trt == 1) & (data$secondary == 1)]),
                      sd(data$Y[(data$trt == 0) & (data$secondary == 1)]))
          
          simtrteff <- (mean(data$Y[(data$trt == 1) & (data$secondary == 0)]) - mean(data$Y[(data$trt == 0) & (data$secondary == 0)])) -
                         (mean(data$Y[(data$trt == 1) & (data$secondary == 1)]) - mean(data$Y[(data$trt == 0) & (data$secondary == 1)]))
        }
        if(num_arms == 1){
          simmeans <- c(mean(data$Y[data$secondary == 0]),
                        mean(data$Y[data$secondary == 1]))
          
          simsds <- c(sd(data$Y[data$secondary == 0]),
                      sd(data$Y[data$secondary == 1])) 
          
          simtrteff <- simmeans[1] - simmeans[2]
        }
        # print(paste("difference in trt eff = ", sum(abs(simtrteff - truetrteff))))
        # print(paste("difference in means = ", sum(abs(simmeans - means))))
        # print(paste("difference in sds = ",sum(abs( simsds - sds))))
    
    }
    
    return(data)
  
}

# #test function
# set.seed(10)
# test1 <- continuousdataset(n = c(100,150,200,250), means = c(10,3,4,5), sds = c(1,1.33,1.4,1.5), num_arms = 2)
# set.seed(11)
# test2 <- continuousdataset(n = c(100,150), means = c(10,3), sds = c(12,1.3), num_arms = 1)


##########################################################
########### SECTION 5: CONTINUOUS CASE STUDIES ########### 
##########################################################

#TRIKAFTA
n <- c(44,38,141,150)
means <- c(11.3,0.0,14.4,-0.2)
sds <- c((13.6-11.3)/qnorm(0.975),2.5/qnorm(0.975), (14.4-13)/qnorm(0.975), (1.1+0.2)/qnorm(0.975))

ediff1 <- means[1] - means[2]
ediff2 <- means[3] - means[4]

casestudy_sim_summary(nsim = 10000, outcome_type = "continuous", num_arms = 2, n = n, means = means, sds = sds, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))

set.seed(15)
RunModels(num_arms = 2,
          outcome_type = "continuous",
          marginal = "BIC",
          data = continuousdataset(n = n, means = means, sds = sds, num_arms = 2))

