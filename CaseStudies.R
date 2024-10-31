#################################################################################
# TITLE: CaseStudies.R
#
# PURPOSE: Run all models with specific binary results for case studies. 
#
# INPUT: data from clinical trials
#
# OUTPUT: p-values and probabilities of exchangeability for each case study
#
# SECTIONS: Section 0 - source functions and load packages
#           Section 1 - function to run simulations given summary stats and   
#                       summarize t1e rates and power
#           Section 2 - function to create binary data set given summary stats
#                       to match the exact counts rather than simulating 
#                       random counts (RunModels() would simulate random counts
#                       if no data was given)
#           Section 3 - run case studies
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
################ SECTION 3:  CASE STUDIES ################ 
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



#TESTING FOR DIFFERENCE IN THOSE WITH NO PA+ CULTURES VS THOSE WITH ANY (significant difference was found by cycled/culture-based groups)
#source/trt = cycled vs culture
n <- c(152,152) #n[1] is samp size for cycled
                #n[2] is samp size for culture-based

events <- c(109,85)

#test one-arm setting
(events[1]/n[1] - events[2]/n[2])  #estimated difference in effect size
casestudy_sim_summary(nsim = 1000, num_arms = 1, outcome_type = "binary", n = n, events = events, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
RunModels(num_arms = 1,
          outcome_type = "binary",
          marginal = "BIC",
          data = binarydataset(n = n, events = events, num_arms = 1))










# #ALBIGLUTIDE (MAY WORK BUT REALLY LOW POWER)
# #TESTING FOR DIFFERENCE IN PROPORTION WITH COMPOSITE ENDPOINT BY SMOKING STATUS
# #source = current smoker vs never/former smoker
# #trt = albiglutide
# n <- c(737,751,2083+1910,1999+1981) #n[1] is samp size for current smoker + albiglutide
#                                     #n[2] is samp size for current smoker + placebo
#                                     #n[3] is samp size for nonsmoker (former or never) + albiglutide
#                                     #n[4] is samp size for nonsmoker (former or never) + placebo
# 
# events <- c(72,63,159+107,210+155)
# 
# 
# effsize <- (events[1]/n[1] - events[2]/n[2]) - (events[3]/n[3] - events[4]/n[4]) #estimated difference in effect size
# effsize
# 
# casestudy_sim_summary(nsim = 1000, num_arms = 2, outcome_type = "binary", n = n, events = events, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
# RunModels(num_arms = 2,
#           outcome_type = "binary",
#           marginal = "BIC",
#           data = binarydataset(n = n, events = events, num_arms = 2))
# 
# #ALBIGLUTIDE (MAY WORK BUT REALLT LOW POWER)
# #TESTING FOR DIFFERENCE IN PROPORTION WITH COMPOSITE ENDPOINT BY NUMBER OF ARTERIAL BEDS
# #source = 1 arterial bed vs 2/3 arterial beds
# #trt = albiglutide
# n <- c(3846,3856,877,865) #n[1] is samp size for pt with 1 affected bed + albiglutide
#                           #n[2] is samp size for pt with 1 affected bed + placebo
#                           #n[3] is samp size for pt with 2 or 3 affected beds + albiglutide
#                           #n[4] is samp size for pt with 2 or 3 affected beds + placebo
# 
# events <- c(209,301,129,127)
# 
# 
# effsize <- (events[1]/n[1] - events[2]/n[2]) - (events[3]/n[3] - events[4]/n[4]) #estimated difference in effect size
# effsize
# 
# casestudy_sim_summary(nsim = 1000, num_arms = 2, outcome_type = "binary", n = n, events = events, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
# RunModels(num_arms = 2,
#           outcome_type = "binary",
#           marginal = "BIC",
#           data = binarydataset(n = n, events = events, num_arms = 2))
# 
# 
# 
# 
# #PAOLA-1 (MAY WORK)
# n <- c(255,132,282,137) #n[1] is samp size for HRD positve + trt
#                         #n[2] is samp size for HRD positive + placebo
#                         #n[3] is samp size for HRD neg/unknown + trt
#                         #n[4] is samp size for HRD neg/unknown + placebo
# 
# events <- c(136,104,230,118)
# 
# 
# effsize <- (events[1]/n[1] - events[2]/n[2]) - (events[3]/n[3] - events[4]/n[4]) #estimated difference in effect size
# effsize
# 
# casestudy_sim_summary(num_arms = 2, outcome_type = "binary", n = n, events = events, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
# RunModels(num_arms = 2,
#           outcome_type = "binary",
#           marginal = "BIC",
#           data = binarydataset(n = n, events = events, num_arms = 2))
# 
# 
# 
# #ISIS-2 (MAY WORK BUT REALLY LOW POWER)
# n <- c(1357,1442,7228,7157) #n[1] is samp size for gemini/libra + aspirin
# #n[2] is samp size for gemini/libra + placebo
# #n[3] is samp size for other + aspirin
# #n[4] is samp size for other + placebo
# 
# events <- c(150,147,654,868)
# 
# 
# effsize <- (events[1]/n[1] - events[2]/n[2]) - (events[3]/n[3] - events[4]/n[4]) #estimated difference in effect size
# effsize
# 
# casestudy_sim_summary(nsim = 1000, num_arms = 2, outcome_type = "binary", n = n, events = events, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
# RunModels(num_arms = 2,
#           outcome_type = "binary",
#           marginal = "BIC",
#           data = binarydataset(n = n, events = events, num_arms = 2))
# 
# 
# #ISIS-2 (MAY WORK BUT REALLY LOW POWER)
# n <- c(1357,1442,7228,7157) #n[1] is samp size for gemini/libra + aspirin
# #n[2] is samp size for gemini/libra + placebo
# #n[3] is samp size for other + aspirin
# #n[4] is samp size for other + placebo
# 
# events <- c(150,147,654,868)
# 
# 
# effsize <- (events[1]/n[1] - events[2]/n[2]) - (events[3]/n[3] - events[4]/n[4]) #estimated difference in effect size
# effsize
# 
# casestudy_sim_summary(nsim = 1000, num_arms = 2, outcome_type = "binary", n = n, events = events, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
# RunModels(num_arms = 2,
#           outcome_type = "binary",
#           marginal = "BIC",
#           data = binarydataset(n = n, events = events, num_arms = 2))
# 
# 
# #WEGOVY (MAY WORK BUT REALLY LOW POWER)
# n <- c(1082,1016,832,813) #n[1] is samp size for BMI < 35 + semaglutide
# #n[2] is samp size for BMI < 35 + placebo
# #n[3] is samp size for BMI >= 35 + semaglutide
# #n[4] is samp size for BMI >= 35 + placebo
# 
# events <- c(60,59,43,79)
# 
# 
# effsize <- (events[1]/n[1] - events[2]/n[2]) - (events[3]/n[3] - events[4]/n[4]) #estimated difference in effect size
# effsize
# 
# casestudy_sim_summary(nsim = 1000, num_arms = 2, outcome_type = "binary", n = n, events = events, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
# RunModels(num_arms = 2,
#           outcome_type = "binary",
#           marginal = "BIC",
#           data = binarydataset(n = n, events = events, num_arms = 2))
# 
# 
# 
# #IMpassion 130 X (NOPE)
# n <- c(185,184,266,267) 
# 
# events <- c(94,110,161,169)
# 
# 
# effsize <- (events[1]/n[1] - events[2]/n[2]) - (events[3]/n[3] - events[4]/n[4]) #estimated difference in effect size
# effsize
# 
# 
# #Keynote 042 (NOPE)
# n <- c(299,300,338,337)
# events <- c(142,101,124,98)
# 
# 
# effsize <- (events[1]/n[1] - events[2]/n[2]) - (events[3]/n[3] - events[4]/n[4]) #estimated difference in effect size
# effsize
# 
# 
# #mammograms (HAVENT RUN YET BUT PROBABLY NOT)
# n <- c(14842,7103,25476,12840)
# events <- c(24,12,42,33)
# 
# 
# effsize <- (events[1]/n[1] - events[2]/n[2]) - (events[3]/n[3] - events[4]/n[4]) #estimated difference in effect size
# effsize
# casestudy_sim_summary(nsim = 1000, num_arms = 2, outcome_type = "binary",n = n, events = events, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
# RunModels(num_arms = 2,
#           outcome_type = "binary",
#           marginal = "BIC",
#           data = binarydataset(n = n, events = events, num_arms = 2))
# 
# 
# #ischemic vs nonischemic amlodipine
# n <- c(362,370,209,212)
# events <- c(123,126,37,66)
# 
# 
# effsize <- (events[1]/n[1] - events[2]/n[2]) - (events[3]/n[3] - events[4]/n[4]) #estimated difference in effect size
# effsize
# casestudy_sim_summary(nsim = 1000, outcome_type = "binary", num_arms = 2, n = n, events = events, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))
# RunModels(num_arms = 2,
#           outcome_type = "binary",
#           marginal = "BIC",
#           data = binarydataset(n = n, events = events, num_arms = 2))
# 


#TRIKAFTA
n <- c(44,38,141,150)
means <- c(11.3,0.0,14.4,-0.2)
sds <- c((13.6-11.3)/qnorm(0.975),2.5/qnorm(0.975), (14.4-13)/qnorm(0.975), (1.1+0.2)/qnorm(0.975))

ediff1 <- means[1] - means[2]
ediff2 <- means[3] - means[4]

casestudy_sim_summary(nsim = 1000, outcome_type = "continuous", num_arms = 2, n = n, means = means, sds = sds, MEMcutoffs = c(0.2,0.8), pvalcutoffs = c(0.05,0.05))

set.seed(15)
RunModels(n = c(n[1]+n[2],n[3]+n[4]), means = c(means[2],means[4]), 
          sds = c(max(sds[1],sds[2]), max(sds[3],sds[4])), 
          trt_effect = c(ediff1, ediff2), num_arms = 2, outcome_type = "continuous", marginal = "BIC")

