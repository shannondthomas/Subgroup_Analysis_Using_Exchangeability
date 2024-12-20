#################################################################################
# TITLE: RecommendedUse.R
#
# PURPOSE: Create Recommended Use Case Tables
#
#
# INPUT: simulation_results_summary.csv - contains information for every simulation
#                                      setting including means, number of arms,
#                                      outcome type, and power for each test,
#                                      and all other relevant information
#
#        simulation_results_summary_08cutoff.csv - summary of power for each 
#                                                    test for a 0.8 MEM cutoff 
#                                                    and 0.05 pval
#        simulation_results_unequaln_summary.csv - same as above with unequal n
#
#        simulation_results_unequaln_summary_08cutoff.csv - same as above with 
#                                                           unequal n
#
# OUTPUT: Recommended Use Case Tables
#
#
#
# AUTHOR: Shannon Thomas
# DATE CREATED: OCT 24, 2024
#################################################################################


#load in packages
#R version 4.3.2 was used for this project.
library(tidyverse) #version 2.0.0
library(flextable) #version 0.9.6


#read in summary data
allresults_summary <- read.csv("simulation_results_summary.csv")
allresults_summary_08 <- read.csv("simulation_results_summary_08cutoff.csv")
allresults_summary_calibrated <- read.csv("simulation_results_summary_calibratedcutoff.csv")
allresults_unequaln_summary <- read.csv("simulation_results_unequaln_summary.csv")
allresults_unequaln_summary_08 <- read.csv("simulation_results_unequaln_summary_08cutoff.csv")
allresults_unequaln_summary_calibrated <- read.csv("simulation_results_unequaln_summary_calibratedcutoff.csv")
allresults_summary <- rbind(allresults_summary, allresults_summary_08, 
                            allresults_summary_calibrated,
                            allresults_unequaln_summary, allresults_unequaln_summary_08, 
                            allresults_unequaln_summary_calibrated) %>% subset(select = c(-X))

#rearrange data to define confirmatory and exploratory use case settings

#pull t1e data
allresults_summary_t1e <- allresults_summary[((allresults_summary$mean2 - allresults_summary$mean1 == 0) & (is.na(allresults_summary$trteff1))) |
                                               ((allresults_summary$trteff1 - allresults_summary$trteff2 == 0) & !is.na(allresults_summary$trteff1)) ,]

colnames(allresults_summary_t1e)[(ncol(allresults_summary_t1e) - 6):(ncol(allresults_summary_t1e)-4)] <- c("MEM_t1e", "MEMr_t1e","pval_t1e")

allresults_summary_t1e <- allresults_summary_t1e %>% subset(select = c(-mean1, -mean2, -trteff1, -trteff2))

#pull power data
allresults_summary_power <- allresults_summary[((allresults_summary$mean2 - allresults_summary$mean1 != 0) & (is.na(allresults_summary$trteff1))) |
                                                 ((allresults_summary$trteff1 - allresults_summary$trteff2 != 0) & (!is.na(allresults_summary$trteff1))) ,]

#combine to have t1e on every row
allresults_summary_wide <- merge(allresults_summary_t1e, allresults_summary_power, 
                                 by = c("MEMcutoff","pvalcutoff","num_arms","outcome_type","n1","n2","sd1","sd2"))

#define confirmatory indicator
#MEM t1e < standard test t1e AND power within 0.05
allresults_summary_wide$MEMconfirmatory <- ifelse((allresults_summary_wide$MEM_t1e < allresults_summary_wide$pval_t1e) &
                                                    (allresults_summary_wide$pval_power - allresults_summary_wide$MEM_power < 0.05), 1, 0)

allresults_summary_wide$MEMrconfirmatory <- ifelse((allresults_summary_wide$MEMr_t1e < allresults_summary_wide$pval_t1e) &
                                                     (allresults_summary_wide$pval_power - allresults_summary_wide$MEMr_power < 0.05), 1, 0)

#define exploratory indicator
#MEM t1e < 0.15 t1e AND power at least 0.1 higher than standard test
allresults_summary_wide$MEMexploratory <- ifelse((allresults_summary_wide$MEM_t1e < 0.15) &
                                                   (allresults_summary_wide$MEM_power - allresults_summary_wide$pval_power >= 0.1), 1, 0)

allresults_summary_wide$MEMrexploratory <- ifelse((allresults_summary_wide$MEMr_t1e < 0.15) &
                                                    (allresults_summary_wide$MEMr_power - allresults_summary_wide$pval_power >= 0.1), 1, 0)


#generate summary by MEMcutoff, num_arms, and outcome_type
allresults_summary_wide$Variance <- allresults_summary_wide$sd1^2
allresults_summary_wide$N <- allresults_summary_wide$n1 + allresults_summary_wide$n2
allresults_summary_wide$effect_size <-ifelse(allresults_summary_wide$num_arms == 1, allresults_summary_wide$mean2 - allresults_summary_wide$mean1,
                                             allresults_summary_wide$trteff2 - allresults_summary_wide$trteff1)

allresults_summary_wide$n2_mult <- allresults_summary_wide$n2/allresults_summary_wide$N

recommendeduse_summary <- function(name, memtype, recommendeduse){
  sum1 <- allresults_summary_wide[allresults_summary_wide[[name]] == 1,] %>% 
            group_by(n2_mult, MEMcutoff, num_arms, outcome_type, Variance, N) %>%
            summarize(min_effsize = min(effect_size),
                      max_effsize = max(effect_size))
  
  sum2 <- sum1 %>%
            group_by(n2_mult, MEMcutoff, num_arms, outcome_type, Variance, min_effsize, max_effsize) %>%
            summarize(min_N = min(N),
                      max_N = max(N))
  
  sum2 <- sum2[complete.cases(sum2 %>% subset(select = c(-Variance))),]
  
  sum2$MEMtype <- memtype
  sum2$recommendeduse <- recommendeduse
  
  return(sum2)
}

MEM_confirmatory <- recommendeduse_summary("MEMconfirmatory", "MEM", "Confirmatory")
MEMr_confirmatory <- recommendeduse_summary("MEMrconfirmatory", "MEMr", "Confirmatory")
MEM_exploratory <- recommendeduse_summary("MEMexploratory", "MEM", "Exploratory")
MEMr_exploratory <- recommendeduse_summary("MEMrexploratory", "MEMr", "Exploratory")

allrecommendeduse <- rbind(MEM_confirmatory, MEMr_confirmatory, MEM_exploratory, MEMr_exploratory) %>%
                        arrange(n2_mult, recommendeduse, MEMcutoff, num_arms, outcome_type, MEMtype, Variance) 

allrecommendeduse$effsize_range <- ifelse(allrecommendeduse$min_effsize != allrecommendeduse$max_effsize,
                                            paste(allrecommendeduse$min_effsize, "-", allrecommendeduse$max_effsize),
                                            allrecommendeduse$min_effsize)
allrecommendeduse$N_range <- ifelse(allrecommendeduse$min_N != allrecommendeduse$max_N,
                                    paste(allrecommendeduse$min_N, "-", allrecommendeduse$max_N),
                                    allrecommendeduse$min_N)

#confirmatory cases only have a minimum effect size requirement because the methods both converge to 1
allrecommendeduse$effsize_range[allrecommendeduse$recommendeduse == "Confirmatory"] <- paste("\u2265", allrecommendeduse$min_effsize[allrecommendeduse$recommendeduse == "Confirmatory"])

allrecommendeduse <- allrecommendeduse %>% 
                      subset(select = c(-min_N, -max_N, -min_effsize, -max_effsize)) %>%
                      select(n2_mult, recommendeduse, MEMcutoff, num_arms, outcome_type, MEMtype, Variance, N_range, effsize_range) %>%
                      arrange(desc(n2_mult))


colnames(allrecommendeduse) <- c("Sample \nSize Dist","Recommended Use","MEM \nCutoff", "Number \nof Arms", "Outcome \nType", 
                                 "Model", "Variance", "Sample Size", "Effect Size")

allrecommendeduse$Variance[is.na(allrecommendeduse$Variance)] <- "--"

ft1 <- flextable(allrecommendeduse)
ft1
ft2 <- merge_v(ft1, j = c("Sample \nSize Dist","Recommended Use","MEM \nCutoff"))
ft2 %>% theme_box()

##NEED TO UPDATE THIS AFTER RERUNNING
ft3 <- ft2 %>% theme_box() %>%
  merge_at(i = 1:10, j = 4) %>% #0.5 n2mult, confirmatory, 0.2 cutoff, 1 arm section
    merge_at(i = 1:4, j = 5) %>% 
      merge_at(i = 1:2, j = 6) %>% 
        merge_at(i = 1:2, j = 7) %>% 
      merge_at(i = 3:4, j = 6) %>% 
        merge_at(i = 3:4, j = 7) %>% 
    merge_at(i = 5:10, j = 5) %>%
      merge_at(i = 5:6, j = 6) %>%
      merge_at(i = 7:10, j = 6) %>%
       merge_at(i = 8:9, j = 7) %>%
  merge_at(i = 11:17, j = 4) %>% #0.5 n2mult, confirmatory, 0.2 cutoff, 2 arm section
    merge_at(i = 11:12, j = 5) %>% 
    merge_at(i = 13:17, j = 5) %>% 
      merge_at(i = 13:17, j = 6) %>%
        merge_at(i = 13:14, j = 7) %>%
        merge_at(i = 15:17, j = 7) %>%
  merge_at(i = 18:21, j = 4) %>% #0.5 n2mult, confirmatory, 0.8 cutoff, 1 arm section
    merge_at(i = 19:21, j = 5) %>%
      merge_at(i = 19:21, j = 6) %>%
  merge_at(i = 22:25, j = 4) %>% #0.5 n2mult, confirmatory, 0.8 cutoff, 2 arm section
    merge_at(i = 23:25, j = 5) %>%
      merge_at(i = 23:25, j = 6) %>%
  merge_at(i = 26:33, j = 4) %>% #0.5 n2mult, confirmatory, calibrated cutoff, 1 arm section
    merge_at(i = 26:27, j = 5) %>%
    merge_at(i = 28:33, j = 5) %>%
      merge_at(i = 28:30, j = 6) %>%
      merge_at(i = 31:33, j = 6) %>%
  merge_at(i = 34:41, j = 4) %>% #0.5 n2mult, confirmatory, calibrated cutoff, 2 arm section 
    merge_at(i = 34:35, j = 5) %>%
    merge_at(i = 36:41, j = 5) %>%
      merge_at(i = 36:38, j = 6) %>%
      merge_at(i = 39:41, j = 6) %>%
  merge_at(i = 43:45, j = 4) %>% #0.5 n2mult, exploratory, 0.2 cutoff, 2 arm section
    merge_at(i = 43:45, j = 5) %>%
      merge_at(i = 43:45, j = 6) %>%
        merge_at(i = 43:45, j = 7) %>%
  merge_at(i = 46:50, j = 4) %>% #0.5 n2mult, exploratory, 0.8 cutoff, 1 arm section
    merge_at(i = 46:47, j = 5) %>%
    merge_at(i = 48:50, j = 5) %>%
      merge_at(i = 48:50, j = 6) %>%
        merge_at(i = 49:50, j = 7) %>%
  merge_at(i = 51:54, j = 4) %>% #0.5 n2mult, exploratory, 0.8 cutoff, 2 arm section
    merge_at(i = 52:54, j = 5) %>%
      merge_at(i = 52:54, j = 6) %>%
        merge_at(i = 53:54, j = 7) %>%
  merge_at(i = 55:58, j = 4) %>% #0.5 n2mult, exploratory, calibrated cutoff, 2 arm section
    merge_at(i = 55:58, j = 5) %>% 
      merge_at(i = 55:58, j = 6) %>% 
        merge_at(i = 55:58, j = 7) %>% 
  
  merge_at(i = 59:64, j = 4) %>% #0.25 n2mult, confirmatory, 0.2 cutoff, 1 arm section
    merge_at(i = 59:62, j = 5) %>%
      merge_at(i = 59:60, j = 6) %>%
        merge_at(i = 59:60, j = 7) %>%
      merge_at(i = 61:62, j = 6) %>%
        merge_at(i = 61:62, j = 7) %>%
    merge_at(i = 63:64, j = 5) %>%
      merge_at(i = 63:64, j = 6) %>%
        merge_at(i = 63:64, j = 7) %>%
  merge_at(i = 65:66, j = 4) %>% #0.25 n2mult, confirmatory, 0.2 cutoff, 2 arm section
    merge_at(i = 65:66, j = 5) %>%
      merge_at(i = 65:66, j = 6) %>%
        merge_at(i = 65:66, j = 7) %>%
  merge_at(i = 67:68, j = 4) %>% #0.25 n2mult, confirmatory, 0.8 cutoff, 1 arm section
  merge_at(i = 69:70, j = 4) %>% #0.25 n2mult, confirmatory, 0.8 cutoff, 2 arm section
  merge_at(i = 71:74, j = 4) %>% #0.25 n2mult, confirmatory, calibrated cutoff, 1 arm section
    merge_at(i = 71:72, j = 5) %>%
    merge_at(i = 73:74, j = 5) %>%
      merge_at(i = 73:74, j = 6) %>%
        merge_at(i = 73:74, j = 7) %>%
  merge_at(i = 75:79, j = 4) %>% #0.25 n2mult, confirmatory, calibrated cutoff, 2 arm section
    merge_at(i = 75:77, j = 5) %>%
      merge_at(i = 76:77, j = 6) %>%
        merge_at(i = 76:77, j = 7) %>%
    merge_at(i = 78:79, j = 5) %>%
  merge_at(i = 81:84, j = 4) %>% #0.25 n2mult, exploratory, 0.2 cutoff, 2 arm section
    merge_at(i = 81:84, j = 5) %>%
      merge_at(i = 81:84, j = 6) %>%
        merge_at(i = 81:84, j = 7) %>%
  merge_at(i = 85:87, j = 4) %>% #0.25 n2mult, exploratory, 0.8 cutoff, 1 arm section
    merge_at(i = 85:86, j = 5) %>%
  merge_at(i = 88:89, j = 4) %>% #0.25 n2mult, exploratory, 0.8 cutoff, 2 arm section
  merge_at(i = 90:93, j = 4) %>% #0.25 n2mult, exploratory, calibrated cutoff, 2 arm section
    merge_at(i = 90:93, j = 5) %>%
      merge_at(i = 90:93, j = 6) %>%
        merge_at(i = 90:93, j = 7) %>%
  
  merge_at(i = 94:102, j = 4) %>% #0.1 n2mult, confirmatory, 0.2 cutoff, 1 arm section
    merge_at(i = 94:99, j = 5) %>%
      merge_at(i = 94:96, j = 6) %>%
        merge_at(i = 94:96, j = 7) %>%
      merge_at(i = 97:99, j = 6) %>%
        merge_at(i = 97:99, j = 7) %>%
    merge_at(i = 100:102, j = 5) %>%
      merge_at(i = 100:102, j = 6) %>%
        merge_at(i = 100:102, j = 7) %>%
  merge_at(i = 104:105, j = 4) %>% #0.1 n2mult, confirmatory, 0.8 cutoff, 1 arm section
  merge_at(i = 106:107, j = 4) %>% #0.1 n2mult, confirmatory, 0.8 cutoff, 2 arm section
  merge_at(i = 108:112, j = 4) %>% #0.1 n2mult, confirmatory, calibrated cutoff, 1 arm section
    merge_at(i = 108:110, j = 5) %>%
      merge_at(i = 109:110, j = 6) %>%
        merge_at(i = 109:110, j = 7) %>%
    merge_at(i = 111:112, j = 5) %>%
  merge_at(i = 113:117, j = 4) %>% #0.1 n2mult, confirmatory, calibrated cutoff, 2 arm section
    merge_at(i = 113:115, j = 5) %>%
      merge_at(i = 114:115, j = 6) %>%
        merge_at(i = 114:115, j = 7) %>%
    merge_at(i = 116:117, j = 5) %>%
  merge_at(i = 119:121, j = 4) %>% #0.1 n2mult, exploratory, 0.2 cutoff, 2 arm section
    merge_at(i = 119:121, j = 5) %>%
      merge_at(i = 119:121, j = 6) %>%
        merge_at(i = 119:121, j = 7) %>%
  merge_at(i = 122:123, j = 4) %>% #0.1 n2mult, exploratory, 0.8 cutoff, 1 arm section
  merge_at(i = 124:126, j = 4) %>% #0.1 n2mult, exploratory, 0.8 cutoff, 2 arm section
    merge_at(i = 125:126, j = 5) %>%
      merge_at(i = 125:126, j = 6) %>%
        merge_at(i = 125:126, j = 7) %>%
  merge_at(i = 127:130, j = 4) %>% #0.1 n2mult, exploratory, calibrated cutoff, 2 arm section
    merge_at(i = 127:130, j = 5) %>%
      merge_at(i = 127:130, j = 6) %>%
        merge_at(i = 127:130, j = 7) 
  
ft3 %>% save_as_docx(path = "RecommendedUse_Master.docx")


















