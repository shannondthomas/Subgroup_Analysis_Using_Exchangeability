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
#         simulation_results_summary_08cutoff.csv  - summary of power for each 
#                                                    test for a 0.8 MEM cutoff 
#                                                    and 0.05 pval
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
allresults_unequaln_summary <- read.csv("simulation_results_unequaln_summary.csv")
allresults_unequaln_summary_08 <- read.csv("simulation_results_unequaln_summary_08cutoff.csv")
allresults_summary <- rbind(allresults_summary, allresults_summary_08,
                            allresults_unequaln_summary, allresults_unequaln_summary_08) %>% subset(select = c(-X))


#rearrange data to define confirmatory and exploratory use case settings

#pull t1e data
allresults_summary_t1e <- allresults_summary[((allresults_summary$mean2 - allresults_summary$mean1 == 0) & (is.na(allresults_summary$trteff1))) |
                                               ((allresults_summary$trteff1 - allresults_summary$trteff2 == 0) & !is.na(allresults_summary$trteff1)) ,]

colnames(allresults_summary_t1e)[(ncol(allresults_summary_t1e) - 4):(ncol(allresults_summary_t1e)-2)] <- c("MEM_t1e", "MEMr_t1e","pval_t1e")

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
allrecommendeduse$effsize_range[allrecommendeduse$recommendeduse == "Confirmatory"] <- paste("\u2265", allrecommendeduse$min_effsize)

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

ft3 <- ft2 %>% theme_box() %>%
  merge_at(i = 1:10, j = 4) %>% 
    merge_at(i = 1:4, j = 5) %>% 
      merge_at(i = 1:2, j = 6) %>% 
        merge_at(i = 1:2, j = 7) %>% 
      merge_at(i = 3:4, j = 6) %>% 
        merge_at(i = 3:4, j = 7) %>% 
    merge_at(i = 5:10, j = 5) %>%
      merge_at(i = 5:6, j = 6) %>%
      merge_at(i = 7:10, j = 6) %>%
        merge_at(i = 7:10, j = 7) %>%
  merge_at(i = 11:16, j = 4) %>% 
    merge_at(i = 12:16, j = 5) %>% 
      merge_at(i = 12:16, j = 6) %>%
        merge_at(i = 12:13, j = 7) %>%
        merge_at(i = 14:16, j = 7) %>%
  merge_at(i = 17:20, j = 4) %>% 
    merge_at(i = 18:20, j = 5) %>%
      merge_at(i = 18:20, j = 6) %>%
  merge_at(i = 21:24, j = 4) %>% 
    merge_at(i = 22:24, j = 5) %>%
      merge_at(i = 22:24, j = 6) %>%
  merge_at(i = 26:30, j = 4) %>% 
    merge_at(i = 26:27, j = 5) %>%
    merge_at(i = 28:30, j = 5) %>%
      merge_at(i = 29:30, j = 6) %>%
        merge_at(i = 29:30, j = 7) %>%
  merge_at(i = 31:35, j = 4) %>% 
    merge_at(i = 31:32, j = 5) %>%
      merge_at(i = 31:32, j = 6) %>%
        merge_at(i = 31:32, j = 7) %>%
    merge_at(i = 33:35, j = 5) %>%
      merge_at(i = 33:35, j = 6) %>%
        merge_at(i = 34:35, j = 7)


ft3 %>% save_as_docx( path = "RecommendedUse_Master.docx")



















