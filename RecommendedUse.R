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

#set default p-value power to be the chisq power in the OA binary setting
allresults_summary$pval_power[!is.na(allresults_summary$pval_chisq_power)] <- allresults_summary$pval_chisq_power[!is.na(allresults_summary$pval_chisq_power)]

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
allresults_summary_wide$effect_size <-ifelse(allresults_summary_wide$num_arms == 1, 
                                             allresults_summary_wide$mean2 - allresults_summary_wide$mean1,
                                             allresults_summary_wide$trteff2 - allresults_summary_wide$trteff1)

allresults_summary_wide$n2_mult <- allresults_summary_wide$n2/allresults_summary_wide$N


recommendeduse_summary <- function(name, memtype, recommendeduse){
  sum1 <- allresults_summary_wide[allresults_summary_wide[[name]] == 1,] %>% 
    group_by(n2_mult, MEMcutoff, num_arms, outcome_type, Variance, N, .data[[name]]) %>%
    summarize(min_effsize = min(effect_size),
              max_effsize = max(effect_size))
  
  
  testsum2 <- sum1 %>% 
    group_by(n2_mult, MEMcutoff, num_arms, outcome_type, Variance, min_effsize, max_effsize) %>% 
    arrange(n2_mult, MEMcutoff, num_arms, outcome_type, Variance, N, min_effsize, max_effsize) %>% 
    mutate(sampgap = ifelse((lag(N) == (N-100)), 0, 1))
  
  testsum2$sampgap[is.na(testsum2$sampgap)] <- 0
  
  testsum3 <- testsum2 %>%
    group_by(n2_mult, MEMcutoff, num_arms, outcome_type, Variance, min_effsize, max_effsize) %>% 
    arrange(n2_mult, MEMcutoff, num_arms, outcome_type, Variance, N, min_effsize, max_effsize) %>% 
    mutate(rangeind = 1+cumsum(sampgap == 1))

  sum2 <- testsum3 %>%
   group_by(n2_mult, MEMcutoff, num_arms, outcome_type, Variance, min_effsize, max_effsize, rangeind) %>%
   summarize(min_N = min(N),
             max_N = max(N))
  
  sum2 <- sum2[complete.cases(sum2 %>% subset(select = c(-Variance))),]
  
  sum2$MEMtype <- memtype
  sum2$recommendeduse <- recommendeduse
  
  return(sum2)
}

################GETTING WEIRD RESULTS -- CHECK MONDAY################
#have unnecessary row in 0.1 confirmatory calibrated TA binary MEMr results 
################CHECK FOR OTHER THINGS LIKE THIS#############################################
MEM_confirmatory <- recommendeduse_summary("MEMconfirmatory", "MEM", "Confirmatory")
MEM_exploratory <- recommendeduse_summary("MEMexploratory", "MEM", "Exploratory")
MEMr_confirmatory <- recommendeduse_summary("MEMrconfirmatory", "MEMr", "Confirmatory")
MEMr_exploratory <- recommendeduse_summary("MEMrexploratory", "MEMr", "Exploratory")


allrecommendeduse <- rbind(MEM_confirmatory, MEMr_confirmatory, MEM_exploratory, MEMr_exploratory) %>%
                        arrange(n2_mult, recommendeduse, MEMcutoff, num_arms, outcome_type, MEMtype, Variance, min_N) 

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

##subset tables
equal_conf_02 <- flextable(allrecommendeduse[1:18,4:9]) %>% theme_box() %>%
  merge_at(i = 1:10, j = 4-3) %>% #0.5 n2mult, confirmatory, 0.2 cutoff, 1 arm section
    merge_at(i = 1:4, j = 5-3) %>% 
      merge_at(i = 1:2, j = 6-3) %>% 
        merge_at(i = 1:2, j = 7-3) %>% 
      merge_at(i = 3:4, j = 6-3) %>% 
        merge_at(i = 3:4, j = 7-3) %>% 
    merge_at(i = 5:10, j = 5-3) %>%
      merge_at(i = 5:6, j = 6-3) %>%
      merge_at(i = 7:10, j = 6-3) %>%
       merge_at(i = 8:9, j = 7-3) %>%
  merge_at(i = 11:18, j = 4-3) %>% #0.5 n2mult, confirmatory, 0.2 cutoff, 2 arm section
    merge_at(i = 11:13, j = 5-3) %>% 
      merge_at(i = 11:12, j = 6-3) %>%
        merge_at(i = 11:12, j = 7-3) %>%
    merge_at(i = 14:18, j = 5-3) %>% 
      merge_at(i = 14:18, j = 6-3) %>%
        merge_at(i = 14:15, j = 7-3) %>%
        merge_at(i = 16:18, j = 7-3) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption(caption = "Table S1: Recommended use cases for the confirmatory setting (type 1 error rate less than standard test and power within 0.05 of the standard test) using the 0.2 cutoff for equal subgroup sample sizes.")

equal_conf_08 <- flextable(allrecommendeduse[19:25,4:9]) %>% theme_box() %>%
  merge_at(i = (19-18):(21-18), j = 4-3) %>% #0.5 n2mult, confirmatory, 0.8 cutoff, 1 arm section
    merge_at(i = (19-18):(21-18), j = 5-3) %>%
      merge_at(i = (19-18):(21-18), j = 6-3) %>%
  merge_at(i = (22-18):(25-18), j = 4-3) %>% #0.5 n2mult, confirmatory, 0.8 cutoff, 2 arm section
    merge_at(i = (23-18):(25-18), j = 5-3) %>%
      merge_at(i = (23-18):(25-18), j = 6-3)%>%
  set_table_properties(layout = "autofit") %>%
  set_caption(caption = "Table S2: Recommended use cases for the confirmatory setting (type 1 error rate less than standard test and power within 0.05 of the standard test) using the 0.8 cutoff for equal subgroup sample sizes.")
  
equal_conf_cal <- flextable(allrecommendeduse[26:53,4:9]) %>% theme_box() %>%
  merge_at(i = (26-25):(37-25), j = 4-3) %>% #0.5 n2mult, confirmatory, calibrated cutoff, 1 arm section
    merge_at(i = (26-25):(37-25), j = 5-3) %>%
      merge_at(i = (26-25):(31-25), j = 6-3) %>%
        merge_at(i = (26-25):(27-25), j = 7-3) %>%
        merge_at(i = (28-25):(29-25), j = 7-3) %>%
        merge_at(i = (30-25):(31-25), j = 7-3) %>%
      merge_at(i = (32-25):(37-25), j = 6-3) %>%
        merge_at(i = (32-25):(33-25), j = 7-3) %>%
        merge_at(i = (34-25):(35-25), j = 7-3) %>%
        merge_at(i = (36-25):(37-25), j = 7-3) %>%
  merge_at(i = (38-25):(53-25), j = 4-3) %>% #0.5 n2mult, confirmatory, calibrated cutoff, 2 arm section 
    merge_at(i = (38-25):(41-25), j = 5-3) %>%
      merge_at(i = (38-25):(39-25), j = 6-3) %>%
        merge_at(i = (38-25):(39-25), j = 7-3) %>%
      merge_at(i = (40-25):(41-25), j = 6-3) %>%
        merge_at(i = (40-25):(41-25), j = 7-3) %>%
    merge_at(i = (42-25):(53-25), j = 5-3) %>%
      merge_at(i = (42-25):(47-25), j = 6-3) %>%
        merge_at(i = (42-25):(43-25), j = 7-3) %>%
        merge_at(i = (44-25):(45-25), j = 7-3) %>%
        merge_at(i = (46-25):(47-25), j = 7-3) %>%
      merge_at(i = (48-25):(53-25), j = 6-3) %>%
        merge_at(i = (48-25):(49-25), j = 7-3) %>%
        merge_at(i = (50-25):(51-25), j = 7-3) %>%
        merge_at(i = (52-25):(53-25), j = 7-3) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption(caption = "Table S3: Recommended use cases for the confirmatory setting (type 1 error rate less than standard test and power within 0.05 of the standard test) using the calibrated cutoff for equal subgroup sample sizes.")


equal_exp_02 <- flextable(allrecommendeduse[54:57,4:9]) %>% theme_box() %>%
  merge_at(i = (55-53):(57-53), j = 4-3) %>% #0.5 n2mult, exploratory, 0.2 cutoff, 2 arm section
    merge_at(i = (55-53):(57-53), j = 5-3) %>%
      merge_at(i = (55-53):(57-53), j = 6-3) %>%
        merge_at(i = (55-53):(57-53), j = 7-3) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption(caption = "Table S4: Recommended use cases for the exploratory setting (type 1 error rate < 0.15 and power at least 0.1 higher than standard test) using the 0.2 cutoff for equal subgroup sample sizes.")

equal_exp_08 <- flextable(allrecommendeduse[58:66,4:9]) %>% theme_box() %>%
  merge_at(i = (58-57):(62-57), j = 4-3) %>% #0.5 n2mult, exploratory, 0.8 cutoff, 1 arm section
    merge_at(i = (58-57):(59-57), j = 5-3) %>%
    merge_at(i = (60-57):(62-57), j = 5-3) %>%
      merge_at(i = (60-57):(62-57), j = 6-3) %>%
        merge_at(i = (61-57):(62-57), j = 7-3) %>%
  merge_at(i = (63-57):(66-57), j = 4-3) %>% #0.5 n2mult, exploratory, 0.8 cutoff, 2 arm section
    merge_at(i = (64-57):(66-57), j = 5-3) %>%
      merge_at(i = (64-57):(66-57), j = 6-3) %>%
        merge_at(i = (65-57):(66-57), j = 7-3) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption(caption = "Table S5: Recommended use cases for the exploratory setting (type 1 error rate < 0.15 and power at least 0.1 higher than standard test) using the 0.8 cutoff for equal subgroup sample sizes.")

equal_exp_cal <- flextable(allrecommendeduse[67:71,4:9]) %>% theme_box() %>%
  merge_at(i = (68-66):(71-66), j = 4-3) %>% #0.5 n2mult, exploratory, calibrated cutoff, 2 arm section
    merge_at(i = (68-66):(71-66), j = 5-3) %>% 
      merge_at(i = (68-66):(71-66), j = 6-3) %>% 
        merge_at(i = (68-66):(71-66), j = 7-3) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption(caption = "Table S6: Recommended use cases for the exploratory setting (type 1 error rate < 0.15 and power at least 0.1 higher than standard test) using the calibrated cutoff for equal subgroup sample sizes.")

    
    
unequal_conf_02 <- flextable(rbind(allrecommendeduse[72:79,c(1,4:9)],
                                   allrecommendeduse[109:118,c(1,4:9)])) %>% theme_box() %>%
  merge_v(j = c("Sample \nSize Dist")) %>%
  merge_at(i = (72-71):(77-71), j = 4-2) %>% #0.25 n2mult, confirmatory, 0.2 cutoff, 1 arm section
    merge_at(i = (72-71):(75-71), j = 5-2) %>%
      merge_at(i = (72-71):(73-71), j = 6-2) %>%
        merge_at(i = (72-71):(73-71), j = 7-2) %>%
      merge_at(i = (74-71):(75-71), j = 6-2) %>%
        merge_at(i = (74-71):(75-71), j = 7-2) %>%
    merge_at(i = (76-71):(77-71), j = 5-2) %>%
      merge_at(i = (76-71):(77-71), j = 6-2) %>%
        merge_at(i = (76-71):(77-71), j = 7-2) %>%
  merge_at(i = (78-71):(79-71), j = 4-2) %>% #0.25 n2mult, confirmatory, 0.2 cutoff, 2 arm section
    merge_at(i = (78-71):(79-71), j = 5-2) %>%
      merge_at(i = (78-71):(79-71), j = 6-2) %>%
        merge_at(i = (78-71):(79-71), j = 7-2) %>%
                   
  merge_at(i = (109-108+8):(117-108+8), j = 4-2) %>% #0.1 n2mult, confirmatory, 0.2 cutoff, 1 arm section
    merge_at(i = (109-108+8):(114-108+8), j = 5-2) %>%
      merge_at(i = (109-108+8):(111-108+8), j = 6-2) %>%
        merge_at(i = (109-108+8):(111-108+8), j = 7-2) %>%
      merge_at(i = (112-108+8):(114-108+8), j = 6-2) %>%
        merge_at(i = (112-108+8):(114-108+8), j = 7-2) %>%
    merge_at(i = (115-108+8):(117-108+8), j = 5-2) %>%
      merge_at(i = (115-108+8):(117-108+8), j = 6-2) %>%
        merge_at(i = (115-108+8):(117-108+8), j = 7-2) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption(caption = "Table S7: Recommended use cases for the confirmatory setting (type 1 error rate less than standard test and power within 0.05 of the standard test) using the 0.2 cutoff for unequal subgroup sample sizes.")
                 
unequal_conf_08 <- flextable(rbind(allrecommendeduse[80:82,c(1,4:9)],
                                   allrecommendeduse[119:121,c(1,4:9)])) %>% theme_box() %>%
  merge_v(j = c("Sample \nSize Dist")) %>%
  merge_at(i = (81-79):(82-79), j = 4-2) %>% #0.25 n2mult, confirmatory, 0.8 cutoff, 2 arm section
  merge_at(i = (120-118+3):(121-118+3), j = 4-2) %>% #0.1 n2mult, confirmatory, 0.8 cutoff, 2 arm section
  set_table_properties(layout = "autofit") %>%
  set_caption(caption = "Table S8: Recommended use cases for the confirmatory setting (type 1 error rate less than standard test and power within 0.05 of the standard test) using the 0.8 cutoff for unequal subgroup sample sizes.")

unequal_conf_cal <- flextable(rbind(allrecommendeduse[83:93,c(1,4:9)],
                                    allrecommendeduse[122:137,c(1,4:9)])) %>% theme_box() %>%
  merge_v(j = c("Sample \nSize Dist")) %>%
  merge_at(i = (83-82):(86-82), j = 4-2) %>% #0.25 n2mult, confirmatory, calibrated cutoff, 1 arm section
    merge_at(i = (83-82):(86-82), j = 5-2) %>%
      merge_at(i = (83-82):(84-82), j = 6-2) %>%
        merge_at(i = (83-82):(84-82), j = 7-2) %>%
      merge_at(i = (85-82):(86-82), j = 6-2) %>%
        merge_at(i = (85-82):(86-82), j = 7-2) %>%
  merge_at(i = (87-82):(93-82), j = 4-2) %>% #0.25 n2mult, confirmatory, calibrated cutoff, 2 arm section
    merge_at(i = (87-82):(89-82), j = 5-2) %>%
      merge_at(i = (88-82):(89-82), j = 6-2) %>%
        merge_at(i = (88-82):(89-82), j = 7-2) %>%
    merge_at(i = (90-82):(93-82), j = 5-2) %>%
      merge_at(i = (90-82):(91-82), j = 6-2) %>%
        merge_at(i = (90-82):(91-82), j = 7-2) %>%
      merge_at(i = (92-82):(93-82), j = 6-2) %>%
        merge_at(i = (92-82):(93-82), j = 7-2) %>%
    
  merge_at(i = (122-121+11):(125-121+11), j = 4-2) %>% #0.1 n2mult, confirmatory, calibrated cutoff, 1 arm section
    merge_at(i = (122-121+11):(125-121+11), j = 5-2) %>%
      merge_at(i = (122-121+11):(123-121+11), j = 6-2) %>%
        merge_at(i = (122-121+11):(123-121+11), j = 7-2) %>%
      merge_at(i = (124-121+11):(125-121+11), j = 6-2) %>%
        merge_at(i = (124-121+11):(125-121+11), j = 7-2) %>%
  merge_at(i = (126-121+11):(137-121+11), j = 4-2) %>% #0.1 n2mult, confirmatory, calibrated cutoff, 2 arm section
    merge_at(i = (126-121+11):(129-121+11), j = 5-2) %>%
      merge_at(i = (127-121+11):(129-121+11), j = 6-2) %>%
        merge_at(i = (127-121+11):(129-121+11), j = 7-2) %>%
    merge_at(i = (130-121+11):(137-121+11), j = 5-2) %>%
      merge_at(i = (130-121+11):(133-121+11), j = 6-2) %>%
        merge_at(i = (130-121+11):(133-121+11), j = 7-2) %>%
      merge_at(i = (134-121+11):(137-121+11), j = 6-2) %>%
        merge_at(i = (134-121+11):(137-121+11), j = 7-2) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption(caption = "Table S9: Recommended use cases for the confirmatory setting (type 1 error rate less than standard test and power within 0.05 of the standard test) using the calibrated cutoff for unequal subgroup sample sizes.")

                   
unequal_exp_02 <- flextable(rbind(allrecommendeduse[94:98,c(1,4:9)],
                                    allrecommendeduse[138:142,c(1,4:9)])) %>% theme_box() %>%
  merge_v(j = c("Sample \nSize Dist")) %>%
  merge_at(i = (95-93):(98-93), j = 4-2) %>% #0.25 n2mult, exploratory, 0.2 cutoff, 2 arm section
    merge_at(i = (95-93):(98-93), j = 5-2) %>%
      merge_at(i = (95-93):(98-93), j = 6-2) %>%
        merge_at(i = (95-93):(98-93), j = 7-2) %>%
  merge_at(i = (139-137+5):(142-137+5), j = 4-2) %>% #0.1 n2mult, exploratory, 0.2 cutoff, 2 arm section
    merge_at(i = (139-137+5):(142-137+5), j = 5-2) %>%
      merge_at(i = (139-137+5):(142-137+5), j = 6-2) %>%
        merge_at(i = (139-137+5):(142-137+5), j = 7-2) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption(caption = "Table S10: Recommended use cases for the exploratory setting (type 1 error rate < 0.15 and power at least 0.1 higher than standard test) using the 0.2 cutoff for unequal subgroup sample sizes.")

unequal_exp_08 <- flextable(rbind(allrecommendeduse[99:103,c(1,4:9)],
                                  allrecommendeduse[143:147,c(1,4:9)])) %>% theme_box() %>%
  merge_v(j = c("Sample \nSize Dist")) %>%
  merge_at(i = (99-98):(101-98), j = 4-2) %>% #0.25 n2mult, exploratory, 0.8 cutoff, 1 arm section
    merge_at(i = (99-98):(100-98), j = 5-2) %>%
  merge_at(i = (102-98):(103-98), j = 4-2) %>% #0.25 n2mult, exploratory, 0.8 cutoff, 2 arm section
             
  merge_at(i = (143-142+5):(144-142+5), j = 4-2) %>% #0.1 n2mult, exploratory, 0.8 cutoff, 1 arm section
  merge_at(i = (145-142+5):(147-142+5), j = 4-2) %>% #0.1 n2mult, exploratory, 0.8 cutoff, 2 arm section
    merge_at(i = (146-142+5):(147-142+5), j = 5-2) %>%
      merge_at(i = (146-142+5):(147-142+5), j = 6-2) %>%
        merge_at(i = (146-142+5):(147-142+5), j = 7-2) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption(caption = "Table S11: Recommended use cases for the exploratory setting (type 1 error rate < 0.15 and power at least 0.1 higher than standard test) using the 0.8 cutoff for unequal subgroup sample sizes.")

unequal_exp_cal <- flextable(rbind(allrecommendeduse[104:108,c(1,4:9)],
                                   allrecommendeduse[148:151,c(1,4:9)])) %>% theme_box() %>%
  merge_v(j = c("Sample \nSize Dist")) %>%
  merge_at(i = (105-103):(108-103), j = 4-2) %>% #0.25 n2mult, exploratory, calibrated cutoff, 2 arm section
    merge_at(i = (105-103):(108-103), j = 5-2) %>%
      merge_at(i = (105-103):(108-103), j = 6-2) %>%
        merge_at(i = (105-103):(108-103), j = 7-2) %>%
  merge_at(i = (148-147+5):(151-147+5), j = 4-2) %>% #0.1 n2mult, exploratory, calibrated cutoff, 2 arm section
    merge_at(i = (148-147+5):(151-147+5), j = 5-2) %>%
      merge_at(i = (148-147+5):(151-147+5), j = 6-2) %>%
        merge_at(i = (148-147+5):(151-147+5), j = 7-2) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption(caption = "Table S12: Recommended use cases for the exploratory setting (type 1 error rate < 0.15 and power at least 0.1 higher than standard test) using the calibrated cutoff for unequal subgroup sample sizes.")

save_as_docx(equal_conf_02,
             equal_conf_08,
             equal_conf_cal,
             equal_exp_02,
             equal_exp_08,
             equal_exp_cal,
             unequal_conf_02,
             unequal_conf_08,
             unequal_conf_cal,
             unequal_exp_02,
             unequal_exp_08,
             unequal_exp_cal,
             path = "RecommendedUseTables.docx")


# ##FULL TABLE
# ft3 <- ft2 %>% theme_box() %>%
#   merge_at(i = 1:10, j = 4) %>% #0.5 n2mult, confirmatory, 0.2 cutoff, 1 arm section
#     merge_at(i = 1:4, j = 5) %>%
#       merge_at(i = 1:2, j = 6) %>%
#         merge_at(i = 1:2, j = 7) %>%
#       merge_at(i = 3:4, j = 6) %>%
#         merge_at(i = 3:4, j = 7) %>%
#     merge_at(i = 5:10, j = 5) %>%
#       merge_at(i = 5:6, j = 6) %>%
#       merge_at(i = 7:10, j = 6) %>%
#        merge_at(i = 8:9, j = 7) %>%
#   merge_at(i = 11:18, j = 4) %>% #0.5 n2mult, confirmatory, 0.2 cutoff, 2 arm section
#     merge_at(i = 11:13, j = 5) %>%
#       merge_at(i = 11:12, j = 6) %>%
#         merge_at(i = 11:12, j = 7) %>%
#     merge_at(i = 14:18, j = 5) %>%
#       merge_at(i = 14:18, j = 6) %>%
#         merge_at(i = 14:15, j = 7) %>%
#         merge_at(i = 16:18, j = 7) %>%
#   merge_at(i = 19:21, j = 4) %>% #0.5 n2mult, confirmatory, 0.8 cutoff, 1 arm section
#     merge_at(i = 19:21, j = 5) %>%
#       merge_at(i = 19:21, j = 6) %>%
#   merge_at(i = 22:25, j = 4) %>% #0.5 n2mult, confirmatory, 0.8 cutoff, 2 arm section
#     merge_at(i = 23:25, j = 5) %>%
#       merge_at(i = 23:25, j = 6) %>%
#   merge_at(i = 26:37, j = 4) %>% #0.5 n2mult, confirmatory, calibrated cutoff, 1 arm section
#     merge_at(i = 26:37, j = 5) %>%
#       merge_at(i = 26:31, j = 6) %>%
#         merge_at(i = 26:27, j = 7) %>%
#         merge_at(i = 28:29, j = 7) %>%
#         merge_at(i = 30:31, j = 7) %>%
#       merge_at(i = 32:37, j = 6) %>%
#         merge_at(i = 32:33, j = 7) %>%
#         merge_at(i = 34:35, j = 7) %>%
#         merge_at(i = 36:37, j = 7) %>%
#   merge_at(i = 38:53, j = 4) %>% #0.5 n2mult, confirmatory, calibrated cutoff, 2 arm section
#     merge_at(i = 38:41, j = 5) %>%
#       merge_at(i = 38:39, j = 6) %>%
#         merge_at(i = 38:39, j = 7) %>%
#       merge_at(i = 40:41, j = 6) %>%
#         merge_at(i = 40:41, j = 7) %>%
#     merge_at(i = 42:53, j = 5) %>%
#       merge_at(i = 42:47, j = 6) %>%
#         merge_at(i = 42:43, j = 7) %>%
#         merge_at(i = 44:45, j = 7) %>%
#         merge_at(i = 46:47, j = 7) %>%
#       merge_at(i = 48:53, j = 6) %>%
#         merge_at(i = 48:49, j = 7) %>%
#         merge_at(i = 50:51, j = 7) %>%
#         merge_at(i = 52:53, j = 7) %>%
#   merge_at(i = 55:57, j = 4) %>% #0.5 n2mult, exploratory, 0.2 cutoff, 2 arm section
#     merge_at(i = 55:57, j = 5) %>%
#       merge_at(i = 55:57, j = 6) %>%
#         merge_at(i = 55:57, j = 7) %>%
#   merge_at(i = 58:62, j = 4) %>% #0.5 n2mult, exploratory, 0.8 cutoff, 1 arm section
#     merge_at(i = 58:59, j = 5) %>%
#     merge_at(i = 60:62, j = 5) %>%
#       merge_at(i = 60:62, j = 6) %>%
#         merge_at(i = 61:62, j = 7) %>%
#   merge_at(i = 63:66, j = 4) %>% #0.5 n2mult, exploratory, 0.8 cutoff, 2 arm section
#     merge_at(i = 64:66, j = 5) %>%
#       merge_at(i = 64:66, j = 6) %>%
#         merge_at(i = 65:66, j = 7) %>%
#   merge_at(i = 68:71, j = 4) %>% #0.5 n2mult, exploratory, calibrated cutoff, 2 arm section
#     merge_at(i = 68:71, j = 5) %>%
#       merge_at(i = 68:71, j = 6) %>%
#         merge_at(i = 68:71, j = 7) %>%
# 
#   merge_at(i = 72:77, j = 4) %>% #0.25 n2mult, confirmatory, 0.2 cutoff, 1 arm section
#     merge_at(i = 72:75, j = 5) %>%
#       merge_at(i = 72:73, j = 6) %>%
#         merge_at(i = 72:73, j = 7) %>%
#       merge_at(i = 74:75, j = 6) %>%
#         merge_at(i = 74:75, j = 7) %>%
#     merge_at(i = 76:77, j = 5) %>%
#       merge_at(i = 76:77, j = 6) %>%
#         merge_at(i = 76:77, j = 7) %>%
#   merge_at(i = 78:79, j = 4) %>% #0.25 n2mult, confirmatory, 0.2 cutoff, 2 arm section
#     merge_at(i = 78:79, j = 5) %>%
#       merge_at(i = 78:79, j = 6) %>%
#         merge_at(i = 78:79, j = 7) %>%
#   merge_at(i = 81:82, j = 4) %>% #0.25 n2mult, confirmatory, 0.8 cutoff, 2 arm section
#   merge_at(i = 83:86, j = 4) %>% #0.25 n2mult, confirmatory, calibrated cutoff, 1 arm section
#     merge_at(i = 83:86, j = 5) %>%
#       merge_at(i = 83:84, j = 6) %>%
#         merge_at(i = 83:84, j = 7) %>%
#       merge_at(i = 85:86, j = 6) %>%
#         merge_at(i = 85:86, j = 7) %>%
#   merge_at(i = 87:93, j = 4) %>% #0.25 n2mult, confirmatory, calibrated cutoff, 2 arm section
#     merge_at(i = 87:89, j = 5) %>%
#       merge_at(i = 88:89, j = 6) %>%
#         merge_at(i = 88:89, j = 7) %>%
#     merge_at(i = 90:93, j = 5) %>%
#       merge_at(i = 90:91, j = 6) %>%
#         merge_at(i = 90:91, j = 7) %>%
#       merge_at(i = 92:93, j = 6) %>%
#         merge_at(i = 92:93, j = 7) %>%
#   merge_at(i = 95:98, j = 4) %>% #0.25 n2mult, exploratory, 0.2 cutoff, 2 arm section
#     merge_at(i = 95:98, j = 5) %>%
#       merge_at(i = 95:98, j = 6) %>%
#         merge_at(i = 95:98, j = 7) %>%
#   merge_at(i = 99:101, j = 4) %>% #0.25 n2mult, exploratory, 0.8 cutoff, 1 arm section
#     merge_at(i = 99:100, j = 5) %>%
#   merge_at(i = 102:103, j = 4) %>% #0.25 n2mult, exploratory, 0.8 cutoff, 2 arm section
#   merge_at(i = 105:108, j = 4) %>% #0.25 n2mult, exploratory, calibrated cutoff, 2 arm section
#     merge_at(i = 105:108, j = 5) %>%
#       merge_at(i = 105:108, j = 6) %>%
#         merge_at(i = 105:108, j = 7) %>%
# 
#   merge_at(i = 109:117, j = 4) %>% #0.1 n2mult, confirmatory, 0.2 cutoff, 1 arm section
#     merge_at(i = 109:114, j = 5) %>%
#       merge_at(i = 109:111, j = 6) %>%
#         merge_at(i = 109:111, j = 7) %>%
#       merge_at(i = 112:114, j = 6) %>%
#         merge_at(i = 112:114, j = 7) %>%
#     merge_at(i = 115:117, j = 5) %>%
#       merge_at(i = 115:117, j = 6) %>%
#         merge_at(i = 115:117, j = 7) %>%
#   merge_at(i = 120:121, j = 4) %>% #0.1 n2mult, confirmatory, 0.8 cutoff, 2 arm section
#   merge_at(i = 122:125, j = 4) %>% #0.1 n2mult, confirmatory, calibrated cutoff, 1 arm section
#     merge_at(i = 122:125, j = 5) %>%
#       merge_at(i = 122:123, j = 6) %>%
#         merge_at(i = 122:123, j = 7) %>%
#       merge_at(i = 124:125, j = 6) %>%
#         merge_at(i = 124:125, j = 7) %>%
#   merge_at(i = 126:137, j = 4) %>% #0.1 n2mult, confirmatory, calibrated cutoff, 2 arm section
#     merge_at(i = 126:129, j = 5) %>%
#       merge_at(i = 127:129, j = 6) %>%
#         merge_at(i = 127:129, j = 7) %>%
#     merge_at(i = 130:137, j = 5) %>%
#       merge_at(i = 130:133, j = 6) %>%
#         merge_at(i = 130:133, j = 7) %>%
#       merge_at(i = 134:137, j = 6) %>%
#         merge_at(i = 134:137, j = 7) %>%
#   merge_at(i = 139:142, j = 4) %>% #0.1 n2mult, exploratory, 0.2 cutoff, 2 arm section
#     merge_at(i = 139:142, j = 5) %>%
#       merge_at(i = 139:142, j = 6) %>%
#         merge_at(i = 139:142, j = 7) %>%
#   merge_at(i = 143:144, j = 4) %>% #0.1 n2mult, exploratory, 0.8 cutoff, 1 arm section
#   merge_at(i = 145:147, j = 4) %>% #0.1 n2mult, exploratory, 0.8 cutoff, 2 arm section
#     merge_at(i = 146:147, j = 5) %>%
#       merge_at(i = 146:147, j = 6) %>%
#         merge_at(i = 146:147, j = 7) %>%
#   merge_at(i = 148:151, j = 4) %>% #0.1 n2mult, exploratory, calibrated cutoff, 2 arm section
#     merge_at(i = 148:151, j = 5) %>%
#       merge_at(i = 148:151, j = 6) %>%
#         merge_at(i = 148:151, j = 7)
  
# ft3 %>% save_as_docx(path = "RecommendedUse_Master.docx")


















