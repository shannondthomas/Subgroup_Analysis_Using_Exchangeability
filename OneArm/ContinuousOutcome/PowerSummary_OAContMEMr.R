library(ggplot2)
library(tidyverse)


#############################################
##############EQUAL SAMPLE###################
#############################################
#results <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/lmresults.csv")
results_raw <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/lmresults_raw.csv")

results_raw <- results_raw[,2:10]
colnames(results_raw)[6] <- 'pnexchgt0798'

#create character variable for nice labels on graphs
results_raw$var_c <- NA
results_raw$var_c[results_raw$var == 1] <- 'Var = 1'
results_raw$var_c[results_raw$var == 10] <- 'Var = 10'
results_raw$var_c[results_raw$var == 100] <- 'Var = 100'


#use raw data to calculate different cutoff performances
results_raw$pnexchgt020 <- results_raw$pnexch > 0.2
results_raw$pnexchgt080 <- results_raw$pnexch > 0.8
results_raw$pnexchgt090 <- results_raw$pnexch > 0.9

results_summary <- results_raw %>% group_by(totsampsize,var_c,efdiff) %>% 
  summarize(plt05 = mean(plt05), 
            pnexchgt090 = mean(pnexchgt090), 
            pnexchgt080 = mean(pnexchgt080), 
            pnexchgt0798 = mean(pnexchgt0798), 
            pnexchgt020 = mean(pnexchgt020),
  )

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEMr/PowerFtest_OAContMEMr.jpeg",
     width=1000, height=750)
ggplot(results_summary, aes(efdiff, totsampsize, fill=plt05)) + 
  geom_tile() + geom_text(aes(label=plt05),size = 5) + facet_grid(var_c~.) +
  labs(title = "Heatmap of F-test Power for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()


jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEMr/PowerPNEXCH_OAContMEMr_02.jpeg",
     width=1000, height=750)
ggplot(results_summary, aes(efdiff, totsampsize, fill=pnexchgt020)) + 
  geom_tile() + geom_text(aes(label=pnexchgt020),size = 5) + facet_grid(var_c~.) +
  labs(title = "Heatmap of  MEMr P(Not Exchangeable) > 0.2 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEMr/PowerPNEXCH_OAContMEMr_0798.jpeg",
     width=1000, height=750)
ggplot(results_summary, aes(efdiff, totsampsize, fill=pnexchgt0798)) + 
  geom_tile() + geom_text(aes(label=pnexchgt0798),size = 5) + facet_grid(var_c~.) +
  labs(title = "Heatmap of  MEMr P(Not Exchangeable) > 0.7976455 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEMr/PowerPNEXCH_OAContMEMr_08.jpeg",
     width=1000, height=750)
ggplot(results_summary, aes(efdiff, totsampsize, fill=pnexchgt080)) + 
  geom_tile() + geom_text(aes(label=pnexchgt080),size = 5) + facet_grid(var_c~.) +
  labs(title = "Heatmap of  MEMr P(Not Exchangeable) > 0.8 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEMr/PowerPNEXCH_OAContMEMr_09.jpeg",
     width=1000, height=750)
ggplot(results_summary, aes(efdiff, totsampsize, fill=pnexchgt090)) + 
  geom_tile() + geom_text(aes(label=pnexchgt090),size = 5) + facet_grid(var_c~.) +
  labs(title = "Heatmap of  MEMr P(Not Exchangeable) > 0.9 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()



######################
#####EXPAND VAR 1#####
######################
results_raw_expandvar1 <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/lmresults_expandvar1_raw.csv")

results_raw_expandvar1 <- results_raw_expandvar1[,2:10]
colnames(results_raw_expandvar1)[6] <- 'pnexchgt0798'

#create character variable for nice labels on graphs
results_raw_expandvar1$var_c <- NA
results_raw_expandvar1$var_c[results_raw_expandvar1$var == 1] <- 'Var = 1'


#use raw data to calculate different cutoff performances
results_raw_expandvar1$pnexchgt020 <- results_raw_expandvar1$pnexch > 0.2
results_raw_expandvar1$pnexchgt080 <- results_raw_expandvar1$pnexch > 0.8
results_raw_expandvar1$pnexchgt090 <- results_raw_expandvar1$pnexch > 0.9

results_summary_expandvar1 <- results_raw_expandvar1 %>% group_by(totsampsize,var_c,efdiff) %>% 
  summarize(plt05 = mean(plt05), 
            pnexchgt090 = mean(pnexchgt090), 
            pnexchgt080 = mean(pnexchgt080), 
            pnexchgt0798 = mean(pnexchgt0798), 
            pnexchgt020 = mean(pnexchgt020),
  )

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEMr/PowerFtest_OAContMEMr_expandvar1.jpeg",
     width=1000, height=750)
ggplot(results_summary_expandvar1, aes(efdiff, totsampsize, fill=plt05)) + 
  geom_tile() + geom_text(aes(label=plt05),size = 7) + facet_grid(var_c~.) +
  labs(title = "Heatmap of F-test Power for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()


jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEMr/PowerPNEXCH_OAContMEMr_02_expandvar1.jpeg",
     width=1000, height=750)
ggplot(results_summary_expandvar1, aes(efdiff, totsampsize, fill=pnexchgt020)) + 
  geom_tile() + geom_text(aes(label=pnexchgt020),size = 7) + facet_grid(var_c~.) +
  labs(title = "Heatmap of  MEMr P(Not Exchangeable) > 0.2 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEMr/PowerPNEXCH_OAContMEMr_0798_expandvar1.jpeg",
     width=1000, height=750)
ggplot(results_summary_expandvar1, aes(efdiff, totsampsize, fill=pnexchgt0798)) + 
  geom_tile() + geom_text(aes(label=pnexchgt0798),size = 7) + facet_grid(var_c~.) +
  labs(title = "Heatmap of  MEMr P(Not Exchangeable) > 0.7976455 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEMr/PowerPNEXCH_OAContMEMr_08_expandvar1.jpeg",
     width=1000, height=750)
ggplot(results_summary_expandvar1, aes(efdiff, totsampsize, fill=pnexchgt080)) + 
  geom_tile() + geom_text(aes(label=pnexchgt080),size = 7) + facet_grid(var_c~.) +
  labs(title = "Heatmap of  MEMr P(Not Exchangeable) > 0.8 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEMr/PowerPNEXCH_OAContMEMr_09_expandvar1.jpeg",
     width=1000, height=750)
ggplot(results_summary_expandvar1, aes(efdiff, totsampsize, fill=pnexchgt090)) + 
  geom_tile() + geom_text(aes(label=pnexchgt090),size = 7) + facet_grid(var_c~.) +
  labs(title = "Heatmap of  MEMr P(Not Exchangeable) > 0.9 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()




#############################################
#############UNEQUAL SAMPLE##################
#############################################

results_raw_unequalsamp <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/OA_cont_MEMr_result_raw.csv")
results_raw_unequalsamp <- results_raw_unequalsamp[,2:10]
colnames(results_raw_unequalsamp)[9] <- 'pnexchgt0798'

results_raw_unequalsamp$n2mult_c <- NA
results_raw_unequalsamp$n2mult_c[results_raw_unequalsamp$n2mult == 0.1] <- 'n2 = 0.1*N'
results_raw_unequalsamp$n2mult_c[results_raw_unequalsamp$n2mult == 0.25] <- 'n2 = 0.25*N'

#use raw data to calculate different cutoff performances
results_raw_unequalsamp$pnexchgt020 <- results_raw_unequalsamp$pnexch > 0.2
results_raw_unequalsamp$pnexchgt080 <- results_raw_unequalsamp$pnexch > 0.8
results_raw_unequalsamp$pnexchgt090 <- results_raw_unequalsamp$pnexch > 0.9

results_summary_unequalsamp <- results_raw_unequalsamp %>% group_by(totsampsize,n2mult_c,efdiff) %>% 
  summarize(plt05 = mean(plt05), 
            pnexchgt090 = mean(pnexchgt090), 
            pnexchgt080 = mean(pnexchgt080), 
            pnexchgt0798 = mean(pnexchgt0798), 
            pnexchgt020 = mean(pnexchgt020),
  )

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-UnequalSamp/MEMr/PowerFtest_OAContMEMr_UnequalSamp.jpeg",
     width=1000, height=750)
ggplot(results_summary_unequalsamp, aes(efdiff, totsampsize, fill=plt05)) + 
  geom_tile() + geom_text(aes(label=plt05),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of F-test Power for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-UnequalSamp/MEMr/PowerPNEXCH_OAContMEMr_UnequalSamp_09.jpeg",
     width=1000, height=750)
ggplot(results_summary_unequalsamp, aes(efdiff, totsampsize, fill=pnexchgt090)) + 
  geom_tile() + geom_text(aes(label=pnexchgt090),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.9 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-UnequalSamp/MEMr/PowerPNEXCH_OAContMEMr_UnequalSamp_08.jpeg",
     width=1000, height=750)
ggplot(results_summary_unequalsamp, aes(efdiff, totsampsize, fill=pnexchgt080)) + 
  geom_tile() + geom_text(aes(label=pnexchgt080),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.8 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-UnequalSamp/MEMr/PowerPNEXCH_OAContMEMr_UnequalSamp_0798.jpeg",
     width=1000, height=750)
ggplot(results_summary_unequalsamp, aes(efdiff, totsampsize, fill=pnexchgt0798)) + 
  geom_tile() + geom_text(aes(label=pnexchgt0798),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.7976455 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-UnequalSamp/MEMr/PowerPNEXCH_OABinMEMr_UnequalSamp_02.jpeg",
     width=1000, height=750)
ggplot(results_summary_unequalsamp, aes(efdiff, totsampsize, fill=pnexchgt020)) + 
  geom_tile() + geom_text(aes(label=pnexchgt020),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.2 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()


