library(tidyverse)
library(ggplot2)


#############################################
###############EQUAL SAMPLE##################
#############################################
OAContMEM <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/OA_cont_MEM_result.csv")

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEM/PowerTtest_OAContMEM.jpeg",
     width=1000, height=750)
ggplot(OAContMEM, aes(es2-es1, N, fill=power)) +
  geom_tile() + geom_text(aes(label=power), size = 7) + facet_grid(var_c~.) +
  labs(title = "Heatmap of One-Arm T-test Power for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = 'Type I Error/Power') +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEM/PowerPNEXCH_OAContMEM_Cutoff02.jpeg",
     width=1000, height=750)
ggplot(OAContMEM, aes(es2-es1, N, fill=pnexch_cutoff2_power)) + 
    geom_tile() + geom_text(aes(label=pnexch_cutoff2_power), size = 7) + facet_grid(var_c~.) + 
    labs(title = "Heatmap of One-Arm MEM P(Not Exchangeable) >= 0.2 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
         fill = "Type I Error/Power") +
    xlab("Difference in Effect Sizes") +
    ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEM/PowerPNEXCH_OAContMEM_Cutoff0798.jpeg",
     width=1000, height=750)
ggplot(OAContMEM, aes(es2-es1, N, fill=pnexch_cutoff797_power)) + 
  geom_tile() + geom_text(aes(label=pnexch_cutoff797_power), size = 7) + facet_grid(var_c~.) + 
  labs(title = "Heatmap of One-Arm MEM P(Not Exchangeable) >= 0.7976455 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEM/PowerPNEXCH_OAContMEM_Cutoff08.jpeg",
     width=1000, height=750)
ggplot(OAContMEM, aes(es2-es1, N, fill=pnexch_cutoff8_power)) + 
  geom_tile() + geom_text(aes(label=pnexch_cutoff8_power), size = 7) + facet_grid(var_c~.) + 
  labs(title = "Heatmap of One-Arm MEM P(Not Exchangeable) >= 0.8 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()


jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEM/PowerPNEXCH_OAContMEM_Cutoff09.jpeg",
     width=1000, height=750)
ggplot(OAContMEM, aes(es2-es1, N, fill=pnexch_cutoff9_power)) + 
  geom_tile() + geom_text(aes(label=pnexch_cutoff9_power), size = 7) + facet_grid(var_c~.) + 
  labs(title = "Heatmap of One-Arm MEM P(Not Exchangeable) >= 0.9 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()


######################
#####expand var=1#####
######################
OAContMEM_expandvar1 <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/OA_cont_MEM_result_ExpandVar1.csv")


jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEM/PowerTtest_OAContMEM_ExpandVar1.jpeg",
     width=1000, height=750)
ggplot(OAContMEM_expandvar1, aes(es2-es1, N, fill=power)) +
  geom_tile() + geom_text(aes(label=power), size = 7) + 
  labs(title = "Heatmap of One-Arm T-test Power for Varying Sample Sizes and Differences in Effect Sizes",
       fill = 'Type I Error/Power') +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEM/PowerPNEXCH_OAContMEM_Cutoff02_ExpandVar1.jpeg",
     width=1000, height=750)
ggplot(OAContMEM_expandvar1, aes(es2-es1, N, fill=pnexch_cutoff2_power)) + 
  geom_tile() + geom_text(aes(label=pnexch_cutoff2_power), size = 7) + 
  labs(title = "Heatmap of One-Arm MEM P(Not Exchangeable) >= 0.2 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEM/PowerPNEXCH_OAContMEM_Cutoff0798_ExpandVar1.jpeg",
     width=1000, height=750)
ggplot(OAContMEM_expandvar1, aes(es2-es1, N, fill=pnexch_cutoff797_power)) + 
  geom_tile() + geom_text(aes(label=pnexch_cutoff797_power), size = 7) + 
  labs(title = "Heatmap of One-Arm MEM P(Not Exchangeable) >= 0.7976455 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEM/PowerPNEXCH_OAContMEM_Cutoff08_ExpandVar1.jpeg",
     width=1000, height=750)
ggplot(OAContMEM_expandvar1, aes(es2-es1, N, fill=pnexch_cutoff8_power)) + 
  geom_tile() + geom_text(aes(label=pnexch_cutoff8_power), size = 7) + 
  labs(title = "Heatmap of One-Arm MEM P(Not Exchangeable) >= 0.8 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()


jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-EqualSamp/MEM/PowerPNEXCH_OAContMEM_Cutoff09_ExpandVar1.jpeg",
     width=1000, height=750)
ggplot(OAContMEM_expandvar1, aes(es2-es1, N, fill=pnexch_cutoff9_power)) + 
  geom_tile() + geom_text(aes(label=pnexch_cutoff9_power), size = 7) + 
  labs(title = "Heatmap of One-Arm MEM P(Not Exchangeable) >= 0.9 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()




#############################################
##############UNEQUAL SAMPLE#################
#############################################
OAContMEM_unequalsamp <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/OA_cont_MEM_diffn_result.csv")

#create character variable for nice labels on plot
OAContMEM_unequalsamp$n2mult_c <- NA
OAContMEM_unequalsamp$n2mult_c[OAContMEM_unequalsamp$n2 == 0.1] <- "n2 = 0.1*N"
OAContMEM_unequalsamp$n2mult_c[OAContMEM_unequalsamp$n2 == 0.25] <- "n2 = 0.25*N"

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-UnequalSamp/MEM/PowerTtest_OAContMEM_UnequalSamp.jpeg",
     width=1050, height=750)
ggplot(OAContMEM_unequalsamp, aes(es2-es1, N, fill=power)) +
  geom_tile() + geom_text(aes(label=power), size = 7) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of One-Arm T-test Power for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = 'Type I Error/Power') +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-UnequalSamp/MEM/PowerPNEXCH_OAContMEM_UnequalSamp_Cutoff02.jpeg",
     width=1050, height=750)
ggplot(OAContMEM_unequalsamp, aes(es2-es1, N, fill=pnexch_cutoff2_power)) + 
  geom_tile() + geom_text(aes(label=pnexch_cutoff2_power), size = 7) + facet_grid(n2mult_c~.) + 
  labs(title = "Heatmap of One-Arm MEM P(Not Exchangeable) >= 0.2 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-UnequalSamp/MEM/PowerPNEXCH_OAContMEM_UnequalSamp_Cutoff0798.jpeg",
     width=1050, height=750)
ggplot(OAContMEM_unequalsamp, aes(es2-es1, N, fill=pnexch_cutoff797_power)) + 
  geom_tile() + geom_text(aes(label=pnexch_cutoff797_power), size = 7) + facet_grid(n2mult_c~.) + 
  labs(title = "Heatmap of One-Arm MEM P(Not Exchangeable) >= 0.7976455 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-UnequalSamp/MEM/PowerPNEXCH_OAContMEM_UnequalSamp_Cutoff08.jpeg",
     width=1050, height=750)
ggplot(OAContMEM_unequalsamp, aes(es2-es1, N, fill=pnexch_cutoff8_power)) + 
  geom_tile() + geom_text(aes(label=pnexch_cutoff8_power), size = 7) + facet_grid(n2mult_c~.) + 
  labs(title = "Heatmap of One-Arm MEM P(Not Exchangeable) >= 0.8 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()


jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/Heatmaps-UnequalSamp/MEM/PowerPNEXCH_OAContMEM_UnequalSamp_Cutoff09.jpeg",
     width=1050, height=750)
ggplot(OAContMEM_unequalsamp, aes(es2-es1, N, fill=pnexch_cutoff9_power)) + 
  geom_tile() + geom_text(aes(label=pnexch_cutoff9_power), size = 7) + facet_grid(n2mult_c~.) + 
  labs(title = "Heatmap of One-Arm MEM P(Not Exchangeable) >= 0.9 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()


