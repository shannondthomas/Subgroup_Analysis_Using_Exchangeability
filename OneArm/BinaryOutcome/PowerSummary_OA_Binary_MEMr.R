library(ggplot2)
library(tidyverse)

#############################################
##############EQUAL SAMPLE###################
#############################################

results <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/OA_binary_MEMr_result_summary.csv")
# results_raw <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/OA_binary_MEMr_result_raw.csv")
# results_raw <- results_raw[,2:ncol(results_raw)]

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/Heatmaps-EqualSamp/MEMr/PowerLRT_OABinMEMr.jpeg",
    width=1000, height=500)
  ggplot(results, aes(pt-pc, totsampsize, fill=plt05)) + 
    geom_tile() + geom_text(aes(label=plt05),size = 5) +
    labs(title = "Heatmap of LRT Power for Varying Sample Sizes and Differences in Effect Sizes",
         fill = "Type I Error/Power") +
    xlab("Difference in Effect Sizes") +
    ylab("Total Sample Size (N)") +
    theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/Heatmaps-EqualSamp/MEMr/PowerPNEXCH_OABinMEMr_Cutoff09.jpeg",
     width=1000, height=500)
ggplot(results, aes(pt-pc, totsampsize, fill=pnexchgt090)) + 
  geom_tile() + geom_text(aes(label=pnexchgt090),size = 5) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.9 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/Heatmaps-EqualSamp/MEMr/PowerPNEXCH_OABinMEMr_Cutoff08.jpeg",
     width=1000, height=500)
ggplot(results, aes(pt-pc, totsampsize, fill=pnexchgt080)) + 
  geom_tile() + geom_text(aes(label=pnexchgt080),size = 5) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.8 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/Heatmaps-EqualSamp/MEMr/PowerPNEXCH_OABinMEMr_Cutoff0798.jpeg",
     width=1000, height=500)
ggplot(results, aes(pt-pc, totsampsize, fill=pnexchgt0798)) + 
  geom_tile() + geom_text(aes(label=pnexchgt0798),size = 5) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.7976455 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/Heatmaps-EqualSamp/MEMr/PowerPNEXCH_OABinMEMr_Cutoff02.jpeg",
     width=1000, height=500)
ggplot(results, aes(pt-pc, totsampsize, fill=pnexchgt020)) + 
  geom_tile() + geom_text(aes(label=pnexchgt020),size = 5) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.2 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()



#############################################
#############UNEQUAL SAMPLE##################
#############################################

results_unequalsamp <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/OA_binary_MEMr_unequalsamp_result_summary.csv")
# results_raw_unequalsamp <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/OA_binary_MEMr_unequalsamp_result_raw.csv")
# results_raw_unequalsamp <- results_raw_unequalsamp[,2:10]


results_unequalsamp$n2mult_c <- NA
results_unequalsamp$n2mult_c[results_unequalsamp$n2mult == 0.1] <- 'n2 = 0.1*N'
results_unequalsamp$n2mult_c[results_unequalsamp$n2mult == 0.25] <- 'n2 = 0.25*N'

# #use raw data to calculate different cutoff performances
# results_unequalsamp$pnexchgt020 <- results_unequalsamp$pnexch > 0.2
# results_unequalsamp$pnexchgt080 <- results_raw_unequalsamp$pnexch > 0.8
# results_raw_unequalsamp$pnexchgt090 <- results_raw_unequalsamp$pnexch > 0.9
# 
# results_summary_unequalsamp <- results_raw_unequalsamp %>% group_by(totsampsize,n2mult_c,s1b0,s2b0) %>% 
#   summarize(plt05 = mean(plt05), 
#             pnexchgt090 = mean(pnexchgt090), 
#             pnexchgt080 = mean(pnexchgt080), 
#             pnexchgt0798 = mean(pnexchgt0798), 
#             pnexchgt020 = mean(pnexchgt020),
#   )

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/Heatmaps-UnequalSamp/MEMr/PowerLRT_OABinMEMr_UnequalSamp.jpeg",
     width=1000, height=750)
ggplot(results_unequalsamp, aes(pt-pc, totsampsize, fill=plt05)) + 
  geom_tile() + geom_text(aes(label=plt05),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of LRT Power for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/Heatmaps-UnequalSamp/MEMr/PowerPNEXCH_OABinMEMr_UnequalSamp_09.jpeg",
     width=1000, height=750)
ggplot(results_unequalsamp, aes(pt-pc, totsampsize, fill=pnexchgt090)) + 
  geom_tile() + geom_text(aes(label=pnexchgt090),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.9 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/Heatmaps-UnequalSamp/MEMr/PowerPNEXCH_OABinMEMr_UnequalSamp_08.jpeg",
     width=1000, height=750)
ggplot(results_unequalsamp, aes(pt-pc, totsampsize, fill=pnexchgt080)) + 
  geom_tile() + geom_text(aes(label=pnexchgt080),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.8 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/Heatmaps-UnequalSamp/MEMr/PowerPNEXCH_OABinMEMr_UnequalSamp_0798.jpeg",
     width=1000, height=750)
ggplot(results_unequalsamp, aes(pt-pc, totsampsize, fill=pnexchgt0798)) + 
  geom_tile() + geom_text(aes(label=pnexchgt0798),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.7976455 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/Heatmaps-UnequalSamp/MEMr/PowerPNEXCH_OABinMEMr_UnequalSamp_02.jpeg",
     width=1000, height=750)
ggplot(results_unequalsamp, aes(pt-pc, totsampsize, fill=pnexchgt020)) + 
  geom_tile() + geom_text(aes(label=pnexchgt020),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.2 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()


