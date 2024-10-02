library(ggplot2)
library(tidyverse)
library(gridExtra)

#########unequal samp
results <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/resultsdata/TA_binary_MEMr_result_summary.csv")
#results_raw <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/resultsdata/TA_binary_MEMr_result_raw.csv")
#results_raw <- results_raw[,2:11]

# #use raw data to calculate different cutoff performances
# results_raw$pnexchgt020 <- results_raw$pnexch > 0.2
# results_raw$pnexchgt080 <- results_raw$pnexch > 0.8
# results_raw$pnexchgt090 <- results_raw$pnexch > 0.9
# 
# results_summary <- results_raw %>% group_by(totsampsize,var,source1te,source2te) %>% 
#   summarize(plt05 = mean(plt05), 
#             pnexchgt090 = mean(pnexchgt090), 
#             pnexchgt080 = mean(pnexchgt080), 
#             pnexchgt0798 = mean(pnexchgt0798), 
#             pnexchgt020 = mean(pnexchgt020),
#   )


#DIFFERENCE PLOTS
results$efdiff <- (results$pt2-results$pc2) - (results$pt1-results$pc1)

p1 <- ggplot(results[results$efdiff != 0,], aes(efdiff, totsampsize, fill=round(pnexchgt020 - plt05,7))) + 
  geom_tile() + 
  geom_text(aes(label=round(pnexchgt020 - plt05,7)),size = 5) + 
  xlab("Difference in Effect Size") +
  ylab("") +
  labs(#title= "Power",
    fill = "Power") +
  theme(text=element_text(size=15),
        plot.title = element_text(size = 10),
        legend.position = "top", 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.6, "cm"),
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        plot.margin = unit(c(5.5,5.5,5.5,-1),"pt")) +
  #scale_fill_gradientn(colors = c("red3","white","green4"), values = scales::rescale(c(-0.028,0,0.187)))
  scale_fill_gradientn(limits = c(-0.2,0.2), colors = c("red4","white", "green4"), breaks = c(-0.2,0,0.2))

p2 <- ggplot(results[results$efdiff == 0,], aes(as.factor(efdiff), totsampsize, fill=round(pnexchgt020 - plt05,7))) + 
  geom_tile() + 
  geom_text(aes(label=round(pnexchgt020 - plt05,7)),size = 5) + 
  labs(#title = "Type I Error",
    fill = "Type I \nError") +
  xlab("") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15),
        plot.title = element_text(size = 10), 
        legend.position = "top", 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.3, "cm"),
        strip.text = element_blank(),
        plot.margin = unit(c(5.5,-1,5.5,5.5),"pt")) +
  #scale_fill_gradientn(colors = c("green4","white","red3"), values = scales::rescale(c(-0.028,0,0.187)))
  scale_fill_gradientn(limits = c(-0.1,0.1), colors = c("green4","white", "red4"), breaks = c(-0.1,0,0.1))


grid.arrange(p2,p1,nrow = 1, widths = c(1,4), top = "Two-Arm Binary Outcome: Difference in Type I Error and Power for Exchangeability Method with 0.2 Cutoff vs LRT")



jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/Heatmaps-EqualSamp/PowerLRT_TABinMEMr.jpeg",
     width=1300, height=750)
ggplot(results, aes((pt2-pc2) - (pt1-pc1), totsampsize, fill=plt05)) + 
  geom_tile() + geom_text(aes(label=plt05),size = 7) +
  labs(title = "Heatmap of LRT Power for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=20))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/Heatmaps-EqualSamp/PowerPNEXCH_TABinMEMr_Cutoff09.jpeg",
     width=1300, height=750)
ggplot(results, aes((pt2-pc2) - (pt1-pc1), totsampsize, fill=pnexchgt090)) + 
  geom_tile() + geom_text(aes(label=pnexchgt090), size = 7) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.9 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=20))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/Heatmaps-EqualSamp/PowerPNEXCH_TABinMEMr_Cutoff08.jpeg",
     width=1300, height=750)
ggplot(results, aes((pt2-pc2) - (pt1-pc1), totsampsize, fill=pnexchgt080)) + 
  geom_tile() + geom_text(aes(label=pnexchgt080), size = 7) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.8 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=20))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/Heatmaps-EqualSamp/PowerPNEXCH_TABinMEMr_Cutoff0798.jpeg",
     width=1300, height=750)
ggplot(results, aes((pt2-pc2) - (pt1-pc1), totsampsize, fill=pnexchgt0798)) + 
  geom_tile() + geom_text(aes(label=pnexchgt0798),size = 7) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.7976455 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=20))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/Heatmaps-EqualSamp/PowerPNEXCH_TABinMEMr_Cutoff02.jpeg",
     width=1300, height=750)
ggplot(results, aes((pt2-pc2) - (pt1-pc1), totsampsize, fill=pnexchgt020)) + 
  geom_tile() + geom_text(aes(label=pnexchgt020), size = 7) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.2 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=20))
dev.off()




#########unequal samp
results_unequalsamp <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/resultsdata/TA_binary_MEMr_result_unequalsamp_summary.csv")
# results_raw_unequalsamp <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/resultsdata/TA_binary_MEMr_result_unequalsamp_raw.csv")

results_unequalsamp$n2mult_c <- NA
results_unequalsamp$n2mult_c[results_unequalsamp$n2mult == 0.1] <- 'n2 = 0.1*N'
results_unequalsamp$n2mult_c[results_unequalsamp$n2mult == 0.25] <- 'n2 = 0.25*N'

# results_raw_unequalsamp <- results_raw_unequalsamp[,2:11]
# colnames(results_raw_unequalsamp)[8] <- 'pnexchgt0798'
# 
# results_raw_unequalsamp$n2mult_c <- NA
# results_raw_unequalsamp$n2mult_c[results_raw_unequalsamp$n2mult == 0.1] <- 'n2 = 0.1*N'
# results_raw_unequalsamp$n2mult_c[results_raw_unequalsamp$n2mult == 0.25] <- 'n2 = 0.25*N'
# 
# #use raw data to calculate different cutoff performances
# results_raw_unequalsamp$pnexchgt020 <- results_raw_unequalsamp$pnexch > 0.2
# results_raw_unequalsamp$pnexchgt080 <- results_raw_unequalsamp$pnexch > 0.8
# results_raw_unequalsamp$pnexchgt090 <- results_raw_unequalsamp$pnexch > 0.9
# 
# results_summary_unequalsamp <- results_raw_unequalsamp %>% group_by(totsampsize,n2mult_c,source1te,source2te) %>% 
#   summarize(plt05 = mean(plt05), 
#             pnexchgt090 = mean(pnexchgt090), 
#             pnexchgt080 = mean(pnexchgt080), 
#             pnexchgt0798 = mean(pnexchgt0798), 
#             pnexchgt020 = mean(pnexchgt020),
#   )

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/Heatmaps-UnequalSamp/PowerLRT_TABinMEMr_UnequalSamp.jpeg",
     width=1300, height=750)
ggplot(results_unequalsamp, aes((pt2-pc2) - (pt1-pc1), totsampsize, fill=plt05)) + 
  geom_tile() + geom_text(aes(label=plt05),size = 7) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of LRT Power for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=20))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/Heatmaps-UnequalSamp/PowerPNEXCH_TABinMEMr_UnequalSamp_Cutoff09.jpeg",
     width=1300, height=750)
ggplot(results_unequalsamp, aes((pt2-pc2) - (pt1-pc1), totsampsize, fill=pnexchgt090)) + 
  geom_tile() + geom_text(aes(label=pnexchgt090), size = 7) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.9 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=20))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/Heatmaps-UnequalSamp/PowerPNEXCH_TABinMEMr_UnequalSamp_Cutoff08.jpeg",
     width=1300, height=750)
ggplot(results_unequalsamp, aes((pt2-pc2) - (pt1-pc1), totsampsize, fill=pnexchgt080)) + 
  geom_tile() + geom_text(aes(label=pnexchgt080), size = 7) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.8 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=20))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/Heatmaps-UnequalSamp/PowerPNEXCH_TABinMEMr_UnequalSamp_Cutoff0798.jpeg",
     width=1300, height=750)
ggplot(results_unequalsamp, aes((pt2-pc2) - (pt1-pc1), totsampsize, fill=pnexchgt0798)) + 
  geom_tile() + geom_text(aes(label=pnexchgt0798),size = 7) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.7976455 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=20))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/BinaryOutcome/Heatmaps-UnequalSamp/PowerPNEXCH_TABinMEMr_UnequalSamp_Cutoff02.jpeg",
     width=1300, height=750)
ggplot(results_unequalsamp, aes((pt2-pc2) - (pt1-pc1), totsampsize, fill=pnexchgt020)) + 
  geom_tile() + geom_text(aes(label=pnexchgt020), size = 7) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of MEMr P(Not Exchangeable) > 0.2 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=20))
dev.off()

