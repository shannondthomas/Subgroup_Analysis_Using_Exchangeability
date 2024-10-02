library(ggplot2)
library(tidyverse)
library(gridExtra)

#############################################
##############EQUAL SAMPLE###################
#############################################
#results <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/resultsdata/results.csv")
results_raw <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/resultsdata/results_raw.csv")

results_raw <- results_raw[,2:11]
colnames(results_raw)[8] <- 'pnexchgt0798'

#create character variable for nice labels on graphs
results_raw$var_c <- NA
results_raw$var_c[results_raw$var == 1] <- 'Var = 1'
results_raw$var_c[results_raw$var == 10] <- 'Var = 10'
results_raw$var_c[results_raw$var == 100] <- 'Var = 100'

#create difference of effect variable
results_raw$efdiff <- results_raw$source2te - results_raw$source1te

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


#DIFFERENCE PLOTS
p1 <- ggplot(results_summary[results_summary$efdiff != 0,], aes(efdiff, totsampsize, fill=round(pnexchgt020 - plt05,7))) + 
        geom_tile() + 
        geom_text(aes(label=round(pnexchgt020 - plt05,7)),size = 5) + facet_grid(var_c~.) +
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

p2 <- ggplot(results_summary[results_summary$efdiff == 0,], aes(as.factor(efdiff), totsampsize, fill=round(pnexchgt020 - plt05,7))) + 
        geom_tile() + 
        geom_text(aes(label=round(pnexchgt020 - plt05,7)),size = 5) + facet_grid(var_c~.) +
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


grid.arrange(p2,p1,nrow = 1, widths = c(1,4), top = "Two-Arm Continuous Outcome: Difference in Type I Error and Power for Exchangeability Method with 0.2 Cutoff vs F-test")


p18 <- ggplot(results_summary[results_summary$efdiff != 0,], aes(efdiff, totsampsize, fill=round(pnexchgt080 - plt05,7))) + 
  geom_tile() + 
  geom_text(aes(label=round(pnexchgt080 - plt05,7)),size = 5) + facet_grid(var_c~.) +
  xlab("Difference in Effect Size") +
  ylab("") +
  labs(#title = "Power", 
        fill = "Power") +
  theme(text=element_text(size=15),
        plot.title = element_text(size = 10),
        legend.position = "top", 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.width = unit(1, "cm"),
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        plot.margin = unit(c(5.5,5.5,5.5,-1),"pt")) +
  #scale_fill_gradientn(colors = c("red3","white"), values = scales::rescale(c(-0.443,0)))
  #scale_fill_gradient2(high="green4", mid = "white", low="red", midpoint=0 ,na.value="gray20")
  scale_fill_gradientn(limits = c(-0.5,0.5), colors = c("red4","white", "green4"), breaks = c(-0.5,0,0.5))

  
p28 <- ggplot(results_summary[results_summary$efdiff == 0,], aes(as.factor(efdiff), totsampsize, fill=round(pnexchgt080 - plt05,7))) + 
  geom_tile() + 
  geom_text(aes(label=round(pnexchgt080 - plt05,7)),size = 5) + facet_grid(var_c~.) +
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
  #scale_fill_gradientn(colors = c("white","lightpink","red3"), values = scales::rescale(c(-0.028,0,0.187)))
  #scale_fill_gradientn(colors = c("green4","green"), values = scales::rescale(c(-0.06,-0.04)))
  #scale_fill_gradient2(low="green4", mid = "white", high="red", midpoint=0 ,na.value="gray20")
  scale_fill_gradientn(limits = c(-0.1,0.1), colors = c("green4","white", "red4"), breaks = c(-0.1,0,0.1))

grid.arrange(p28,p18,nrow = 1, widths = c(1,4), top = "Two-Arm Continuous Outcome: Difference in Type I Error and Power for Exchangeability Method with 0.8 Cutoff vs F-test")







+jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/Heatmaps-EqualSamp/PowerFtest_TAContMEMr.jpeg",
     width=1000, height=750)
ggplot(results_summary, aes(efdiff, totsampsize, fill=plt05)) + 
  geom_tile() + geom_text(aes(label=plt05),size = 5) + facet_grid(var_c~.) +
  labs(title = "Heatmap of Two-Arm F-test Power for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()


jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/Heatmaps-EqualSamp/PowerPNEXCH_TAContMEMr_02.jpeg",
     width=1000, height=750)
ggplot(results_summary, aes(efdiff, totsampsize, fill=pnexchgt020)) + 
  geom_tile() + geom_text(aes(label=pnexchgt020),size = 5) + facet_grid(var_c~.) +
  labs(title = "Heatmap of Two-Arm MEMr P(Not Exchangeable) > 0.2 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/Heatmaps-EqualSamp/PowerPNEXCH_TAContMEMr_0798.jpeg",
     width=1000, height=750)
ggplot(results_summary, aes(efdiff, totsampsize, fill=pnexchgt0798)) + 
  geom_tile() + geom_text(aes(label=pnexchgt0798),size = 5) + facet_grid(var_c~.) +
  labs(title = "Heatmap of Two-Arm MEMr P(Not Exchangeable) > 0.7976455 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/Heatmaps-EqualSamp/PowerPNEXCH_TAContMEMr_08.jpeg",
     width=1000, height=750)
ggplot(results_summary, aes(efdiff, totsampsize, fill=pnexchgt080)) + 
  geom_tile() + geom_text(aes(label=pnexchgt080),size = 5) + facet_grid(var_c~.) +
  labs(title = "Heatmap of Two-Arm MEMr P(Not Exchangeable) > 0.8 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/Heatmaps-EqualSamp/PowerPNEXCH_TAContMEMr_09.jpeg",
     width=1000, height=750)
ggplot(results_summary, aes(efdiff, totsampsize, fill=pnexchgt090)) + 
  geom_tile() + geom_text(aes(label=pnexchgt090),size = 5) + facet_grid(var_c~.) +
  labs(title = "Heatmap of Two-Arm MEMr P(Not Exchangeable) > 0.9 for Varying Sample Sizes, Variances, and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()




#############################################
#############UNEQUAL SAMPLE##################
#############################################

results_raw_unequalsamp <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/resultsdata/TA_cont_MEMr_result_unequalsamp_raw.csv")

results_raw_unequalsamp <- results_raw_unequalsamp[,2:12]
colnames(results_raw_unequalsamp)[8] <- 'pnexchgt0798'

results_raw_unequalsamp$n2mult_c <- NA
results_raw_unequalsamp$n2mult_c[results_raw_unequalsamp$n2mult == 0.1] <- 'n2 = 0.1*N'
results_raw_unequalsamp$n2mult_c[results_raw_unequalsamp$n2mult == 0.25] <- 'n2 = 0.25*N'

#use raw data to calculate different cutoff performances
results_raw_unequalsamp$pnexchgt020 <- results_raw_unequalsamp$pnexch > 0.2
results_raw_unequalsamp$pnexchgt080 <- results_raw_unequalsamp$pnexch > 0.8
results_raw_unequalsamp$pnexchgt090 <- results_raw_unequalsamp$pnexch > 0.9

#create difference of effect variable
results_raw_unequalsamp$efdiff <- results_raw_unequalsamp$source2te - results_raw_unequalsamp$source1te

results_summary_unequalsamp <- results_raw_unequalsamp %>% group_by(totsampsize,n2mult_c,efdiff) %>% 
  summarize(plt05 = mean(plt05), 
            pnexchgt090 = mean(pnexchgt090), 
            pnexchgt080 = mean(pnexchgt080), 
            pnexchgt0798 = mean(pnexchgt0798), 
            pnexchgt020 = mean(pnexchgt020),
  )

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/Heatmaps-UnequalSamp/PowerFtest_TAContMEMr_UnequalSamp.jpeg",
     width=1000, height=750)
ggplot(results_summary_unequalsamp, aes(efdiff, totsampsize, fill=plt05)) + 
  geom_tile() + geom_text(aes(label=plt05),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of Two-Arm F-test Power for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/Heatmaps-UnequalSamp/PowerPNEXCH_TAContMEMr_UnequalSamp_09.jpeg",
     width=1000, height=750)
ggplot(results_summary_unequalsamp, aes(efdiff, totsampsize, fill=pnexchgt090)) + 
  geom_tile() + geom_text(aes(label=pnexchgt090),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of Two-Arm MEMr P(Not Exchangeable) > 0.9 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/Heatmaps-UnequalSamp/PowerPNEXCH_TAContMEMr_UnequalSamp_08.jpeg",
     width=1000, height=750)
ggplot(results_summary_unequalsamp, aes(efdiff, totsampsize, fill=pnexchgt080)) + 
  geom_tile() + geom_text(aes(label=pnexchgt080),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of Two-Arm MEMr P(Not Exchangeable) > 0.8 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/Heatmaps-UnequalSamp/PowerPNEXCH_TAContMEMr_UnequalSamp_0798.jpeg",
     width=1000, height=750)
ggplot(results_summary_unequalsamp, aes(efdiff, totsampsize, fill=pnexchgt0798)) + 
  geom_tile() + geom_text(aes(label=pnexchgt0798),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of Two-Arm MEMr P(Not Exchangeable) > 0.7976455 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()

jpeg(file="C:/Users/mushanno/Desktop/Work/Dissertation1/TwoArm/ContinuousOutcome/Heatmaps-UnequalSamp/PowerPNEXCH_TABinMEMr_UnequalSamp_02.jpeg",
     width=1000, height=750)
ggplot(results_summary_unequalsamp, aes(efdiff, totsampsize, fill=pnexchgt020)) + 
  geom_tile() + geom_text(aes(label=pnexchgt020),size = 5) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of Two-Arm MEMr P(Not Exchangeable) > 0.2 for Varying Sample Sizes and Differences in Effect Sizes",
       fill = "Type I Error/Power") +
  xlab("Difference in Effect Size") +
  ylab("Total Sample Size (N)") +
  theme(text=element_text(size=15))
dev.off()


