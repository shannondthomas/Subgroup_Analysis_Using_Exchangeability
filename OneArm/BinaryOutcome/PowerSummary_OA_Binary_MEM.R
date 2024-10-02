library(ggplot2)
library(tidyverse)

OA_binary_MEM_result <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/OA_binary_MEM_result.csv")
OA_binary_MEM_result <- OA_binary_MEM_result[,2:13] #remove column of row names
OA_binary_MEM_diffn_result <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/BinaryOutcome/resultsdata/OA_binary_MEM_diffn_result.csv")
OA_binary_MEM_diffn_result <- OA_binary_MEM_diffn_result[,2:14] #remove column of row names


cutoffs <- c(0.2,0.7976455,0.8,0.9) #from PowerSimulation_OA_Binary_MEM.R



## DIFFERENCE PLOTS

OA_binary_MEM_result$efdiff <- OA_binary_MEM_result$es2 - OA_binary_MEM_result$es1
results_summary <- OA_binary_MEM_result
p1 <- ggplot(results_summary[results_summary$efdiff != 0,], aes(efdiff, N, fill=round(pnexch_cutoff2_power - power,7))) + 
  geom_tile() + 
  geom_text(aes(label=round(pnexch_cutoff2_power - power,7)),size = 5) + 
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
  scale_fill_gradientn(limits = c(-0.5,0.5), colors = c("red4","white", "green4"), breaks = c(-0.5,0,0.5))

p2 <- ggplot(results_summary[results_summary$efdiff == 0,], aes(as.factor(efdiff), N, fill=round(pnexch_cutoff2_power - power,7))) + 
  geom_tile() + 
  geom_text(aes(label=round(pnexch_cutoff2_power - power,7)),size = 5) + 
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
  scale_fill_gradientn(limits = c(-0.2,0.2), colors = c("green4","white", "red4"), breaks = c(-0.2,0,0.2))


grid.arrange(p2,p1,nrow = 1, widths = c(1,4), top = "One-Arm Binary Outcome: Difference in Type I Error and Power for Exchangeability Method with 0.2 Cutoff vs Prop.Test")


p18 <- ggplot(results_summary[results_summary$efdiff != 0,], aes(efdiff, N, fill=round(pnexch_cutoff8_power - power,7))) + 
  geom_tile() + 
  geom_text(aes(label=round(pnexch_cutoff8_power - power,7)),size = 5) + 
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
  scale_fill_gradientn(limits = c(-0.5,0.5), colors = c("red4","white", "green4"), breaks = c(-0.5,0,0.5))

p28 <- ggplot(results_summary[results_summary$efdiff == 0,], aes(as.factor(efdiff), N, fill=round(pnexch_cutoff8_power - power,7))) + 
  geom_tile() + 
  geom_text(aes(label=round(pnexch_cutoff8_power - power,7)),size = 5) + 
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
  scale_fill_gradientn(limits = c(-0.2,0.2), colors = c("green4","white", "red4"), breaks = c(-0.2,0,0.2))


grid.arrange(p28,p18,nrow = 1, widths = c(1,4), top = "One-Arm Binary Outcome: Difference in Type I Error and Power for Exchangeability Method with 0.8 Cutoff vs Prop.Test")


## Plot results
## equal sample
ggplot(OA_binary_MEM_result, aes(es2-es1, N, fill=power)) +
  geom_tile() + geom_text(aes(label=power)) +
  labs(title = "Heatmap of prop.test Power for Varying Sample Sizes and Differences in Effect Sizes",
       fill = 'Type I Error/Power') +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)")

for(i in 1:length(cutoffs)){
  fillvar <- colnames(OA_binary_MEM_result)[8+i]
  
  print(
    ggplot(OA_binary_MEM_result, aes(es2-es1, N, fill=.data[[fillvar]])) + 
      geom_tile() + geom_text(aes(label=.data[[fillvar]])) + 
      labs(title = paste("Heatmap of MEM P(Not Exchangeable) >=", cutoffs[i], 
                         "for Varying Sample Sizes and Differences in Effect Sizes",
                         sep = ' '),
           fill = "Type I Error/Power") +
      xlab("Difference in Effect Sizes") +
      ylab("Total Sample Size (N)")
  )
}

#unequal sample
ggplot(OA_binary_MEM_diffn_result, aes(es2-es1, N, fill=power)) +
  geom_tile() + geom_text(aes(label=power)) + facet_grid(n2mult_c~.) +
  labs(title = "Heatmap of prop.test Power for Varying Sample Sizes per Group and Differences in Effect Sizes",
       fill = 'Type I Error/Power') +
  xlab("Difference in Effect Sizes") +
  ylab("Total Sample Size (N)")

for(i in 1:length(cutoffs)){
  fillvar <- colnames(OA_binary_MEM_diffn_result)[8+i]
  
  print(
    ggplot(OA_binary_MEM_diffn_result, aes(es2-es1, N, fill=.data[[fillvar]])) + 
      geom_tile() + geom_text(aes(label=.data[[fillvar]])) + facet_grid(n2mult_c~.) +
      labs(title = paste("Heatmap of MEM P(Not Exchangeable) >=", cutoffs[i], 
                         "for Varying Sample Sizes per Group and Differences in Effect Sizes",
                         sep = ' '),
           fill = "Type I Error/Power") +
      xlab("Difference in Effect Sizes") +
      ylab("Total Sample Size (N)")
  )
}


