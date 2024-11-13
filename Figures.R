#################################################################################
# TITLE: Figures.R
#
# PURPOSE: Run simulations in parallel using the snowfall package and create
#          a master data output file. 
#
#
# INPUT: simulation_results_summary.csv - contains information for every simulation
#                                      setting including means, number of arms,
#                                      outcome type, and power for each test,
#                                      and all other relevant information
#        simulation_unequaln_results_summary.csv - same as above for unequal n
#
# OUTPUT: lots of heatmaps, 
#         histograms of probability of exchangeability and corresponding p-vals
#
#
# SECTIONS: Section 1
#              heatmap_func() - function to make heatmaps and save to appropriate
#                               figure folder
#           Section 2
#              run heatmap function for all scenarios
#              2.1 One-Arm Binary Outcome
#              2.2 One-Arm Continuous Outcome
#              2.3 Two-Arm Binary Outcome
#              2.4 Two-Arm Continuous Outcome
#           Section 3
#              create histograms 
#
# AUTHOR: Shannon Thomas
# DATE CREATED: OCT 22, 2024
#################################################################################



###########################################################
####### SECTION 0: LOAD DEPENDENCIES & READ IN DATA #######
###########################################################

# alldat <- read.csv("simulation_results_raw.csv")
alldat <- read.csv("simulation_results_summary.csv")
alldat_unequaln <- read.csv("simulation_results_unequaln_summary.csv")
alldat_08 <- read.csv("simulation_results_summary_08cutoff.csv")
alldat_08_unequaln <- read.csv("simulation_results_unequaln_summary_08cutoff.csv")

#R version 4.3.2 was used for this project.
library(tidyverse) #version 2.0.0
library(gridExtra) #version 2.3
library(grid)      #base package


###########################################################
############### SECTION 1: HEATMAP FUNCTION ###############
###########################################################

nsim <- 10000

heatmap_func <- function(dat, num_arms, outcome_type, ht, width, MEM_cutoff, p_cutoff, unequal_n = FALSE) {
  
  
  if(outcome_type == "binary"){
    standard_test_name <- "LRT"
    outcome_type_c <- "Binary"
  }
  if(outcome_type == "continuous"){
    standard_test_name <- "Partial F-test"
    outcome_type_c <- "Continuous"
  }
  
  #subset data to OA/TA and Binary/Continuous setting
  subdat_summary <- dat[((dat$num_arms == num_arms) & (dat$outcome_type == outcome_type)),]
  
  # subdat_summary <- subdat_summary %>% 
  #                     group_by(n1, n2, mean1, mean2, sd1, sd2, trteff1, trteff2) %>%
  #                     summarize(MEM_power = sum(MEMpexch < MEM_cutoff)/nsim,
  #                               MEMr_power = sum(MEMrpexch < MEM_cutoff)/nsim,
  #                               pval_power = sum(pval < p_cutoff)/nsim)
  
  subdat_summary$totsampsize <- subdat_summary$n1 + subdat_summary$n2
  subdat_summary$var_c <- paste("Var = ", subdat_summary$sd1^2)
  
  #define n2mult
  if(unequal_n != FALSE){
    if(MEM_cutoff != 0.2){
      foldername <- paste("Figures/UnequalSamp/", MEM_cutoff, "Cutoff/", sep = "")
    }
    else{
      foldername <- "Figures/UnequalSamp/"
    }
    subdat_summary$n2_prop <- paste("n2 =",round(subdat_summary$n2/(subdat_summary$totsampsize),2),"N")
  }
  else{
    if(MEM_cutoff != 0.2){
      foldername <- paste("Figures/EqualSamp/", MEM_cutoff, "Cutoff/", sep = "")
    }
    else{
      foldername <- "Figures/EqualSamp/"
    }
  }
  
  #define effect size difference for plots
  if(num_arms == 1){
    subdat_summary$efdiff <- subdat_summary$mean2 - subdat_summary$mean1
    num_arm_c <- "One"
  }
  
  if(num_arms == 2){
    subdat_summary$efdiff <- subdat_summary$trteff2 - subdat_summary$trteff1
    num_arm_c <- "Two"
  }

    
  
  
  #First plot is the difference plot for MEMr vs standard test.
  #This plot has two pieces to allow for opposite coloring of power and t1e
    p1 <- ggplot(subdat_summary[subdat_summary$efdiff != 0,], aes(efdiff, totsampsize, fill=round(MEMr_power - pval_power,7))) + 
      geom_tile() + 
      geom_text(aes(label=round(MEMr_power - pval_power,digits = 4)),size = 5) + 
      {if((outcome_type == "continuous") & (unequal_n == FALSE))facet_grid(var_c~.)} +
      {if(unequal_n != FALSE)facet_grid(n2_prop~.)} +
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
      {if(MEM_cutoff == 0.2)scale_fill_gradientn(limits = c(-0.61,0.61), colors = c("purple4","white", "green4"), breaks = c(-0.61,0,0.61))} +
      {if(MEM_cutoff == 0.8)scale_fill_gradientn(limits = c(-1,1), colors = c("purple4","white", "green4"), breaks = c(-1,0,1))} 
    
    
    p2 <- ggplot(subdat_summary[subdat_summary$efdiff == 0,], aes(as.factor(efdiff), totsampsize, fill=round(MEMr_power - pval_power,7))) + 
      geom_tile() + 
      geom_text(aes(label=round(MEMr_power - pval_power,digits = 4)),size = 5) + 
      {if((outcome_type == "continuous") & (unequal_n == FALSE))facet_grid(var_c~.)} +
      {if(unequal_n != FALSE)facet_grid(n2_prop~.)} +
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
      {if((outcome_type == "continuous")&(MEM_cutoff == 0.2))scale_fill_gradientn(limits = c(-0.6,0.6), colors = c("green4","white", "purple4"), breaks = c(-0.6,0,0.6))} +
      {if((outcome_type == "binary")&(MEM_cutoff == 0.2))scale_fill_gradientn(limits = c(-0.1,0.1), colors = c("green4","white", "purple4"), breaks = c(-0.1,0,0.1))} +
      {if((outcome_type == "continuous")&(MEM_cutoff == 0.8))scale_fill_gradientn(limits = c(-1,1), colors = c("green4","white", "purple4"), breaks = c(-1,0,1))} +
      {if((outcome_type == "binary")&(MEM_cutoff == 0.8))scale_fill_gradientn(limits = c(-0.3,0.3), colors = c("green4","white", "purple4"), breaks = c(-0.3,0,0.3))}
    
    jpeg(paste(foldername,"DiffMEMr_",num_arm_c,"Arm",outcome_type_c,"_cutoff",MEM_cutoff,".jpg",sep = ""), 
         height = 0.25*ht, width = 0.25*width)
      grid.arrange(p2,p1,nrow = 1, widths = c(1,4), 
                   top = textGrob(paste(num_arm_c, "Arm", outcome_type_c,
                                        "Outcome: Difference in Type I Error and Power for \nMEMr Exchangeability Method with", 
                                         MEM_cutoff, "Cutoff vs", standard_test_name), 
                                  gp = gpar(fontsize = 14)))
    dev.off()
    
    #Second plot is the difference plot for MEM vs standard test.
    #This plot has two pieces to allow for opposite coloring of power and t1e
    if(num_arms == 1){
      p1_MEM <- ggplot(subdat_summary[subdat_summary$efdiff != 0,], aes(efdiff, totsampsize, fill=round(MEM_power - pval_power,7))) + 
        geom_tile() + 
        geom_text(aes(label=round(MEM_power - pval_power,4)),size = 5) + 
        {if((outcome_type == "continuous") & (unequal_n == FALSE))facet_grid(var_c~.)} +
        {if(unequal_n != FALSE)facet_grid(n2_prop~.)} +
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
        {if(MEM_cutoff == 0.2)scale_fill_gradientn(limits = c(-0.61,0.61), colors = c("purple4","white", "green4"), breaks = c(-0.61,0,0.61))} + 
        {if(MEM_cutoff == 0.8)scale_fill_gradientn(limits = c(-1,1), colors = c("purple4","white", "green4"), breaks = c(-1,0,1))} 
        
      p2_MEM <- ggplot(subdat_summary[subdat_summary$efdiff == 0,], aes(as.factor(efdiff), totsampsize, fill=round(MEM_power - pval_power,7))) + 
        geom_tile() + 
        geom_text(aes(label=round(MEM_power - pval_power,4)),size = 5) + 
        {if((outcome_type == "continuous") & (unequal_n == FALSE))facet_grid(var_c~.)} +
        {if(unequal_n != FALSE)facet_grid(n2_prop~.)} +
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
        {if((outcome_type == "continuous")&(MEM_cutoff == 0.2))scale_fill_gradientn(limits = c(-0.6,0.6), colors = c("green4","white", "purple4"), breaks = c(-0.6,0,0.6))} +
        {if((outcome_type == "binary")&(MEM_cutoff == 0.2))scale_fill_gradientn(limits = c(-0.1,0.1), colors = c("green4","white", "purple4"), breaks = c(-0.1,0,0.1))} +
        {if((outcome_type == "continuous")&(MEM_cutoff == 0.8))scale_fill_gradientn(limits = c(-1,1), colors = c("green4","white", "purple4"), breaks = c(-1,0,1))} +
        {if((outcome_type == "binary")&(MEM_cutoff == 0.8))scale_fill_gradientn(limits = c(-0.3,0.3), colors = c("green4","white", "purple4"), breaks = c(-0.3,0,0.3))}
      
      jpeg(paste(foldername, "DiffMEM_",num_arm_c,"Arm",outcome_type_c,"_cutoff",MEM_cutoff,".jpg",sep = ""), 
           height = 0.25*ht, width = 0.25*width)
        grid.arrange(p2_MEM,p1_MEM,nrow = 1, widths = c(1,4), 
                   top = textGrob(paste(num_arm_c, "Arm", outcome_type_c,
                                        "Outcome: Difference in Type I Error and Power for \nMEM Exchangeability Method with", 
                                        MEM_cutoff, "Cutoff vs", standard_test_name), 
                                  gp = gpar(fontsize = 14)))
      dev.off()
      }
    
    
    
    
    #The next three plots are the heatmaps with t1e rates and power 
    #for MEM (if applicable), MEMr, and standard test. 
    
    if(num_arms == 1){
            ggplot(subdat_summary, aes(efdiff, totsampsize, fill=MEM_power)) + 
              geom_tile() + geom_text(aes(label=round(MEM_power, digits = 4)),size = 7) + 
              {if((outcome_type == "continuous") & (unequal_n == FALSE))facet_grid(var_c~.)} +
              {if(unequal_n != FALSE)facet_grid(n2_prop~.)} +
              labs(title = paste("Heatmap of MEM Power for Varying Sample Sizes and Differences \nin Effect Sizes in", 
                                 num_arm_c, "Arm", outcome_type_c, "Outcome Scenario"),
                   fill = "Type I Error/\nPower") +
              xlab("Difference in Effect Sizes") +
              ylab("Total Sample Size (N)") +
              theme(text=element_text(size=12), plot.title = element_text(hjust = 0.5), legend.position = "top") +
              scale_fill_gradientn(limits = c(0,1), colors = c("white","dodgerblue3"), breaks = c(0,1))
      ggsave(paste(foldername,"PowerMEM_",num_arm_c,"Arm",outcome_type_c,"_Cutoff", MEM_cutoff,".jpg",sep = ""), 
             height = ht, width = width, units = "px")
    }
    
    #standard test
          ggplot(subdat_summary, aes(efdiff, totsampsize, fill=pval_power)) + 
            #geom_tile() + geom_text(aes(label=pval_power),size = 7) + 
            geom_tile() + geom_text(aes(label=round(pval_power, digits = 4)),size = 7) + 
            {if((outcome_type == "continuous") & (unequal_n == FALSE))facet_grid(var_c~.)} +
            {if(unequal_n != FALSE)facet_grid(n2_prop~.)} +
            labs(title = paste("Heatmap of", standard_test_name,"Power for Varying Sample Sizes and Differences \nin Effect Sizes in", 
                               num_arm_c, "Arm", outcome_type_c, "Outcome Scenario"),
                 fill = "Type I Error/\nPower") +
            xlab("Difference in Effect Sizes") +
            ylab("Total Sample Size (N)") +
            theme(text=element_text(size=12), plot.title = element_text(hjust = 0.5), legend.position = "top") +
            scale_fill_gradientn(limits = c(0,1), colors = c("white","dodgerblue3"), breaks = c(0,1))
    ggsave(paste(foldername,"Power",standard_test_name,"_",num_arm_c,"Arm",outcome_type_c,".jpg",sep = ""), 
           height = ht, width = width, units = "px")
    
    #MEMr
          ggplot(subdat_summary, aes(efdiff, totsampsize, fill=MEMr_power)) + 
            geom_tile() + geom_text(aes(label=round(MEMr_power, digits = 4)),size = 7) + 
            {if((outcome_type == "continuous") & (unequal_n == FALSE))facet_grid(var_c~.)} +
            {if(unequal_n != FALSE)facet_grid(n2_prop~.)} +
            labs(title = paste("Heatmap of MEMr Power for Varying Sample Sizes and Differences \nin Effect Sizes in", 
                               num_arm_c, "Arm", outcome_type_c, "Outcome Scenario"),
                 fill = "Type I Error/\nPower") +
            xlab("Difference in Effect Sizes") +
            ylab("Total Sample Size (N)") +
            theme(text=element_text(size=12), plot.title = element_text(hjust = 0.5), legend.position = "top") +
            scale_fill_gradientn(limits = c(0,1), colors = c("white","dodgerblue3"), breaks = c(0,1))
    ggsave(paste(foldername, "PowerMEMr_",num_arm_c,"Arm",outcome_type_c,"_Cutoff", MEM_cutoff,".jpg",sep = ""), 
           height = ht, width = width, units = "px")
}







############################################################
################# SECTION 3: MAKE HEATMAPs #################
############################################################

#0.2 cutoff
heatmap_func(alldat, num_arms = 1, outcome_type = "binary", ht = 2000, width = 3000, MEM_cutoff = 0.2, p_cutoff = 0.05)
heatmap_func(alldat, num_arms = 1, outcome_type = "continuous", ht = 3000, width = 3200, MEM_cutoff = 0.2, p_cutoff = 0.05)
heatmap_func(alldat, num_arms = 2, outcome_type = "binary", ht = 2000, width = 3000, MEM_cutoff = 0.2, p_cutoff = 0.05)
heatmap_func(alldat, num_arms = 2, outcome_type = "continuous", ht = 3000, width = 3200, MEM_cutoff = 0.2, p_cutoff = 0.05)


heatmap_func(alldat_unequaln, num_arms = 1, outcome_type = "binary", ht = 2000, width = 3000, MEM_cutoff = 0.2, p_cutoff = 0.05, 
             unequal_n = TRUE)
heatmap_func(alldat_unequaln, num_arms = 1, outcome_type = "continuous", ht = 2000, width = 3200, MEM_cutoff = 0.2, p_cutoff = 0.05, 
             unequal_n = TRUE)
heatmap_func(alldat_unequaln, num_arms = 2, outcome_type = "binary", ht = 2000, width = 3000, MEM_cutoff = 0.2, p_cutoff = 0.05, 
             unequal_n = TRUE)
heatmap_func(alldat_unequaln, num_arms = 2, outcome_type = "continuous", ht = 2000, width = 3200, MEM_cutoff = 0.2, p_cutoff = 0.05, 
             unequal_n = TRUE)


#0.8 cutoff
heatmap_func(alldat_08, num_arms = 1, outcome_type = "binary", ht = 2000, width = 3000, MEM_cutoff = 0.8, p_cutoff = 0.05)
heatmap_func(alldat_08, num_arms = 1, outcome_type = "continuous", ht = 3000, width = 3200, MEM_cutoff = 0.8, p_cutoff = 0.05)
heatmap_func(alldat_08, num_arms = 2, outcome_type = "binary", ht = 2000, width = 3000, MEM_cutoff = 0.8, p_cutoff = 0.05)
heatmap_func(alldat_08, num_arms = 2, outcome_type = "continuous", ht = 3000, width = 3200, MEM_cutoff = 0.8, p_cutoff = 0.05)


heatmap_func(alldat_08_unequaln, num_arms = 1, outcome_type = "binary", ht = 2000, width = 3000, MEM_cutoff = 0.8, p_cutoff = 0.05, 
             unequal_n = TRUE)
heatmap_func(alldat_08_unequaln, num_arms = 1, outcome_type = "continuous", ht = 2000, width = 3200, MEM_cutoff = 0.8, p_cutoff = 0.05, 
             unequal_n = TRUE)
heatmap_func(alldat_08_unequaln, num_arms = 2, outcome_type = "binary", ht = 2000, width = 3000, MEM_cutoff = 0.8, p_cutoff = 0.05, 
             unequal_n = TRUE)
heatmap_func(alldat_08_unequaln, num_arms = 2, outcome_type = "continuous", ht = 2000, width = 3200, MEM_cutoff = 0.8, p_cutoff = 0.05, 
             unequal_n = TRUE)















