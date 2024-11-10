# Supplementary materials analysis and figure plotting

rm(list=ls())

#load packages
{
  library(boot)
  library(tidyverse)
  library(multipanelfigure)
  library(ggpubr)
}

#######################################################################################################################################

#load data - requires download from repository at https://zenodo.org/records/12820416
load("All_model_outputs_20240531.RData")

#-------------------------------------------------------------------------------
#Plot posterior distributions for all parameters comparing between models with different data inclusions

{
  MR3_all_gamma<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$gamma, xmax=prod.mcmc_A_base$q97.5$gamma, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$gamma, xmax=prod.mcmc_A_base_long$q97.5$gamma, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$gamma, xmax=prod.mcmc_A_base_met$q97.5$gamma, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$gamma, xmax=prod.mcmc_A_base_long_met$q97.5$gamma, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$gamma), aes(x=prod.mcmc_A_base$sims.list$gamma), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$gamma), aes(x=prod.mcmc_A_base_long$sims.list$gamma), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$gamma), aes(x=prod.mcmc_A_base_met$sims.list$gamma), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$gamma), aes(x=prod.mcmc_A_base_long_met$sims.list$gamma), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$gamma), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$gamma), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$gamma), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$gamma), col=5, linewidth=1)+
    labs(x=expression(paste("Immune period ", gamma, " (Days)")), y="Density")+xlim(100,400)+
    theme_bw()
  
  MR3_all_alpha<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$alpha, xmax=prod.mcmc_A_base$q97.5$alpha, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$alpha, xmax=prod.mcmc_A_base_long$q97.5$alpha, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$alpha, xmax=prod.mcmc_A_base_met$q97.5$alpha, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$alpha, xmax=prod.mcmc_A_base_long_met$q97.5$alpha, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$alpha), aes(x=prod.mcmc_A_base$sims.list$alpha), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$alpha), aes(x=prod.mcmc_A_base_long$sims.list$alpha), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$alpha), aes(x=prod.mcmc_A_base_met$sims.list$alpha), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$alpha), aes(x=prod.mcmc_A_base_long_met$sims.list$alpha), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$alpha), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$alpha), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$alpha), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$alpha), col=5, linewidth=1)+
    labs(x=expression(paste("Infectious period ", alpha, " (Days)")), y="Density")+xlim(3,12)+
    theme_bw()
  
  
  
  
  MR3_all_delta<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$delta, xmax=prod.mcmc_A_base$q97.5$delta, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$delta, xmax=prod.mcmc_A_base_long$q97.5$delta, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$delta, xmax=prod.mcmc_A_base_met$q97.5$delta, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$delta, xmax=prod.mcmc_A_base_long_met$q97.5$delta, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$delta), aes(x=prod.mcmc_A_base$sims.list$delta), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$delta), aes(x=prod.mcmc_A_base_long$sims.list$delta), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$delta), aes(x=prod.mcmc_A_base_met$sims.list$delta), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$delta), aes(x=prod.mcmc_A_base_long_met$sims.list$delta), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$delta), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$delta), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$delta), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$delta), col=5, linewidth=1)+
    labs(x=expression(paste("Incubation period ", delta, " (Days)")), y="Density")+xlim(0,7)+
    theme_bw()
  
  MR3_all_cull_int<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$cull_int, xmax=prod.mcmc_A_base$q97.5$cull_int, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$cull_int, xmax=prod.mcmc_A_base_long$q97.5$cull_int, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$cull_int, xmax=prod.mcmc_A_base_met$q97.5$cull_int, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$cull_int, xmax=prod.mcmc_A_base_long_met$q97.5$cull_int, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$cull_int), aes(x=prod.mcmc_A_base$sims.list$cull_int), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$cull_int), aes(x=prod.mcmc_A_base_long$sims.list$cull_int), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$cull_int), aes(x=prod.mcmc_A_base_met$sims.list$cull_int), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$cull_int), aes(x=prod.mcmc_A_base_long_met$sims.list$cull_int), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$cull_int), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$cull_int), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$cull_int), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$cull_int), col=5, linewidth=1)+
    labs(x="Culling scaling factor", y="Density")+xlim(0,0.006)+
    theme_bw()
  
  MR3_all_beta0_south<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$beta0_south, xmax=prod.mcmc_A_base$q97.5$beta0_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$beta0_south, xmax=prod.mcmc_A_base_long$q97.5$beta0_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$beta0_south, xmax=prod.mcmc_A_base_met$q97.5$beta0_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta0_south, xmax=prod.mcmc_A_base_long_met$q97.5$beta0_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$beta0_south), aes(x=prod.mcmc_A_base$sims.list$beta0_south), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$beta0_south), aes(x=prod.mcmc_A_base_long$sims.list$beta0_south), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$beta0_south), aes(x=prod.mcmc_A_base_met$sims.list$beta0_south), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta0_south), aes(x=prod.mcmc_A_base_long_met$sims.list$beta0_south), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$beta0_south), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$beta0_south), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$beta0_south), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$beta0_south), col=5, linewidth=1)+
    labs(x="beta0 south", y="Density")+xlim(2,16)+
    theme_bw()
  
  MR3_all_beta1_south<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$beta1_south, xmax=prod.mcmc_A_base$q97.5$beta1_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$beta1_south, xmax=prod.mcmc_A_base_long$q97.5$beta1_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$beta1_south, xmax=prod.mcmc_A_base_met$q97.5$beta1_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta1_south, xmax=prod.mcmc_A_base_long_met$q97.5$beta1_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$beta1_south), aes(x=prod.mcmc_A_base$sims.list$beta1_south), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$beta1_south), aes(x=prod.mcmc_A_base_long$sims.list$beta1_south), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$beta1_south), aes(x=prod.mcmc_A_base_met$sims.list$beta1_south), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta1_south), aes(x=prod.mcmc_A_base_long_met$sims.list$beta1_south), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$beta1_south), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$beta1_south), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$beta1_south), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$beta1_south), col=5, linewidth=1)+
    labs(x="beta1 south", y="Density")+xlim(0,1)+
    theme_bw()
  
  MR3_all_beta2_south<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$beta2_south, xmax=prod.mcmc_A_base$q97.5$beta2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$beta2_south, xmax=prod.mcmc_A_base_long$q97.5$beta2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$beta2_south, xmax=prod.mcmc_A_base_met$q97.5$beta2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta2_south, xmax=prod.mcmc_A_base_long_met$q97.5$beta2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$beta2_south), aes(x=prod.mcmc_A_base$sims.list$beta2_south), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$beta2_south), aes(x=prod.mcmc_A_base_long$sims.list$beta2_south), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$beta2_south), aes(x=prod.mcmc_A_base_met$sims.list$beta2_south), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta2_south), aes(x=prod.mcmc_A_base_long_met$sims.list$beta2_south), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$beta2_south), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$beta2_south), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$beta2_south), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$beta2_south), col=5, linewidth=1)+
    labs(x="beta2 south", y="Density")+xlim(0,1)+
    theme_bw()
  
  MR3_all_period2_south<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$period2_south, xmax=prod.mcmc_A_base$q97.5$period2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$period2_south, xmax=prod.mcmc_A_base_long$q97.5$period2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$period2_south, xmax=prod.mcmc_A_base_met$q97.5$period2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$period2_south, xmax=prod.mcmc_A_base_long_met$q97.5$period2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$period2_south), aes(x=prod.mcmc_A_base$sims.list$period2_south), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$period2_south), aes(x=prod.mcmc_A_base_long$sims.list$period2_south), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$period2_south), aes(x=prod.mcmc_A_base_met$sims.list$period2_south), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$period2_south), aes(x=prod.mcmc_A_base_long_met$sims.list$period2_south), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$period2_south), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$period2_south), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$period2_south), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$period2_south), col=5, linewidth=1)+
    labs(x="Period 2 (Days) south", y="Density")+xlim(1000,2000)+
    theme_bw()
  
  MR3_all_x1_south<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$x1_south, xmax=prod.mcmc_A_base$q97.5$x1_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$x1_south, xmax=prod.mcmc_A_base_long$q97.5$x1_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$x1_south, xmax=prod.mcmc_A_base_met$q97.5$x1_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$x1_south, xmax=prod.mcmc_A_base_long_met$q97.5$x1_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$x1_south), aes(x=prod.mcmc_A_base$sims.list$x1_south), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$x1_south), aes(x=prod.mcmc_A_base_long$sims.list$x1_south), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$x1_south), aes(x=prod.mcmc_A_base_met$sims.list$x1_south), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$x1_south), aes(x=prod.mcmc_A_base_long_met$sims.list$x1_south), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$x1_south), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$x1_south), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$x1_south), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$x1_south), col=5, linewidth=1)+
    labs(x="x1 south", y="Density")+xlim(0,2)+
    theme_bw()
  
  MR3_all_x2_south<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$x2_south, xmax=prod.mcmc_A_base$q97.5$x2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$x2_south, xmax=prod.mcmc_A_base_long$q97.5$x2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$x2_south, xmax=prod.mcmc_A_base_met$q97.5$x2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$x2_south, xmax=prod.mcmc_A_base_long_met$q97.5$x2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$x2_south), aes(x=prod.mcmc_A_base$sims.list$x2_south), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$x2_south), aes(x=prod.mcmc_A_base_long$sims.list$x2_south), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$x2_south), aes(x=prod.mcmc_A_base_met$sims.list$x2_south), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$x2_south), aes(x=prod.mcmc_A_base_long_met$sims.list$x2_south), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$x2_south), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$x2_south), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$x2_south), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$x2_south), col=5, linewidth=1)+
    labs(x="x2 south", y="Density")+xlim(0,2)+
    theme_bw()
  
  
  MR3_all_beta0_north<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$beta0_north, xmax=prod.mcmc_A_base$q97.5$beta0_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$beta0_north, xmax=prod.mcmc_A_base_long$q97.5$beta0_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$beta0_north, xmax=prod.mcmc_A_base_met$q97.5$beta0_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta0_north, xmax=prod.mcmc_A_base_long_met$q97.5$beta0_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$beta0_north), aes(x=prod.mcmc_A_base$sims.list$beta0_north), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$beta0_north), aes(x=prod.mcmc_A_base_long$sims.list$beta0_north), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$beta0_north), aes(x=prod.mcmc_A_base_met$sims.list$beta0_north), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta0_north), aes(x=prod.mcmc_A_base_long_met$sims.list$beta0_north), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$beta0_north), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$beta0_north), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$beta0_north), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$beta0_north), col=5, linewidth=1)+
    labs(x="beta0 north", y="Density")+xlim(2,16)+
    theme_bw()
  
  MR3_all_beta1_north<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$beta1_north, xmax=prod.mcmc_A_base$q97.5$beta1_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$beta1_north, xmax=prod.mcmc_A_base_long$q97.5$beta1_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$beta1_north, xmax=prod.mcmc_A_base_met$q97.5$beta1_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta1_north, xmax=prod.mcmc_A_base_long_met$q97.5$beta1_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$beta1_north), aes(x=prod.mcmc_A_base$sims.list$beta1_north), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$beta1_north), aes(x=prod.mcmc_A_base_long$sims.list$beta1_north), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$beta1_north), aes(x=prod.mcmc_A_base_met$sims.list$beta1_north), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta1_north), aes(x=prod.mcmc_A_base_long_met$sims.list$beta1_north), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$beta1_north), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$beta1_north), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$beta1_north), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$beta1_north), col=5, linewidth=1)+
    labs(x="beta1 north", y="Density")+xlim(0,1)+
    theme_bw()
  
  MR3_all_beta2_north<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$beta2_north, xmax=prod.mcmc_A_base$q97.5$beta2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$beta2_north, xmax=prod.mcmc_A_base_long$q97.5$beta2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$beta2_north, xmax=prod.mcmc_A_base_met$q97.5$beta2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta2_north, xmax=prod.mcmc_A_base_long_met$q97.5$beta2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$beta2_north), aes(x=prod.mcmc_A_base$sims.list$beta2_north), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$beta2_north), aes(x=prod.mcmc_A_base_long$sims.list$beta2_north), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$beta2_north), aes(x=prod.mcmc_A_base_met$sims.list$beta2_north), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta2_north), aes(x=prod.mcmc_A_base_long_met$sims.list$beta2_north), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$beta2_north), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$beta2_north), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$beta2_north), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$beta2_north), col=5, linewidth=1)+
    labs(x="beta2 north", y="Density")+xlim(0,1)+
    theme_bw()
  
  MR3_all_period2_north<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$period2_north, xmax=prod.mcmc_A_base$q97.5$period2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$period2_north, xmax=prod.mcmc_A_base_long$q97.5$period2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$period2_north, xmax=prod.mcmc_A_base_met$q97.5$period2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$period2_north, xmax=prod.mcmc_A_base_long_met$q97.5$period2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$period2_north), aes(x=prod.mcmc_A_base$sims.list$period2_north), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$period2_north), aes(x=prod.mcmc_A_base_long$sims.list$period2_north), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$period2_north), aes(x=prod.mcmc_A_base_met$sims.list$period2_north), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$period2_north), aes(x=prod.mcmc_A_base_long_met$sims.list$period2_north), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$period2_north), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$period2_north), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$period2_north), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$period2_north), col=5, linewidth=1)+
    labs(x="Period 2 (Days) north", y="Density")+xlim(1000,2000)+
    theme_bw()
  
  MR3_all_x1_north<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$x1_north, xmax=prod.mcmc_A_base$q97.5$x1_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$x1_north, xmax=prod.mcmc_A_base_long$q97.5$x1_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$x1_north, xmax=prod.mcmc_A_base_met$q97.5$x1_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$x1_north, xmax=prod.mcmc_A_base_long_met$q97.5$x1_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$x1_north), aes(x=prod.mcmc_A_base$sims.list$x1_north), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$x1_north), aes(x=prod.mcmc_A_base_long$sims.list$x1_north), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$x1_north), aes(x=prod.mcmc_A_base_met$sims.list$x1_north), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$x1_north), aes(x=prod.mcmc_A_base_long_met$sims.list$x1_north), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$x1_north), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$x1_north), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$x1_north), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$x1_north), col=5, linewidth=1)+
    labs(x="x1 north", y="Density")+xlim(0,2)+
    theme_bw()
  
  MR3_all_x2_north<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$x2_north, xmax=prod.mcmc_A_base$q97.5$x2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$x2_north, xmax=prod.mcmc_A_base_long$q97.5$x2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$x2_north, xmax=prod.mcmc_A_base_met$q97.5$x2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$x2_north, xmax=prod.mcmc_A_base_long_met$q97.5$x2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$x2_north), aes(x=prod.mcmc_A_base$sims.list$x2_north), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$x2_north), aes(x=prod.mcmc_A_base_long$sims.list$x2_north), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$x2_north), aes(x=prod.mcmc_A_base_met$sims.list$x2_north), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$x2_north), aes(x=prod.mcmc_A_base_long_met$sims.list$x2_north), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$x2_north), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$x2_north), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$x2_north), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$x2_north), col=5, linewidth=1)+
    labs(x="x2 north", y="Density")+xlim(0,2)+
    theme_bw()
  
  MR3_all_beta0_lima<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$beta0_lima, xmax=prod.mcmc_A_base$q97.5$beta0_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$beta0_lima, xmax=prod.mcmc_A_base_long$q97.5$beta0_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$beta0_lima, xmax=prod.mcmc_A_base_met$q97.5$beta0_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta0_lima, xmax=prod.mcmc_A_base_long_met$q97.5$beta0_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$beta0_lima), aes(x=prod.mcmc_A_base$sims.list$beta0_lima), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$beta0_lima), aes(x=prod.mcmc_A_base_long$sims.list$beta0_lima), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$beta0_lima), aes(x=prod.mcmc_A_base_met$sims.list$beta0_lima), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta0_lima), aes(x=prod.mcmc_A_base_long_met$sims.list$beta0_lima), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$beta0_lima), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$beta0_lima), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$beta0_lima), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$beta0_lima), col=5, linewidth=1)+
    labs(x="beta0 lima", y="Density")+xlim(2,16)+
    theme_bw()
  
  MR3_all_beta1_lima<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$beta1_lima, xmax=prod.mcmc_A_base$q97.5$beta1_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$beta1_lima, xmax=prod.mcmc_A_base_long$q97.5$beta1_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$beta1_lima, xmax=prod.mcmc_A_base_met$q97.5$beta1_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta1_lima, xmax=prod.mcmc_A_base_long_met$q97.5$beta1_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$beta1_lima), aes(x=prod.mcmc_A_base$sims.list$beta1_lima), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$beta1_lima), aes(x=prod.mcmc_A_base_long$sims.list$beta1_lima), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$beta1_lima), aes(x=prod.mcmc_A_base_met$sims.list$beta1_lima), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta1_lima), aes(x=prod.mcmc_A_base_long_met$sims.list$beta1_lima), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$beta1_lima), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$beta1_lima), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$beta1_lima), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$beta1_lima), col=5, linewidth=1)+
    labs(x="beta1 lima", y="Density")+xlim(0,1)+
    theme_bw()
  
  MR3_all_beta2_lima<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$beta2_lima, xmax=prod.mcmc_A_base$q97.5$beta2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$beta2_lima, xmax=prod.mcmc_A_base_long$q97.5$beta2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$beta2_lima, xmax=prod.mcmc_A_base_met$q97.5$beta2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta2_lima, xmax=prod.mcmc_A_base_long_met$q97.5$beta2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$beta2_lima), aes(x=prod.mcmc_A_base$sims.list$beta2_lima), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$beta2_lima), aes(x=prod.mcmc_A_base_long$sims.list$beta2_lima), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$beta2_lima), aes(x=prod.mcmc_A_base_met$sims.list$beta2_lima), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta2_lima), aes(x=prod.mcmc_A_base_long_met$sims.list$beta2_lima), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$beta2_lima), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$beta2_lima), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$beta2_lima), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$beta2_lima), col=5, linewidth=1)+
    labs(x="beta2 lima", y="Density")+xlim(0,1)+
    theme_bw()
  
  MR3_all_period2_lima<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$period2_lima, xmax=prod.mcmc_A_base$q97.5$period2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$period2_lima, xmax=prod.mcmc_A_base_long$q97.5$period2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$period2_lima, xmax=prod.mcmc_A_base_met$q97.5$period2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$period2_lima, xmax=prod.mcmc_A_base_long_met$q97.5$period2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$period2_lima), aes(x=prod.mcmc_A_base$sims.list$period2_lima), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$period2_lima), aes(x=prod.mcmc_A_base_long$sims.list$period2_lima), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$period2_lima), aes(x=prod.mcmc_A_base_met$sims.list$period2_lima), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$period2_lima), aes(x=prod.mcmc_A_base_long_met$sims.list$period2_lima), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$period2_lima), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$period2_lima), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$period2_lima), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$period2_lima), col=5, linewidth=1)+
    labs(x="Period 2 (Days) lima", y="Density")+xlim(1000,2000)+
    theme_bw()
  
  MR3_all_x1_lima<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$x1_lima, xmax=prod.mcmc_A_base$q97.5$x1_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$x1_lima, xmax=prod.mcmc_A_base_long$q97.5$x1_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$x1_lima, xmax=prod.mcmc_A_base_met$q97.5$x1_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$x1_lima, xmax=prod.mcmc_A_base_long_met$q97.5$x1_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$x1_lima), aes(x=prod.mcmc_A_base$sims.list$x1_lima), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$x1_lima), aes(x=prod.mcmc_A_base_long$sims.list$x1_lima), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$x1_lima), aes(x=prod.mcmc_A_base_met$sims.list$x1_lima), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$x1_lima), aes(x=prod.mcmc_A_base_long_met$sims.list$x1_lima), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$x1_lima), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$x1_lima), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$x1_lima), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$x1_lima), col=5, linewidth=1)+
    labs(x="x1 lima", y="Density")+xlim(-1,1)+
    theme_bw()
  
  MR3_all_x2_lima<-ggplot()+
    annotate("rect", xmin=prod.mcmc_A_base$q2.5$x2_lima, xmax=prod.mcmc_A_base$q97.5$x2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=1)+
    annotate("rect", xmin=prod.mcmc_A_base_long$q2.5$x2_lima, xmax=prod.mcmc_A_base_long$q97.5$x2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=2)+
    annotate("rect", xmin=prod.mcmc_A_base_met$q2.5$x2_lima, xmax=prod.mcmc_A_base_met$q97.5$x2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=4)+
    annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$x2_lima, xmax=prod.mcmc_A_base_long_met$q97.5$x2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill=5)+
    geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$x2_lima), aes(x=prod.mcmc_A_base$sims.list$x2_lima), col=1, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long$sims.list$x2_lima), aes(x=prod.mcmc_A_base_long$sims.list$x2_lima), col=2, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_met$sims.list$x2_lima), aes(x=prod.mcmc_A_base_met$sims.list$x2_lima), col=4, linewidth=1)+
    geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$x2_lima), aes(x=prod.mcmc_A_base_long_met$sims.list$x2_lima), col=5, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base$q50$x2_lima), col=1, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long$q50$x2_lima), col=2, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_met$q50$x2_lima), col=4, linewidth=1)+
    geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$q50$x2_lima), col=5, linewidth=1)+
    labs(x="x2 lima", y="Density")+xlim(0,2)+
    theme_bw()
  
}

#plot the figure legend
key<-
ggplot()+
  geom_line(aes(x=c(2, 3), y=c(4, 4)), col=1, linewidth=1.5)+
  geom_line(aes(x=c(2, 3), y=c(3, 3)), col=2, linewidth=1.5)+
  geom_line(aes(x=c(2, 3), y=c(1.8, 1.8)), col=4, linewidth=1.5)+
  geom_line(aes(x=c(2, 3), y=c(0.65, 0.65)), col=5, linewidth=1.5)+
  annotate("text", x = 3.5, y = 4, label = "Population-level serology data",
           hjust = 0, vjust = 0, angle = 0, size = 3)+
  annotate("text", x = 3.5, y = 2.7, label = "Population-level \n+ individual-level serology data",
           hjust = 0, vjust = 0, angle = 0, size = 3)+
  annotate("text", x = 3.5, y = 1.5, label = "Population-level serology \n+ metagenomic data",
           hjust = 0, vjust = 0, angle = 0, size = 3)+
  annotate("text", x = 3.5, y = 0.0, label = "All data: population-level \n+ individual-level serology \n+ metagenomic data",
           hjust = 0, vjust = 0, angle = 0, size = 3)+
  theme_void()+ylim(0,4.5)+xlim(1.5,6)


#plot all parameters not included in the rate of transmission term
comp_figure<-multi_panel_figure(width=170, height=120, rows=2, columns=2)
comp_figure

fill_panel(comp_figure, MR3_all_gamma, col=1, row=1)%>%
  fill_panel(MR3_all_alpha, col=2, row=1)%>%
  fill_panel(MR3_all_delta, col=1, row=2)%>%
  fill_panel(key, col=2, row=2)

####################### these figures are not actually included in the supplement, 
####################### but are here comparing the posteriors for all region-specific
####################### parameters with different data inclusion.

comp_figure_N<-multi_panel_figure(width=170, height=180, rows=3, columns=2)
comp_figure_N

comp_figure_L<-multi_panel_figure(width=170, height=180, rows=3, columns=2)
comp_figure_L

comp_figure_S<-multi_panel_figure(width=170, height=240, rows=4, columns=2)
comp_figure_S

fill_panel(comp_figure_N, MR3_all_beta0_north, col=1, row=1, label = "A")%>%
  fill_panel(MR3_all_beta1_north, col=1, row=2, label = "C")%>%
  fill_panel(MR3_all_beta2_north, col=1, row=3, label = "E")%>%
  fill_panel(MR3_all_period2_north, col=2, row=1, label = "B")%>%
  fill_panel(MR3_all_x1_north, col=2, row=2, label = "D")%>%
  fill_panel(MR3_all_x2_north, col=2, row=3, label = "F")

fill_panel(comp_figure_L, MR3_all_beta0_lima, col=1, row=1, label = "G")%>%
  fill_panel(MR3_all_beta1_lima, col=1, row=2, label = "I")%>%
  fill_panel(MR3_all_beta2_lima, col=1, row=3, label = "K")%>%
  fill_panel(MR3_all_period2_lima, col=2, row=1, label = "H")%>%
  fill_panel(MR3_all_x1_lima, col=2, row=2, label = "J")%>%
  fill_panel(MR3_all_x2_lima, col=2, row=3, label = "L")

fill_panel(comp_figure_S, MR3_all_beta0_south, col=1, row=1, label = "M")%>%
  fill_panel(MR3_all_beta1_south, col=1, row=2, label = "O")%>%
  fill_panel(MR3_all_beta2_south, col=1, row=3, label = "Q")%>%
  fill_panel(MR3_all_cull_int, col=1, row=4, label = "S")%>%
  fill_panel(MR3_all_period2_south, col=2, row=1, label = "N")%>%
  fill_panel(MR3_all_x1_south, col=2, row=2, label = "P")%>%
  fill_panel(MR3_all_x2_south, col=2, row=3, label = "R")%>%
  fill_panel(key, col=2, row=4, label = "T")

#-------------------------------------------------------------------------------
#Plot posterior distributions for all parameters from the final model inclusive of all data types

{
beta0_plot<-
  ggplot()+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta0_north, xmax=prod.mcmc_A_base_long_met$q97.5$beta0_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#99ddff")+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta0_south, xmax=prod.mcmc_A_base_long_met$q97.5$beta0_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#ee8866")+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta0_lima, xmax=prod.mcmc_A_base_long_met$q97.5$beta0_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#44BB99")+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta0_north), aes(x=prod.mcmc_A_base_long_met$sims.list$beta0_north), col="#99ddff", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta0_south), aes(x=prod.mcmc_A_base_long_met$sims.list$beta0_south), col="#ee8866", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta0_lima), aes(x=prod.mcmc_A_base_long_met$sims.list$beta0_lima), col="#44BB99", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$beta0_north), col="#99ddff", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$beta0_south), col="#ee8866", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$beta0_lima), col="#44BB99", linewidth=1)+
  labs(x="Beta0", y="Density")+xlim(0,15)+ylim(0,1)+
  theme_bw()
stat_function(fun = dbeta, args = list(shape=1.5, rate=400), col = "black", linetype="dashed")


beta1_plot<-
  ggplot()+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta1_north, xmax=prod.mcmc_A_base_long_met$q97.5$beta1_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#99ddff")+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta1_south, xmax=prod.mcmc_A_base_long_met$q97.5$beta1_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#ee8866")+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta1_lima, xmax=prod.mcmc_A_base_long_met$q97.5$beta1_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#44BB99")+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta1_north), aes(x=prod.mcmc_A_base_long_met$sims.list$beta1_north), col="#99ddff", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta1_south), aes(x=prod.mcmc_A_base_long_met$sims.list$beta1_south), col="#ee8866", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta1_lima), aes(x=prod.mcmc_A_base_long_met$sims.list$beta1_lima), col="#44BB99", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$beta1_north), col="#99ddff", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$beta1_south), col="#ee8866", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$beta1_lima), col="#44BB99", linewidth=1)+
  labs(x="Beta1", y="Density")+xlim(0,1)+ylim(0,10)+
  theme_bw()
stat_function(fun = dbeta, args = list(shape=1.5, rate=400), col = "black", linetype="dashed")


beta2_plot<-
  ggplot()+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta2_north, xmax=prod.mcmc_A_base_long_met$q97.5$beta2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#99ddff")+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta2_south, xmax=prod.mcmc_A_base_long_met$q97.5$beta2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#ee8866")+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$beta2_lima, xmax=prod.mcmc_A_base_long_met$q97.5$beta2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#44BB99")+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta2_north), aes(x=prod.mcmc_A_base_long_met$sims.list$beta2_north), col="#99ddff", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta2_south), aes(x=prod.mcmc_A_base_long_met$sims.list$beta2_south), col="#ee8866", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$beta2_lima), aes(x=prod.mcmc_A_base_long_met$sims.list$beta2_lima), col="#44BB99", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$beta2_north), col="#99ddff", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$beta2_south), col="#ee8866", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$beta2_lima), col="#44BB99", linewidth=1)+
  labs(x="Beta2", y="Density")+xlim(0,1)+ylim(0,10)+
  theme_bw()
stat_function(fun = dbeta, args = list(shape=1.5, rate=400), col = "black", linetype="dashed")


period2_plot<-
  ggplot()+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$period2_north, xmax=prod.mcmc_A_base_long_met$q97.5$period2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#99ddff")+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$period2_south, xmax=prod.mcmc_A_base_long_met$q97.5$period2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#ee8866")+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$period2_lima, xmax=prod.mcmc_A_base_long_met$q97.5$period2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#44BB99")+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$period2_north), aes(x=prod.mcmc_A_base_long_met$sims.list$period2_north), col="#99ddff", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$period2_south), aes(x=prod.mcmc_A_base_long_met$sims.list$period2_south), col="#ee8866", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$period2_lima), aes(x=prod.mcmc_A_base_long_met$sims.list$period2_lima), col="#44BB99", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$period2_north), col="#99ddff", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$period2_south), col="#ee8866", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$period2_lima), col="#44BB99", linewidth=1)+
  labs(x="Period 2 (Days)", y="Density")+xlim(1000,1899)+ylim(0,.02)+
  theme_bw()
stat_function(fun = dbeta, args = list(shape=1.5, rate=400), col = "black", linetype="dashed")


x1_plot<-
  ggplot()+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$x1_north, xmax=prod.mcmc_A_base_long_met$q97.5$x1_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#99ddff")+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$x1_south, xmax=prod.mcmc_A_base_long_met$q97.5$x1_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#ee8866")+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$x1_lima, xmax=prod.mcmc_A_base_long_met$q97.5$x1_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#44BB99")+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$x1_north), aes(x=prod.mcmc_A_base_long_met$sims.list$x1_north), col="#99ddff", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$x1_south), aes(x=prod.mcmc_A_base_long_met$sims.list$x1_south), col="#ee8866", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$x1_lima), aes(x=prod.mcmc_A_base_long_met$sims.list$x1_lima), col="#44BB99", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$x1_north), col="#99ddff", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$x1_south), col="#ee8866", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$x1_lima), col="#44BB99", linewidth=1)+
  labs(x="x1", y="Density")+xlim(-1,2)+ylim(0,10)+
  theme_bw()
stat_function(fun = dbeta, args = list(shape=1.5, rate=400), col = "black", linetype="dashed")



x2_plot<-
  ggplot()+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$x2_north, xmax=prod.mcmc_A_base_long_met$q97.5$x2_north, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#99ddff")+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$x2_south, xmax=prod.mcmc_A_base_long_met$q97.5$x2_south, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#ee8866")+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$x2_lima, xmax=prod.mcmc_A_base_long_met$q97.5$x2_lima, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#44BB99")+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$x2_north), aes(x=prod.mcmc_A_base_long_met$sims.list$x2_north), col="#99ddff", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$x2_south), aes(x=prod.mcmc_A_base_long_met$sims.list$x2_south), col="#ee8866", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$x2_lima), aes(x=prod.mcmc_A_base_long_met$sims.list$x2_lima), col="#44BB99", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$x2_north), col="#99ddff", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$x2_south), col="#ee8866", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$x2_lima), col="#44BB99", linewidth=1)+
  labs(x="x2", y="Density")+xlim(0,2)+ylim(0,10)+
  theme_bw()
stat_function(fun = dbeta, args = list(shape=1.5, rate=400), col = "black", linetype="dashed")
}

figure_s2_plot<-multi_panel_figure(width=c(85,85), height=c(60,60,60), unit="mm")
figure_s2_plot

fill_panel(figure_s2_plot, beta0_plot, row=1, col=1)%>%
  fill_panel(period2_plot, row=1, col=2)%>%
  fill_panel(beta1_plot, row=2, col=1)%>%
  fill_panel(beta2_plot, row=2, col=2)%>%
  fill_panel(x1_plot, row=3, col=1)%>%
  fill_panel(x2_plot, row=3, col=2)
#-------------------------------------------------------------------------------





#######################################################################################################################################
 #plot the results of model fitting with test data

#load data
setwd("C:/Users/mg321b/OneDrive - University of Glasgow/RA_projects/Git/Bat_flu_modelling")
model_long_test_output_both<-read.csv("model_output_data/test_data_JAGS_output.csv")

parameter_plot<-function(parameter_in, parameter_m, parameter_sd){
  
  output<-model_long_test_output_both[,c(parameter_in, parameter_m, parameter_sd, "model")]
  output$upp<-output[, parameter_m] + (output[, parameter_sd] * 2.5)
  output$lwr<-output[, parameter_m] - (output[, parameter_sd] * 2.5)
  
  nms <- names(output)
  x <- nms[1]
  y <- nms[2]
  p <- ggplot(output, aes(x = !!ensym(x), y = !!ensym(y), group=as.factor(model), col=as.factor(model)))+ 
    geom_point(alpha=0)+
    geom_smooth(method='lm',formula=y~x, se=FALSE)+
    geom_smooth(method="lm", aes(x=!!ensym(x), y=upp), formula=y~x, se=FALSE, alpha=0.1)+
    geom_smooth(method="lm", aes(x=!!ensym(x), y=lwr), formula=y~x, se=FALSE, alpha=0.1)+
    geom_abline()+
    theme_bw()
  
  lines<-ggplot_build(p)$data[[2]]
  upp_lines<-ggplot_build(p)$data[[3]]
  lwr_lines<-ggplot_build(p)$data[[4]]
  
  ribbon<-cbind(lines[,c(2,3,5)], upp_lines[,c(3)], lwr_lines[,c(3)])
  names(ribbon)[4]<-"yu"
  names(ribbon)[5]<-"yl"
  
  q<- ggplot(data=subset(ribbon))+
    geom_ribbon(aes(ymax=yu, ymin=yl, x=x, group=as.factor(group), fill=as.factor(group)), alpha=0.1)+
    geom_line(aes(x=x, y=y, group=as.factor(group), col=as.factor(group)), linewidth=1)+
    geom_abline(linetype=2)+
    theme_bw()+
    theme(legend.position = "none")
  labs(x="Input", y="Output")+
    # labs(fill="Data", col="Data")+
    scale_color_manual(labels = c("Population-level \ndata only", "All data"), values=c(2,4))+
    scale_fill_manual(labels = c("Population-level \ndata only", "All data"), values=c(2,4))
  
  return(q)
}



cull_int_plot<-parameter_plot("cull_int", "cull_int_m", "cull_int_sd")
cull_int_plot<-
  cull_int_plot+labs(x="Input csf", y="Output csf")

gamma_plot<-parameter_plot("gamma", "gamma_m", "gamma_sd")
gamma_plot<-
  gamma_plot+labs(x=expression(paste("Input ", gamma)), y=expression(paste("Output ", gamma)))

alpha_plot<-parameter_plot("alpha", "alpha_m", "alpha_sd")
alpha_plot<-
  alpha_plot+labs(x=expression(paste("Input ", alpha)), y=expression(paste("Output ", alpha)))

delta_plot<-parameter_plot("delta", "delta_m", "delta_sd")
delta_plot<-
  delta_plot+labs(x=expression(paste("Input ", delta)), y=expression(paste("Output ", delta)))

beta0_plot<-parameter_plot("beta0", "beta0_m", "beta0_sd")
beta0_plot<-
  beta0_plot+labs(x=expression(paste("Input ", beta, "0")), y=expression(paste("Output ", beta, "0")))

beta1_plot<-parameter_plot("beta1", "beta1_m", "beta1_sd")
beta1_plot<-
  beta1_plot+labs(x=expression(paste("Input ", beta, "1")), y=expression(paste("Output ", beta, "1")))

beta2_plot<-parameter_plot("beta2", "beta2_m", "beta2_sd")
beta2_plot<-
  beta2_plot+labs(x=expression(paste("Input ", beta, "2")), y=expression(paste("Output ", beta, "2")))

period2_plot<-parameter_plot("P2", "period2_n_m", "period2_n_sd")
period2_plot<-
  period2_plot+labs(x="Input period2", y="Output period2")



x1_data<-data.frame(model=c(1:192), mean=c(1:192))

x1_data$model<-model_long_test_output_both$model
x1_data$mean<-(model_long_test_output_both$x1_m)
x1_data$x<-"x1"

x2_data<-data.frame(model=c(1:192), mean=c(1:192))

x2_data$model<-model_long_test_output_both$model
x2_data$mean<-(model_long_test_output_both$x2_m)
x2_data$x<-"x2"

x_data<-rbind(x1_data, x2_data)


x_plot<-
  ggplot()+
  geom_hline(aes(yintercept=1), linetype="dashed")+
  geom_boxplot(data=x_data, aes(x=x, y=mean, fill=as.factor(model)))+
  theme_bw()+ylim(0,2)+
  labs(y="posterior mean", x="parameter")+
  scale_fill_discrete(name = "Data set", labels = c("\nPopulation-level \ndata only\n", "\nAll data: population \nand individual-level\n"))+
  theme(legend.position="none")



#plot a figure with the legend to get the legend to plot separately

legend_plot<-parameter_plot("P2", "period2_n_m", "period2_n_sd")
legend_plot<-
  legend_plot+labs(x="Input period2", y="Output period2")+
  scale_color_discrete(name = "Data set", labels = c("\nPopulation-level \nserology data only\n", "\nAll data"))+
  scale_fill_discrete(name = "Data set", labels = c("\nPopulation-level \nserology data only\n", "\nAll data"))+
  
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=11),
        legend.key.size = unit(1.3, "cm"),
        legend.position = "right")

legend_plot<-get_legend(legend_plot)
legend_plot<-as_ggplot(legend_plot)



test_data_plot<-multi_panel_figure(width=c(85,85), height=c(50,50,50,50,40), unit="mm")
test_data_plot

fill_panel(test_data_plot, alpha_plot, row=3, col=1, label = "C")%>%
  fill_panel(gamma_plot, row=2, col=1, label = "B")%>%
  fill_panel(delta_plot, row=4, col=1, label = "D")%>%
  fill_panel(cull_int_plot, row=1, col=1, label = "A")%>%
  fill_panel(beta0_plot, row=1, col=2, label = "F")%>%
  fill_panel(beta1_plot, row=2, col=2, label = "G")%>%
  fill_panel(beta2_plot, row=3, col=2, label = "H")%>%
  fill_panel(period2_plot, row=4, col=2, label = "I")%>%
  fill_panel(x_plot, row=5, col=1, label = "E")%>%
  fill_panel(legend_plot, row=5, col=2)







