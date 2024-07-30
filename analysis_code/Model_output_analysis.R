#############Model analysis

rm(list=ls())

#load packages
library(tidyverse)
library(dplyr)
library(loo)
library(Hmisc)


#load data - requires download from repository at https://zenodo.org/records/12820416
load("All_model_outputs_20240531.RData")

############################################################################################################

################## Reformat input population-level data for plotting purposes

{
  #reformat the input data for plotting purposes
  
  sample_y_lima[is.na(sample_y_lima)] <- 0
  sample_y_north[is.na(sample_y_north)] <- 0
  sample_y_south[is.na(sample_y_south)] <- 0
  
  data_points_north<-as.data.frame(binconf(x=sample_y_north, n=sample_n_north))
  data_points_north<-cbind(data_points_north, time=(seq(1:144)/12)+2007)
  data_points_north<-cbind(data_points_north, array=data_points_north$Upper-data_points_north$PointEst)
  data_points_north<-cbind(data_points_north, arrayminus=data_points_north$PointEst-data_points_north$Lower)
  data_points_north<-cbind(data_points_north, sample_y_north)%>% drop_na()
  
  data_points_south<-as.data.frame(binconf(x=sample_y_south, n=sample_n_south))
  data_points_south<-cbind(data_points_south, time=(seq(1:144)/12)+2007)
  data_points_south<-cbind(data_points_south, array=data_points_south$Upper-data_points_south$PointEst)
  data_points_south<-cbind(data_points_south, arrayminus=data_points_south$PointEst-data_points_south$Lower)
  data_points_south<-cbind(data_points_south, sample_y_south)%>% drop_na()
  
  data_points_lima<-as.data.frame(binconf(x=sample_y_lima, n=sample_n_lima))
  data_points_lima<-cbind(data_points_lima, time=(seq(1:144)/12)+2007)
  data_points_lima<-cbind(data_points_lima, array=data_points_lima$Upper-data_points_lima$PointEst)
  data_points_lima<-cbind(data_points_lima, arrayminus=data_points_lima$PointEst-data_points_lima$Lower)
  data_points_lima<-cbind(data_points_lima, sample_y_lima)%>% drop_na()
  
}

############################################################################################################

################## Leave-one-out analysis for models 1 and 2

{

# prod.mcmc_A_M1_LOO contains all of the data required for LOO of model 1
# prod.mcmc_A_M2_LOO contains all of the data required for LOO of model 2

##################################################################### model 1

p_south<-as.data.frame(prod.mcmc_A_M1_LOO$sims.list$p_south)
p_north<-as.data.frame(prod.mcmc_A_M1_LOO$sims.list$p_north)
p_lima<-as.data.frame(prod.mcmc_A_M1_LOO$sims.list$p_lima)

p_south<-rbind(p_south, sample_n_south, sample_y_south)
p_north<-rbind(p_north, sample_n_north, sample_y_north)
p_lima<-rbind(p_lima, sample_n_lima, sample_y_lima)

p_all<-cbind(p_north, p_south, p_lima)

p_all[is.na(p_all)] <- 0

loglik_all<-matrix(nrow=4000, ncol=432)

for(j in 1:432){
  
  X<-p_all[4002,j]
  n<-p_all[4001,j]
  p<-p_all[(1:4000),j]
  
  loglike <- dbinom(X, size = n, prob = p, log = TRUE)  # log-likelihood of p  
  
  loglik_all[1:4000,j]<-loglike
  
}



loglik_all<-loglik_all[,loglik_all[1,]!=0]


#get the metagenomic data

#average the I for each region for the entire time period for each iteration
infecteds_means<-matrix(0, nrow=4002, ncol=3)

for(i in 1:4000){
  
  infecteds_means[i,1]<-mean(prod.mcmc_A_M1_LOO$sims.list$i_north[i,1:4320])
  infecteds_means[i,2]<-mean(prod.mcmc_A_M1_LOO$sims.list$i_south[i,1:4320])
  infecteds_means[i,3]<-mean(prod.mcmc_A_M1_LOO$sims.list$i_lima[i,1:4320])
  
}

infecteds_means[4001,1]<-50
infecteds_means[4001,2]<-107
infecteds_means[4001,3]<-15


loglik_all<-cbind(loglik_all, NA, NA, NA)

for(j in 1:3){
  
  X<-infecteds_means[4002,j]
  n<-infecteds_means[4001,j]
  p<-infecteds_means[(1:4000),j]
  
  loglike <- dbinom(X, size = n, prob = p, log = TRUE)  # log-likelihood of p  
  
  loglik_all[1:4000,(j+54)]<-loglike
  
}

#get the longitudinal likelihood data

loglik_long_north<-matrix(0,nrow=4000, ncol=nrow(B_sim_north))
loglik_long_south<-matrix(0,nrow=4000, ncol=nrow(B_sim_south))
loglik_long_lima<-matrix(0,nrow=4000, ncol=nrow(B_sim_lima))


for (i in 1:nrow(B_sim_north)){
  for (j in 1:4000){
    
    p<-prod.mcmc_A_M1_LOO$sims.list$theta_north[j,i,Ts_north[i,2]]
    n<-1
    X<-B_sim_north[i,2]
    
    loglike <- dbinom(X, size = n, prob = p, log = TRUE)  # log-likelihood of p 
    loglik_long_north[j,i]<-loglike
    
  }
}

for (i in 1:nrow(B_sim_south)){
  for (j in 1:4000){
    
    p<-prod.mcmc_A_M1_LOO$sims.list$theta_south[j,i,Ts_south[i,2]]
    n<-1
    X<-B_sim_south[i,2]
    
    loglike <- dbinom(X, size = n, prob = p, log = TRUE)  # log-likelihood of p 
    loglik_long_south[j,i]<-loglike
    
  }
}

for (i in 1:nrow(B_sim_lima)){
  for (j in 1:4000){
    
    p<-prod.mcmc_A_M1_LOO$sims.list$theta_lima[j,i,Ts_lima[i,2]]
    n<-1
    X<-B_sim_lima[i,2]
    
    loglike <- dbinom(X, size = n, prob = p, log = TRUE)  # log-likelihood of p 
    loglik_long_lima[j,i]<-loglike
    
  }
}

loglik_all<-cbind(loglik_all, loglik_long_north, loglik_long_south, loglik_long_lima)


##############################

loo_n<-loo(loglik_all)
print(loo_n)

plot(loo_n)

##################################################################### model 2


p2_north<-data.frame(prod.mcmc_A_M2_LOO$sims.list$p_north)
p2_south<-data.frame(prod.mcmc_A_M2_LOO$sims.list$p_south)
p2_lima<-data.frame(prod.mcmc_A_M2_LOO$sims.list$p_lima)

p2_north[4001,1:144]<-p_north[4001,1:144]
p2_north[4002,1:144]<-p_north[4002,1:144]

p2_south[4001,1:144]<-p_south[4001,1:144]
p2_south[4002,1:144]<-p_south[4002,1:144]

p2_lima[4001,1:144]<-p_lima[4001,1:144]
p2_lima[4002,1:144]<-p_lima[4002,1:144]

p2_all<-cbind(p2_north, p2_south, p2_lima)

p2_all[is.na(p2_all)] <- 0

loglik2_all<-matrix(nrow=4000, ncol=432)

for(j in 1:432){
  
  X<-p2_all[4002,j]
  n<-p2_all[4001,j]
  p<-p2_all[(1:4000),j]
  
  loglike <- dbinom(X, size = n, prob = p, log = TRUE)  # log-likelihood of p  
  
  loglik2_all[1:4000,j]<-loglike
  
}


loglik2_all<-loglik2_all[,loglik2_all[1,]!=0]

#get the metagenomic data

#average the I for each region for the entire time period for each iteration
infecteds_means<-matrix(0, nrow=4002, ncol=3)

for(i in 1:4000){
  
  infecteds_means[i,1]<-mean(prod.mcmc_A_M2_LOO$sims.list$i_north[i,1:4320])
  infecteds_means[i,2]<-mean(prod.mcmc_A_M2_LOO$sims.list$i_south[i,1:4320])
  infecteds_means[i,3]<-mean(prod.mcmc_A_M2_LOO$sims.list$i_lima[i,1:4320])
  
}

infecteds_means[4001,1]<-50
infecteds_means[4001,2]<-107
infecteds_means[4001,3]<-15


loglik2_all<-cbind(loglik2_all, NA, NA, NA)

for(j in 1:3){
  
  X<-infecteds_means[4002,j]
  n<-infecteds_means[4001,j]
  p<-infecteds_means[(1:4000),j]
  
  loglike <- dbinom(X, size = n, prob = p, log = TRUE)  # log-likelihood of p  
  
  loglik2_all[1:4000,(j+54)]<-loglike
  
}


#get the longitudinal likelihood data

loglik2_long_north<-matrix(0,nrow=4000, ncol=nrow(B_sim_north))
loglik2_long_south<-matrix(0,nrow=4000, ncol=nrow(B_sim_south))
loglik2_long_lima<-matrix(0,nrow=4000, ncol=nrow(B_sim_lima))


for (i in 1:nrow(B_sim_north)){
  for (j in 1:4000){
    
    p<-prod.mcmc_A_M2_LOO$sims.list$theta_north[j,i,Ts_north[i,2]]
    n<-1
    X<-B_sim_north[i,2]
    
    loglike <- dbinom(X, size = n, prob = p, log = TRUE)  # log-likelihood of p 
    loglik2_long_north[j,i]<-loglike
    
  }
}

for (i in 1:nrow(B_sim_south)){
  for (j in 1:4000){
    
    p<-prod.mcmc_A_M2_LOO$sims.list$theta_south[j,i,Ts_south[i,2]]
    n<-1
    X<-B_sim_south[i,2]
    
    loglike <- dbinom(X, size = n, prob = p, log = TRUE)  # log-likelihood of p 
    loglik2_long_south[j,i]<-loglike
    
  }
}

for (i in 1:nrow(B_sim_lima)){
  for (j in 1:4000){
    
    p<-prod.mcmc_A_M2_LOO$sims.list$theta_lima[j,i,Ts_lima[i,2]]
    n<-1
    X<-B_sim_lima[i,2]
    
    loglike <- dbinom(X, size = n, prob = p, log = TRUE)  # log-likelihood of p 
    loglik2_long_lima[j,i]<-loglike
    
  }
}

loglik2_all<-cbind(loglik2_all, loglik2_long_north, loglik2_long_south, loglik2_long_lima)


loo_2n<-loo(loglik2_all)
print(loo_2n)

plot(loo_2n)

loo_diff<-loo_compare(loo_n, loo_2n)

print(loo_diff)

}

#Comparing the models on PSIS-LOO reveals an estimated difference in elpd of 73.0 
#(with a standard error of 23.8) in favor of Model 1.

############################################################################################################

################# Model 1 - change in total population size due to culling
################# and formatting compartment predictions for plotting

{

  #The data set you just loaded already contains the data frames generated below 
  #because it takes a while, but feel free to run them yourself as well.
  
  
  #the model fitting has been set up to provide us with the number of individuals in each compartment 
  #for each time point, with each parameter set tested.
  #However, since the population size in the South meta-region is not constant, this code converts
  #these raw numbers into the proportion of the total population size.
  
  #Make data frames
  N_all<-data.frame(time=c(1:4320))
  r_p_all<-data.frame(time=c(1:4320))
  i_p_all<-data.frame(time=c(1:4320))
  e_p_all<-data.frame(time=c(1:4320))
  
  loop_length<-(length(prod.mcmc_A_base_long_met$sims.list$beta0_south))
  
  #for each parameter set
  for(x in 1:loop_length){
    
    #sum s+e+i+r south at each time point to get the total population size  
    N_run<-c(prod.mcmc_A_base_long_met$sims.list$s_south[x,]+ prod.mcmc_A_base_long_met$sims.list$e_south[x,]+ prod.mcmc_A_base_long_met$sims.list$i_south[x,]+ prod.mcmc_A_base_long_met$sims.list$r_south[x,] )
    
    #divide the relevant compartment by the total population size at each time point
    r_p_run<- c(prod.mcmc_A_base_long_met$sims.list$r_south[x,]/(prod.mcmc_A_base_long_met$sims.list$s_south[x,]+ prod.mcmc_A_base_long_met$sims.list$e_south[x,]+ prod.mcmc_A_base_long_met$sims.list$i_south[x,]+ prod.mcmc_A_base_long_met$sims.list$r_south[x,]))
    i_p_run<- c(prod.mcmc_A_base_long_met$sims.list$i_south[x,]/(prod.mcmc_A_base_long_met$sims.list$s_south[x,]+ prod.mcmc_A_base_long_met$sims.list$e_south[x,]+ prod.mcmc_A_base_long_met$sims.list$i_south[x,]+ prod.mcmc_A_base_long_met$sims.list$r_south[x,]))
    e_p_run<- c(prod.mcmc_A_base_long_met$sims.list$e_south[x,]/(prod.mcmc_A_base_long_met$sims.list$s_south[x,]+ prod.mcmc_A_base_long_met$sims.list$e_south[x,]+ prod.mcmc_A_base_long_met$sims.list$i_south[x,]+ prod.mcmc_A_base_long_met$sims.list$r_south[x,]))
    
    #add to data frame
    N_all<-cbind(N_all, N_run)
    r_p_all<-cbind(r_p_all, r_p_run)
    i_p_all<-cbind(i_p_all, i_p_run)
    e_p_all<-cbind(e_p_all, e_p_run)
    
  }
  
  #For each time point in the data frames created above, get the mean and 95% CI.
  for (row in 1:nrow(N_all)){
    N_all$mean[row]<-as.numeric((sum(N_all[row,(2:(loop_length+1))])/loop_length))
    N_all$q2.5[row]<-as.numeric(quantile(N_all[row,(2:(loop_length+1))], probs=0.025))
    N_all$q97.5[row]<-as.numeric(quantile(N_all[row,(2:(loop_length+1))], probs=0.975))
    
    r_p_all$mean[row]<-as.numeric((sum(r_p_all[row,(2:(loop_length+1))])/loop_length))
    r_p_all$q2.5[row]<-as.numeric(quantile(r_p_all[row,(2:(loop_length+1))], probs=0.025))
    r_p_all$q97.5[row]<-as.numeric(quantile(r_p_all[row,(2:(loop_length+1))], probs=0.975))
    
    i_p_all$mean[row]<-as.numeric((sum(i_p_all[row,(2:(loop_length+1))])/loop_length))
    i_p_all$q2.5[row]<-as.numeric(quantile(i_p_all[row,(2:(loop_length+1))], probs=0.025))
    i_p_all$q97.5[row]<-as.numeric(quantile(i_p_all[row,(2:(loop_length+1))], probs=0.975))
    
    e_p_all$mean[row]<-as.numeric((sum(e_p_all[row,(2:(loop_length+1))])/loop_length))
    e_p_all$q2.5[row]<-as.numeric(quantile(e_p_all[row,(2:(loop_length+1))], probs=0.025))
    e_p_all$q97.5[row]<-as.numeric(quantile(e_p_all[row,(2:(loop_length+1))], probs=0.975))
  }
  
  
  
  {
    N_plot<-as.data.frame(N_all[,((loop_length+2):(loop_length+4))])
    N_plot$time<-c(1:4320)
    
    r_p_plot<-as.data.frame(r_p_all[,((loop_length+2):(loop_length+4))])
    r_p_plot$time<-c(1:4320)
    
    i_p_plot<-as.data.frame(i_p_all[,((loop_length+2):(loop_length+4))])
    i_p_plot$time<-c(1:4320)
    
    e_p_plot<-as.data.frame(e_p_all[,((loop_length+2):(loop_length+4))])
    e_p_plot$time<-c(1:4320)
  }
  
  #Get the mean estimate + 95% CI for the lowest population size to occur due to culling
  
  (1-min(N_plot$mean))*100
  (1-min(N_plot$q2.5))*100
  (1-min(N_plot$q97.5))*100
  
  
  
  
  #get the mean values and upper and lower bounds for all time varying outputs
  
  {
    time <- c(1:((144)*30))
    time2 <- 30*seq(1,144) 
    
    s_north.mean<- apply(prod.mcmc_A_base_long_met$sims.list$s_north,2,mean)
    s_north.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$s_north,2,"quantile",probs=0.025)
    s_north.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$s_north,2,"quantile",probs=0.975)
    e_north.mean<- apply(prod.mcmc_A_base_long_met$sims.list$e_north,2,mean)
    e_north.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$e_north,2,"quantile",probs=0.025)
    e_north.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$e_north,2,"quantile",probs=0.975)
    i_north.mean<- apply(prod.mcmc_A_base_long_met$sims.list$i_north,2,mean)
    i_north.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$i_north,2,"quantile",probs=0.025)
    i_north.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$i_north,2,"quantile",probs=0.975)
    r_north.mean<- apply(prod.mcmc_A_base_long_met$sims.list$r_north,2,mean)
    r_north.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$r_north,2,"quantile",probs=0.025)
    r_north.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$r_north,2,"quantile",probs=0.975)
    
    s_lima.mean<- apply(prod.mcmc_A_base_long_met$sims.list$s_lima,2,mean)
    s_lima.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$s_lima,2,"quantile",probs=0.025)
    s_lima.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$s_lima,2,"quantile",probs=0.975)
    e_lima.mean<- apply(prod.mcmc_A_base_long_met$sims.list$e_lima,2,mean)
    e_lima.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$e_lima,2,"quantile",probs=0.025)
    e_lima.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$e_lima,2,"quantile",probs=0.975)
    i_lima.mean<- apply(prod.mcmc_A_base_long_met$sims.list$i_lima,2,mean)
    i_lima.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$i_lima,2,"quantile",probs=0.025)
    i_lima.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$i_lima,2,"quantile",probs=0.975)
    r_lima.mean<- apply(prod.mcmc_A_base_long_met$sims.list$r_lima,2,mean)
    r_lima.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$r_lima,2,"quantile",probs=0.025)
    r_lima.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$r_lima,2,"quantile",probs=0.975)
    
    s_south.mean<- apply(prod.mcmc_A_base_long_met$sims.list$s_south,2,mean)
    s_south.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$s_south,2,"quantile",probs=0.025)
    s_south.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$s_south,2,"quantile",probs=0.975)
    e_south.mean<- apply(prod.mcmc_A_base_long_met$sims.list$e_south,2,mean)
    e_south.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$e_south,2,"quantile",probs=0.025)
    e_south.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$e_south,2,"quantile",probs=0.975)
    i_south.mean<- apply(prod.mcmc_A_base_long_met$sims.list$i_south,2,mean)
    i_south.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$i_south,2,"quantile",probs=0.025)
    i_south.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$i_south,2,"quantile",probs=0.975)
    r_south.mean<- apply(prod.mcmc_A_base_long_met$sims.list$r_south,2,mean)
    r_south.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$r_south,2,"quantile",probs=0.025)
    r_south.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$r_south,2,"quantile",probs=0.975)
    
    s_p2_south.mean<- apply(prod.mcmc_A_base_long_met$sims.list$s_p2_south,2,mean)
    s_p2_south.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$s_p2_south,2,"quantile",probs=0.025)
    s_p2_south.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$s_p2_south,2,"quantile",probs=0.975)
    e_p2_south.mean<- apply(prod.mcmc_A_base_long_met$sims.list$e_p2_south,2,mean)
    e_p2_south.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$e_p2_south,2,"quantile",probs=0.025)
    e_p2_south.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$e_p2_south,2,"quantile",probs=0.975)
    i_p2_south.mean<- apply(prod.mcmc_A_base_long_met$sims.list$i_p2_south,2,mean)
    i_p2_south.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$i_p2_south,2,"quantile",probs=0.025)
    i_p2_south.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$i_p2_south,2,"quantile",probs=0.975)
    r_p2_south.mean<- apply(prod.mcmc_A_base_long_met$sims.list$r_p2_south,2,mean)
    r_p2_south.lwr<-  apply(prod.mcmc_A_base_long_met$sims.list$r_p2_south,2,"quantile",probs=0.025)
    r_p2_south.upp<-  apply(prod.mcmc_A_base_long_met$sims.list$r_p2_south,2,"quantile",probs=0.975)
    
    
    
    beta_north.mean<- c(NA,apply(prod.mcmc_A_base_long_met$sims.list$beta_north[(1:4000),(2:4320)],2,mean))
    beta_south.mean<- c(NA,apply(prod.mcmc_A_base_long_met$sims.list$beta_south[(1:4000),(2:4320)],2,mean))
    beta_lima.mean<- c(NA,apply(prod.mcmc_A_base_long_met$sims.list$beta_lima[(1:4000),(2:4320)],2,mean))
    
    beta_north.upp<- c(NA,apply(prod.mcmc_A_base_long_met$sims.list$beta_north[(1:4000),(2:4320)],2,"quantile", probs=0.975))
    beta_south.upp<- c(NA,apply(prod.mcmc_A_base_long_met$sims.list$beta_south[(1:4000),(2:4320)],2,"quantile", probs=0.975))
    beta_lima.upp<- c(NA,apply(prod.mcmc_A_base_long_met$sims.list$beta_lima[(1:4000),(2:4320)],2,"quantile", probs=0.975))
    
    beta_north.lwr<- c(NA,apply(prod.mcmc_A_base_long_met$sims.list$beta_north[(1:4000),(2:4320)],2,"quantile", probs=0.025))
    beta_south.lwr<- c(NA,apply(prod.mcmc_A_base_long_met$sims.list$beta_south[(1:4000),(2:4320)],2,"quantile", probs=0.025))
    beta_lima.lwr<- c(NA,apply(prod.mcmc_A_base_long_met$sims.list$beta_lima[(1:4000),(2:4320)],2,"quantile", probs=0.025))
    
    
    
    time_series_data_lima<-data.frame(time=(time/360)+2007, 
                                      s.mean=s_lima.mean, s.lwr=s_lima.lwr, s.upp=s_lima.upp,
                                      e.mean=e_lima.mean, e.lwr=e_lima.lwr, e.upp=e_lima.upp,
                                      i.mean=i_lima.mean, i.lwr=i_lima.lwr, i.upp=i_lima.upp,
                                      r.mean=r_lima.mean, r.lwr=r_lima.lwr, r.upp=r_lima.upp
                                      ,beta.mean=beta_lima.mean, beta.lwr=beta_lima.lwr, beta.upp=beta_lima.upp
    )
    
    time_series_data_north<-data.frame(time=(time/360)+2007, 
                                       s.mean=s_north.mean, s.lwr=s_north.lwr, s.upp=s_north.upp,
                                       e.mean=e_north.mean, e.lwr=e_north.lwr, e.upp=e_north.upp,
                                       i.mean=i_north.mean, i.lwr=i_north.lwr, i.upp=i_north.upp,
                                       r.mean=r_north.mean, r.lwr=r_north.lwr, r.upp=r_north.upp
                                       ,beta.mean=beta_north.mean, beta.lwr=beta_north.lwr, beta.upp=beta_north.upp
    )
    
    time_series_data_south<-data.frame(time=(time/360)+2007, 
                                       s.mean=s_south.mean, s.lwr=s_south.lwr, s.upp=s_south.upp,
                                       e.mean=e_south.mean, e.lwr=e_south.lwr, e.upp=e_south.upp,
                                       i.mean=i_south.mean, i.lwr=i_south.lwr, i.upp=i_south.upp,
                                       r.mean=r_south.mean, r.lwr=r_south.lwr, r.upp=r_south.upp
                                       , beta.mean=beta_south.mean, beta.lwr=beta_south.lwr, beta.upp=beta_south.upp
    )
    
    time_series_data_south_p2<-data.frame(time=(time/360)+2007, 
                                          s.mean=s_p2_south.mean, s.lwr=s_p2_south.lwr, s.upp=s_p2_south.upp,
                                          e.mean=e_p2_south.mean, e.lwr=e_p2_south.lwr, e.upp=e_p2_south.upp,
                                          i.mean=i_p2_south.mean, i.lwr=i_p2_south.lwr, i.upp=i_p2_south.upp,
                                          r.mean=r_p2_south.mean, r.lwr=r_p2_south.lwr, r.upp=r_p2_south.upp
                                          
                                          
    )
  }   
}

############################################################################################################

################# Model 1 - model fit analysis

{
  
  p_south<-as.data.frame(prod.mcmc_A_base_long_met$sims.list$p_south)
  p_north<-as.data.frame(prod.mcmc_A_base_long_met$sims.list$p_north)
  p_lima<-as.data.frame(prod.mcmc_A_base_long_met$sims.list$p_lima)
  
  p_south<-rbind(p_south, sample_n_south, sample_y_south)
  p_north<-rbind(p_north, sample_n_north, sample_y_north)
  p_lima<-rbind(p_lima, sample_n_lima, sample_y_lima)
  
  
  P_south_y_dist<-as.data.frame(c(1:((nrow(p_south)-2)*1000)))
  for (x in 1:144){
    cum_dens<-c()
    if(p_south[(nrow(p_south)-1),x]==0){
      P_south_y_dist<-P_south_y_dist
    } else{ 
      for (y in 1:(nrow(p_south)-2)){
        cum_dens<-append(cum_dens, rbinom(n=100, size=p_south[(nrow(p_south)-1), x], prob=p_south[y,x]))
      }
      P_south_y_dist<-cbind(P_south_y_dist, cum_dens)
    }
  }
  
  P_north_y_dist<-as.data.frame(c(1:((nrow(p_north)-2)*1000)))
  for (x in 1:144){
    cum_dens<-c()
    if(p_north[(nrow(p_north)-1),x]==0){
      P_north_y_dist<-P_north_y_dist
    } else{ 
      for (y in 1:(nrow(p_north)-2)){
        cum_dens<-append(cum_dens, rbinom(n=1000, size=p_north[(nrow(p_north)-1), x], prob=p_north[y,x]))
      }
      P_north_y_dist<-cbind(P_north_y_dist, cum_dens)
    }
  }
  
  P_lima_y_dist<-as.data.frame(c(1:((nrow(p_lima)-2)*1000)))
  for (x in 1:144){
    cum_dens<-c()
    if(p_lima[(nrow(p_lima)-1),x]==0){
      P_lima_y_dist<-P_lima_y_dist
    } else{ 
      for (y in 1:(nrow(p_lima)-2)){
        cum_dens<-append(cum_dens, rbinom(n=1000, size=p_lima[(nrow(p_lima)-1), x], prob=p_lima[y,x]))
      }
      P_lima_y_dist<-cbind(P_lima_y_dist, cum_dens)
    }
  }
  
  
  
  percentiles_south<-c()
  for (i in 1:nrow(data_points_south)){
    cdf<-ecdf(P_south_y_dist[,i+1])
    percentiles_south<-append(percentiles_south,(cdf(data_points_south[i,7])))
    
  }
  
  inside_q_south<-c()
  for (i in 1:nrow(data_points_south)){
    lower_south<-quantile(P_south_y_dist[,i+1], 0.025)
    higher_south<-quantile(P_south_y_dist[,i+1], 0.975)
    test<-data_points_south[i,7] >= lower_south & data_points_south[i,7] <= higher_south
    inside_q_south<-append(inside_q_south, test)
    
  }
  
  
  percentiles_north<-c()
  for (i in 1:nrow(data_points_north)){
    cdf<-ecdf(P_north_y_dist[,i+1])
    percentiles_north<-append(percentiles_north,(cdf(data_points_north[i,7])))
    
  }
  
  inside_q_north<-c()
  for (i in 1:nrow(data_points_north)){
    lower_north<-quantile(P_north_y_dist[,i+1], 0.025)
    higher_north<-quantile(P_north_y_dist[,i+1], 0.975)
    test<-data_points_north[i,7] >= lower_north & data_points_north[i,7] <= higher_north
    inside_q_north<-append(inside_q_north, test)
    
  }
  
  percentiles_lima<-c()
  for (i in 1:nrow(data_points_lima)){
    cdf<-ecdf(P_lima_y_dist[,i+1])
    percentiles_lima<-append(percentiles_lima,(cdf(data_points_lima[i,7])))
    
  }
  
  inside_q_lima<-c()
  for (i in 1:nrow(data_points_lima)){
    lower_lima<-quantile(P_lima_y_dist[,i+1], 0.025)
    higher_lima<-quantile(P_lima_y_dist[,i+1], 0.975)
    test<-data_points_lima[i,7] >= lower_lima & data_points_lima[i,7] <= higher_lima
    inside_q_lima<-append(inside_q_lima, test)
    
  }
  
  
  percentiles_all<-append(percentiles_south, c(percentiles_north, percentiles_lima))
  
  
  rmse<-function(percentiles){
    
    sqrt((sum((0.5-percentiles)^2)/length(percentiles)))
    
  }
  
  rmse(percentiles_south)
  rmse(percentiles_north)
  rmse(percentiles_lima)
  
  
  rmse(percentiles_all)
  
  
  data_points_south<-cbind(data_points_south, percentiles_south)
  data_points_north<-cbind(data_points_north, percentiles_north)
  data_points_lima<-cbind(data_points_lima, percentiles_lima)
  
  data_points_south<-cbind(data_points_south, inside_q_south)
  data_points_north<-cbind(data_points_north, inside_q_north)
  data_points_lima<-cbind(data_points_lima, inside_q_lima)
  
  
}

################# Model 2 - model fit analysis

{
  #at present these will overwite the dataframes from the model 1 fit analysis
  
  p_south<-as.data.frame(prod.mcmc_A_M2$sims.list$p_south)
  p_north<-as.data.frame(prod.mcmc_A_M2$sims.list$p_north)
  p_lima<-as.data.frame(prod.mcmc_A_M2$sims.list$p_lima)
  
  p_south<-rbind(p_south, sample_n_south, sample_y_south)
  p_north<-rbind(p_north, sample_n_north, sample_y_north)
  p_lima<-rbind(p_lima, sample_n_lima, sample_y_lima)
  
  
  P_south_y_dist<-as.data.frame(c(1:((nrow(p_south)-2)*1000)))
  for (x in 1:144){
    cum_dens<-c()
    if(p_south[(nrow(p_south)-1),x]==0){
      P_south_y_dist<-P_south_y_dist
    } else{ 
      for (y in 1:(nrow(p_south)-2)){
        cum_dens<-append(cum_dens, rbinom(n=100, size=p_south[(nrow(p_south)-1), x], prob=p_south[y,x]))
      }
      P_south_y_dist<-cbind(P_south_y_dist, cum_dens)
    }
  }
  
  P_north_y_dist<-as.data.frame(c(1:((nrow(p_north)-2)*1000)))
  for (x in 1:144){
    cum_dens<-c()
    if(p_north[(nrow(p_north)-1),x]==0){
      P_north_y_dist<-P_north_y_dist
    } else{ 
      for (y in 1:(nrow(p_north)-2)){
        cum_dens<-append(cum_dens, rbinom(n=1000, size=p_north[(nrow(p_north)-1), x], prob=p_north[y,x]))
      }
      P_north_y_dist<-cbind(P_north_y_dist, cum_dens)
    }
  }
  
  P_lima_y_dist<-as.data.frame(c(1:((nrow(p_lima)-2)*1000)))
  for (x in 1:144){
    cum_dens<-c()
    if(p_lima[(nrow(p_lima)-1),x]==0){
      P_lima_y_dist<-P_lima_y_dist
    } else{ 
      for (y in 1:(nrow(p_lima)-2)){
        cum_dens<-append(cum_dens, rbinom(n=1000, size=p_lima[(nrow(p_lima)-1), x], prob=p_lima[y,x]))
      }
      P_lima_y_dist<-cbind(P_lima_y_dist, cum_dens)
    }
  }
  
  
  
  percentiles_south<-c()
  for (i in 1:nrow(data_points_south)){
    cdf<-ecdf(P_south_y_dist[,i+1])
    percentiles_south<-append(percentiles_south,(cdf(data_points_south[i,7])))
    
  }
  
  inside_q_south<-c()
  for (i in 1:nrow(data_points_south)){
    lower_south<-quantile(P_south_y_dist[,i+1], 0.025)
    higher_south<-quantile(P_south_y_dist[,i+1], 0.975)
    test<-data_points_south[i,7] >= lower_south & data_points_south[i,7] <= higher_south
    inside_q_south<-append(inside_q_south, test)
    
  }
  
  
  percentiles_north<-c()
  for (i in 1:nrow(data_points_north)){
    cdf<-ecdf(P_north_y_dist[,i+1])
    percentiles_north<-append(percentiles_north,(cdf(data_points_north[i,7])))
    
  }
  
  inside_q_north<-c()
  for (i in 1:nrow(data_points_north)){
    lower_north<-quantile(P_north_y_dist[,i+1], 0.025)
    higher_north<-quantile(P_north_y_dist[,i+1], 0.975)
    test<-data_points_north[i,7] >= lower_north & data_points_north[i,7] <= higher_north
    inside_q_north<-append(inside_q_north, test)
    
  }
  
  percentiles_lima<-c()
  for (i in 1:nrow(data_points_lima)){
    cdf<-ecdf(P_lima_y_dist[,i+1])
    percentiles_lima<-append(percentiles_lima,(cdf(data_points_lima[i,7])))
    
  }
  
  inside_q_lima<-c()
  for (i in 1:nrow(data_points_lima)){
    lower_lima<-quantile(P_lima_y_dist[,i+1], 0.025)
    higher_lima<-quantile(P_lima_y_dist[,i+1], 0.975)
    test<-data_points_lima[i,7] >= lower_lima & data_points_lima[i,7] <= higher_lima
    inside_q_lima<-append(inside_q_lima, test)
    
  }
  
  
  percentiles_all<-append(percentiles_south, c(percentiles_north, percentiles_lima))
  
  
  rmse<-function(percentiles){
    
    sqrt((sum((0.5-percentiles)^2)/length(percentiles)))
    
  }
  
  rmse(percentiles_south)
  rmse(percentiles_north)
  rmse(percentiles_lima)
  
  
  rmse(percentiles_all)
  
  
  data_points_south<-cbind(data_points_south, percentiles_south)
  data_points_north<-cbind(data_points_north, percentiles_north)
  data_points_lima<-cbind(data_points_lima, percentiles_lima)
  
  data_points_south<-cbind(data_points_south, inside_q_south)
  data_points_north<-cbind(data_points_north, inside_q_north)
  data_points_lima<-cbind(data_points_lima, inside_q_lima)
  
  
}

############################################################################################################

