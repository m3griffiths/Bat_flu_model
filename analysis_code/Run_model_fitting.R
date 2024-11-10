#This code collates the serology data and uses this data to run model fitting using JAGS.

#clear the environment
rm(list=ls())

#load packages
{
  library(boot)
  library(R2jags)
  library(coda)
  library(runjags)
  library(tidyverse)
  library(deSolve)
}



#load data
flu_data<-read.csv("raw_data/bat_final_unique_MG_day.csv")
bats_flu_od<-read.csv("raw_data/bat_final_unique_MG_day.csv")
bat_cull <- read.csv("raw_data/Bat_cull.csv", header=T)
Laura_metagenomic_samples <- read.csv("raw_data/Laura_metagenomic_samples.csv", header=T)


#convert the culling data into the format required by the model

##by month
culling_vector<-c(rep(0, times=144))

for(i in 1:nrow(bat_cull)){
  n<-bat_cull[i,8]
  culling_vector[n-1]<-bat_cull[i,6]
}

##by day
culling_vector_day<-rep(culling_vector, each=30)



# function for organising the population-level serology data
##counts the number of x in l
count <- function(x,l){
  sum(l[which(!(is.na(l)))] == x)
}



###################################
###################################
##        Format datasets        ##
###################################
###################################

################
##format serology data sets
{
  #Lima
  vect_prev_lima = c()
  vect_tot_lima = c()
  for (i in c(min(bats_flu_od$Year):max(bats_flu_od$Year))){
    for (j in c(1:12)){
      n0 = count(0,bats_flu_od[which(bats_flu_od$Year == i & bats_flu_od$Month == j & bats_flu_od$Site %in% c("LMA4", "LMA5", "LMA6", "LMA10")),]$Serol)
      n1 = count(1,bats_flu_od[which(bats_flu_od$Year == i & bats_flu_od$Month == j & bats_flu_od$Site %in% c("LMA4", "LMA5", "LMA6", "LMA10")),]$Serol)
      vect_tot_lima <- c(vect_tot_lima,n0+n1)
      if (n0+n1 == 0){
        vect_prev_lima <- c(vect_prev_lima,NA)}
      else{
        vect_prev_lima <- c(vect_prev_lima,n1)
      }
    }
  }
  
  sample_y_lima <- vect_prev_lima
  sample_n_lima <- vect_tot_lima
  
  t_max <- length(vect_tot_lima)
  
  
  #North
  vect_prev_north = c()
  vect_tot_north = c()
  for (i in c(min(bats_flu_od$Year):max(bats_flu_od$Year))){
    for (j in c(1:12)){
      n0 = count(0,bats_flu_od[which(bats_flu_od$Year == i & bats_flu_od$Month == j & bats_flu_od$MetaRegionB == "North"),]$Serol)
      n1 = count(1,bats_flu_od[which(bats_flu_od$Year == i & bats_flu_od$Month == j & bats_flu_od$MetaRegionB == "North"),]$Serol)
      vect_tot_north <- c(vect_tot_north,n0+n1)
      if (n0+n1 == 0){
        vect_prev_north <- c(vect_prev_north,NA)}
      else{
        vect_prev_north <- c(vect_prev_north,n1)
      }
    }
  }
  
  sample_y_north <- vect_prev_north
  sample_n_north <- vect_tot_north
  
 
  #South 
  vect_prev = c()
  vect_tot = c()
  for (i in c(min(bats_flu_od$Year):max(bats_flu_od$Year))){
    for (j in c(1:12)){
      n0 = count(0,bats_flu_od[which(bats_flu_od$Year == i & bats_flu_od$Month == j & bats_flu_od$MetaRegionB %in% c("South_b", "South_a")),]$Serol)
      n1 = count(1,bats_flu_od[which(bats_flu_od$Year == i & bats_flu_od$Month == j & bats_flu_od$MetaRegionB %in% c("South_b", "South_a")),]$Serol)
      vect_tot <- c(vect_tot,n0+n1)
      if (n0+n1 == 0){
        vect_prev <- c(vect_prev,NA)}
      else{
        vect_prev <- c(vect_prev,n1)
      }
    }
  }
  
  sample_y_south <- vect_prev
  sample_n_south <- vect_tot
}


################
##format metagenomic data sets

{
#North 
vect_prev = c()
vect_tot = c()
for (i in c(min(bats_flu_od$Year):max(bats_flu_od$Year))){
  for (j in c(1:12)){
    n0 = count(0,Laura_metagenomic_samples[which(Laura_metagenomic_samples$Year == i & Laura_metagenomic_samples$Month == j & Laura_metagenomic_samples$MR %in% c("North")),]$Test)
    n1 = count(1,Laura_metagenomic_samples[which(Laura_metagenomic_samples$Year == i & Laura_metagenomic_samples$Month == j & Laura_metagenomic_samples$MR %in% c("North")),]$Test)
    vect_tot <- c(vect_tot,n0+n1)
    if (n0+n1 == 0){
      vect_prev <- c(vect_prev,NA)}
    else{
      vect_prev <- c(vect_prev,n1)
    }
  }
}

sample_i_y_north <- vect_prev
sample_i_n_north <- vect_tot


#South 
vect_prev = c()
vect_tot = c()
for (i in c(min(bats_flu_od$Year):max(bats_flu_od$Year))){
  for (j in c(1:12)){
    n0 = count(0,Laura_metagenomic_samples[which(Laura_metagenomic_samples$Year == i & Laura_metagenomic_samples$Month == j & Laura_metagenomic_samples$MR %in% c("Lima")),]$Test)
    n1 = count(1,Laura_metagenomic_samples[which(Laura_metagenomic_samples$Year == i & Laura_metagenomic_samples$Month == j & Laura_metagenomic_samples$MR %in% c("Lima")),]$Test)
    vect_tot <- c(vect_tot,n0+n1)
    if (n0+n1 == 0){
      vect_prev <- c(vect_prev,NA)}
    else{
      vect_prev <- c(vect_prev,n1)
    }
  }
}

sample_i_y_lima <- vect_prev
sample_i_n_lima <- vect_tot



#South 
vect_prev = c()
vect_tot = c()
for (i in c(min(bats_flu_od$Year):max(bats_flu_od$Year))){
  for (j in c(1:12)){
    n0 = count(0,Laura_metagenomic_samples[which(Laura_metagenomic_samples$Year == i & Laura_metagenomic_samples$Month == j & Laura_metagenomic_samples$MR %in% c("South")),]$Test)
    n1 = count(1,Laura_metagenomic_samples[which(Laura_metagenomic_samples$Year == i & Laura_metagenomic_samples$Month == j & Laura_metagenomic_samples$MR %in% c("South")),]$Test)
    vect_tot <- c(vect_tot,n0+n1)
    if (n0+n1 == 0){
      vect_prev <- c(vect_prev,NA)}
    else{
      vect_prev <- c(vect_prev,n1)
    }
  }
}

sample_i_y_south <- vect_prev
sample_i_n_south <- vect_tot

sum(sample_i_n_north)
sum(sample_i_n_south)
sum(sample_i_n_lima)

}

################
##get multiply sampled bats for each region

##extract bat which were sampled more than once (whose Id appears multiple times) from the sample data frame
bats_flu_multi <- flu_data[which(flu_data$Id %in% flu_data[duplicated(flu_data$Id),]$Id),]

##make column for the number of the month which aligns to the other data
bats_flu_multi$Month_number <- (bats_flu_multi$Year-2007)*12+bats_flu_multi$Month


{
  bats_flu_multi_North<-filter(bats_flu_multi, MetaRegionB=="North")
  bats_flu_multi_South<-filter(bats_flu_multi, MetaRegionB %in% c("South_a", "South_b"))
  bats_flu_multi_Lima<-filter(bats_flu_multi, MetaRegionB=="Lima")
  
  samples_ID_North<-c(names(which(table(bats_flu_multi_North$Id) > 1)))
  samples_ID_South<-c(names(which(table(bats_flu_multi_South$Id) > 1)))
  samples_ID_Lima<-c(names(which(table(bats_flu_multi_Lima$Id) > 1)))
  
  repeat_samples_df_North<-setNames(data.frame(matrix(ncol = length(samples_ID_North)+1, nrow = 144)), c("Time", samples_ID_North))
  repeat_samples_df_South<-setNames(data.frame(matrix(ncol = length(samples_ID_South)+1, nrow = 144)), c("Time", samples_ID_South))
  repeat_samples_df_Lima<-setNames(data.frame(matrix(ncol = length(samples_ID_Lima)+1, nrow = 144)), c("Time", samples_ID_Lima))
  repeat_samples_df_North$Time<-c(1:144)
  repeat_samples_df_South$Time<-c(1:144)
  repeat_samples_df_Lima$Time<-c(1:144)
  
  n<-2
  for(x in samples_ID_North){
    bat_single<-filter(bats_flu_multi_North, Id==x)
    for(row in 1:nrow(bat_single)){
      repeat_samples_df_North[bat_single$Month_number[row],n]<-bat_single$Serol[row]
      
    }
    n<-n+1
  }
  
  
  n<-2
  for(x in samples_ID_South){
    bat_single<-filter(bats_flu_multi_South, Id==x)
    for(row in 1:nrow(bat_single)){
      repeat_samples_df_South[bat_single$Month_number[row],n]<-bat_single$Serol[row]
      
    }
    n<-n+1
  }
  
  n<-2
  for(x in samples_ID_Lima){
    bat_single<-filter(bats_flu_multi_Lima, Id==x)
    for(row in 1:nrow(bat_single)){
      repeat_samples_df_Lima[bat_single$Month_number[row],n]<-bat_single$Serol[row]
      
    }
    n<-n+1
  }
  
  
}

#get longitudinal samples in the correct format
{
  
  long_bats_subset1<-repeat_samples_df_Lima
  long_bats_subset2<-repeat_samples_df_North
  long_bats_subset3<-repeat_samples_df_South
  
  
  #Lima
  Ts_lima<-data.frame(T1=NA, T2=NA)
  B_sim_lima<-data.frame(S1=NA, S2=NA)
  n<-1
  
  for(i in 2:(ncol(long_bats_subset1))){
    Ts_interm<-filter(long_bats_subset1, is.na(long_bats_subset1[,i])==FALSE)
    for(x in 1:(nrow(Ts_interm)-1)){
      Ts_lima[n,]<-Ts_interm$Time[x:(x+1)] 
      B_sim_lima[n,1]<-Ts_interm[x,i]
      B_sim_lima[n,2]<-Ts_interm[(x+1),i]
      n<-n+1
    }
    
  }
  
  #North
  Ts_north<-data.frame(T1=NA, T2=NA)
  B_sim_north<-data.frame(S1=NA, S2=NA)
  n<-1
  
  for(i in 2:(ncol(long_bats_subset2))){
    Ts_interm<-filter(long_bats_subset2, is.na(long_bats_subset2[,i])==FALSE)
    for(x in 1:(nrow(Ts_interm)-1)){
      Ts_north[n,]<-Ts_interm$Time[x:(x+1)] 
      B_sim_north[n,1]<-Ts_interm[x,i]
      B_sim_north[n,2]<-Ts_interm[(x+1),i]
      n<-n+1
    }
    
  }
  
  #South
  Ts_south<-data.frame(T1=NA, T2=NA)
  B_sim_south<-data.frame(S1=NA, S2=NA)
  n<-1
  
  for(i in 2:(ncol(long_bats_subset3))){
    Ts_interm<-filter(long_bats_subset3, is.na(long_bats_subset3[,i])==FALSE)
    for(x in 1:(nrow(Ts_interm)-1)){
      Ts_south[n,]<-Ts_interm$Time[x:(x+1)] 
      B_sim_south[n,1]<-Ts_interm[x,i]
      B_sim_south[n,2]<-Ts_interm[(x+1),i]
      n<-n+1
    }
    
  }
  
  B_sim_lima<-as.matrix(B_sim_lima)
  Ts_lima<-as.matrix(Ts_lima)
  B_sim_north<-as.matrix(B_sim_north)
  Ts_north<-as.matrix(Ts_north)
  B_sim_south<-as.matrix(B_sim_south)
  Ts_south<-as.matrix(Ts_south)
  
}
##


###################################
###################################
##          Set up JAGS          ##
###################################
###################################


##list of parameters we want JAGS to save the output for


#Model 1

pars_A1<-c("beta1_lima","beta2_lima","period2_lima","beta0_lima",
          "beta1_north","beta2_north","period2_north","beta0_north",
          "beta1_south",
          "beta_north", 
          "beta_south", 
          "beta_lima",
          "x1_north", "x2_north", 
          "x1_south", "x2_south", 
          "x1_lima", "x2_lima",
          "beta2_south","period2_south","beta0_south",
          "gamma", "alpha", "delta",
          "cull_int",
          "s_lima",
          "i_lima","r_lima","e_lima",
          "s_north",
          "i_north","r_north","e_north",
          "s_south","i_south","r_south","e_south",
          "R0_lima","I0_lima","E0_lima",
          "R0_north","I0_north","E0_north",
          "R0_south","I0_south","E0_south",
          "s_p_south","i_p_south","r_p_south","e_p_south",
          "s_p2_south","i_p2_south","r_p2_south","e_p2_south",
          "theta_north", "theta_south", "theta_lima",
          "p_met_north", "p_met_south", "p_met_lima"
)


#Model 2

pars_A2<-c("beta1_lima","beta2_lima","period2_lima","beta0_lima",
          "beta1_north","beta2_north","period2_north","beta0_north",
          "beta1_south",
          "beta_north", 
          "beta_south", 
          "beta_lima",
          "x1_north", "x2_north", 
          "x1_south", "x2_south", 
          "x1_lima", "x2_lima",
          "beta2_south","period2_south","beta0_south",
          "gamma", "alpha", "delta",
          "cull_int",
          "s_lima","i_lima","r1_lima","r2_lima","e_lima",
          "s_north","i_north","r1_north","r2_north","e_north",
          "s_south","i_south","r1_south","r2_south","e_south",
          "R10_lima","R20_lima","I0_lima","E0_lima",
          "R10_north","R20_north","I0_north","E0_north",
          "R10_south","R20_south","I0_south","E0_south",
          "theta_north", "theta_south", "theta_lima",
          "p_met_north", "p_met_south", "p_met_lima"
)


##list of inputs to JAGS

mydata_model<-list(t_max=144,
                   #population level serology data
                   sample_y_lima = sample_y_lima, sample_n_lima = sample_n_lima, 
                   sample_y_north = sample_y_north, sample_n_north = sample_n_north, 
                   sample_y_south = sample_y_south, sample_n_south = sample_n_south, 
                   #add metagenomic samples here
                   met_north_n = sum(sample_i_n_north), met_north_y = 0,
                   met_south_n = sum(sample_i_n_south), met_south_y = 0,
                   met_lima_n = sum(sample_i_n_lima), met_lima_y = 0,
                   #add longitudinal data here
                   nbats_lima = nrow(Ts_lima), B_lima = B_sim_lima, Ts_lima = Ts_lima,
                   nbats_north = nrow(Ts_north), B_north = B_sim_north, Ts_north = Ts_north,
                   nbats_south = nrow(Ts_south), B_south = B_sim_south, Ts_south = Ts_south,
                   pi = pi, 
                   cull_time = culling_vector_day,
                   # add any fixed parameters here
                   r_bd = (3.5*365),
                   period1_north = 365,
                   period1_lima = 365
)


##select appropriate model file

##Model 1 - different levels of data inclusion

#modelfile <- "analysis_code/Base_model_final.txt"
#modelfile <- "analysis_code/Base_long_model_final.txt"
#modelfile <- "analysis_code/Base_met_model_final.txt"
modelfile <- "analysis_code/Base_long_met_model_final.txt"


##Model 2 - full dataset only

#modelfile <- "analysis_code/Base_long_met_model_final_R1R2.txt"


###################################
###################################
##     Run JAGS model fitting    ##
###################################
###################################


myinitials <- NULL
time1 <- Sys.time()


#For testing the impact of data inclusion, the model fitting was run with n.iter=5000.
#For the final model including all of the data, n.iter=100,000.

prod.mcmc_A <- do.call(jagsUI::jags,
                       list(data=mydata_model, inits=myinitials, parameters.to.save=pars_A1, 
                            model.file = modelfile,
                            n.chains=4, n.iter=100000, n.burnin=100000/2,
                            n.thin=max(c(((100000/2)/1000),1)),  DIC = TRUE, parallel = TRUE,
                            set.seed(123), bugs.format=FALSE))

time2<-Sys.time()
working_time<-time2-time1
print(working_time)


#save output of the model run
save.image("model_output_data/model_output_n_internations_yyyymmdd.RData")

