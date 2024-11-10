#JAGS - testing with and without longitudinal data - model 4, fixed

rm(list=ls())


{
  library(boot)
  library(R2jags)
  library(coda)
  library(runjags)
  library(tidyverse)
  library(deSolve)
  library(jagsUI)
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

met_samples_n<-as.numeric(nrow(filter(Laura_metagenomic_samples, Test==0)))
met_samples_y<-as.numeric(nrow(filter(Laura_metagenomic_samples, Test==1)))

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



###############################################################################
#simulated ODE to run with test parameters in order to simulate data


###################################
###################################
##     Deterministic model       ##
###################################
###################################

#function to control when culling takes place
intro_cull <- function(time, cull_time_list){
  if(is.na(cull_time_list[round((time/30)+1)])){
    return(0)
  }else{
    return(cull_time_list[round((time/30)+1)])
  }
}



si <- function(time, state, parameters) {
  with(as.list(c(state, parameters,time)), {
    
    N <- S+E+I+R
    
    b <- b
    
    m <- b
    
    
    gamma <- gamma
    delta <- delta
    alpha <- alpha
    
    beta_0 <- beta0
    
    
    beta_1 <- beta1
    beta_2 <- beta2
    
    P1 <- P1
    P2 <- P2
    
    x1 <- x1
    x2 <- x2
    
    cull_int <- cull_int
    cull_time_list <- c(0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.28,
                        1,0.95,0.81,0,0,0,0,0.37,0.58,0.52,0.51,0.42,
                        0.51,0.39,0.27,0.22,0,0.28,0.1,0,0,0,0.25,0.24,
                        0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0
    )
    
    
    cull_time <- intro_cull(time, cull_time_list)
    
    beta<-(1+((beta_0/4)*(1+beta_1*(sin((pi*x1)+(2*pi*time/P1))))*(1+beta_2*(sin((pi*x2)+(2*pi*time/P2))))))
    
    dS <- -(1/b)*S + (1/m)*N+((1/m)*N*(1-(N/K))) + (1/gamma)*R - (1/beta)*I*S   -cull_int*cull_time*S
    dE <- -(1/b)*E - (1/delta)*E +         (1/beta)*I*S -cull_int*cull_time*E
    dI <- -(1/b)*I + (1/delta)*E - (1/alpha)*I          -cull_int*cull_time*I
    dR <- -(1/b)*R + (1/alpha)*I - (1/gamma)*R     -cull_int*cull_time*R
    dI_c <- + (1/delta)*E 
    
    
    
    
    return(list(c( dS, dE, dI, dR, dI_c)))
  })
}



###################################
###################################
##           Test data           ##
###################################
###################################

#-------------------------------------------------------------------------------

#make dataframe for output of the model fitting
model_long_test_output<-data.frame(i=numeric(0), model=numeric(0),
                                   alpha_m=numeric(0), alpha_sd=numeric(0),
                                   delta_m=numeric(0), delta_sd=numeric(0),
                                   gamma_m=numeric(0), gamma_sd=numeric(0),
                                   beta0_m=numeric(0), beta0_sd=numeric(0),
                                   beta1_m=numeric(0), beta1_sd=numeric(0),
                                   beta2_m=numeric(0), beta2_sd=numeric(0),
                                   x1_m=numeric(0), x1_sd=numeric(0),
                                   x2_m=numeric(0), x2_sd=numeric(0),
                                   cull_int_m=numeric(0), cull_int_sd=numeric(0),
                                   period2_n_m=numeric(0), period2_n_sd=numeric(0)
)





#get random combinations of parameters for the test data

library(pomp)

sims<-100


{
  set.seed(41)
  runif_design( lower=c(beta0=5,  delta=2, alpha=4,  gamma=30*3, beta1=0.15, beta2=0.15, P2=1300, x1=1, x2=1, cull_int=0), 
                upper=c(beta0=12, delta=4, alpha=10, gamma=800,  beta1=0.85, beta2=0.85, P2=1700, x1=1, x2=1, cull_int=0.005), 
                nseq=sims) -> M4_params
}

M4_params<-M4_params[-c(47,68,94,97), ]  # remove parameter sets in which the bat population is driven to 0 by culling.


#x1 and x2 are not varied here, as the cyclical nature of this parameter means that prior bound adjustment would be needed for each run.

#for each parameter set
for (i in 1:sims){
  
  pset<-i
  
  
  #parameters for north and lima i.e., cull_int==0
  parameters<-c(beta0=M4_params[i,1], delta=M4_params[i,2], alpha=M4_params[i,3], gamma=M4_params[i,4], 
                beta1=M4_params[i,5], beta2=M4_params[i,6], P2=M4_params[i,7], x1=M4_params[i,8], x2=M4_params[i,9],
                cull_int=0, b=1277.50, P1=365, K=1)
  
  #parameters for south i.e., with cull_int taken from M4_params
  parameters_c<-c(beta0=M4_params[i,1], delta=M4_params[i,2], alpha=M4_params[i,3], gamma=M4_params[i,4], 
                  beta1=M4_params[i,5], beta2=M4_params[i,6], P2=M4_params[i,7], x1=M4_params[i,8], x2=M4_params[i,9],
                  cull_int=M4_params[i,10], b=1277.50, P1=365, K=1)
  
  
  #choose starting values for each compartment
  init <- c(S = 0.49, E=0.005, I=.01, R=0.495, I_c=0)
  #choose length of simulation in days
  times <- seq(0, (144*30), by = 1) 
  
  
  #run the deterministic model
  #south
  out_c <- as.data.frame(ode(y = init, times = times, func = si, parms = parameters_c))
  #north and lima
  out <- as.data.frame(ode(y = init, times = times, func = si, parms = parameters))
  
  
  
  I_new_month_c<-c()
  I_new_month<-c()
  
  for (t in 1:144){
    I_new_month_c[t]<-(out_c$I_c[30*t]-out_c$I_c[30*(t-1)+1])/(out_c$S[30*t]+out_c$E[30*t]+out_c$I[30*t]+out_c$R[30*t])
    I_new_month[t]<-(out$I_c[30*t]-out$I_c[30*(t-1)+1])/(out$S[30*t]+out$E[30*t]+out$I[30*t]+out$R[30*t])
    
  }
  
  
  #Get the test data from this model run
  #Month:
  R_test_m_c<-c()
  R_test_m<-c()
  
  
  for(t in 1:144){
    R_test_m_c[t]<- mean(out_c$R[(30*(t-1)+1):(30*t)])/(mean(out_c$S[(30*(t-1)+1):(30*t)])+mean(out_c$E[(30*(t-1)+1):(30*t)])+mean(out_c$I[(30*(t-1)+1):(30*t)])+mean(out_c$R[(30*(t-1)+1):(30*t)]))
    R_test_m[t]<- mean(out$R[(30*(t-1)+1):(30*t)])/(mean(out$S[(30*(t-1)+1):(30*t)])+mean(out$E[(30*(t-1)+1):(30*t)])+mean(out$I[(30*(t-1)+1):(30*t)])+mean(out$R[(30*(t-1)+1):(30*t)]))
    
  }
  
  
  
  
  sample_y_south <- round(R_test_m_c*sample_n_south,0)
  sample_n_south <- sample_n_south
  
  sample_y_north <- round(R_test_m*sample_n_north,0)
  sample_n_north <- sample_n_north
  
  sample_y_lima <- round(R_test_m*sample_n_lima,0)
  sample_n_lima <- sample_n_lima
  
  #-------------------------------------------------------------------------------
  
  #Generate longitudinal samples
  
  {
    nmonths <- 144 #12*12
    nbats <- 1000
    
    B_c <- matrix(NA, nbats, nmonths)
    B <- matrix(NA, nbats, nmonths)
    
    gamma_month<-parameters_c["gamma"]/30
    pe <- 1/gamma_month # probability of "extinction"
    
    # initial state
    
    B_c[ , 1] <- rbinom(nbats, size = 1, prob = out_c[1,5])
    B[ , 1] <- rbinom(nbats, size = 1, prob = out[1,5])
    
    for(t in 2:nmonths){
      
      
      pc_c <- I_new_month_c[t]
      
      pc <- I_new_month[t]
      
      # infected can become healthy with probability pe
      
      
      B_c[ , t] <- rbinom(nbats, size = B_c[ , t-1], prob = 1 - pe) + 
        rbinom(nbats, size = 1 - B_c[, t-1], prob = pc_c)
      
      B[ , t] <- rbinom(nbats, size = B[ , t-1], prob = 1 - pe) + 
        rbinom(nbats, size = 1 - B[, t-1], prob = pc)
    }
    
  }
  
  #################################################################################
  
  
  #################################################################################
  #From the generated data above (both population level and individual), generate
  #test data i.e., use actual data number of sample times/points
  
  #get the correct number of longitudinally sampled pairs as per the data
  
  {
    nbats_south<-59
    nbats_north<-62
    nbats_lima<-148
    
    
    {
      set.seed(55)
      
      long_bats_subset_south<-repeat_samples_df_South
      long_bats_subset_north<-repeat_samples_df_North
      long_bats_subset_lima<-repeat_samples_df_Lima
      
    }
    
    
    Ts_south<-data.frame(T1=NA, T2=NA)
    Ts_north<-data.frame(T1=NA, T2=NA)
    Ts_lima<-data.frame(T1=NA, T2=NA)
    
    
    
    for(i in 2:(ncol(long_bats_subset_south))){
      Ts_interm<-filter(long_bats_subset_south, is.na(long_bats_subset_south[,i])==FALSE)
      Ts_south[i-1,]<-Ts_interm$Time[1:2]
    }
    
    for(i in 2:(ncol(long_bats_subset_north))){
      Ts_interm<-filter(long_bats_subset_north, is.na(long_bats_subset_north[,i])==FALSE)
      Ts_north[i-1,]<-Ts_interm$Time[1:2]
    }
    
    for(i in 2:(ncol(long_bats_subset_lima))){
      Ts_interm<-filter(long_bats_subset_lima, is.na(long_bats_subset_lima[,i])==FALSE)
      Ts_lima[i-1,]<-Ts_interm$Time[1:2]
    }
    
    #We have the times at which the samples were taken, now we need to get the data at these points from the 
    #simulated longitudinal bats in B_df
    
    {
      set.seed(55)
      
      long_bats_sim_subset_c<-B_c[c(1,sample(1:nrow(B_c),nbats_south)),]
      
      long_bats_sim_subset<-B[c(1,sample(1:nrow(B),nbats_lima)),]
      
    }
    
    
    
    B_sim_south<-matrix(nrow=nbats_south, ncol=2)
    B_sim_north<-matrix(nrow=nbats_north, ncol=2)
    B_sim_lima<-matrix(nrow=nbats_lima, ncol=2)
    
    
    for(r in 1:nrow(Ts_south)){
      B_sim_south[r,1]<-long_bats_sim_subset_c[r,Ts_south[r,1]]
      B_sim_south[r,2]<-long_bats_sim_subset_c[r,Ts_south[r,2]]
    }
    
    for(r in 1:nrow(Ts_north)){
      B_sim_north[r,1]<-long_bats_sim_subset[r,Ts_north[r,1]]
      B_sim_north[r,2]<-long_bats_sim_subset[r,Ts_north[r,2]]
    }
    
    for(r in 1:nrow(Ts_lima)){
      B_sim_lima[r,1]<-long_bats_sim_subset[r,Ts_lima[r,1]]
      B_sim_lima[r,2]<-long_bats_sim_subset[r,Ts_lima[r,2]]
    }
    
    
    Ts_south<-as.matrix(Ts_south)
    Ts_north<-as.matrix(Ts_north)
    Ts_lima<-as.matrix(Ts_lima)
    
    
  }
  #################################################################################
  
  
  #################################################################################
  
  #parameters to save the output of - just recording the joint and South parameters here for computational time/disc space
  
  pars_A<-c("beta1_south",
            "x1_south", "x2_south", 
            "beta2_south","period2_south","beta0_south",
            "gamma", "alpha", "delta",
            "cull_int"
  )
  
  mydata_model<-list(t_max=144,
                     sample_y_south = sample_y_south, sample_n_south = sample_n_south, 
                     sample_y_north = sample_y_north, sample_n_north = sample_n_north, 
                     sample_y_lima = sample_y_lima, sample_n_lima = sample_n_lima, 
                     #add metagenomic samples here
                     met_all_n = 172, met_all_y = 0,
                     #add longitudinal data here
                     nbats_south = nrow(Ts_south), B_south = B_sim_south, Ts_south = Ts_south,
                     nbats_north = nrow(Ts_north), B_north = B_sim_north, Ts_north = Ts_north,
                     nbats_lima = nrow(Ts_lima), B_lima = B_sim_lima, Ts_lima = Ts_lima,
                     pi = pi, 
                     cull_time = culling_vector_day,
                     # add any fixed parameters here
                     r_bd = (3.5*365)
                     , period1_north = 365
                     , period1_lima = 365

  )
  
  #setwd to JAGS script location 
  setwd("C:/Users/mg321b/OneDrive - University of Glasgow/RA_projects/Git/Bat_flu_modelling")
  
  models<-c("analysis_code/Base_model_final.txt","analysis_code/Base_long_model_final.txt")
  
  for(m in 1:2){
    
    modelfile<-models[m]
    
    #set up
    myinitials <- NULL
    time1 <- Sys.time()
    
    #Here, the model was run using only 2 chains for 1000 iterations to keep processing time down.
    prod.mcmc_A <- do.call(jagsUI::jags,
                           list(data=mydata_model, inits=myinitials, parameters.to.save=pars_A, 
                                model.file = modelfile,
                                n.chains=2, n.iter=1000, n.burnin=1000/2,
                                n.thin=max(c(((1000/2)/1000),1)),  DIC = TRUE, parallel = TRUE,
                                set.seed(123), bugs.format=FALSE))
    
    time2<-Sys.time()
    working_time<-time2-time1
    print(pset)
    print(working_time)
    
    
    
    #update the output dataframe based on the results of model fitting to this parameter set.
    model_long_test_output[nrow(model_long_test_output)+1,]<-c(pset, m,
                                                               (prod.mcmc_A$mean$alpha),   (prod.mcmc_A$sd$alpha),
                                                               (prod.mcmc_A$mean$delta),   (prod.mcmc_A$sd$delta),
                                                               (prod.mcmc_A$mean$gamma),   (prod.mcmc_A$sd$gamma),
                                                               (prod.mcmc_A$mean$beta0_south),   (prod.mcmc_A$sd$beta0_south),
                                                               (prod.mcmc_A$mean$beta1_south),   (prod.mcmc_A$sd$beta1_south),
                                                               (prod.mcmc_A$mean$beta2_south),   (prod.mcmc_A$sd$beta2_south),
                                                               (prod.mcmc_A$mean$x1_south),      (prod.mcmc_A$sd$x1_south),
                                                               (prod.mcmc_A$mean$x2_south),      (prod.mcmc_A$sd$x2_south),
                                                               (prod.mcmc_A$mean$cull_int),(prod.mcmc_A$sd$cull_int),
                                                               (prod.mcmc_A$mean$period2_south), (prod.mcmc_A$sd$period2_south)
    )
    
  }
}



model_long_test_output_base<-filter(model_long_test_output_all, model==1)
model_long_test_output_full<-filter(model_long_test_output_all, model==2)

model_long_test_output_base<-cbind(M4_params, model_long_test_output_base)
model_long_test_output_full<-cbind(M4_params, model_long_test_output_full)

model_long_test_output_both<-rbind(model_long_test_output_base, model_long_test_output_full)


write.csv(model_long_test_output_both, "model_output_data/test_data_JAGS_output.csv")

#################################################################################

