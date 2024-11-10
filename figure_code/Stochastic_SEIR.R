#BIV stochastic model

#load packages
library(adaptivetau)
library(ggplot2)
library(dplyr)

#get test simulated data

#load data
flu_data<-read.csv("raw_data/bat_final_unique_MG_day.csv")
bats_flu_od<-read.csv("raw_data/bat_final_unique_MG_day.csv")
bat_cull <- read.csv("raw_data/Bat_cull.csv", header=T)

culling_vector<-c(rep(0, times=144))

for(i in 1:nrow(bat_cull)){
  n<-bat_cull[i,8]
  culling_vector[n-1]<-bat_cull[i,6]
}

culling_vector_day<-rep(culling_vector, each=30)

# counts the number of x in l
count <- function(x,l){
  sum(l[which(!(is.na(l)))] == x)
}



###################################
###################################
##        Format datasets        ##
###################################
###################################

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
  
  #South A 
  vect_prev = c()
  vect_tot = c()
  for (i in c(min(bats_flu_od$Year):max(bats_flu_od$Year))){
    for (j in c(1:12)){
      n0 = count(0,bats_flu_od[which(bats_flu_od$Year == i & bats_flu_od$Month == j & bats_flu_od$MetaRegionB == "South_a"),]$Serol)
      n1 = count(1,bats_flu_od[which(bats_flu_od$Year == i & bats_flu_od$Month == j & bats_flu_od$MetaRegionB == "South_a"),]$Serol)
      vect_tot <- c(vect_tot,n0+n1)
      if (n0+n1 == 0){
        vect_prev <- c(vect_prev,NA)}
      else{
        vect_prev <- c(vect_prev,n1)
      }
    }
  }
  
  sample_y_south_a <- vect_prev
  sample_n_south_a <- vect_tot
  
  
  
  #South B 
  vect_prev = c()
  vect_tot = c()
  for (i in c(min(bats_flu_od$Year):max(bats_flu_od$Year))){
    for (j in c(1:12)){
      n0 = count(0,bats_flu_od[which(bats_flu_od$Year == i & bats_flu_od$Month == j & bats_flu_od$MetaRegionB == "South_b"),]$Serol)
      n1 = count(1,bats_flu_od[which(bats_flu_od$Year == i & bats_flu_od$Month == j & bats_flu_od$MetaRegionB == "South_b"),]$Serol)
      vect_tot <- c(vect_tot,n0+n1)
      if (n0+n1 == 0){
        vect_prev <- c(vect_prev,NA)}
      else{
        vect_prev <- c(vect_prev,n1)
      }
    }
  }
  
  sample_y_south_b <- vect_prev
  sample_n_south_b <- vect_tot
  
  
  
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




library(Hmisc)

##########################
##      By Month        ##
##########################

Seropos_df<-data.frame(matrix(ncol=9,nrow=0, dimnames=list(NULL, c("Year", "Month", "MetaRegion", "total", "positive", "negative", "PointEst", "Lower", "Upper"))))

for (year in c(seq(2007, 2019, by=1))){
  for (month in c(1:12)){
    for (MR in c("Lima", "North", "South_a", "South_b")){
      df_sub<-subset(flu_data, Year == year & MetaRegionB == MR & Month == month)
      total<-nrow(df_sub)
      pos<-nrow(subset(df_sub, Serol == 1))
      bin_sero<-data.frame(binconf(pos, total))
      row<-data.frame(Year=year, Month=month, MetaRegion=MR, total=total, positive=pos, negative=total-pos)
      row<-cbind(row, bin_sero)
      Seropos_df<-rbind(Seropos_df, row)
    }
  }
}

Seropos_df$Date_continuous<-Seropos_df$Year+(Seropos_df$Month/12)









###################################################################################################

# function to give simulation output at specified time steps only
ssa.sptime<-function(sim,time.steps){
  i <- 2
  t(sapply(time.steps,function(t){
    while(sim[i,1]<t) i <<- i + 1
    c(time=t,sim[i-1,-1])
  }))
}


transitions_SEIR=list(c(S=+1), 
                      c(S=-1),
                      c(E=-1),
                      c(I=-1),
                      c(R1=-1),
                      c(R2=-1),
                      
                      c(S=-1, E=+1),
                      
                      c(S=-1, E=+1),
                      c(E=-1, I=+1),
                      c(I=-1, R1=+1),
                      c(R1=-1, R2=+1),
                      c(R2=-1, R1=+1),
                      c(R2=-1, E=+1)


)                                       # transitions between compartments, each corresponding to a rate specified below


SIRrateF_SEIR <- function(x, p, t){
  
  beta0 <- p["beta0"]
  beta1 <- p["beta1"]
  beta2 <- p["beta2"]
  P1 <- p["period1"]
  P2 <- p["period2"]
  x1 <- p["x1"]
  x2 <- p["x2"]
  
  phi <- p["phi"]
  theta <- p["theta"]
  
  beta <- (1/(1+((beta0/4)*(1+beta1*(sin((pi*x1)+((2*pi*t)/P1))))*(1+beta2*(sin((pi*x2)+((2*pi*t)/P2)))))))
  
  delta <- 1/p["D_inc"]                 # incubation rate = 1/incubation period
  alpha <- 1/p["D_inf"]                 # recovery rate = 1/infectious period
  gamma <- 1/p["D_imm"]                 # rate of immune waning = 1/immune period
  b <- 1/p["b"]                           # rate of natural births/deaths = 1/natural lifespan


  

  S <- x["S"]                           # susceptibles
  E <- x["E"]                          # Exposed to BIV
  I <- x["I"]                          # Infected with BIV and shedding
  R1 <- x["R1"]                          # Recovered and immune to BIV
  R2 <- x["R2"]                          # antibody negative but still immune (mostly, theta)


  N <- S + E + I + R1 + R2   # total number of bats 
  
  return(c(                             # the rates of each transition between compartments as specified above
    b * N,                              # births into the susceptible compartment
    b * S,                              # natural death
    b * E,                              # natural death
    b * I,                              # natural death
    b * R1,                              # natural death
    b * R2,

    phi * S,
    
    if(I>0){
    (beta-(phi/I)) * (S * I)
    }else{
      beta*S*I
    },
    delta * E,
    alpha * I,
    gamma * R1,
    beta * I * R2 * theta,
    beta * I *R2 * (1-theta)

    
    
    
  ))
}

total_pop<-10000

init.values=c(S=total_pop*0.295, E=total_pop*0.005, I=total_pop*0.005, R1=total_pop*0.495, R2=total_pop*0.20
              #, v=0
              ) # initial population values
sum(init.values)==total_pop

params<-c(b=365*3.5,
          D_inf=6.16,
          D_inc=2.03,
          D_imm=234, 
          beta0=5.08*total_pop,
          beta1=0.54,
          beta2=0.72,
          period1=365,
          period2=1423,
          x1=0.36,
          x2=0.74,
          phi=0.00000,
          theta=0
          # input parameters
          
)   


lima_params= c(b=365*3.5,
          D_inf=6.16,
          D_inc=2.03,
          D_imm=234, 
          beta0=5.91,
          beta1=0.31,
          beta2=0.72,
          period1=365,
          period2=1457,
          x1=-0.24,
          x2=0.73,
          phi=0.00000,
          theta=0
                             # input parameters

)                             # all time period values are in days.

north_params= c(b=365*3.5,
               D_inf=6.16,
               D_inc=2.03,
               D_imm=234, 
               beta0=5.08,
               beta1=0.54,
               beta2=0.72,
               period1=365,
               period2=1423,
               x1=0.36,
               x2=0.74,
               phi=0.00000,
               theta=0
               # input parameters
               
)                             # all time period values are in days.

south_params= c(b=365*3.5,
               D_inf=6.16,
               D_inc=2.03,
               D_imm=234, 
               beta0=8.01,
               beta1=0.31,
               beta2=0.70,
               period1=365,
               period2=1581,
               x1=0.81,
               x2=0.84,
               phi=0.00000,
               theta=0
               # input parameters
               
)                             # all time period values are in days.


x1<-runif(1, 1, 10000) # randomly generate number for the stochastic simulation
set.seed(x1)

time_run<-12*365 # length of time in days to run the simulation

# run the simulation

r<-data.frame(ssa.sptime(ssa.adaptivetau(init.values = init.values,transitions = transitions_SEIR,
                                         SIRrateF_SEIR, params = params,tf=time_run,
                                         tl.params =  list(maxtau = 1,extraChecks = F),
                                         set.seed(x1)),
                         time.steps = seq(0,time_run,30))); # set the time steps that will be output

(mean(r$I)/total_pop)*100

tail(r$I,1)

ggplot()+
  geom_line(data=r,aes(x=(time/365)+2007, y=E/(S+E+I+R1+R2)), col="red")+                                                   # plot number of active DrBHV infections
  geom_line(data=r,aes(x=(time/365)+2007, y=I/(S+E+I+R1+R2)), col="blue")+                                                  # plot number of latent DrBHV infections
  geom_line(data=r,aes(x=(time/365)+2007, y=R1/(S+E+I+R1+R2)), col="green")+               
  #geom_point(aes(x=Date_continuous, y=PointEst))+
  #geom_errorbar(data=filter(Seropos_df, MetaRegion=="North"), aes(x=Date_continuous, ymin=Lower, ymax=Upper), width=1/12)+
  ylim(0,1)+theme_bw()


################################
################################
#                              #
#                              #
#   Stochastic model loop      #
#                              #
#                              #
################################
################################


#explore the proportions of extinctions in the stochastic model at different population sizes

stoch_sim<-function(MR, sims, low, high, stepby){
  
  
  population_size_test<-data.frame( 
    PopSize=rep(seq(from=low, to=high, by=stepby)),
    Count=NA)
  n<-1  
  
  for (pop in seq(from=low, to=high, by=stepby)){
    
    count<-c()
    total_pop<-pop
    
    init.values=c(S=round(total_pop*0.295), E=round(total_pop*0.005), I=round(total_pop*0.005), R1=round(total_pop*0.495), R2=round(total_pop*0.20)
    ) 
    
    for (x1 in c(1:sims)){
      
      if (MR == "North"){
        params<-c(b=365*3.5,
                  D_inf=sample(prod.mcmc_A_base_long_met$sims.list$alpha, size=1)+0.25,
                  D_inc=sample(prod.mcmc_A_base_long_met$sims.list$delta, size=1),
                  D_imm=sample(prod.mcmc_A_base_long_met$sims.list$gamma, size=1), 
                  beta0=(sample(prod.mcmc_A_base_long_met$sims.list$beta0_north, size=1)-0.5)*total_pop,
                  beta1=sample(prod.mcmc_A_base_long_met$sims.list$beta1_north, size=1),
                  beta2=sample(prod.mcmc_A_base_long_met$sims.list$beta2_north, size=1),
                  period1=365,
                  period2=sample(prod.mcmc_A_base_long_met$sims.list$period2_north, size=1),
                  x1=sample(prod.mcmc_A_base_long_met$sims.list$x1_north, size=1),
                  x2=sample(prod.mcmc_A_base_long_met$sims.list$x2_north, size=1),
                  phi=0,
                  theta=0
                  # input parameters
                  
        )   
      }else if (MR == "South"){
        params<-c(b=365*3.5,
                  D_inf=sample(prod.mcmc_A_base_long_met$sims.list$alpha, size=1)+0.25,
                  D_inc=sample(prod.mcmc_A_base_long_met$sims.list$delta, size=1),
                  D_imm=sample(prod.mcmc_A_base_long_met$sims.list$gamma, size=1), 
                  beta0=(sample(prod.mcmc_A_base_long_met$sims.list$beta0_south, size=1)-0.5)*total_pop,
                  beta1=sample(prod.mcmc_A_base_long_met$sims.list$beta1_south, size=1),
                  beta2=sample(prod.mcmc_A_base_long_met$sims.list$beta2_south, size=1),
                  period1=365,
                  period2=sample(prod.mcmc_A_base_long_met$sims.list$period2_south, size=1),
                  x1=sample(prod.mcmc_A_base_long_met$sims.list$x1_south, size=1),
                  x2=sample(prod.mcmc_A_base_long_met$sims.list$x2_south, size=1),
                  phi=0,
                  theta=0
                  # input parameters
                  
        )   
      }else if (MR == "Lima"){
        params<-c(b=365*3.5,
                  D_inf=sample(prod.mcmc_A_base_long_met$sims.list$alpha, size=1)+0.25,
                  D_inc=sample(prod.mcmc_A_base_long_met$sims.list$delta, size=1),
                  D_imm=sample(prod.mcmc_A_base_long_met$sims.list$gamma, size=1), 
                  beta0=(sample(prod.mcmc_A_base_long_met$sims.list$beta0_lima, size=1)-0.5)*total_pop,
                  beta1=sample(prod.mcmc_A_base_long_met$sims.list$beta1_lima, size=1),
                  beta2=sample(prod.mcmc_A_base_long_met$sims.list$beta2_lima, size=1),
                  period1=365,
                  period2=sample(prod.mcmc_A_base_long_met$sims.list$period2_lima, size=1),
                  x1=sample(prod.mcmc_A_base_long_met$sims.list$x1_lima, size=1),
                  x2=sample(prod.mcmc_A_base_long_met$sims.list$x2_lima, size=1),
                  phi=0,
                  theta=0
                  # input parameters
                  
        )   
      }
      
      r<-data.frame(ssa.sptime(ssa.adaptivetau(init.values = init.values,transitions = transitions_SEIR,
                                               SIRrateF_SEIR, params = params,tf=time_run,
                                               tl.params =  list(maxtau = 1,extraChecks = F),
                                               set.seed(x1)),
                               time.steps = seq(0,time_run,30))); # set the time steps that will be output
      
      count<-append(count, tail(r$I,1))
      
    }
    
    Count<-(sims-length(count[count==0]))
    
    population_size_test$Count[n]<-Count/sims
    
    n<-n+1
    
  }
  
  return(population_size_test)
  
}

start_time<-Sys.time()

stoch_sim_lima<-stoch_sim(MR="Lima", 250, 2500, 250000, 2500)
stoch_sim_north<-stoch_sim(MR="North", 250, 2500, 250000, 2500)
stoch_sim_south<-stoch_sim(MR="South", 250, 2500, 250000, 2500)

end_time<-Sys.time()
print(end_time-start_time)

load("model_output_data/stochastic_model_pop_size_output_20241021.RData")

stoch_sim_lima<-rbind(data.frame(PopSize=0, Count=0.00), stoch_sim_lima)
stoch_sim_north<-rbind(data.frame(PopSize=0, Count=0.00), stoch_sim_north)
stoch_sim_south<-rbind(data.frame(PopSize=0, Count=0.00), stoch_sim_south)

ggplot()+
  geom_point(data=stoch_sim_lima, aes(x=PopSize/1000, y=Count), col="#44BB99")+
  #geom_smooth(data=stoch_sim_lima, aes(x=PopSize/1000, y=Count), method='lm', formula=y~poly(x,15), se=FALSE, col="#44BB99", linewidth=1)+
  geom_point(data=stoch_sim_north, aes(x=PopSize/1000, y=Count), col="#99DDFF")+
  #geom_smooth(data=stoch_sim_north, aes(x=PopSize/1000, y=Count), method='lm', formula=y~poly(x,14), se=FALSE, col="#99DDFF", linewidth=1)+
  geom_point(data=stoch_sim_south, aes(x=PopSize/1000, y=Count), col="#EE8866")+
  #geom_smooth(data=stoch_sim_south, aes(x=PopSize/1000, y=Count), method='lm', formula=y~poly(x,14), se=FALSE, col="#EE8866", linewidth=1)+
  theme_bw()+ylim(-0.01,1.01)+
  geom_hline(aes(yintercept=0.90), linetype="dashed")+
  labs(x="Population size (K)", y="BIV survival rate")


