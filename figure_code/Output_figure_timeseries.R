rm(list=ls())

#load packages
{
  library(boot)
  library(tidyverse)
  library(multipanelfigure)
  library(Hmisc)
}


#######################################################################################################################################
#load data - requires download from repository at https://zenodo.org/records/12820416
load("All_model_outputs_20240531.RData")

##############################
##############################
##                          ##
##     FIGURE PLOTTING      ##
##                          ##
##############################
##############################

############ Figure 1

#packages
{
  library(sf)
  library(raster)
  library(ggiraph)
  library(ggrepel)
  library(geodata)
  library(tidyverse)
  library(geos)
  library(RColorBrewer)
  library(gridExtra)
  library(ggmap)
  library(ggpubr)
  library(ggspatial)
}


#register_google(key = "xxx", write = TRUE)

flu_data<-read.csv("raw_data/bat_final_unique_MG_day.csv")
site_coords<-read.csv("raw_data/Sample_sites_coords_all.csv")
bat_cull <- read.csv("raw_data/Bat_cull.csv", header=T)
Laura_metagenomic_samples <- read.csv("raw_data/Laura_metagenomic_samples.csv", header=T)



culling_vector<-c(rep(0, times=144))

for(i in 1:nrow(bat_cull)){
  n<-bat_cull[i,8]
  culling_vector[n-1]<-bat_cull[i,5]
}


culling_combined<-as.data.frame(culling_vector)

# counts the number of x in l
count <- function(x,l){
  sum(l[which(!(is.na(l)))] == x)
}

culling_combined$time<-c(1:144)
culling_combined$transformed<-(-culling_combined$culling_vector/3000)+1


positions <- data.frame(
  x = c(rep(0:143, each=2),144),
  y = c(1,rep(culling_combined$transformed, each=2))
)


{
  sample_y_north[is.na(sample_y_north)]<-0
  sample_y_south[is.na(sample_y_south)]<-0
  sample_y_lima[is.na(sample_y_lima)]<-0
  
  sample_i_y_north[is.na(sample_i_y_north)]<-0
  sample_i_y_south[is.na(sample_i_y_south)]<-0
  sample_i_y_lima[is.na(sample_i_y_lima)]<-0
  
  north_bin<-data.frame(time=c(1:144), binconf(x=sample_y_north, n=sample_n_north))
  south_bin<-data.frame(time=c(1:144), binconf(x=sample_y_south, n=sample_n_south))
  lima_bin<-data.frame(time=c(1:144), binconf(x=sample_y_lima, n=sample_n_lima))
  
  north_i_bin<-data.frame(time=c(1:144), binconf(x=sample_i_y_north, n=sample_i_n_north), sample_i_n_north)
  south_i_bin<-data.frame(time=c(1:144), binconf(x=sample_i_y_south, n=sample_i_n_south), sample_i_n_south)
  lima_i_bin<-data.frame(time=c(1:144), binconf(x=sample_i_y_lima, n=sample_i_n_lima), sample_i_n_lima)
  
}

###################################
###################################
##      Plot seroprevalence      ##
###################################
###################################

#Figure 1A - recovereds in the North+input population-level serology data
north_R<-ggplot()+
  geom_line(data=time_series_data_north, aes(x=time, y=r.mean),col="#99ddff", size=1)+
  geom_ribbon(data=time_series_data_north, aes(x=time, ymax=r.upp, ymin=r.lwr), fill="#99ddff", alpha=0.3 )+
  theme_bw()+labs(x="Time (years)", y="Proportion \nSeropositive")+
  geom_point(data=data_points_north, aes(x=time, y=PointEst, col=inside_q_north), size=1)+
  geom_errorbar(data=data_points_north, aes(x=time, ymax=Upper, ymin=Lower, col=inside_q_north))+
  scale_x_continuous(limits=c(2007.2,2019), breaks=seq(2008,2018,2))+
  scale_y_continuous( breaks=seq(0,1,0.2))+
  scale_colour_manual(values=c("red", "black"))+theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 1), clip="off")+
  scale_shape_manual(values=c(4,19))+
  geom_point(data=filter(north_i_bin, PointEst==0), aes(x=(time/12)+2007, y=1.1))+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(2,1,1,0), "lines"))

#Figure 1B - recovereds in Lima+input population-level serology data
lima_R<-ggplot()+
  geom_line(data=time_series_data_lima, aes(x=time, y=r.mean),col="#44BB99", size=1)+
  geom_ribbon(data=time_series_data_lima, aes(x=time, ymax=r.upp, ymin=r.lwr), fill="#44BB99", alpha=0.2 )+
  theme_bw()+labs(x="Time (years)", y="Proportion \nSeropositive")+
  geom_point(data=data_points_lima, aes(x=time, y=PointEst, col=inside_q_lima), size=1)+
  geom_errorbar(data=data_points_lima, aes(x=time, ymax=Upper, ymin=Lower, col=inside_q_lima))+
  scale_x_continuous(limits=c(2007.2,2019), breaks=seq(2008,2018,2))+
  scale_y_continuous( breaks=seq(0,1,0.2))+
  scale_colour_manual(values=c("black", "red"))+theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 1), clip="off")+
  scale_shape_manual(values=c(4,19))+
  geom_point(data=filter(lima_i_bin, PointEst==0), aes(x=(time/12)+2007, y=1.1))+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(2,1,1,0), "lines"))

#Figure 1C - recovereds in the South+input population-level serology data
south_R<-ggplot()+
  geom_bar(data = culling_combined, aes(x=c(1:144)/12+2007, y=culling_vector/3000), stat="identity", alpha=0.2)+
  geom_line(data=r_p_plot, aes(x=(time/360)+2007, y=mean),col="#ee8866", size=1)+
  geom_ribbon(data=r_p_plot, aes(x=(time/360)+2007, ymax=q2.5, ymin=q97.5), fill="#ee8866", alpha=0.2 )+
  theme_bw()+labs(x="Time (years)", y="Proportion \nSeropositive")+
  geom_point(data=data_points_south, aes(x=time, y=PointEst, col=inside_q_south), size=1)+
  geom_errorbar(data=data_points_south, aes(x=time, ymax=Upper, ymin=Lower, col=inside_q_south))+
  scale_x_continuous(limits=c(2007.2,2019), breaks=seq(2008,2018,2))+
  scale_colour_manual(values=c("black", "red"))+theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 1), clip="off")+
  scale_shape_manual(values=c(4,19))+
  geom_point(data=filter(south_i_bin, PointEst==0), aes(x=(time/12)+2007, y=1.1))+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(2,1,1,0), "lines"))+
  coord_cartesian(ylim = c(0, 1), clip="off")+
  scale_shape_manual(values=c(4,19))+
  scale_y_continuous(breaks=seq(0,1.15,0.2),
                     
                     # Features of the first axis
                     
                     
                     # Add a second axis and specify its features
                     sec.axis = sec_axis(~.*3000, name="Culling effort"))



###################################
###################################
##  Plot longitudinal summary    ##
###################################
###################################



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




longitudinal_summary<-data.frame("MR"=rep(c("North", "Central", "South"), each=4), "State"=rep(c("0-0", "0-1", "1-0", "1-1"), times=3), "Count"=NA)
#1=00, 2=01, 3=10, 4=11

count_transitions<-function(df, s1, s2){
  mid<-filter(df, df[,1]==s1 & df[,2]==s2)
  return(nrow(mid))
}

x<-1

for (n in c(1,2,3)){
  for (i in c(0,1)){
    for (j in c(0,1)){
      
      if(n==1){
        df<-as.data.frame(B_sim_north)
      }else if(n==2){
        df<-as.data.frame(B_sim_lima)
      }  else{
        df<-as.data.frame(B_sim_south)
      }
      
      Count <- count_transitions(df, i, j)
      
      longitudinal_summary[x,3]<-Count
      
      x<-x+1
      
    }
  }
}

longitudinal_summary

long_plot<- ggplot()+
  geom_bar(data=longitudinal_summary, aes(x=State, y=Count, fill=factor(MR, levels = c("North", "Central", "South"))), stat="identity")+
  theme_bw()+scale_fill_manual(values=c( "#0077BB", "#44BB99","#EE8866"), labels=c("N", "C", "S"))+
  labs(fill=" ", x="Transition", y="Longitudinal sample pairs")+theme(legend.position = "top")

long_plot 





###################################
###################################
##         Map plotting          ##
###################################
###################################

newd <-  flu_data %>% group_by(Site) %>% filter(mean(Date_continuous)/max(Date_continuous)!=0) #
newd <- filter(newd, MetaRegionB %in% c("South_a", "South_b", "North", "Lima"))
sites<-as.vector(unique(newd$Site))
plot_sites<-filter(site_coords, Code %in% sites)
plot_sites$MR<-c("N", "N", "N", "N", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "N", "N", "N", "N", "S", "S", "L", "L", "L", "L")

PER_G <- gadm(country="PE", level=1, path="path_to_your_map_data")
PERU_G_df <-sf::st_as_sf(PER_G, region = "NAME_1")

PERU_G_df_South <- filter((PERU_G_df), NAME_1 %in% c("Apurímac", "Ayacucho", "Cusco"))
PERU_G_df_North <- filter((PERU_G_df), NAME_1 %in% c("Cajamarca", "Amazonas"))
PERU_G_df_Lima <- filter((PERU_G_df), NAME_1 %in% c("Lima", "Lima Province"))

map<-
  ggplot(data=PERU_G_df)+geom_sf()+theme_void()+
  geom_sf(data=PERU_G_df_North, fill="#0077BB", alpha=0.25)+
  geom_sf(data=PERU_G_df_South, fill="#EE8866", alpha=0.25)+
  geom_sf(data=PERU_G_df_Lima, fill="#44BB99", alpha=0.25)+
  geom_point(data=plot_sites, aes(x=longitude, y=latitude, col=MR), size=2)+
  scale_colour_manual(values = c("#44BB99", "#0077BB", "#EE8866"))+
  theme(legend.position = "none")+annotation_scale()


map


###################################
###################################
##       Plot all together       ##
###################################
###################################

library(multipanelfigure)



fig_1_plot<-multi_panel_figure(width = c(102,8,45), height=c(58,40,20,58), unit="mm", row_spacing = 0, column_spacing = c(0,5,0,0))
fig_1_plot


fig_1<-
  fill_panel(fig_1_plot, north_R, column = 1, row=1, label = "A")%>%
  fill_panel(lima_R, column=1, row=c(2:3))%>%
  fill_panel(south_R, column=c(1:2), row=4)%>%
  fill_panel(map, column=c(2:3), row=c(1:2))%>%
  fill_panel(long_plot, column=3, row=c(3:4))

####################################################################################################


####################################################################################################

############ Figure 2

model_long_test_output_both<-read.csv("model_output_data/test_data_JAGS_output.csv")

#Figure 2A - alpha and beta parameter posterior distribution and comparison with using the population-level data only
MR3_all_alpha_delta<-ggplot()+
  annotate("rect", xmin=prod.mcmc_A_base$q2.5$alpha, xmax=prod.mcmc_A_base$q97.5$alpha, ymin=-Inf, ymax=Inf, alpha=0.2, fill="grey80")+
  geom_vline(aes(xintercept=prod.mcmc_A_base$mean$alpha), col="grey80", linewidth=1)+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$alpha, xmax=prod.mcmc_A_base_long_met$q97.5$alpha, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#bb5566")+
  geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$alpha), aes(x=prod.mcmc_A_base$sims.list$alpha), col="grey80", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$alpha), aes(x=prod.mcmc_A_base_long_met$sims.list$alpha), col="#bb5566", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$alpha), col="#bb5566", linewidth=1)+
  annotate("rect", xmin=prod.mcmc_A_base$q2.5$delta, xmax=prod.mcmc_A_base$q97.5$delta, ymin=-Inf, ymax=Inf, alpha=0.2, fill="grey80")+
  geom_vline(aes(xintercept=prod.mcmc_A_base$mean$delta), col="grey80", linewidth=1)+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$delta, xmax=prod.mcmc_A_base_long_met$q97.5$delta, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#ddaa33")+
  geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$delta), aes(x=prod.mcmc_A_base$sims.list$delta), col="grey80", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$delta), aes(x=prod.mcmc_A_base_long_met$sims.list$delta), col="#ddaa33", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$delta), col="#ddaa33", linewidth=1)+
  labs(x=expression(paste("Incubation ", delta,  " and Infectious ", alpha, " Periods (Days)",)), y="Density")+ylim(0,0.8)+
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,2))+
  stat_function(fun = dgamma, args = list(shape=3.5, rate=0.5), col = "#bb5566", linetype="dashed")+
  stat_function(fun = dgamma, args = list(shape=6, rate=2), col = "#ddaa33", linetype="dashed")+
  theme_bw()

#Figure 2B - gamma parameter posterior distribution and comparison with using the population-level data only
MR3_all_gamma<-ggplot()+
  annotate("rect", xmin=prod.mcmc_A_base$q2.5$gamma, xmax=prod.mcmc_A_base$q97.5$gamma, ymin=-Inf, ymax=Inf, alpha=0.15, fill="grey70")+
  annotate("rect", xmin=quantile(prod.mcmc_A_base_long_met$sims.list$gamma, prob=0.025), xmax=quantile(prod.mcmc_A_base_long_met$sims.list$gamma, prob=0.975), ymin=-Inf, ymax=Inf, alpha=0.1, fill="#0097D5")+
  geom_vline(aes(xintercept=prod.mcmc_A_base$q50$gamma), col="grey70", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base$sims.list$gamma), aes(x=prod.mcmc_A_base$sims.list$gamma), col="grey70", linewidth=1)+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$gamma), aes(x=prod.mcmc_A_base_long_met$sims.list$gamma), col="#0097D5", linewidth=1)+
  geom_vline(aes(xintercept=median(prod.mcmc_A_base_long_met$sims.list$gamma)), col="#0097D5", linewidth=1)+
  labs(x=expression(paste("Immune Period " , gamma, " (Days)  ")), y="Density")+xlim(0,500)+ylim(0,0.025)+
  theme_bw()+
  stat_function(fun = dgamma, args = list(shape=13.3225, rate=0.0365), col = "#0097D5", linetype="dashed") 


#Figure 2C - gamma test data comparison

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
    scale_color_manual(labels = c("Population-level \ndata only", "All data"), values=c("grey80","#0097D5"))+
    scale_fill_manual(labels = c("Population-level \ndata only", "All data"), values=c("grey80","#0097D5"))
  
  return(q)
}




gamma_plot<-parameter_plot("gamma", "gamma_m", "gamma_sd")
gamma_plot<-
  gamma_plot+labs(x=expression(paste("Input: Immune Period ",  gamma, " (Days)")), y=expression(paste("Output: Immune Period ",  gamma, " (Days)")))+
  scale_color_manual(labels = c("Population-level \ndata only", "All data"), values=c("grey70","#0097D5"))+
  scale_fill_manual(labels = c("Population-level \ndata only", "All data"), values=c("grey70","#0097D5"))+
  
  annotate("text", x = 390, y = 1100, label = "Full dataset:     RMSE = 57.53",
           col="#0097D5") +
  annotate("text", x = 390, y = 1030, label = "Base dataset:   RMSE = 74.10",
           col="grey70") +
  
  coord_cartesian(ylim = c(100,900), clip = "off") +
  theme(plot.margin = unit(c(3,1,1,0), "lines"))


gamma_plot



#plot parts A, B, and C together
figure_2_plot<-multi_panel_figure(width=c(90,70), height=c(45,45), unit = "mm")
figure_2_plot

fig_2<-
  fill_panel(figure_2_plot, MR3_all_alpha_delta, column=1, row=1)%>%
  fill_panel(MR3_all_gamma, column=1, row=2)%>%
  fill_panel(gamma_plot, column=2, row=c(1:2))


####################################################################################################

#--------------------------------------------------------------------------------

############ Figure 3

#Figure 3A - The transmission rate beta over time for each meta-region
all_beta<-ggplot()+
  geom_line(data=time_series_data_north, aes(x=time, y=beta.mean),col="#99ddff", size=1)+
  geom_ribbon(data=time_series_data_north, aes(x=time, ymax=beta.upp, ymin=beta.lwr), fill="#99ddff", alpha=0.3 )+
  geom_line(data=time_series_data_south, aes(x=time, y=beta.mean),col="#ee8866", size=1)+
  geom_ribbon(data=time_series_data_south, aes(x=time, ymax=beta.upp, ymin=beta.lwr), fill="#ee8866", alpha=0.2 )+
  geom_line(data=time_series_data_lima, aes(x=time, y=beta.mean),col="#44BB99", size=1)+
  geom_ribbon(data=time_series_data_lima, aes(x=time, ymax=beta.upp, ymin=beta.lwr), fill="#44BB99", alpha=0.2 )+
  theme_bw()+labs(x="Time (years)", y="Transmission Rate")+
  scale_x_continuous(limits=c(2007.2,2019), breaks=seq(2008,2018,2))+
  scale_y_continuous(limits=c(0,1))

#Figure 3B - the proportion of infected bats in the population over time for each meta-region
all_I<-ggplot()+
  geom_line(data=time_series_data_north, aes(x=time, y=i.mean),col="#99ddff", size=1)+
  geom_ribbon(data=time_series_data_north, aes(x=time, ymax=i.upp, ymin=i.lwr), fill="#99ddff", alpha=0.3 )+
  geom_line(data=i_p_plot, aes(x=(time/365)+2007, y=mean),col="#ee8866", size=1)+
  geom_ribbon(data=i_p_plot, aes(x=(time/365)+2007, ymax=q2.5, ymin=q97.5), fill="#ee8866", alpha=0.2 )+
  geom_line(data=time_series_data_lima, aes(x=time, y=i.mean),col="#44BB99", size=1)+
  geom_ribbon(data=time_series_data_lima, aes(x=time, ymax=i.upp, ymin=i.lwr), fill="#44BB99", alpha=0.2 )+
  theme_bw()+labs(x="Time (years)", y="Proportion Infected")+
  scale_x_continuous(limits=c(2007.2,2019), breaks=seq(2008,2018,2))+
  scale_y_continuous(limits=c(0,0.15))

#get the mean beta and infecteds over time for each meta-region as quoted in the manuscript
mean(time_series_data_north$i.mean[72:4320])*100
mean(i_p_plot$mean[72:4320])*100
mean(time_series_data_lima$i.mean[72:4320])*100

range(time_series_data_north$i.mean[72:4320])*100
range(i_p_plot$mean[72:4320])*100
range(time_series_data_lima$i.mean[72:4320])*100


mean(time_series_data_north$beta.mean[72:4320])
mean(time_series_data_south$beta.mean[72:4320])
mean(time_series_data_lima$beta.mean[72:4320])

range(time_series_data_north$beta.mean[72:4320])
range(time_series_data_south$beta.mean[72:4320])
range(time_series_data_lima$beta.mean[72:4320])


#Find at what time in the year peaks occur in each meta-region

findPeaks<-function (x, thresh = 0) 
{
  pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0) + 2
  if (!missing(thresh)) {
    pks[x[pks - 1] - x[pks] > thresh]
  }
  else pks
}


for (i in c(2007.2, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018)){
  
  year_lima<-filter(time_series_data_lima,  time >= i & time < round(i)+1 )
  year_north<-filter(time_series_data_north,  time >= i & time < round(i)+1 )
  year_south<-filter(time_series_data_south,  time >= i & time < round(i)+1 )
  
  print((year_lima$time[findPeaks(year_lima$i.mean)]-(round(i)))*12)
  print((year_north$time[findPeaks(year_north$i.mean)]-(round(i)))*12)
  print((year_south$time[findPeaks(year_south$i.mean)]-(round(i)))*12)
  
}

peak_times<-read.csv("model_output_data/Peak_times.csv")

#get mean month of peak occurence for each meta-region
mean(subset(peak_times, MR=="Central")$PeakTime)
mean(subset(peak_times, MR=="North")$PeakTime)
mean(subset(peak_times, MR=="South" & is.na(PeakTime)==FALSE)$PeakTime)


peak_times$MR <- factor(peak_times$MR , levels=c("North", "Central", "South"))

peak_times_plot<-
  ggplot(peak_times, aes(x=MR, y=PeakTimeNew+0.5))+geom_boxplot(aes(fill=MR), show.legend = FALSE)+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12), limits = c(1,12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  labs(x="Zone", y="Month")+scale_fill_manual(values = c( "#99ddff","#44BB99", "#ee8866"))+
  theme_bw()+
  annotate("rect", ymin=2.05+0.5, ymax=4.9+0.5, xmin=-Inf,xmax=Inf, alpha=0.2, fill="#ee8866")+
  annotate("rect", ymin=4.7+0.5, ymax=6.15+0.5, xmin=-Inf,xmax=Inf, alpha=0.2, fill="#99ddff")+
  annotate("rect", ymin=8.07+0.5, ymax=9.65+0.5, xmin=-Inf,xmax=Inf, alpha=0.2, fill="#44bb99")+
  geom_boxplot(aes(fill=MR), show.legend = FALSE)

peak_times_plot


#Plot parts A, B, and C together
fig_3_plot<-multi_panel_figure(width = c(100,60), height=c(50,50), unit="mm")
fig_3_plot

fig_3<-
fill_panel(fig_3_plot, all_beta, column = 1, row=1)%>%
  fill_panel(all_I, column=1, row=2)%>%
  fill_panel(peak_times_plot, column=2, row=c(1,2))


####################################################################################################


####################################################################################################

############ Figure 4

#Plot the posterior distribution for the csf parameter (here, "cull_int") (Figure 5A)
cull_scale_plot<-
  ggplot()+
  annotate("rect", xmin=prod.mcmc_A_base_long_met$q2.5$cull_int, xmax=prod.mcmc_A_base_long_met$q97.5$cull_int, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#0097D5")+
  geom_density(data=as.data.frame(prod.mcmc_A_base_long_met$sims.list$cull_int), aes(x=prod.mcmc_A_base_long_met$sims.list$cull_int), col="#0097D5", linewidth=1)+
  geom_vline(aes(xintercept=prod.mcmc_A_base_long_met$mean$cull_int), col="#0097D5", linewidth=1)+
  labs(x="Culling scaling factor", y="Density")+xlim(0,0.006)+ylim(0,600)+
  theme_bw()+
  stat_function(fun = dgamma, args = list(shape=1.5, rate=400), col = "black", linetype="dashed")


#Plot the change in the total population size for the South metaregion (Figure 5B)
Pop_size<-
  ggplot()+
  #geom_line(data=r_p_plot, aes(x=(time/365)+2007, y=mean),col="#ee8866", size=1, alpha=0.5)+
  #geom_ribbon(data=r_p_plot, aes(x=(time/365)+2007, ymax=q2.5, ymin=q97.5), fill="#ee8866", alpha=0.1 )+
  geom_line(data=N_plot, aes(x=(time/360)+2007, y=mean), col="black", size=1)+
  geom_ribbon(data=N_plot, aes(x=(time/360)+2007, ymax=q97.5, ymin=q2.5), fill="black", alpha=0.2 )+
  geom_line(data=time_series_data_south_p2, aes(x=time, y=i.mean*5),col="grey50", size=1)+
  geom_ribbon(data=time_series_data_south_p2, aes(x=time, ymax=i.upp*5, ymin=i.lwr*5), fill="grey50", alpha=0.2 )+
  geom_line(data=time_series_data_south, aes(x=time, y=i.mean*5),col="#ee8866", size=1)+
  geom_ribbon(data=time_series_data_south, aes(x=time, ymax=i.upp*5, ymin=i.lwr*5), fill="#ee8866", alpha=0.2 )+
  geom_vline(aes(xintercept=2014.5), linetype="dashed", size=0.60)+
  geom_vline(aes(xintercept=2016.5), linetype="dashed", size=0.60)+
 # annotate("text", x = 2012.6, y = 0.15, label = "Culling \nStart",
 #          hjust = 0, vjust = 1, angle = 90, size = 4)+
 # annotate("text", x = 2016.8, y = 0.15, label = "Culling \nEnd",
 #          hjust = 0, vjust = 1, angle = 90, size = 4)+
  theme_bw()+labs(x="Time (years)", y="Total Population")+
  scale_x_continuous(limits=c(2007.2,2019), breaks=seq(2008,2018,2))+
  scale_y_continuous(limits=c(0,1.01), breaks=seq(0,1.01,0.25),
                     
                     # Features of the first axis
                     
                     
                     # Add a second axis and specify its features
                     sec.axis = sec_axis(~./5, name="Proportion Infected"))




fig_4_plot<-multi_panel_figure(width = c(60,100), height=c(70,1), unit="mm", row_spacing = 0)

fig_4<-
fill_panel(fig_4_plot, cull_scale_plot, row=1, column = 1)%>%
  fill_panel(Pop_size, row=1, column=2)


a1<-mean(i_p_plot$mean[2738:4320])*100
b1<-mean(time_series_data_south_p2$i.mean[2738:4320])*100

a2<-mean(i_p_plot$mean[2738:3468])*100
b2<-mean(time_series_data_south_p2$i.mean[2738:3468])*100

#percentage decrease in the mean proportion of the population Infected due to culling compared to predicted had culling not taken place.
100-(a1/b1)*100
100-(a2/b2)*100
