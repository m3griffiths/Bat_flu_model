
library(tidyverse)
library(multipanelfigure)

setwd("Bat_flu_modelling")


flu_data<-read.csv("raw_data/bat_final_unique_MG_day.csv")


#######################################################################################################################################


#######################################################################################################################################
#multi-sampled bats

bats_flu_multi <- flu_data[which(flu_data$Id %in% flu_data[duplicated(flu_data$Id),]$Id),]


bats_flu_multi$state <- NA

id_names<-unique(bats_flu_multi$Id)

for (i in c(unique(bats_flu_multi$Id))){
  
  id_int<-filter(bats_flu_multi, Id == i)
  
  sero_c<-c()
  
  for (n in 1:nrow(id_int)){
    
    sero_c<-append(sero_c, id_int[n,4])
    
  }
  
  if (sum(sero_c)==0){
    bats_flu_multi$state[bats_flu_multi$Id == i]<-0
  }else if (sum(sero_c)==length(sero_c)){
    bats_flu_multi$state[bats_flu_multi$Id == i]<-0
  }else if (sero_c[1]>sero_c[2] & length(sero_c)==2){
    bats_flu_multi$state[bats_flu_multi$Id == i]<-1
  }else if (sero_c[1]<sero_c[2] & length(sero_c)==2){
    bats_flu_multi$state[bats_flu_multi$Id == i]<-2
  }else{
    bats_flu_multi$state[bats_flu_multi$Id == i]<-3
  }
  
  
  
  
}

north_multi<-
ggplot(data = filter(bats_flu_multi, MetaRegion=="North"), aes(x=Date_continuous, y=Id))+
  geom_line(aes(col=as.factor(state)))+
  geom_point(aes(shape=as.factor(Serol), col=as.factor(state)), size=1.5, fill="white", stroke=1)+theme_bw()+
  labs(x="Time (years)", y= "Bat ID", shape="Seropositivity", title="North")+
  scale_shape_manual(values=c(4,19), labels=c("Seronegative (-)", "Seropositive (+)"))+
  scale_color_manual(values = c("grey80", "royalblue", "red", "orange"), labels=c("Unchanged", "Antibody waning", "Seroconversion", "Other"))+
  scale_x_continuous(breaks=seq(2008,2018,2))+theme(legend.position = "none", axis.text.y = element_text(size=5))

south_multi<-
ggplot(data = filter(bats_flu_multi, MetaRegion=="South"), aes(x=Date_continuous, y=Id))+
  geom_line(aes(col=as.factor(state)))+
  geom_point(aes(shape=as.factor(Serol), col=as.factor(state)), size=1.5, fill="white", stroke=1)+theme_bw()+
  labs(x="Time (years)", y= "Bat ID", shape="Seropositivity", title="South")+
  scale_shape_manual(values=c(4,19), labels=c("Seronegative (-)", "Seropositive (+)"))+
  scale_color_manual(values = c("grey80", "royalblue", "red", "orange"), labels=c("Unchanged", "Antibody waning", "Seroconversion", "Other"))+
  scale_x_continuous(breaks=seq(2008,2018,2))+theme(legend.position = "none", axis.text.y = element_text(size=5))

lima1_multi<-
ggplot(data = filter(bats_flu_multi, MetaRegion=="Lima" & Id <= 4500), aes(x=Date_continuous, y=Id))+
  geom_line(aes(col=as.factor(state)))+
  geom_point(aes(shape=as.factor(Serol), col=as.factor(state)), size=1.5, fill="white", stroke=1)+theme_bw()+
  labs(x="Time (years)", y= "Bat ID", shape="Seropositivity", title="Coast")+
  scale_shape_manual(values=c(4,19), labels=c("Seronegative (-)", "Seropositive (+)"))+
  scale_color_manual(values = c("grey80", "royalblue", "red", "orange"), labels=c("Unchanged", "Antibody waning", "Seroconversion", "Other"))+
  scale_x_continuous(breaks=seq(2008,2018,2))+theme(legend.position = "none", axis.text.y = element_text(size=5))

lima2_multi<-
  ggplot(data = filter(bats_flu_multi, MetaRegion=="Lima" & Id > 4500), aes(x=Date_continuous, y=Id))+
  geom_line(aes(col=as.factor(state)))+
  geom_point(aes(shape=as.factor(Serol), col=as.factor(state)), size=1.5, fill="white", stroke=1)+theme_bw()+
  labs(x="Time (years)", y= "Bat ID", shape="Seropositivity", title="")+
  scale_shape_manual(values=c(4,19), labels=c("Seronegative (-)", "Seropositive (+)"))+
  scale_color_manual(values = c("grey80", "royalblue", "red", "orange"), labels=c("Unchanged", "Antibody waning", "Seroconversion", "Other"))+
  scale_x_continuous(breaks=seq(2008,2018,2))+theme(legend.position = "none", axis.text.y = element_text(size=5))



legend1_plot<-ggplot(data = filter(bats_flu_multi, MetaRegion=="Lima"), aes(x=Date_continuous, y=Id))+
  geom_line(aes())+
  geom_point(aes(shape=as.factor(Serol)), size=1.5, fill="white", stroke=1)+theme_bw()+
  labs(x="Time (years)", y= "Bat ID", shape="Seropositivity", title="Coast")+
  scale_shape_manual(values=c(4,19), labels=c("Seronegative (-)", "Seropositive (+)"))+
  scale_x_continuous(breaks=seq(2008,2018,2))+theme(legend.position = "left", axis.text.y = element_text(size=5))


legend2_plot<-ggplot(data = filter(bats_flu_multi, MetaRegion=="Lima"), aes(x=Date_continuous, y=Id))+
  geom_line(aes(col=as.factor(state)))+
  geom_point(aes( col=as.factor(state)), size=1.5, fill="white", stroke=1)+theme_bw()+
  labs(x="Time (years)", y= "Bat ID", shape="Seropositivity", title="Coast", colour="Change in serostatus")+
  scale_color_manual(values = c("grey80", "royalblue", "red", "orange"), labels=c("Serostatus Unchanged", "Antibody waning (1 -> 0)", "Seroconversion (0 -> 1)", "Combination"))+
  scale_x_continuous(breaks=seq(2008,2018,2))+theme(legend.position = "left", axis.text.y = element_text(size=5))



legend1<-as_ggplot(get_legend(legend1_plot))
legend2<-as_ggplot(get_legend(legend2_plot))



fig_S1a_plot<-multi_panel_figure(width = c(80,80), height=c(150,15,15), unit="mm", column_spacing = c(0,0,0,0))
fig_S1a_plot

fig_S1b_plot<-multi_panel_figure(width = c(80,80), height=c(200,20), unit="mm", column_spacing = c(0,0,0,0))
fig_S1b_plot


fill_panel(fig_S1a_plot, north_multi, column = 1, row=1, label = "")%>%
  fill_panel(south_multi, column=2, row=1, label = "")%>%
  fill_panel(legend2, column=c(1), row=c(2:3), label = "")%>%
  fill_panel(legend1, column=c(2), row=c(2:3), label = "")

fill_panel(fig_S1b_plot, lima1_multi, column = 1, row=1, label = "")%>%
  fill_panel(lima2_multi, column=2, row=1, label = "")
