#load packages
library(tidyverse)
library(ggplot2)
library(plotly)

#load data
flu_data<-read.csv("raw_data/bat_final_unique_MG_day.csv")
peak_times<-read.csv("model_output_data/Peak_times.csv")

flu_data2<-filter(flu_data, MetaRegionB %in% c("North", "Lima", "South_a", "South_b"))
flu_data2$MetaRegion<-factor(flu_data2$MetaRegion, levels=c("North", "Lima", "South"),
                             labels = c("North", "Central", "South"))

peak_times$PeakTimeNew_day<-peak_times$PeakTimeNew*30

peak_times_means<-data.frame(MetaRegion=c("North", "South", "Central"), 
                             mean=c(mean(filter(peak_times, MR == "North")$PeakTimeNew_day),
                                    mean(filter(peak_times, MR == "South" & PeakTime != "NA")$PeakTimeNew_day),
                                    mean(filter(peak_times, MR == "Central")$PeakTimeNew_day)), 
                             upper=c(quantile(filter(peak_times, MR == "North")$PeakTimeNew_day, probs = c(0.75)),
                                     quantile(filter(peak_times, MR == "South" & PeakTime != "NA")$PeakTimeNew_day, probs = c(0.75)),
                                     quantile(filter(peak_times, MR == "Central")$PeakTimeNew_day, probs = c(0.75))),
                             lower=c(quantile(filter(peak_times, MR == "North")$PeakTimeNew_day, probs = c(0.25)),
                                     quantile(filter(peak_times, MR == "South" & PeakTime != "NA")$PeakTimeNew_day, probs = c(0.25)),
                                     quantile(filter(peak_times, MR == "Central")$PeakTimeNew_day, probs = c(0.25))))

ggplot(data=flu_data2)+
  geom_vline(aes(xintercept = 1), col="black", linetype="dashed")+
  geom_hline(aes(yintercept = 2350))+
  geom_rect(data=peak_times_means, aes(xmin=lower, xmax=upper, ymin=0, ymax=Inf), fill = "grey50", alpha=0.5)+
  geom_vline(data=peak_times_means, aes(xintercept=mean))+
  geom_point(aes(x=((Month-1)/12+Day/365)*365, y=as.factor(Id), col=Year), alpha = 1)+
  coord_polar()+theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(2.5, "cm"),
        legend.title.position = "top",
        strip.background = element_rect(
          color="white", fill="white", size=1.5, linetype="solid"))+
  labs(x="Collection Day")+
  scale_x_continuous(breaks = seq(from=0, to=365, by=30), limits = c(1, 365))+
  scale_color_stepsn(n.breaks = 12, colours = viridis::viridis(12))+
  facet_wrap(~factor(MetaRegion, levels = c("North", "Central", "South")), scales = "free")
  

