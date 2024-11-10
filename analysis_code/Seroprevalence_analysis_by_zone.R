#Find seroprevalence of BIV across
##(a) all the bat colonies used in this model
##(b) each Zone individually

library(Hmisc)
library(tidyverse)

#load full dataset
bat_flu<-read.csv("raw_data/bat_final_unique_MG_day.csv")

#filter the data to only look at sites used in the model
bat_flu_all_sites<-filter(bat_flu, Site %in% c("AMA1",
                                               "AMA2",
                                               "AMA3",
                                               "API1",
                                               "API13",
                                               "API140",
                                               "API141",
                                               "API15",
                                               "API16",
                                               "API17",
                                               "API18",
                                               "API3",
                                               "API9",
                                               "AYA1",
                                               "AYA11",
                                               "AYA13",
                                               "AYA14",
                                               "AYA15",
                                               "AYA4",
                                               "AYA7",
                                               "CAJ1",
                                               "CAJ2",
                                               "CAJ3",
                                               "CAJ4",
                                               "CUS6",
                                               "CUS8",
                                               "LMA10",
                                               "LMA4",
                                               "LMA5",
                                               "LMA6"
))

#get binomial confidence interval for seroprevalence
binconf(x=sum(bat_flu_all_sites$Serol), n=nrow(bat_flu_all_sites))

#filter for each Zone and get the relevant binconfs
bat_flu_north<-filter(bat_flu_all_sites, MetaRegion == "North")
bat_flu_south<-filter(bat_flu_all_sites, MetaRegion == "South")
bat_flu_central<-filter(bat_flu_all_sites, MetaRegion == "Lima")

binconf(x=sum(bat_flu_north$Serol), n=nrow(bat_flu_north))
binconf(x=sum(bat_flu_south$Serol), n=nrow(bat_flu_south))
binconf(x=sum(bat_flu_central$Serol), n=nrow(bat_flu_central))

#looking for a difference in seroprevalence in colonies that are
## DR only vs co-roosting with Artibeus/Sturnira species
## In the pre-culling period only

a_s_sites<-filter(bat_flu, Site %in% c("AMA1", 
                                       "API13",
                                       "AYA11", 
                                       "LMA4"
                                       ) & Date_continuous < 2014.5
                  )

non_a_s_sites<-filter(bat_flu_all_sites, !Site %in% c("AMA1", 
                                                      "API13", 
                                                      "AYA11", 
                                                      "LMA4"
                                                      ) & Date_continuous < 2014.5
                      )

binconf(x=sum(a_s_sites$Serol), n=nrow(a_s_sites))
binconf(x=sum(non_a_s_sites$Serol), n=nrow(non_a_s_sites))

