setwd()

library(tidyverse)

cetacean_raw_data<-read.csv("cleaned_raw_data_CFD.csv")

ggplot(cetacean_raw_data,aes(x=max_length_m))+
  geom_histogram(bins=11)

ggplot(cetacean_raw_data,aes(x=max_mass_kg))+
  geom_histogram(bins=11)

ggplot(cetacean_raw_data,aes(x=max_mass_max_length_ratio))+
  geom_histogram(bins=11)

ggplot(cetacean_raw_data,aes(x=log(max_length_m, 10)))+
  geom_histogram(bins=22)

ggplot(cetacean_raw_data,aes(x=log(max_mass_kg, 10)))+
  geom_histogram(bins=11)

ggplot(cetacean_raw_data,aes(x=log(max_mass_max_length_ratio, 10)))+
  geom_histogram(bins=12)

cetacean_dataset<-read.csv("mtuck_honours_Rdocument.csv")

ggplot(cetacean_dataset, aes(x=dentition))+
  geom_bar(stat="count")
  
ggplot(cetacean_dataset, aes(x=migratory_behaviour))+
  geom_bar(stat="count")

ggplot(cetacean_dataset, aes(x=max_diving_depth))+
  geom_bar(stat="count")

ggplot(cetacean_dataset, aes(x=average_group_size))+
  geom_bar(stat="count")

ggplot(cetacean_dataset, aes(x=prey_choice))+
  geom_bar(stat="count")


