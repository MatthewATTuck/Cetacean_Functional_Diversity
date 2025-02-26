####Family Breakdown Figure####
setwd()
library(patchwork)
library(tidyverse)
library(ggpubr)

Family_data<-read.csv("Family_proportions_data.csv")
Location<-c(rep("Global", 11), rep("OPB", 11))
Family_data$Location<-Location
Family_data

family_proportions_plot<-ggplot(Family_data, aes(x=Location, y=Family_Proportions, fill=Common_Names))+
  geom_bar(stat = "identity", width=0.5, colour="black")+
  scale_fill_manual(values = c("#C15CCB", "orange", "blue", "springgreen", "grey23", "darksalmon", "#CC0099", "yellow", "midnightblue", "#FF0000", "#00CCFF"))+
  theme_bw(base_size = 10)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title = element_blank(), legend.key.size = unit(0.8, "cm"))+
  scale_y_continuous(limits=c(0,1.1), expand = c(0, 0))+
  theme(axis.title.y = element_text(size=20, vjust = 2)) +
  theme(axis.title.x = element_text(size=20, vjust = 0), axis.text.y =element_text(size=15, colour = "black"), axis.text.x =element_text(size=15, colour="black"))+
  labs(x = "Species Pool", y = "Proportion of Species")

family_count_plot<-ggplot(Family_data, aes(x=Location, y=Family_Counts, fill=Common_Names))+
  geom_bar(stat = "identity", width=0.5, colour="black", show.legend=FALSE)+
  scale_fill_manual(values = c("#C15CCB", "orange", "blue", "springgreen", "grey23", "darksalmon", "#CC0099", "yellow", "midnightblue", "#FF0000", "#00CCFF"))+
  theme_bw(base_size = 10)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())+
  scale_y_continuous(limits=c(0,96.8), expand = c(0, 0))+
  theme(axis.text.y=element_text(size=15, colour="black"))+
  theme(axis.title.y = element_text(size=20, vjust = 2)) +
  theme(axis.title.x = element_text(size=20, vjust = 0), axis.text.y =element_text(size=15, colour = "black"), axis.text.x =element_text(size=15, colour="black"))+
  labs(x = "Species Pool", y = "Number of Species")  

ggarrange(family_count_plot, family_proportions_plot, ncol=2, nrow=1, align = "v")


####PCEEZ Map####
library(rnaturalearthdata)
library(rnaturalearth)
library(sf)
library(ggplot2)

EEZ <- read_sf(dsn = "~/GitHub/Matthew-Tuck-honours/MatthewTuck_Honours/clean_data/eez_iho/eez_iho.shp")
EEZ
#Flanders Marine Institute (2023). Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), version 12. Available online at http://www.marineregions.org/. https://doi.org/10.14284/632

OPB <- read_sf(dsn = "~/GitHub/Matthew-Tuck-honours/MatthewTuck_Honours/clean_data/OPB_data/FederalMarineBioregions_SHP/FederalMarineBioregions.shp")
OPB <-OPB%>%
  filter(NAME_E=="Offshore Pacific")
OPB
## CANADA MAP
## Set the map theme
maptheme_blank <-  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
## Load in the world Data
world <- ne_countries(scale = "medium", returnclass = "sf")
# Create Large Map of Canada
ggplot() +
  xlab("") +
  ylab("") +
  geom_sf(data = world)+
  coord_sf(xlim = c(-170, -100), ylim = c(30, 75))+ # define the plot limits
  maptheme_blank+
  geom_point(aes(y = 51, x = -128), colour = "red", cex = 30, shape = 0, stroke = 2) + maptheme_blank

# Create Map of EEZ
## Set the map theme
maptheme <-  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
## Create the plot of the EEZ
ggplot() +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_sf(data = EEZ, fill = "lightblue", color = "lightblue") +
  geom_sf(data = world)+
  coord_sf(xlim = c(-140, -120), ylim = c(46.5, 56.5))+ # define the plot limits
  maptheme

# Create Map of OPB
## Set the map theme
maptheme <-  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
## Create the plot of the OPB
ggplot() +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_sf(data = OPB, fill = "lightblue", color = "lightblue") +
  maptheme
ggplot() +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_sf(data = world)+
  theme(axis.text.x =element_text(size=25, colour="black"))+
  coord_sf(xlim = c(-140, -120), ylim = c(46.5, 56.5))+ # define the plot limits
  maptheme

####Unique Species BarPlot####
Location<-c("Global", "Global", "PCEEZ", "PCEEZ")
unique_barplot<-as.data.frame(Location, )
unique_barplot$Location<-(unique_barplot)  
unique_barplot$counts<-c(43,45,15,10)
unique_barplot$proportions<-c(0.488636363636364, 0.511363636363636, 0.6, 0.4)
unique_barplot$Category<-c("Unique", "Not Unique", "Unique", "Not Unique")
unique_barplot

ggplot(unique_barplot, aes(x=Location$Location, y=proportions, fill=Category))+
  geom_bar(stat = "identity", width=0.5, colour="black")+
  scale_fill_manual(values = c( "orange",  "blue"), 
                    labels = c("Not Unique", "Unique"))+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(size=15, colour = "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right", legend.title = element_text(size=20, colour="black"), legend.text=element_text(size=18, colour="black")) +
  scale_y_continuous(limits=c(0, 1.1), expand = c(0, 0))+
  theme(axis.text.y=element_text(size=15, colour="black"), axis.text.x=element_text(size=18, colour="black"))+
  theme(axis.title.y = element_text(size=24, vjust = 2)) +
  theme(axis.title.x = element_text(size=24, vjust = 0))+
  labs(x = "Species Pool", y = "Proportion of Total Species")


####Unique_Species_Histogram####
unique_sp_null_plot<-as.data.frame(unique_sp_perm_null)
unique_sp_null_plot
unique_sp_mean_plot<-as.data.frame(mean_unique_sp_perm_nm)

ggplot(unique_sp_null_plot, aes(x=unique_sp_perm_null))+
  geom_histogram(bins=12, colour="#000000", fill="#999999")+
  geom_vline(aes(xintercept = mean(unique_sp_perm_null)), color = "red", linetype = "dashed", linewidth = 1)+
  geom_vline(aes(xintercept = c(mean(unique_sp_perm_null)+(1.96*sd(unique_sp_perm_null)))), color = "blue", linetype = "dashed", linewidth = 1)+
  geom_vline(aes(xintercept = c(mean(unique_sp_perm_null)-(1.96*sd(unique_sp_perm_null)))), color = "blue", linetype = "dashed", linewidth = 1)+
  geom_vline(aes(xintercept = 15), linewidth=1)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x="Number of Unique Species", y="Frequency")

ggplot(unique_sp_null_plot, aes(x=unique_sp_perm_null))+
  geom_histogram(bins=12, colour="#000000", fill="#999999")+
  geom_vline(aes(xintercept = mean(unique_sp_perm_null)), color = "red", linetype = "dashed", linewidth = 1)+
  geom_vline(aes(xintercept = 15), linewidth=1)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x="Number of Unique Species", y="Frequency")


####PCoA Plots####
PCoA_1<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(2)]))+
  geom_point(size=10, aes(colour = ind_pceez)) +
  scale_colour_manual(values=c("orange", "blue"), labels=c("Not Found in OPB", "Found in OPB"))+
  labs(x="PCoA 1 (62.76% of variation)", y="PCoA 2 (31.32% of variation)", colour="Presence in PCEEZ")+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  scale_y_continuous(limits=c(-0.55, 0.45), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.55, 0.55), expand=c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(size=15, colour = "black"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))
PCoA_1

PCoA_2<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(3)]))+
  geom_point(size=10, aes(colour = ind_pceez)) +
  scale_colour_manual(values=c("orange", "blue"), labels=c("Not Found in PCEEZ", "Found in PCEEZ"))+
  labs(x="PCoA 1 (62.76% of variation)", y="PCoA 3 (8.85% of variation)", colour="Presence in PCEEZ")+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right", legend.title = element_text(size=28), legend.text = element_text(size=24)) +
  theme(axis.text=element_text(size=15, colour = "black"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))

PCoA_3<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(2)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(3)]))+
  geom_point(size=10, aes(colour = ind_pceez), show.legend = FALSE) +
  scale_colour_manual(values=c("orange", "blue"), labels=c("Not Found in PCEEZ", "Found in PCEEZ"))+
  labs(x="PCoA 2 (31.32% of variation)", y="PCoA 3 (8.85% of variation)", colour="Presence in PCEEZ")+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.text=element_text(size=15, colour = "black"))+
  theme(axis.title.y = element_text(size= 20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))

ggarrange(PCoA_1, PCoA_2, PCoA_3, ncol=1, labels=c("a","b","c"),font.label=list(color="black",size=20), align = "v")

####Plots for Functional Indicies (NEED TO UPDATE AFTER SETTING SEEDS OR ANY OTHER CHANGE TO FUNCTIONAL INDICIES)####
multif_null_plot<-c(multif_perm_nm)
multif_null_plot<-as.data.frame(multif_null_plot)
multif_null_plot<-t(multif_null_plot)

ggplot(multif_null_plot, aes(x=PCEEZ))+
  geom_histogram(bins=20, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=fdis_mean_nm), colour="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(fdis_mean_nm+(1.96*fdis_sd_nm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(fdis_mean_nm-(1.96*fdis_sd_nm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=fdis_mean_PCEEZ), colour="black", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=15, color="black")) +
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))+
  labs(x="Functional Dispersion", y="Frequency")

ggplot(multif_null_plot, aes(x=PCEEZ))+
  geom_histogram(bins=20, colour="black", fill="grey")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=fdis_mean_nm), colour="black", linetype="dashed", size=2)+
  geom_vline(aes(xintercept=fdis_mean_PCEEZ), colour="blue",linetype="dashed", size=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=15, color="black")) +
  scale_y_continuous(limits=c(0,175), expand=c(0,0))+
  theme(axis.title.y = element_text(size=20, vjust = 3))+
  theme(axis.title.x = element_text(size=20, vjust = 0))+
  labs(x="Functional Dispersion", y="Frequency")

fe_null_plot<-c(fe_perm_nm)
fe_null_plot<-as.data.frame(fe_null_plot)
fe_null_plot<-t(fe_null_plot)
filter<-(rep(c(1,0), 1000))
filter<-as.factor(filter)
fe_null_plot<-as.data.frame(fe_null_plot)
fe_null_plot$filter<-filter

fe_null_plot<-fe_null_plot%>%
  filter(filter=="1")

ggplot(fe_null_plot, aes(x=PCEEZ))+
  geom_histogram(bins=12, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=fred_mean_nm), colour="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(fred_mean_nm+(1.96*fred_sd_nm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(fred_mean_nm-(1.96*fred_sd_nm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=fred_mean_PCEEZ), colour="black", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=15, color="black")) +
  scale_y_continuous(expand=c(0,1))+
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x="Functional Redundancy", y="Frequency")

ggplot(fe_null_plot, aes(x=PCEEZ))+
  geom_histogram(bins=12, colour="black", fill="grey")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=fred_mean_nm), colour="black", linetype="dashed", size=2)+
  geom_vline(aes(xintercept=fred_mean_PCEEZ), colour="blue", linetype="dashed", size=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=15, color="black")) +
  scale_y_continuous(limits=c(0,250), expand=c(0,0))+
  theme(axis.title.y = element_text(size=20, vjust = 3))+
  theme(axis.title.x = element_text(size=20, vjust = 0))+
  labs(x="Functional Redundancy", y="Frequency")



####Histograms for Continuous Traits (UPDATE AFTER SETTING SEED OR MAKING CHANGES TO BASE CODE)####
length_histogram<-ggplot(continuous_trait_plots, aes(x=max_length_perm))+
  geom_histogram(bins=25, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=mean(max_length_perm)), colour="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(mean_max_length_perm+(1.96*sd_max_length_perm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(mean_max_length_perm-(1.96*sd_max_length_perm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=morphometric_means$mean_length), colour="black", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x="Mean Maximum Length (m)", y="Frequency")

mass_histogram<-ggplot(continuous_trait_plots, aes(x=max_mass_perm))+
  geom_histogram(bins=25, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=mean(max_mass_perm)), colour="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(mean_max_mass_perm+(1.96*sd_max_mass_perm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(mean_max_mass_perm-(1.96*sd_max_mass_perm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=morphometric_means$mean_mass), colour="black", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x="Mean Maximum Mass (kg)", y="Frequency")

ratio_histogram<-ggplot(continuous_trait_plots, aes(x=max_ratio_perm))+
  geom_histogram(bins=25, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=mean(max_ratio_perm)), colour="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(mean_max_ratio_perm+(1.96*sd_max_ratio_perm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(mean_max_ratio_perm-(1.96*sd_max_ratio_perm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=morphometric_means$mean_ratio), colour="black", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x="Mean Maximum Mass/Length Ratio (kg/m)", y="Frequency")

ggarrange(length_histogram, mass_histogram, ratio_histogram, ncol=1, labels=c("a","b","c"))


####Boxplots for individual morphometric traits####
box1<-ggplot(comparitive_boxplots, aes(x=Location, y=`Max Length (m)`, fill = Location))+
  geom_boxplot(outliers = FALSE, width=0.5, show.legend = FALSE)+
  scale_fill_manual(values=c("orange", "blue"))+
  labs(x="Species Pool", y="Maximum Length (m)")+
  theme_bw(base_size = 15)+
  theme(text = element_text(size = 10))+
  theme(axis.text=element_text(size=24, colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y = element_text(size=36, vjust = 1)) +
  theme(axis.title.x = element_text(size=36, vjust = 0))
box1


box2<-ggplot(comparitive_boxplots, aes(x=Location, y=`Max Mass (kg)`, fill=Location))+
  geom_boxplot(outliers = FALSE, width=0.5, show.legend = FALSE)+
  scale_fill_manual(values=c("orange", "blue"))+
  labs(x="Species Pool", y="Maximum Mass (kg)")+
  theme_bw(base_size = 15)+
  theme(axis.text=element_text(size=24, colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y = element_text(size=36, vjust = 1)) +
  theme(axis.title.x = element_text(size=36, vjust = 0))
box2

box3<-ggplot(comparitive_boxplots, aes(x=Location, y=`Max Mass/Max Length Ratio (kg/m)`, fill=Location))+
  geom_boxplot(outliers = FALSE, width=0.5, show.legend = FALSE)+
  scale_fill_manual(values=c("orange", "blue"))+
  labs(x="Species Pool", y="Mass/Length (kg/m)")+
  theme_bw(base_size = 15)+
  theme(axis.text=element_text(size=24, colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y = element_text(size=36, vjust = 1)) +
  theme(axis.title.x = element_text(size=36, vjust = 0))
box3

####CWM plots####
CWMs_plots_split<-split(CWMs_plots, f=CWMs_plots$trait_types)
CMW_plot_length<-ggplot(CWMs_plots_split$`Max Length`, aes(fill=Categories, x=location, y=CWMs))+
  geom_bar(position="stack", stat = "identity", width=0.5, colour="black")+
  labs(x="Species Pool", y="Community Weighted Means")+
  scale_fill_manual(values=c("#C15CCB", "orange", "blue", "springgreen", "grey23", "darksalmon"))+
  facet_wrap(~trait_types, ncol=1)+
  theme_bw(base_size = 15)+
  theme(axis.text=element_text(size=24, colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(limits=c(0,1.1), expand=c(0,0))+ 
  theme(axis.title.y = element_text(size=36, vjust = 1)) +
  theme(axis.title.x = element_text(size=36, vjust = 1))+
  theme(strip.text=element_text(size=40, colour="black"))+
  theme(strip.background = element_rect(fill = "white"))+
  theme(legend.position = "right",legend.justification = "left", legend.margin = margin(l=1, unit="cm"), legend.text=element_text(size=28), legend.title = element_blank())
  


CWM_plot_mass<-CMW_plot_length%+%CWMs_plots_split$`Max Mass`
CWM_plot_ratio<-CMW_plot_length%+%CWMs_plots_split$`Max Mass/Length Ratio`
CWM_plot_dentition<-CMW_plot_length%+%CWMs_plots_split$Dentition
CWM_plot_diving<-CMW_plot_length%+%CWMs_plots_split$`Max Diving Depth`
CWM_plot_group_size<-CMW_plot_length%+%CWMs_plots_split$`Average Group Size`
CWM_plot_prey<-CMW_plot_length%+%CWMs_plots_split$`Prey Choice`

CWM_plot_dentition
CWM_plot_diving
CWM_plot_group_size
CWM_plot_prey

####Completed individual trait plot####
box_plots_combined<-ggarrange(box1,box2,box3,ncol=1,labels=c("a", "b", "c"), font.label=list(color="black",size=28), align = "v")
box_plots_combined
CWM_plots_combined<-ggarrange(CWM_plot_dentition, CWM_plot_diving, CWM_plot_group_size, CWM_plot_prey, ncol=2, nrow=2, labels=c("d", "e", "f", "g"), font.label=list(color="black",size=28), align="v")
CWM_plots_combined
Individual_traits_plot<-(((box_plots_combined)|  (CWM_plots_combined)) + plot_layout(widths=c(1,4)))
Individual_traits_plot
