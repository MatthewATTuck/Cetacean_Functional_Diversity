####Getting silhouettes for desired species####
library(rphylopic)
get_uuid(name= "Tursiops truncatus", n=4, filter="sa")
####Family Breakdown Figure####
setwd()
library(patchwork)
library(tidyverse)
library(ggpubr)

Family_data<-read.csv("Family_proportions_data.csv")
Location<-c(rep("Global", 11), rep("OPB", 11))
Family_data$Location<-Location
Family_data

family_proportions_plot<-ggplot(Family_data, aes(x=Location, y=Family_Proportions, fill=Family))+
  geom_bar(stat = "identity", width=0.5, colour="black")+
  scale_fill_manual(values = c("#C15CCB", "#FF9900", "#0000FF","#00ff66","#333333", "#e9967a", "#ff3399", "yellow", "midnightblue", "#FF0000", "#00CCFF"))+
  theme_bw(base_size = 10)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title = element_blank(), legend.key.size = unit(0.8, "cm"))+
  scale_y_continuous(limits=c(0,1.1), expand = c(0, 0))+
  theme(axis.title.y = element_text(size=20, vjust = 2)) +
  theme(axis.title.x = element_text(size=20, vjust = 0), axis.text.y =element_text(size=15, colour = "black"), axis.text.x =element_text(size=15, colour="black"))+
  labs(x = "Species Pool", y = "Proportion of Species")

family_count_plot<-ggplot(Family_data, aes(x=Location, y=Family_Counts, fill=Family))+
  geom_bar(stat = "identity", width=0.5, colour="black", show.legend=FALSE)+
  scale_fill_manual(values = c("#C15CCB", "#FF9900", "#0000FF","#00ff66","#333333", "#e9967a", "#ff3399", "yellow", "midnightblue", "#FF0000", "#00CCFF"))+
  theme_bw(base_size = 10)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())+
  scale_y_continuous(limits=c(0,96.8), expand = c(0, 0))+
  theme(axis.text.y=element_text(size=15, colour="black"))+
  theme(axis.title.y = element_text(size=20, vjust = 2)) +
  theme(axis.title.x = element_text(size=20, vjust = 0), axis.text.y =element_text(size=15, colour = "black"), axis.text.x =element_text(size=15, colour="black"))+
  labs(x = "Species Pool", y = "Number of Species")  

ggarrange(family_count_plot, family_proportions_plot, ncol=2, nrow=1, labels=c("A", "B"), font.label=list(color="black",size=20), align = "v")


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

####Histogram for Number of Baleen####
baleen_sp_null<-as.vector(baleen_sp_perm_null)
view(baleen_sp_null)
baleen_sp_null<-data.frame(baleen_sp_null)

ggplot(baleen_sp_null, aes(x=baleen_sp_null))+
  geom_histogram(bins=10, colour="#000000", fill="#999999")+
  geom_vline(aes(xintercept = mean(baleen_sp_null)), color = "black", linetype = "dashed", linewidth = 1)+
  geom_vline(aes(xintercept = 7), colour="blue", linetype="dashed", linewidth=1)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  theme(axis.text=element_text(size=15, color="black")) +
  scale_y_continuous(limits=c(0,300), expand=c(0,0))+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9))+
  theme(axis.title.y = element_text(vjust = 3, size=20)) +
  theme(axis.title.x = element_text(size=20))+
  labs(x="Number of Mysticete Species", y="Frequency")


####Unique Species BarPlot####
Location<-c("Global", "Global", "OPB", "OPB")
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
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=18, colour="black")) +
  scale_y_continuous(limits=c(0, 1.1), expand = c(0, 0))+
  theme(axis.text.y=element_text(size=15, colour="black"), axis.text.x=element_text(size=18, colour="black"))+
  theme(axis.title.y = element_text(size=20, vjust = 2)) +
  theme(axis.title.x = element_text(size=20, vjust = 0))+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
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

PCoA_1_families<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(2)]))+
  geom_point(size=7, aes(colour = cetacean_dataset$Family)) +
  add_phylopic(name = c("Balaenoptera", "Ziphiidae"), filter="sa", x=c(-0.25, -0.45), y=c(-0.4, 0.35), height=0.1)+
  add_phylopic(uuid = c("388e792c-f8fd-4bd8-ad38-e251b5244ac0", "6e892ffc-4443-4213-8c97-bb9c0950f0b7"), filter="sa", x=c(0.5, -0.5), y=c(0.3, 0.3), height=0.1)+
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
PCoA_1_families

PCoA_length<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(2)]))+
  geom_point(size=7, aes(colour = sp_traits_nm$max_length_m)) +
  labs(x="PCoA 1 (62.76% of variation)", y="PCoA 2 (31.32% of variation)", colour="Trait Category")+
  scale_colour_manual(values = c("#C15CCB", "#FF9900", "#0000FF","#00ff66"), labels=c("Small", "Intermediate", "Large", "Very Large"))+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  scale_y_continuous(limits=c(-0.55, 0.45), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.55, 0.55), expand=c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(size=15, colour = "black"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))
PCoA_length

PCoA_mass<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(2)]))+
  geom_point(size=7, aes(colour = sp_traits_nm$max_mass_kg)) +
  scale_colour_manual(values = c("#C15CCB", "#FF9900", "#0000FF","#00ff66"), labels=c("Small", "Intermediate", "Large", "Very Large"))+
  labs(x="PCoA 1 (62.76% of variation)", y="PCoA 2 (31.32% of variation)")+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  scale_y_continuous(limits=c(-0.55, 0.45), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.55, 0.55), expand=c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(size=15, colour = "black"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))
PCoA_mass

PCoA_ratio<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(2)]))+
  geom_point(size=7, aes(colour = sp_traits_nm$max_mass_max_length_ratio)) +
  scale_colour_manual(values = c("#C15CCB", "#FF9900", "#0000FF"), labels=c("Low", "Medium", "High"))+
  labs(x="PCoA 1 (62.76% of variation)", y="PCoA 2 (31.32% of variation)")+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  scale_y_continuous(limits=c(-0.55, 0.45), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.55, 0.55), expand=c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(size=15, colour = "black"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))
PCoA_ratio

PCoA_dentition<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(2)]))+
  geom_point(size=7, aes(colour = sp_traits_nm$dentition)) +
  scale_colour_manual(values = c("#C15CCB", "#FF9900", "#0000FF","#00ff66"), labels=c("Long Baleen", "Short Baleen", "Functional Calcareous Teeth", "Non-Functional Calcareous Teeth"))+
  labs(x="PCoA 1 (62.76% of variation)", y="PCoA 2 (31.32% of variation)")+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  scale_y_continuous(limits=c(-0.55, 0.45), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.55, 0.55), expand=c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(size=15, colour = "black"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))
PCoA_dentition

PCoA_diving<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(2)]))+
  geom_point(size=7, aes(colour = sp_traits_nm$max_diving_depth)) +
  scale_colour_manual(values = c("#C15CCB", "#FF9900", "#0000FF","#00ff66"), labels=c("Epipelagic", "Upper Mesopelagic", "Lower Mesopelagic", "Bathypelagic"))+
  labs(x="PCoA 1 (62.76% of variation)", y="PCoA 2 (31.32% of variation)")+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  scale_y_continuous(limits=c(-0.55, 0.45), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.55, 0.55), expand=c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(size=15, colour = "black"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))
PCoA_diving

PCoA_groups<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(2)]))+
  geom_point(size=7, aes(colour = sp_traits_nm$average_group_size)) +
  scale_colour_manual(values = c("#C15CCB", "#FF9900", "#0000FF","#00ff66","#333333", "#e9967a"), labels=c("Solitary", "Small", "Small-Medium", "Medium", "Medium-Large", "Large"))+
  labs(x="PCoA 1 (62.76% of variation)", y="PCoA 2 (31.32% of variation)")+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  scale_y_continuous(limits=c(-0.55, 0.45), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.55, 0.55), expand=c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(size=15, colour = "black"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))
PCoA_groups

PCoA_prey<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(2)]))+
  geom_point(size=7, aes(colour = sp_traits_nm$prey_choice)) +
  scale_colour_manual(values = c("#C15CCB", "#FF9900", "#0000FF","#00ff66","#333333", "#e9967a"), labels=c("Zooplankton", "Non-Cephalopod Invertebrates", "Cephalopods", "Fish", "High Vertebrates"))+
  labs(x="PCoA 1 (62.76% of variation)", y="PCoA 2 (31.32% of variation)")+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  scale_y_continuous(limits=c(-0.55, 0.45), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.55, 0.55), expand=c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(size=15, colour = "black"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))
PCoA_prey


PCoA_2<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(3)]))+
  geom_point(size=10, aes(colour = ind_pceez), show.legend=FALSE) +
  scale_colour_manual(values=c("blue", "orange"), labels=c("Not Found in OPB", "Found in OPB"))+
  labs(x="PCoA 1 (62.76% of variation)", y="PCoA 3 (8.85% of variation)", colour="Presence in PCEEZ")+
  theme_bw(base_size = 15)+
  scale_y_continuous(limits=c(-0.45, 0.3), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.55, 0.55), expand=c(0,0))+
  theme(text=element_text(size = 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_blank()) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  theme(axis.text=element_text(size=15, colour = "black"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))
PCoA_2

PCoA_3<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(2)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(3)]))+
  geom_point(size=10, aes(colour = ind_pceez), show.legend = FALSE) +
  scale_colour_manual(values=c("blue", "orange"), labels=c("Not Found in PCEEZ", "Found in PCEEZ"))+
  labs(x="PCoA 2 (31.32% of variation)", y="PCoA 3 (8.85% of variation)", colour="Presence in PCEEZ")+
  theme_bw(base_size = 15)+
  scale_y_continuous(limits=c(-0.45, 0.3), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.55, 0.45), expand=c(0,0))+
  theme(text=element_text(size = 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  theme(axis.text=element_text(size=15, colour = "black"))+
  theme(axis.title.y = element_text(size= 20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))
PCoA_3

ggarrange(PCoA_2, PCoA_3, ncol=2, nrow=1, labels=c("A", "B"),font.label=list(color="black",size=20), align = "v")
ggarrange(PCoA_length, PCoA_mass, PCoA_ratio, PCoA_dentition, PCoA_diving, PCoA_groups, PCoA_prey, ncol=4, nrow=2, labels = c("Maximum Length (m)", "Maximum Mass (kg)", "Maximum Mass/Length Ratio (kg/m)", "Dentition", "Maximum Diving Depth", "Average Group Size", "Prey Choice"), font.label = list(colour="black", size=28), align = "v")

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
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
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
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
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
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  theme(axis.title.y = element_text(size=20, vjust = 3))+
  theme(axis.title.x = element_text(size=20, vjust = 0))+
  labs(x="Functional Redundancy", y="Frequency")



####Histograms for Continuous Traits (UPDATE AFTER SETTING SEED OR MAKING CHANGES TO BASE CODE)####
length_histogram<-ggplot(continuous_trait_plots, aes(x=max_length_perm))+
  geom_histogram(bins=20, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=mean(max_length_perm)), colour="black", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(mean_max_length_perm+(1.96*sd_max_length_perm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(mean_max_length_perm-(1.96*sd_max_length_perm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=morphometric_means$mean_length), colour="black", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x="Mean Maximum Length (m)", y="Frequency")

length_histogram<-ggplot(continuous_trait_plots, aes(x=max_length_perm))+
  geom_histogram(bins=20, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=mean(max_length_perm)), colour="black", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=morphometric_means$mean_length), colour="blue", linetype="dashed", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.text=element_text(size=15, color="black")) +
  scale_y_continuous(limits=c(0,175), expand=c(0,0))+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))+
  labs(x="Mean Maximum Length (m)", y="Frequency")

mass_histogram<-ggplot(continuous_trait_plots, aes(x=max_mass_perm))+
  geom_histogram(bins=20, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=mean(max_mass_perm)), colour="black", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=morphometric_means$mean_mass), colour="blue", linetype="dashed", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.text=element_text(size=15, color="black")) +
  scale_y_continuous(limits=c(0,175), expand=c(0,0))+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))+
  labs(x="Mean Maximum Mass (kg)", y="Frequency")

ratio_histogram<-ggplot(continuous_trait_plots, aes(x=max_ratio_perm))+
  geom_histogram(bins=20, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=mean(max_ratio_perm)), colour="black", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=morphometric_means$mean_ratio), colour="blue", linetype="dashed", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.text=element_text(size=15, color="black")) +
  scale_y_continuous(limits=c(0,175), expand=c(0,0))+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  theme(axis.title.y = element_text(size=20, vjust = 3)) +
  theme(axis.title.x = element_text(size=20, vjust = -1))+
  labs(x="Mean Maximum Mass/Length Ratio (kg/m)", y="Frequency")

ggarrange(length_histogram, mass_histogram, ratio_histogram, ncol=1, labels=c("A","B","C"), font.label=list(color="black",size=20), align = "v")


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
  scale_fill_manual(values=c("#C15CCB", "#FF9900", "#0000FF","#00ff66","#333333", "#e9967a"))+
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
CWM_plots_combined<-ggarrange(CWM_plot_dentition, CWM_plot_diving, CWM_plot_group_size, CWM_plot_prey, ncol=2, nrow=2, labels=c("A", "B", "C", "D"), font.label=list(color="black",size=28), align="v")
CWM_plots_combined
Individual_traits_plot<-(((box_plots_combined)|  (CWM_plots_combined)) + plot_layout(widths=c(1,4)))
Individual_traits_plot




####Best Figures for PCoA Plots####

# Create a FE coordinate matrix in the functional space from species dissimilarities
fe_faxes <- data.frame(PCoA_1 = sp_faxes_nm[, "PC1"], PCoA_2 = sp_faxes_nm[,"PC2"], PCoA_3 = sp_faxes_nm[,"PC3"]) # select species 4 first axes
fe_faxes$fe <- func_ent_nm$sp_fe # add FE to each species
fe_faxes <- fe_faxes[c(which(duplicated(fe_faxes$fe) == FALSE)),] # select unique FE 
rownames(fe_faxes) <- fe_faxes$fe # change rownames
fe_faxes <- as.matrix(fe_faxes[,-4]) # delete FE column

# community matrix of FE per province
fe_comm_biog <- func_ind_fe_nm$details_fdfe$asb_fe_nbsp 

# Trait distribution along the Global Functional Space
######################################################

fe_tr <- as.data.frame(func_ent_nm$fe_tr) # Functional traits per FE
colnames(fe_tr) <- c("max_length", "max_mass", "ratio", "dentition", "max_diving", "group_size", "prey_choice")
fe_color <- data.frame(max_length = rep(NA, 56), max_mass = rep(NA, 56), ratio= rep(NA, 56), dentition = rep(NA, 56), max_diving= rep(NA, 56), group_size = rep(NA, 56), prey_choice = rep(NA, 56)) # dataframe for colors
rownames(fe_color) <- rownames(fe_faxes)
legend_info <- as.list(NA) # list for legend colors

color_gr_length <- c("#C15CCB", "#FF9900", "#0000FF","#00ff66") # color palette as matrix
ord_gr <- c("small", "intermediate", "large", "very_large") # ordinal categories for max length
legend_info[["max_length"]] <- matrix(color_gr_length, nrow = 1, ncol = 4, dimnames = list(1, c(ord_gr)))
for (i in 1:4) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$max_length == ord_gr[i])),])), "max_length"] <- color_gr_length[i]
}

ord_gr <- c("small","intermediate", "large", "very_large") # ordinal categories of max mass
legend_info[["max_mass"]] <- matrix(color_gr_length, nrow = 1, ncol = 4, dimnames = list(1, c(ord_gr)))

for (i in 1:4) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$max_mass == ord_gr[i])),])), "max_mass"] <- color_gr_length[i]
}

color_gr_ratio <- c("#C15CCB", "#FF9900", "#0000FF")
ord_gr <- c("low", "medium", "high") # ordinal categories for mass/length ratios
legend_info[["ratio"]] <- matrix(color_gr_ratio, nrow = 1, ncol = 3, dimnames = list(1, c(ord_gr)))

for (i in 1:3) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$ratio == ord_gr[i])),])), "ratio"] <- color_gr_ratio[i]
}

ord_gr <- c("LB", "SB", "FCT", "NFCT") # ordinal categories for dentition
legend_info[["dentition"]] <- matrix(color_gr_length, nrow = 1, ncol = 4, dimnames = list(1, c(ord_gr)))

for (i in 1:4) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$dentition == ord_gr[i])),])), "dentition"] <- color_gr_length[i]
}


ord_gr <- c("epipelagic", "upper_mesopelagic", "lower_mesopelagic", "bathypelagic") # ordinal categories for max diving depth
legend_info[["max_diving"]] <- matrix(color_gr_length, nrow = 1, ncol = 4, dimnames = list(1, c(ord_gr)))

for (i in 1:4) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$max_diving == ord_gr[i])),])), "max_diving"] <- color_gr_length[i]
}

color_gr_group <- c("#C15CCB", "#FF9900", "#0000FF","#00ff66","#333333", "#e9967a")
ord_gr <- c("Solitary", "Small", "Medium_Small", "Medium", "Medium_Large", "Large") # ordinal categories of average group size
legend_info[["group_size"]] <- matrix(color_gr_group, nrow = 1, ncol = 6, dimnames = list(1, c(ord_gr)))

for (i in 1:6) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$group_size == ord_gr[i])),])), "group_size"] <- color_gr_group[i]
}

color_gr_prey <- c("#C15CCB", "#FF9900", "#0000FF","#00ff66","#333333")
ord_gr <- c("zooplankton", "non_cephalopod_invertebrates", "cephalopods", "fish", "high_vertebrates") # ordinal categories of average group size
legend_info[["prey_choice"]] <- matrix(color_gr_prey, nrow = 1, ncol = 5, dimnames = list(1, c(ord_gr)))

for (i in 1:5) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$prey_choice == ord_gr[i])),])), "prey_choice"] <- color_gr_prey[i]
}


#Setting up colour_fig (just in case)
color_fig = data.frame(color_p = c("#FFA500","#3A5FCD"),
color_ch = c("#FFA50070","#3A5FCD70")) # color points and convex hulls of figures
row.names(color_fig) <- rownames(fe_comm_biog)
write.csv(color_fig, "color_fig.csv")

install.packages("tripack")
library(tripack)
glob_tr <- tri.mesh(fe_faxes[,1],fe_faxes[,2]) # trigonometry Global Functional Space 1st and 2nd PCoA axes
glob_ch <- convex.hull(glob_tr) # convex hull Global Functional Space

# Plot province and Global functional spaces (1st & 2nd PCoA axes)
par(mfrow = c(1,2), mai=c(rep(0.2, 4))) # set in between plots sizes

for (i in rownames(fe_comm_biog)) {
  
  plot(fe_faxes[,1],fe_faxes[,2], type="n",
       xlab = NA, ylab = NA, xaxt = "n", yaxt = "n",
       xlim = c(-0.55, 0.55), ylim = c(-0.55, 0.45), main = i) # plot 
  axis(1, labels = FALSE)
  axis(2, labels = FALSE)
  
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey95") # grey background
  polygon(glob_ch, col="white", border=NA) # Global Functional Space
  
  fe_prov <- colnames(fe_comm_biog)[which(fe_comm_biog[i,] > 0)] # select FE present at each province
  fe_coord <- fe_faxes[rownames(fe_faxes) %in% fe_prov, c(1,2)] # select coordinates for present FE
  tr <-tri.mesh(fe_coord[,1],fe_coord[,2]) # compute trigonometry
  ch <- convex.hull(tr) # convex hull
  
  polygon(ch, col = color_fig[i,2], border = NA) # plot convex hull 
  points(fe_coord[,1], fe_coord[,2], col = color_fig[i,1], pch = 19, cex = 1) # plot FE
  points(func_ind_sp_nm$details$asb_G_coord[[i]][1], func_ind_sp_nm$details$asb_G_coord[[i]][2], col = "black", bg= color_fig[i,1], pch = 23, cex= 1.2) # plot center of gravity
  
} # Plot the functioanl sapces of the PCEEZ and the Global pool


# Create a FE coordinate matrix in the functional space from species dissimilarities
plot(fe_faxes[,1],fe_faxes[,2], type="n", xlab = "PCoA 1", ylab = "PCoA 2", xlim = c(-0.55, 0.55),
     ylim = c(-0.55, 0.45), main = "Global Functional Space") # plot Global Functional space
polygon(glob_ch, col="grey95", border=NA) # Global Functional Space

fe_faxes <- as.data.frame(fe_faxes)
fe_faxes$n_sp <- NA
for (i in rownames(fe_faxes)) {
  fe_faxes[i,"n_sp"] <- as.numeric(func_ent_nm$fe_nb_sp[i])
}

fe_faxes$n_sp <- (fe_faxes$n_sp/88)*100
fe_faxes$col <- NA
fe_faxes[c(which(fe_faxes$n_sp != min(fe_faxes$n_sp))),"col"] = "#3A5FCD70"
fe_faxes[c(which(fe_faxes$n_sp == min(fe_faxes$n_sp))),"col"] = "#CD262670"

points(fe_faxes[,1], fe_faxes[,2], col = "#3A5FCD70", pch = 19, cex = c(fe_faxes$n_sp)) # plot FE
legend(-0.5, 0.3, legend = c(9,6,5, 4, 3, 2, 1), 
       col = "#3A5FCD70", pch = 19, bty = "n", 
       pt.cex = c(10.227273, 6.818182, 5.681818, 4.545455, 3.409091, 2.272727, 1.136364))


for (i in rownames(fe_comm_biog)) {
  points(func_ind_sp_nm$details$asb_G_coord[[i]][1], func_ind_sp_nm$details$asb_G_coord[[i]][2], col = "black", bg= color_fig[i,1], pch = 23, cex= 1.2) # plot center of gravity of each province
}


##
glob_tr2 <- tri.mesh(fe_faxes[,1],fe_faxes[,3]) # trigonometry Global Functional Space 1st and 2nd PCoA axes
glob_ch2 <- convex.hull(glob_tr) # convex hull Global Functional Space

# Plot province and Global functional spaces (1st & 2nd PCoA axes)
par(mfrow = c(1,2), mai=c(rep(0.2, 4))) # set in between plots sizes

for (i in rownames(fe_comm_biog)) {
  
  plot(fe_faxes[,1],fe_faxes[,3], type="n",
       xlab = NA, ylab = NA, xaxt = "n", yaxt = "n",
       xlim = c(-0.55, 0.55), ylim = c(-0.45, 0.3), main = i) # plot 
  axis(1, labels = FALSE)
  axis(2, labels = FALSE)
  
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey95") # grey background
  polygon(glob_ch2, col="white", border=NA) # Global Functional Space
  
  fe_prov2 <- colnames(fe_comm_biog)[which(fe_comm_biog[i,] > 0)] # select FE present at each province
  fe_coord2 <- fe_faxes[rownames(fe_faxes) %in% fe_prov2, c(1,2)] # select coordinates for present FE
  tr2 <-tri.mesh(fe_coord2[,1],fe_coord2[,3]) # compute trigonometry
  ch2 <- convex.hull(tr) # convex hull
  
  polygon(ch, col = color_fig[i,2], border = NA) # plot convex hull 
  points(fe_coord[,1], fe_coord[,3], col = color_fig[i,1], pch = 19, cex = 1) # plot FE
  points(func_ind_sp_nm$details$asb_G_coord[[i]][1], func_ind_sp_nm$details$asb_G_coord[[i]][3], col = "black", bg= color_fig[i,1], pch = 23, cex= 1.2) # plot center of gravity
  
} # Plot the functioanl sapces of the PCEEZ and the Global pool


# Create a FE coordinate matrix in the functional space from species dissimilarities
plot(fe_faxes[,1],fe_faxes[,3], type="n", xlab = "PCoA 1", ylab = "PCoA 2", xlim = c(-0.55, 0.55),
     ylim = c(-0.45, 0.3), main = "Global Functional Space") # plot Global Functional space
polygon(glob_ch, col="grey95", border=NA) # Global Functional Space


points(fe_faxes[,1], fe_faxes[,3], col = "#3A5FCD70", pch = 19, cex = c(fe_faxes$n_sp)) # plot FE
legend(-0.5, 0.3, legend = c(9,6,5, 4, 3, 2, 1), 
       col = "#3A5FCD70", pch = 19, bty = "n", 
       pt.cex = c(10.227273, 6.818182, 5.681818, 4.545455, 3.409091, 2.272727, 1.136364))


for (i in rownames(fe_comm_biog)) {
  points(func_ind_sp_nm$details$asb_G_coord[[i]][1], func_ind_sp_nm$details$asb_G_coord[[i]][3], col = "black", bg= color_fig[i,1], pch = 23, cex= 1.2) # plot center of gravity of each province
}


##
glob_tr <- tri.mesh(fe_faxes[,1],fe_faxes[,2]) # trigonometry Global Functional Space 1st and 2nd PCoA axes
glob_ch <- convex.hull(glob_tr) # convex hull Global Functional Space

# Plot province and Global functional spaces (1st & 2nd PCoA axes)
par(mfrow = c(1,2), mai=c(rep(0.2, 4))) # set in between plots sizes

for (i in rownames(fe_comm_biog)) {
  
  plot(fe_faxes[,1],fe_faxes[,2], type="n",
       xlab = NA, ylab = NA, xaxt = "n", yaxt = "n",
       xlim = c(-0.55, 0.55), ylim = c(-0.55, 0.45), main = i) # plot 
  axis(1, labels = FALSE)
  axis(2, labels = FALSE)
  
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey95") # grey background
  polygon(glob_ch, col="white", border=NA) # Global Functional Space
  
  fe_prov <- colnames(fe_comm_biog)[which(fe_comm_biog[i,] > 0)] # select FE present at each province
  fe_coord <- fe_faxes[rownames(fe_faxes) %in% fe_prov, c(1,2)] # select coordinates for present FE
  tr <-tri.mesh(fe_coord[,1],fe_coord[,2]) # compute trigonometry
  ch <- convex.hull(tr) # convex hull
  
  polygon(ch, col = color_fig[i,2], border = NA) # plot convex hull 
  points(fe_coord[,1], fe_coord[,2], col = color_fig[i,1], pch = 19, cex = 1) # plot FE
  points(func_ind_sp_nm$details$asb_G_coord[[i]][1], func_ind_sp_nm$details$asb_G_coord[[i]][2], col = "black", bg= color_fig[i,1], pch = 23, cex= 1.2) # plot center of gravity
  
} # Plot the functioanl sapces of the PCEEZ and the Global pool


# Create a FE coordinate matrix in the functional space from species dissimilarities
plot(fe_faxes[,1],fe_faxes[,2], type="n", xlab = "PCoA 1", ylab = "PCoA 2", xlim = c(-0.55, 0.55),
     ylim = c(-0.55, 0.45), main = "Global Functional Space") # plot Global Functional space
polygon(glob_ch, col="grey95", border=NA) # Global Functional Space

fe_faxes <- as.data.frame(fe_faxes)
fe_faxes$n_sp <- NA
for (i in rownames(fe_faxes)) {
  fe_faxes[i,"n_sp"] <- as.numeric(func_ent_nm$fe_nb_sp[i])
}

fe_faxes$n_sp <- (fe_faxes$n_sp/88)*100
fe_faxes$col <- NA
fe_faxes[c(which(fe_faxes$n_sp != min(fe_faxes$n_sp))),"col"] = "#3A5FCD70"
fe_faxes[c(which(fe_faxes$n_sp == min(fe_faxes$n_sp))),"col"] = "#CD262670"

points(fe_faxes[,1], fe_faxes[,2], col = "#3A5FCD70", pch = 19, cex = c(fe_faxes$n_sp)) # plot FE
legend(-0.5, 0.3, legend = c(9,6,5, 4, 3, 2, 1), 
       col = "#3A5FCD70", pch = 19, bty = "n", 
       pt.cex = c(10.227273, 6.818182, 5.681818, 4.545455, 3.409091, 2.272727, 1.136364))


for (i in rownames(fe_comm_biog)) {
  points(func_ind_sp_nm$details$asb_G_coord[[i]][1], func_ind_sp_nm$details$asb_G_coord[[i]][2], col = "black", bg= color_fig[i,1], pch = 23, cex= 1.2) # plot center of gravity of each province
}


# Kernel density of species distribution along the Functional space axes PCoA 1 & 2
###################################################################################

par(mfrow = c(2,4), mai=c(rep(0.2, 4))) # set in between plots sizes
for (i in rownames(fe_comm_biog)) {
  d <- density(func_ind_sp_nm$details$asb_sp_faxes_coord[[i]][,1], adjust = 0.8, from = -0.5, to = 0.5)
  plot(d, type = "n", ylab = "", xlab= "", main = "", xlim = c(-0.5, 0.5), ylim = c(0, 8))
  polygon(d, col = color_fig[i,2], border = color_fig[i,1])
} # plot Kernel density of species along PCoA 1

par(mfrow = c(2,4), mai=c(rep(0.2, 4))) # set in between plots sizes
for (i in rownames(fe_comm_biog)) {
  d <- density(func_ind_sp_nm$details$asb_sp_faxes_coord[[i]][,2], adjust = 0.8, from = -0.5, to = 0.5)
  plot(d, type = "n", ylab = "", xlab= "", main = "", xlim = c(-0.5, 0.5), ylim = c(0, 8))
  polygon(d, col = color_fig[i,2], border = color_fig[i,1])
} # plot Kernel density of species along PCoA 2


# PCoA 1

density1<-plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[1], xlim = c(-0.55, 0.55), ylim = c(0, 10))
for (i in levels(sp_traits_nm$max_length_m)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$max_length_m == i)), 1],
               adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$max_length[,i],"70"), border = legend_info$max_length[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}


plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[2], xlim = c(-0.55, 0.55), ylim = c(0, 10))
for (i in levels(sp_traits_nm$max_mass_kg)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$max_mass_kg == i)), 1],
               adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$max_mass[,i],"70"), border = legend_info$max_mass[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[3], xlim = c(-0.55, 0.55), ylim = c(0, 10))
for (i in levels(sp_traits_nm$max_mass_max_length_ratio)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$max_mass_max_length_ratio == i)), 1],
               adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$ratio[,i],"70"), border = legend_info$ratio[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[4], xlim = c(-0.55, 0.55), ylim = c(0, 15))
for (i in levels(sp_traits_nm$dentition)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$dentition == i)), 1],
               adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$dentition[,i],"70"), border = legend_info$dentition[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[5], xlim = c(-0.55, 0.55), ylim = c(0, 10))
for (i in levels(sp_traits_nm$max_diving_depth)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$max_diving_depth == i)), 1],
               adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$max_diving[,i],"70"), border = legend_info$max_diving[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[6], xlim = c(-0.55, 0.55), ylim = c(0, 6))
for (i in levels(sp_traits_nm$average_group_size)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$average_group_size == i)), 1],
               adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$group_size[,i],"70"), border = legend_info$group_size[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[7], xlim = c(-0.55, 0.55), ylim = c(0, 10))
for (i in levels(sp_traits_nm$prey_choice)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$prey_choice == i)), 1],
               bw=0.08, adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$prey_choice[,i],"70"), border = legend_info$prey_choice[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

#PCoA2

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[1], xlim = c(-0.55, 0.45), ylim = c(0, 10))
for (i in levels(sp_traits_nm$max_length_m)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$max_length_m == i)), 2],
               bw=0.05, adjust = 1, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$max_length[,i],"70"), border = legend_info$max_length[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[2], xlim = c(-0.55, 0.45), ylim = c(0, 12))
for (i in levels(sp_traits_nm$max_mass_kg)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$max_mass_kg == i)), 2],
               bw=0.05, adjust = 1, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$max_mass[,i],"70"), border = legend_info$max_mass[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[3], xlim = c(-0.55, 0.45), ylim = c(0, 10))
for (i in levels(sp_traits_nm$max_mass_max_length_ratio)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$max_mass_max_length_ratio == i)), 2],
               bw=0.05, adjust = 1, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$ratio[,i],"70"), border = legend_info$ratio[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[4], xlim = c(-0.55, 0.45), ylim = c(0, 15))
for (i in levels(sp_traits_nm$dentition)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$dentition == i)), 2],
               bw=0.05, adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$dentition[,i],"70"), border = legend_info$dentition[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[5], xlim = c(-0.55, 0.45), ylim = c(0, 10))
for (i in levels(sp_traits_nm$max_diving_depth)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$max_diving_depth == i)), 2],
               bw=0.05, adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$max_diving[,i],"70"), border = legend_info$max_diving[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[6], xlim = c(-0.55, 0.45), ylim = c(0, 6))
for (i in levels(sp_traits_nm$average_group_size)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$average_group_size == i)), 2],
               bw=0.05, adjust = 1, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$group_size[,i],"70"), border = legend_info$group_size[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[7], xlim = c(-0.55, 0.45), ylim = c(0, 10))
for (i in levels(sp_traits_nm$prey_choice)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$prey_choice == i)), 2],
              bw=0.05,  adjust = 1, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$prey_choice[,i],"70"), border = legend_info$prey_choice[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

#PCoA3

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[1], xlim = c(-0.6, 0.5), ylim = c(0, 10))
for (i in levels(sp_traits_nm$max_length_m)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$max_length_m == i)), 3],
               adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$max_length[,i],"70"), border = legend_info$max_length[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[2], xlim = c(-0.6, 0.5), ylim = c(0, 10))
for (i in levels(sp_traits_nm$max_mass_kg)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$max_mass_kg == i)), 3],
               adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$max_mass[,i],"70"), border = legend_info$max_mass[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[3], xlim = c(-0.6, 0.5), ylim = c(0, 10))
for (i in levels(sp_traits_nm$max_mass_max_length_ratio)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$max_mass_max_length_ratio == i)), 3],
               adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$ratio[,i],"70"), border = legend_info$ratio[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[4], xlim = c(-0.6, 0.5), ylim = c(0, 15))
for (i in levels(sp_traits_nm$dentition)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$dentition == i)), 3],
               adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$dentition[,i],"70"), border = legend_info$dentition[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[5], xlim = c(-0.6, 0.5), ylim = c(0, 10))
for (i in levels(sp_traits_nm$max_diving_depth)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$max_diving_depth == i)), 3],
               adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$max_diving[,i],"70"), border = legend_info$max_diving[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[6], xlim = c(-0.6, 0.5), ylim = c(0, 10))
for (i in levels(sp_traits_nm$average_group_size)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$average_group_size == i)), 3],
               bw=0.05,adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$group_size[,i],"70"), border = legend_info$group_size[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

plot(d, type = "n", ylab = "", xlab= "", main = colnames(sp_traits_nm)[7], xlim = c(-0.6, 0.5), ylim = c(0, 10))
for (i in levels(sp_traits_nm$prey_choice)) {
  d <- density(func_ind_sp_nm$details$sp_faxes_coord[c(which(sp_traits_nm$prey_choice == i)), 3],
               bw=0.05,  adjust = 0.8, from = -1, to = 1)
  polygon(d, col = paste0(legend_info$prey_choice[,i],"70"), border = legend_info$prey_choice[,i]) # plot Kernel density of mobility trait along PCoA 1 of the Global Functional Space
}

####END OF Code####
