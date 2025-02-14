####Family Breakdown Figure####
setwd()
library(patchwork)
library(tidyverse)

Family_data<-read.csv("Family_proportions_data.csv")
Location<-c(rep("Global", 11), rep("PCEEZ", 11))
Family_data$Location<-Location


ggplot(Family_data, aes(x=Location, y=Family_Proportions, fill=Family))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#C15CCB", "#00868B", "#FF6A00", "#FFFF00", "#330000", "#0000FF", "#CC0099", "#00CCFF", "#FF99CC", "#FF0000", "#666666"), 
                      labels = c("Balaenidae", "Balaenopteridae", "Delphinidae", "Eschrichtiidae", "Kogiidae", "Monodontidae", "Neobalaenidae", "Phocoenidae", "Physeteridae", "Pontoporiidae", "Ziphiidae"))+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  scale_y_continuous(limits = c(0, 1))+
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x = "Species Pool", y = "Proportion")
  

####PCEEZ Map####
library(rnaturalearthdata)
library(rnaturalearth)
library(sf)
library(ggplot2)

EEZ <- read_sf(dsn = "~/GitHub/Matthew-Tuck-honours/MatthewTuck_Honours/clean_data/eez_iho/eez_iho.shp")
#Flanders Marine Institute (2023). Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), version 12. Available online at http://www.marineregions.org/. https://doi.org/10.14284/632
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

####Unique Species BarPlot####
Location<-c("Global", "Global", "PCEEZ", "PCEEZ")
unique_barplot<-as.data.frame(Location, )
unique_barplot$Location<-(unique_barplot)  
unique_barplot$counts<-c(43,45,15,10)
unique_barplot$proportions<-c(0.488636363636364, 0.511363636363636, 0.6, 0.4)
unique_barplot$category<-c("Unique", "Not Unique", "Unique", "Not Unique")
unique_barplot

ggplot(unique_barplot, aes(x=Location$Location, y=proportions, fill=category))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c( "#FF9900",  "#0000FF"), 
                    labels = c("Unique", "Not Unique"))+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x = "Species Pool", y = "Proportion")


####Unique_Species_Histogram####
unique_sp_null_plot<-as.data.frame(unique_sp_perm_null)
unique_sp_null_plot
unique_sp_mean_plot<-as.data.frame(mean_unique_sp_perm_nm)

ggplot(unique_sp_null_plot, aes(x=unique_sp_perm_null))+
  geom_histogram(bins=15, colour="#000000", fill="#999999")+
  geom_vline(aes(xintercept = mean(unique_sp_perm_null)), color = "red", linetype = "dashed", size = 1)+
  geom_vline(aes(xintercept = c(mean(unique_sp_perm_null)+(1.96*sd(unique_sp_perm_null)))), color = "blue", linetype = "dashed", size = 1)+
  geom_vline(aes(xintercept = c(mean(unique_sp_perm_null)-(1.96*sd(unique_sp_perm_null)))), color = "blue", linetype = "dashed", size = 1)+
  geom_vline(aes(xintercept = 15), size=1)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x="Number of Unique Species", y="Frequency")


####PCoA Plots####
PCoA_1<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(2)]))+
  geom_point(size=3, aes(shape=ind_pceez), alpha=0.5) +
  scale_shape_manual(values=c(17,16), labels=c("Not Found in PCEEZ", "Found in PCEEZ"))+
  labs(x="PCoA 1 (62.76% of variation)", y="PCoA 2 (31.32% of variation)", shape="Presence in PCEEZ")+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))

PCoA_2<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(3)]))+
  geom_point(size=3, aes(shape=ind_pceez), alpha=0.5) +
  scale_shape_manual(values=c(17,16), labels=c("Not Found in PCEEZ", "Found in PCEEZ"))+
  labs(x="PCoA 1 (62.76% of variation)", y="PCoA 3 (8.85% of variation)", shape="Presence in PCEEZ")+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))

PCoA_3<-ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(2)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(3)]))+
  geom_point(size=3, aes(shape=ind_pceez), alpha=0.5) +
  scale_shape_manual(values=c(17,16), labels=c("Not Found in PCEEZ", "Found in PCEEZ"))+ 
  labs(x="PCoA 2 (31.32% of variation)", y="PCoA 3 (8.85% of variation)", shape="Presence in PCEEZ")+
  theme_bw(base_size = 15)+
  theme(text=element_text(size = 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))

ggarrange(PCoA_1, PCoA_2, PCoA_3, nrow=3, ncol=1, labels=c("a","b","c"))

####Plots for Functional Indicies (NEED TO UPDATE AFTER SETTING SEEDS OR ANY OTHER CHANGE TO FUNCTIONAL INDICIES)####
multif_null_plot<-c(multif_perm_nm)
multif_null_plot<-as.data.frame(multif_null_plot)
multif_null_plot<-t(multif_null_plot)

ggplot(multif_null_plot, aes(x=PCEEZ))+
  geom_histogram(bins=50, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=fdis_mean_nm), colour="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(fdis_mean_nm+(1.96*fdis_sd_nm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(fdis_mean_nm-(1.96*fdis_sd_nm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=fdis_mean_PCEEZ), colour="black", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
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
  geom_histogram(bins=11, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=fred_mean_nm), colour="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(fred_mean_nm+(1.96*fred_sd_nm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(fred_mean_nm-(1.96*fred_sd_nm))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=fred_mean_PCEEZ), colour="black", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
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

ggarrange(length_histogram, mass_histogram, ratio_histogram, ncol=3, labels=c("a","b","c"))


####Boxplots for individual morphometric traits####
box1<-ggplot(comparitive_boxplots, aes(x=Location, y=`Max Length (m)`))+
  geom_boxplot(outliers = FALSE)+
  labs(x="Species Pool", y="Maximum Length (m)")+
  theme_bw(base_size = 20)+
  theme(text = element_text(size = 10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))
box1


box2<-ggplot(comparitive_boxplots, aes(x=Location, y=`Max Mass (kg)`))+
  geom_boxplot(outliers = FALSE)+
  labs(x="Species Pool", y="Maximum Mass (kg)")+
  theme_bw(base_size = 20)+
  theme(text = element_text(size = 10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))

box3<-ggplot(comparitive_boxplots, aes(x=Location, y=`Max Mass/Max Length Ratio (kg/m)`))+
  geom_boxplot(outliers = FALSE)+
  labs(x="Species Pool", y="Mass/Length (kg/m)")+
  theme_bw(base_size = 20)+
  theme(text = element_text(size = 10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))


####CWM plots####
CWMs_plots_split<-split(CWMs_plots, f=CWMs_plots$trait_types)
CMW_plot_length<-ggplot(CWMs_plots_split$`Max Length`, aes(fill=Categories, x=location, y=CWMs))+
  geom_bar(position="stack", stat = "identity")+
  labs(x="Species Pool", y="Community Weighted Means")+
  facet_wrap(~trait_types, ncol=1)+
  theme_bw(base_size = 15)+
  theme(text = element_text(size = 10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  theme(legend.position = "right", legend.text=element_text(size=10))+
  scale_fill_brewer(palette = "Set2")


CWM_plot_mass<-CMW_plot_length%+%CWMs_plots_split$`Max Mass`
CWM_plot_ratio<-CMW_plot_length%+%CWMs_plots_split$`Max Mass/Length Ratio`
CWM_plot_dentition<-CMW_plot_length%+%CWMs_plots_split$Dentition
CWM_plot_diving<-CMW_plot_length%+%CWMs_plots_split$`Max Diving Depth`
CWM_plot_group_size<-CMW_plot_length%+%CWMs_plots_split$`Average Group Size`
CWM_plot_prey<-CMW_plot_length%+%CWMs_plots_split$`Prey Choice`

####Completed individual trait plot####
Individual_traits_plot<-(((box1 + labs(tag = 'a'))/(box2 + labs(tag='b'))/ (box3 + labs(tag='c')))|(((CWM_plot_dentition+ labs(tag = 'd')) + (CWM_plot_diving + labs(tag = 'e'))) / ((CWM_plot_group_size + labs(tag = 'f'))+(CWM_plot_prey+ labs(tag = 'g'))))) + plot_layout(widths=c(1,2))
Individual_traits_plot
