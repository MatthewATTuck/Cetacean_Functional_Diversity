
#################Matthew Tuck - Cetacean Functional Diversity OPB###############
#                  (Below this is the important stuff)                         #
################################################################################
#### Set directory ####
setwd()

#### Load R packages ####
library("mFD") # functional trait analyses
library("vegan") # species diversity analyses
library("betapart") # functional B-diversity
library("picante") # null models
library("abind") # arrays
library("tidyverse")

#### CSV loading and setting up data frames####
cetacean_dataset<- read.csv("mtuck_honours_Rdocument.csv") # Import cleaned csv dataset for Cetacean functional diversity

sp_comm_nm <- as.data.frame(cetacean_dataset[,c(19:20)], row.names=cetacean_dataset$species_name)# species occurrence data for PCEEZ without migration
sp_traits_nm <- as.data.frame(cetacean_dataset[,c(3:6, 8:10)], row.names = cetacean_dataset$species_name) #species level trait data without migration
trait_info_nm <- data.frame(trait_name = colnames(sp_traits_nm[,c(1:7)]), trait_type = c("O","O","O","N","O","O","N")) #trait information data frame without migration


sp_traits_nm$max_length_m <- ordered(sp_traits_nm$max_length_m, levels=c("small", "intermediate", "large", "very_large")) 
sp_traits_nm$max_mass_kg <- ordered(sp_traits_nm$max_mass_kg, levels=c("small","intermediate","large","very_large")) 
sp_traits_nm$max_mass_max_length_ratio <- ordered(sp_traits_nm$max_mass_max_length_ratio, levels=c("low", "medium", "high")) 
sp_traits_nm$dentition <- factor(sp_traits_nm$dentition, levels = c("LB", "SB", "FCT", "NFCT")) 
sp_traits_nm$max_diving_depth <- ordered(sp_traits_nm$max_diving_depth, levels=c("epipelagic", "upper_mesopelagic", "lower_mesopelagic", "bathypelagic")) 
sp_traits_nm$average_group_size <- ordered(sp_traits_nm$average_group_size, levels = c("Solitary", "Small", "Medium_Small", "Medium", "Medium_Large", "Large")) 
sp_traits_nm$prey_choice <- factor(sp_traits_nm$prey_choice, levels = c("zooplankton", "non_cephalopod_invertebrates", "cephalopods", "fish", "high_vertebrates")) #setting up species traits as ordinal and nominal variables

####Baleen whale count randimization test####
n_perm=1000
baleen_sp_perm<-as.list(rep(NA, n_perm))

for(i in seq(n_perm)){
  
  baleen_sp_perm_nm<-sample(cetacean_dataset$baleen_pa, size=25, replace=FALSE, prob=NULL)
  baleen_sp_perm_nm<-(sum(baleen_sp_perm_nm==1))
  baleen_sp_perm[[i]]<-baleen_sp_perm_nm}
#for loop which randomly selects 25 cetaceans from the global pool and calculates the 
#number of Mysticete species for each of the 1000 permutations


baleen_sp_perm
baleen_sp_perm_null<-unlist(baleen_sp_perm, use.names = FALSE)
#convert the data into a usable format

mean_baleen_sp_perm_nm<-mean(baleen_sp_perm_null)
mean_baleen_sp_perm_nm
sd_baleen_sp_perm_nm<-sd(baleen_sp_perm_null)
sd_baleen_sp_perm_nm

#SES
ses_baleen_sp_perm_nm<-(7 - mean_baleen_sp_perm_nm)/sd_baleen_sp_perm_nm
ses_baleen_sp_perm_nm
2*pnorm(ses_baleen_sp_perm_nm, mean=0, sd=1, lower.tail=FALSE)

#### Establishing Gower's Functional Distance between species ####
func_dist_sp_nm <- funct.dist(sp_traits_nm, 
                           tr_cat = trait_info_nm, 
                           metric = "gower", 
                           ordinal_var = "metric",
                           weight_type = "equal") # Functional dissimilarity (Gower) distance between species
view(func_dist_sp_nm)

#### Adjustments to species location matrix for later use ####
sp_comm_biog_nm = as.matrix(sp_comm_nm) # create matrix with species classified based on occurrence in the PCEEZ
sp_comm_biog_nm <- decostand(sp_comm_biog_nm, "pa") # presence/abscense transformation (1/0)
sp_comm_biog_nm<-(t(sp_comm_biog_nm)) #transpose matrix for use in the alpha.fd.multidim function


####Unique Species Null Distribution and SES####
library(funrar)
sp_traits_unique_nm<-as.matrix(func_dist_sp_nm) #setting the functional dissimilarity data as a matrix
unique_sp_nm<-funrar(pres_matrix = sp_comm_biog_nm, dist_matrix = sp_traits_unique_nm, rel_abund = FALSE) #calculating the number of unique trait assemblages in the data set
uniqueness_nm<-uniqueness(sp_comm_biog_nm, sp_traits_unique_nm) #calculating "uniqueness" for each species in the data set
uniqueness_nm
number_of_uniques<-(88-colSums(uniqueness_nm==0)) #this means that 43 species are unique (~49% of the total species richness)
number_of_uniques

uniqueness_plus_location_nm<-as.data.frame(uniqueness_nm)
uniqueness_plus_location_nm$Location<-as.factor(cetacean_dataset$PCEEZ)
#Constructs a data frame which includes the uniqueness data of the global pool and presence/absence data for the PCEEZ

uniqueness_PCEEZ<-uniqueness_plus_location_nm%>%
  filter(Location=="1")
#Filtering the data for 25 species in PCEEZ


number_of_uniques_PCEEZ<-(25-colSums(uniqueness_PCEEZ==0))
number_of_uniques_PCEEZ
#Calculating the number of unique species in the PCEEZ species pool (15)

n_perm=1000
unique_sp_perm<-as.list(rep(NA, n_perm))
#setting up randomization test

for(i in seq(n_perm)){
  
  uniqueness_sp_perm_nm<-sample(uniqueness_nm$Ui, size=25, replace=FALSE, prob=NULL)
  unique_sp_perm_nm<-(25-sum(uniqueness_sp_perm_nm==0))
  unique_sp_perm[[i]]<-unique_sp_perm_nm}
#for loop which randomly selects 25 cetaceans from the global pool and calculates the 
#number of unique species for each of the 1000 permutations


unique_sp_perm
unique_sp_perm_null<-unlist(unique_sp_perm, use.names = FALSE)
#convert the data into a usable format

mean_unique_sp_perm_nm<-mean(unique_sp_perm_null)
mean_unique_sp_perm_nm
sd_unique_sp_perm_nm<-sd(unique_sp_perm_null)
sd_unique_sp_perm_nm
#calculate the mean and standard deviation of the uniqueness null distribution


ses_unique_sp_perm_nm<-(number_of_uniques_PCEEZ - mean_unique_sp_perm_nm)/sd_unique_sp_perm_nm
ses_unique_sp_perm_nm
2*pnorm(ses_unique_sp_perm_nm, mean=0, sd=1, lower.tail=FALSE)
#calculating the SES of the PCEEZ region and p-value

#### PCoA coordinate construction ####
func_qual_nm <- quality.fspaces(func_dist_sp_nm, 
                             deviation_weighting = c("absolute", "squared"),
                             maxdim_pcoa = 9,
                             fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

func_qual_nm$details_fspaces

apply(func_qual_nm$quality_fspaces, 2, which.min) # best number of dimesions calculated with MAD index (3D) (Maire et al. 2015)
sp_faxes_nm <- func_qual_nm$details_fspaces$sp_pc_coord # species coordinates in the functional space

sp_tr_faxes <-traits.faxes.cor(
  sp_tr          = sp_traits_nm, 
  sp_faxes_coord = sp_faxes_nm[ , c("PC1", "PC2", "PC3")], 
  plot           = TRUE)
sp_tr_faxes
#This creates a data frame displaying which traits have a significant effect on data position for each PCoA axis.
#It also creates a plot which indicates highlights traits significantly impacting position in blue. 


####Calculation of PCEEZ Functional Dispersion (and other indices)####
func_ind_sp_nm <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_nm[,c(1:3)],
                                 asb_sp_w = sp_comm_biog_nm,
                                 scaling = TRUE,
                                 details_returned = T,
                                 verbose = T) # mutidimensinal functional diversity indices in a 3D functional space

func_ind_fdis_nm <- func_ind_sp_nm$functional_diversity_indices$fdis # functional dispersion of PCEEZ region and global pool without migration trait 
func_ind_fdis_nm

#### Unique Trait Combinations ####
func_ent_nm <- sp.to.fe(sp_traits_nm, trait_info_nm, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe_nm <- alpha.fd.fe(sp_comm_biog_nm, sp_to_fe = func_ent_nm) # compute UTC derived indices
func_ind_fe_nm$asb_fdfe # UTC Functional indexes
  
view(func_ind_fe_nm$asb_fdfe)
####Functional Indicies NULL Distribution Setup ####
n_perm = 1000 # number of permutations
multif_perm_nm <- as.list(rep(NA, n_perm)) # list to save permutations results of multidimensional functional alpha-diversity indices (FDis)
fe_perm_nm <- as.list(rep(NA, n_perm)) # list to save permutations results of FE indices (FRed)


for(i in seq(n_perm)){
  
  sp_comm_n_nm = randomizeMatrix(sp_comm_biog_nm, null.model = "richness") # Randomizations maintaining species richness in the PCEEZ
 
   multif_n_nm <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_nm[,c(1:3)],
                                asb_sp_w = sp_comm_n_nm,
                                ind_vect= "fdis",
                                scaling = T,
                                details_returned = T,
                                verbose = F) #  Null Distribution for FDis
  
  multif_perm_nm[[i]] <- multif_n_nm$functional_diversity_indices["fdis"] # store null FDis results 
  
  func_ind_fe_n_nm <- alpha.fd.fe(sp_comm_n_nm, sp_to_fe = func_ent_nm) # Null FE indices
  fe_perm_nm[[i]] <- func_ind_fe_n_nm$asb_fdfe[,c("fred", "fvuln")]} # store null FRed indices result


#### Wrap Null Distribution Results ####
multif_perm_a_nm <- array(NA, dim = c(2,1,0)) # create an array with 1 col & 2 rows for null FDis results
fe_perm_a_nm <- array(NA, dim = c(2,2,0)) # create an array with 1 col & 2 rows for FRed indices results

for (i in c(1:n_perm)) { 
  multif_perm_a_nm <- abind(multif_perm_a_nm, multif_perm_nm[[i]]) # FDis
  fe_perm_a_nm <- abind(fe_perm_a_nm, fe_perm_nm[[i]])} #FRed


multif_perm_mean_nm <- matrix(NA, ncol = 1, nrow= 2, dimnames = list(rownames(multif_perm_a_nm), colnames(multif_perm_a_nm))) # matrix to save mean values of FDis permutations
multif_perm_sd_nm <- matrix(NA, ncol = 1, nrow= 2, dimnames = list(rownames(multif_perm_a_nm), colnames(multif_perm_a_nm))) # matrix to save sd values of FDis permutations 
for (i in 1:2) {
  for (j in "fdis") {
    multif_perm_mean_nm[i,j] = mean(multif_perm_a_nm[i,j,])
    multif_perm_sd_nm[i,j] = sd(multif_perm_a_nm[i,j,])
  }
} # mean and sd values of FDis permutations

fe_perm_mean_nm <- matrix(NA, ncol = 2, nrow= 2, dimnames = list(rownames(fe_perm_a_nm), colnames(fe_perm_a_nm))) # matrix to save mean values of FRed and FVuln indices permutations
fe_perm_sd_nm <- matrix(NA, ncol = 2, nrow= 2, dimnames = list(rownames(fe_perm_a_nm), colnames(fe_perm_a_nm))) # matrix to save sd values of FRed and FVuln indices permutations 
for (i in 1:2) {
  for (j in 1:2) {
    fe_perm_mean_nm[i,j] = mean(fe_perm_a_nm[i,j,])
    fe_perm_sd_nm[i,j] = sd(fe_perm_a_nm[i,j,]) }
  
}#mean and sd for the FRed indicies

fdis_mean_nm<-as.data.frame(multif_perm_mean_nm)
fdis_mean_nm<-fdis_mean_nm[1,1]
fdis_sd_nm<-as.data.frame(multif_perm_sd_nm)
fdis_sd_nm<-fdis_sd_nm[1,1]
fdis_mean_PCEEZ<-func_ind_fdis_nm[1]
fdis_mean_nm
fdis_sd_nm
fdis_mean_PCEEZ


fred_mean_nm<-as.data.frame(fe_perm_mean_nm)
fred_mean_nm<-fred_mean_nm[1,1]
fred_sd_nm<-as.data.frame(fe_perm_sd_nm)
fred_sd_nm<-fred_sd_nm[1,1]
fred_mean_PCEEZ<-as.data.frame(func_ind_fe_nm$asb_fdfe)
fred_mean_PCEEZ<-fred_mean_PCEEZ[1,3]
fred_mean_nm
fred_sd_nm
fred_mean_PCEEZ
#making it so these values will update automatically in the ggplot code later after analysis is rerun


#### SES and p-values ####
ses_multif_nm <- (func_ind_sp_nm$functional_diversity_indices[,c("fdis")] - multif_perm_mean_nm) / multif_perm_sd_nm # SES functional multidimensional alpha diversity indices
ses_multif_nm
2*pnorm(ses_multif_nm, mean=0, sd=1, lower.tail = FALSE)

ses_fe_nm <- (func_ind_fe_nm$asb_fdfe[,c("fred")] - fe_perm_mean_nm) / fe_perm_sd_nm # SES FE indices FRed
ses_fe_nm
2*pnorm(ses_fe_nm, mean=0, sd=1, lower.tail = TRUE)



####Individual trait analyses setup####
cetacean_ind_traits<-read.csv("cleaned_raw_data_CFD.csv")

cetacean_ind_traits$max_length_m <- cetacean_ind_traits$max_length_m
cetacean_ind_traits$max_mass_kg <- cetacean_ind_traits$max_mass_kg
cetacean_ind_traits$max_mass_max_length_ratio <- cetacean_ind_traits$max_mass_max_length_ratio
ind_dentition <- cetacean_ind_traits$dentition
ind_migratory_behaviour <- cetacean_ind_traits$migratory_behaviour 
ind_max_diving_depth <- cetacean_ind_traits$max_diving_depth
ind_average_group_size <- cetacean_ind_traits$average_group_size 
ind_prey_choice <- cetacean_ind_traits$prey_choice 
ind_pceez<- as.factor(cetacean_ind_traits$Present_in_PCEEZ_Region) #setting up species traits for individual analyses
cetacean_ind_traits


cetacean_ind_traits_pceez <- cetacean_ind_traits %>%
  filter(ind_pceez=="Yes")
#Filtering for the PCEEZ species pool

morphometric_means<-cetacean_ind_traits_pceez %>%
  summarize(mean_length=mean(max_length_m), mean_mass=mean(max_mass_kg), mean_ratio=mean(max_mass_max_length_ratio))
morphometric_means
#means of continuous traits from OPb


morphometric_sd<-cetacean_ind_traits_pceez%>%
  summarise(sd_length=sd(max_length_m), sd_mass=sd(max_mass_kg), sd_ratio=sd(max_mass_max_length_ratio))
morphometric_sd
#sd of continuous traits from OPB

morphometric_se<-cetacean_ind_traits_pceez%>%
  summarise(se_length=(sd(max_length_m)/5), sd_mass=(sd(max_mass_kg)/5), sd_ratio=(sd(max_mass_max_length_ratio)/5))
morphometric_se

#se of continuous traits from the OPB

####Null Distributions for the continuous traits####
n_perm=1000
max_length_perm <- as.list(rep(NA, n_perm))
max_mass_perm<-as.list(rep(NA, n_perm))
max_ratio_perm<-as.list(rep(NA, n_perm))

max_length_se_perm <- as.list(rep(NA, n_perm))
max_mass_se_perm<-as.list(rep(NA, n_perm))
max_ratio_se_perm<-as.list(rep(NA, n_perm))

for(i in seq(n_perm)){
  
  max_length_perm_n<-sample(cetacean_ind_traits$max_length_m, size=25, replace=FALSE, prob=NULL)
  max_length_perm_means_n<-mean(max_length_perm_n)
  max_length_perm[[i]]<-max_length_perm_means_n
  max_length_perm_se_n<-(sd(max_length_perm_n)/5)
  max_length_se_perm[[i]]<-max_length_perm_se_n
  
  max_mass_perm_n<-sample(cetacean_ind_traits$max_mass_kg, size=25, replace=FALSE, prob=NULL)
  max_mass_perm_means_n<-mean( max_mass_perm_n)
  max_mass_perm[[i]]<-max_mass_perm_means_n
  max_mass_perm_se_n<-(sd(max_mass_perm_n)/5)
  max_mass_se_perm[[i]]<-max_mass_perm_se_n
  
  max_ratio_perm_n<-sample(cetacean_ind_traits$max_mass_max_length_ratio, size=25, replace=FALSE, prob=NULL)
  max_ratio_perm_means_n<-mean( max_ratio_perm_n)
  max_ratio_perm[[i]]<-max_ratio_perm_means_n
  max_ratio_perm_se_n<-(sd(max_ratio_perm_n)/5)
  max_ratio_se_perm[[i]]<-max_ratio_perm_se_n
  }
#For loop calculating the mean of each trait for each random selection of 25 species

####Means and SD's of continuous trait null distributions####
view(max_length_perm)
max_length_perm<-unlist(max_length_perm)
max_length_se_perm<-unlist(max_length_se_perm)


mean_max_length_perm<-mean(max_length_perm)
mean_max_length_perm
sd_max_length_perm<-sd(max_length_perm)
sd_max_length_perm
mean_max_length_se_perm<-mean(max_length_se_perm)
mean_max_length_se_perm
sd_max_length_se_perm<-sd(max_length_se_perm)
sd_max_length_se_perm

view(max_mass_perm) 
max_mass_perm<-unlist(max_mass_perm)
max_mass_se_perm<-unlist(max_mass_se_perm)

mean_max_mass_perm<-mean(max_mass_perm)
mean_max_mass_perm
sd_max_mass_perm<-sd(max_mass_perm)
sd_max_mass_perm
mean_max_mass_se_perm<-mean(max_mass_se_perm)
mean_max_mass_se_perm
sd_max_mass_se_perm<-sd(max_mass_se_perm)
sd_max_mass_se_perm

view(max_ratio_perm) 
max_ratio_perm<-unlist(max_ratio_perm)
max_ratio_se_perm<-unlist(max_ratio_se_perm)

mean_max_ratio_perm<-mean(max_ratio_perm)
mean_max_ratio_perm
sd_max_ratio_perm<-sd(max_ratio_perm)
sd_max_ratio_perm
mean_max_ratio_se_perm<-mean(max_ratio_se_perm)
mean_max_ratio_se_perm
sd_max_ratio_se_perm<-sd(max_ratio_se_perm)
sd_max_ratio_se_perm
#Unlist functions convert the list into a usable format for analyses
#calculated means and standard deviations for null distribution of each variable


####SES and p-values for continuous traits in PCEEZ####
ses_max_length_pceez<-(morphometric_means$mean_length - mean_max_length_perm)/sd_max_length_perm
ses_max_length_pceez
2*pnorm(ses_max_length_pceez, mean=0, sd=1, lower.tail = FALSE)

ses_max_mass_pceez<-(morphometric_means$mean_mass - mean_max_mass_perm)/sd_max_mass_perm
ses_max_mass_pceez
2*pnorm(ses_max_mass_pceez, mean=0, sd=1, lower.tail = FALSE)

ses_max_ratio_pceez<-(morphometric_means$mean_ratio - mean_max_ratio_perm)/sd_max_ratio_perm
ses_max_ratio_pceez
2*pnorm(ses_max_ratio_pceez, mean=0, sd=1, lower.tail = FALSE)

ses_max_length_se_pceez<-(morphometric_se$se_length - mean_max_length_se_perm)/sd_max_length_se_perm
ses_max_length_pceez
2*pnorm(ses_max_length_se_pceez, mean=0, sd=1, lower.tail = FALSE)

ses_max_mass_se_pceez<-(morphometric_se$sd_mass - mean_max_mass_se_perm)/sd_max_mass_se_perm
ses_max_mass_se_pceez
2*pnorm(ses_max_mass_se_pceez, mean=0, sd=1, lower.tail = FALSE)

ses_max_ratio_se_pceez<-(morphometric_se$sd_ratio - mean_max_ratio_se_perm)/sd_max_ratio_se_perm
ses_max_ratio_se_pceez
2*pnorm(ses_max_ratio_se_pceez, mean=0, sd=1, lower.tail = FALSE)



continuous_trait_plots<-data.frame(max_length_perm, max_mass_perm, max_ratio_perm, max_length_se_perm, max_mass_se_perm, max_ratio_se_perm)
continuous_trait_plots
#Setting up plot creation in Figures_and_tables script


####CWM plots set up####
library(gridExtra)
library(FD)
cetacean_CWMs<-data.frame(functcomp(sp_traits_nm, sp_comm_biog_nm, CWM.type = c("all"))) #calculating the community weighted means for the global pool and the PCEEZ. 
cetacean_CWMs<-data.frame(t(cetacean_CWMs))
cetacean_CWMs<-rownames_to_column(cetacean_CWMs, "trait_designation")
CWMs<-c(cetacean_CWMs$PCEEZ, cetacean_CWMs$Global)
CWMs_plots<-data.frame(CWMs)
CWMs_plots$Categories<-factor(c("Small", "Medium", "Large", "Very Large", "Small", "Medium", "Large", "Very Large", "Low", "Medium", "High", "Functional Calcareous Teeth", "Long Baleen", "Non Functional Calcareous Teeth", "Short Baleen", "Epipelagic", "Upper Mesopelagic", "Lower Mesopelagic", "Bathypelagic", "Solitary", "Small", "Small-Medium", "Medium", "Medium-Large", "Large", "Cephalopods", "Fish", "High Vertebrates", "Non-Cephalopod Invertebrates", "Zooplankton"), levels = c("Solitary", "Low", "Small", "Small-Medium", "Medium", "Medium-Large", "Large", "Very Large", "High", "Long Baleen", "Short Baleen", "Functional Calcareous Teeth", "Non Functional Calcareous Teeth", "Epipelagic", "Upper Mesopelagic", "Lower Mesopelagic", "Bathypelagic", "Zooplankton", "Non-Cephalopod Invertebrates", "Cephalopods", "Fish", "High Vertebrates"))
CWMs_plots$title_labels<-c("Max Length", "Max Length", "Max Length", "Max Length", "Max Mass", "Max Mass", "Max Mass", "Max Mass", "Max Mass/Length Ratio", "Max Mass/Length Ratio", "Max Mass/Length Ratio", "Dentition", "Dentition", "Dentition", "Dentition", "Max Diving Depth", "Max Diving Depth", "Max Diving Depth", "Max Diving Depth", "Average Group Size", "Average Group Size", "Average Group Size", "Average Group Size", "Average Group Size", "Average Group Size", "Prey Choice", "Prey Choice", "Prey Choice", "Prey Choice", "Prey Choice", "Max Length", "Max Length", "Max Length", "Max Length", "Max Mass", "Max Mass", "Max Mass", "Max Mass", "Max Mass/Length Ratio", "Max Mass/Length Ratio", "Max Mass/Length Ratio", "Dentition", "Dentition", "Dentition", "Dentition", "Max Diving Depth", "Max Diving Depth", "Max Diving Depth", "Max Diving Depth", "Average Group Size", "Average Group Size", "Average Group Size", "Average Group Size", "Average Group Size", "Average Group Size", "Prey Choice", "Prey Choice", "Prey Choice", "Prey Choice", "Prey Choice")
CWMs_plots$trait_types<-c("Max Length", "Max Length", "Max Length", "Max Length", "Max Mass", "Max Mass", "Max Mass", "Max Mass", "Max Mass/Length Ratio", "Max Mass/Length Ratio", "Max Mass/Length Ratio", "Dentition", "Dentition", "Dentition", "Dentition", "Max Diving Depth", "Max Diving Depth", "Max Diving Depth", "Max Diving Depth", "Average Group Size", "Average Group Size", "Average Group Size", "Average Group Size", "Average Group Size", "Average Group Size", "Prey Choice", "Prey Choice", "Prey Choice", "Prey Choice", "Prey Choice", "Max Length", "Max Length", "Max Length", "Max Length", "Max Mass", "Max Mass", "Max Mass", "Max Mass", "Max Mass/Length Ratio", "Max Mass/Length Ratio", "Max Mass/Length Ratio", "Dentition", "Dentition", "Dentition", "Dentition", "Max Diving Depth", "Max Diving Depth", "Max Diving Depth", "Max Diving Depth", "Average Group Size", "Average Group Size", "Average Group Size", "Average Group Size", "Average Group Size", "Average Group Size", "Prey Choice", "Prey Choice", "Prey Choice", "Prey Choice", "Prey Choice")
CWMs_plots$location<-c("OPB", "OPB","OPB","OPB","OPB","OPB","OPB","OPB", "OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","OPB","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global")
CWMs_plots$trait_labels<-c("A", "B", "C", "D", "A", "B", "C", "D", "A", "B", "C", "A", "B", "C", "D", "A", "B", "C", "D", "A", "B", "C", "D", "E", "F", "A", "B", "C", "D", "E", "A", "B", "C", "D", "A", "B", "C", "D", "A", "B", "C", "A", "B", "C", "D", "A", "B", "C", "D", "A", "B", "C", "D", "E", "F", "A", "B", "C", "D", "E")
CWMs_plots #Creating a dataframe which contains relevant location and community weighted means for every trait assessed.

####Setting up data frames for Fisher's exact tests####

fisher_dentition<-data.frame(LB=(c(0.04*25, 0.04545455*88)), SB=(c(0.24*25, 0.12500000*88)), FCT=(c(0.56*25, 0.59090909*88)), NFCT=(c(0.16*25, 0.23863636*88)))
row.names(fisher_dentition)<-c("PCEEZ", "Global")
fisher_dentition
fisher_test_dentition<-fisher.test(fisher_dentition)
fisher_test_dentition


fisher_diving<-data.frame(epipelagic=(c(0.20*25, 0.38636364*88)), upper_mesopleagic=(c(0.36*25, 0.18181818*88)), lower_mesopelagic=(c(0.16*25, 0.11363636*88)), bathypelagic=(c(0.28*25, 0.31818182*88)))
row.names(fisher_diving)<-c("PCEEZ", "Global")
fisher_diving
fisher_test_diving<-fisher.test(fisher_diving)
fisher_test_diving

fisher_group_size<-data.frame(Solitary=(c(0.16*25, 0.12500000*88)), Small=(c(0.32*25, 0.44318182*88)), Small_Medium=(c(0.12*25, 0.17045455*88)), Medium=(c(0.16*25, 0.09090909*88)), Medium_Large=(c(0.08*25, 0.05681818*88)), Large=(c(0.16*25, 0.11363636*88)))
row.names(fisher_group_size)<-c("PCEEZ", "Global")
fisher_group_size
fisher_test_group_size<-fisher.test(fisher_group_size)
fisher_test_group_size

fisher_prey<-data.frame(Cephalopods=(c(0.52*25, 0.45454545*88)), Fish=(c(0.16*25, 0.37500000*88)), High_Vertebrates=(c(0.04*25, 0.01136364*88)), Non_Cephalopod_Invertebrates=(c(0.04*25, 0.01136364*88)), Zooplankton=(c(0.24*25, 0.14772727*88)))
row.names(fisher_prey)<-c("PCEEZ", "Global")
fisher_prey
fisher_test_prey<-fisher.test(fisher_prey)
fisher_test_prey

fisher_dentition_trial<-data.frame(LB=(c(0.04*100, 0.05*100)), SB=(c(0.24*100, 0.13*100)), FCT=(c(0.56*100, 0.59*100)), NFCT=(c(0.16*100, 0.24*100)))
row.names(fisher_dentition_trial)<-c("PCEEZ", "Global")
fisher_test_dentition_trial<-fisher.test(fisher_dentition_trial)
fisher_test_dentition_trial

fisher_diving_trial<-data.frame(epipelagic=(c(0.20*100, 0.39*100)), upper_mesopleagic=(c(0.36*100, 0.18*100)), lower_mesopelagic=(c(0.16*100, 0.11*100)), bathypelagic=(c(0.28*100, 0.32*100)))
row.names(fisher_diving_trial)<-c("PCEEZ", "Global")
fisher_test_diving_trial<-fisher.test(fisher_diving_trial)
fisher_test_diving_trial

fisher_group_size_trial<-data.frame(Solitary=(c(0.16*100, 0.13*100)), Small=(c(0.32*100, 0.44*100)), Small_Medium=(c(0.12*100, 0.17*100)), Medium=(c(0.16*100, 0.09*100)), Medium_Large=(c(0.08*100, 0.06*100)), Large=(c(0.16*100, 0.11*100)))
row.names(fisher_group_size_trial)<-c("PCEEZ", "Global")
fisher_test_group_size_trial<-fisher.test(fisher_group_size_trial)
fisher_test_group_size_trial

fisher_prey_trial<-data.frame(Cephalopods=(c(0.52*100, 0.45*100)), Fish=(c(0.16*100, 0.38*100)), High_Vertebrates=(c(0.04*100, 0.01*100)), Non_Cephalopod_Invertebrates=(c(0.04*100, 0.01*100)), Zooplankton=(c(0.24*100, 0.15*100)))
row.names(fisher_prey_trial)<-c("PCEEZ", "Global")
fisher_test_prey_trial<-fisher.test(fisher_prey_trial)
fisher_test_prey_trial

####box plots for continuous variables dataframe setup####
comparitive_boxplots_PCEEZ<-data.frame(cetacean_ind_traits_pceez$max_length_m, cetacean_ind_traits_pceez$max_mass_kg, cetacean_ind_traits_pceez$max_mass_max_length_ratio)
comparitive_boxplots_PCEEZ$Location<-rep(c("OPB"), 25)
comparitive_boxplots_PCEEZ<-comparitive_boxplots_PCEEZ%>%
  rename("Max Length (m)" = cetacean_ind_traits_pceez.max_length_m, "Max Mass (kg)"= cetacean_ind_traits_pceez.max_mass_kg, "Max Mass/Max Length Ratio (kg/m)"=cetacean_ind_traits_pceez.max_mass_max_length_ratio)
comparitive_boxplots_PCEEZ


comparitive_boxplots_Global<-data.frame(cetacean_ind_traits$max_length_m, cetacean_ind_traits$max_mass_kg, cetacean_ind_traits$max_mass_max_length_ratio)
comparitive_boxplots_Global$Location<-rep(c("Global"), 88)
comparitive_boxplots_Global<-comparitive_boxplots_Global%>%
  rename("Max Length (m)" = cetacean_ind_traits.max_length_m, "Max Mass (kg)"= cetacean_ind_traits.max_mass_kg , "Max Mass/Max Length Ratio (kg/m)"=cetacean_ind_traits.max_mass_max_length_ratio)
comparitive_boxplots_Global

comparitive_boxplots<-rbind(comparitive_boxplots_PCEEZ, comparitive_boxplots_Global)
comparitive_boxplots


citation()
citation(package="mFD")
citation(package="funrar")
citation(package="FD")

####END OF CODE####