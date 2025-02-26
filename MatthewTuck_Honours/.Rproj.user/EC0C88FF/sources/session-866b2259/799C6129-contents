###############################################################################################
### Script for Cetacean Functional Diversity                                                ###
### Matthew Tuck (mtuck@uvic.ca)                                                            ###
###############################################################################################

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

sp_comm <- as.data.frame(cetacean_dataset[,c(19:20)], row.names=cetacean_dataset$species_name)# species occurrence data for PCEEZ region 
sp_traits <- as.data.frame(cetacean_dataset[,c(3:10)], row.names = cetacean_dataset$species_name)# species level trait data
trait_info <- data.frame(trait_name = colnames(sp_traits[,c(1:8)]), trait_type = c("O","O","O","N","N","O","O","N"))# trait information dataframe


sp_traits$max_length_m <- ordered(sp_traits$max_length_m, levels=c("small", "intermediate", "large", "very_large")) 
sp_traits$max_mass_kg <- ordered(sp_traits$max_mass_kg, levels=c("small","intermediate","large","very_large")) 
sp_traits$max_mass_max_length_ratio <- ordered(sp_traits$max_mass_max_length_ratio, levels=c("low", "medium", "high")) 
sp_traits$dentition <- as.factor(sp_traits$dentition) 
sp_traits$migratory_behaviour <- as.factor(sp_traits$migratory_behaviour) 
sp_traits$max_diving_depth <- ordered(sp_traits$max_diving_depth, levels=c("epipelagic", "upper_mesopelagic", "lower_mesopelagic", "bathypelagic")) 
sp_traits$average_group_size <- ordered(sp_traits$average_group_size, levels = c("Solitary", "Small", "Medium_Small", "Medium", "Medium_Large", "Large")) 
sp_traits$prey_choice <- as.factor(sp_traits$prey_choice) #setting up species traits as ordinal and nominal variables

#### Establishing Gower's Functional Distance between species ####
func_dist_sp <- funct.dist(sp_traits, 
                           tr_cat = trait_info, 
                           metric = "gower", 
                           ordinal_var = "metric",
                           weight_type = "equal") # computes the functional dissimilarity (Gower) distance between species

#### PCoA construction ####
func_qual <- quality.fspaces(func_dist_sp, 
                             deviation_weighting = c("absolute", "squared"),
                             maxdim_pcoa = 9,
                             fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

apply(func_qual$quality_fspaces, 2, which.min) # best number of dimesions calculated with MAD index (4D) (Maire et al. 2015)
sp_faxes <- func_qual$details_fspaces$sp_pc_coord # species coordinates in the functional space
func_qual$details_fspaces

####Rough PCoA plots####
ggplot(func_qual$details_fspaces$sp_pc_coord, aes(x=func_qual$details_fspaces$sp_pc_coord[,c(1)], y= func_qual$details_fspaces$sp_pc_coord[, c(2)], colour=cetacean_dataset$Family))+
  geom_point() +
  labs(x="PCoA 1", y="PCoA 2")+
  theme_light()

ggplot(func_qual$details_fspaces$sp_pc_coord, aes(x=func_qual$details_fspaces$sp_pc_coord[,c(3)], y= func_qual$details_fspaces$sp_pc_coord[, c(4)], colour=cetacean_dataset$Family))+
  geom_point() +
  labs(x="PCoA 3", y="PCoA 4")+
  theme_light()  

#### Adjustments to species location matrix for later use ####
sp_comm_biog = as.matrix(sp_comm) # creating a matrix from the occurrence data 
sp_comm_biog <- decostand(sp_comm_biog, "pa") # presence/abscense transformation (1/0)
sp_comm_biog<-(t(sp_comm_biog)) #transposing the matrix for use in the alpha.fd.multidim function

####Calculation of PCEEZ Functional Dispersion (and other indices)####
func_ind_sp <- alpha.fd.multidim(sp_faxes_coord = sp_faxes[,c(1:4)],
                                 asb_sp_w = sp_comm_biog,
                                 scaling = TRUE,
                                 details_returned = T,
                                 verbose = T) # mutidimensinal functional diversity indices in a 4D functional space



func_ind_fdis <- func_ind_sp$functional_diversity_indices$fdis # functional dispersion calculation for PCEEZ 
func_ind_fdis


#### Unique Trait Combinations ####
func_ent <- sp.to.fe(sp_traits, trait_info, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe <- alpha.fd.fe(sp_comm_biog, sp_to_fe = func_ent) # compute UTC derived indices
func_ind_fe$asb_fdfe # UTC Functional indexes

func_ind_fred <- func_ind_fe$asb_fdfe[,"fred"] # functional redundancy

library(funrar)

sp_traits_unique<-as.matrix(func_dist_sp) #setting the functional dissimilarity data as a matrix
unique_sp<-funrar(pres_matrix = sp_comm_biog, dist_matrix = sp_traits_unique, rel_abund = FALSE) #calculating the number of unique trait assemblages in the data set
uniqueness<-uniqueness(sp_comm_biog, sp_traits_unique) #calculating "uniqueness" for each species in the data set
uniqueness
colSums(uniqueness==0) ## 41 species identified which had uniqueness=0, which means they shared trait assemblage with at least one other species 
#this means that 47 species are unique (>50% of the total species richness)

#### NULL Distribution Setup ####
n_perm = 1000 # number of permutations
multif_perm <- as.list(rep(NA, n_perm)) # list to save permutations results of multidimensional functional alpha-diversity indices (FDis)
fe_perm <- as.list(rep(NA, n_perm)) # list to save permutations results of FE indices (FRed)

for(i in seq(n_perm)){
  sp_comm_n = randomizeMatrix(sp_comm_biog, null.model = "richness") # Randomizations maintaining species richness in the PCEEZ
  
  multif_n <- alpha.fd.multidim(sp_faxes_coord = sp_faxes[,c(1:4)],
                                asb_sp_w = sp_comm_n,
                                ind_vect= "fdis",
                                scaling = T,
                                details_returned = T,
                                verbose = F) #  Null Distribution for FDis
  
  multif_perm[[i]] <- multif_n$functional_diversity_indices["fdis"] # store null FDis results
  
  func_ind_fe_n <- alpha.fd.fe(sp_comm_n, sp_to_fe = func_ent) # Null FE indices
  fe_perm[[i]] <- func_ind_fe_n$asb_fdfe[,c("fred", "fvuln")]} # store null FRed indices result

#### Wrap Null Distribution Results ####
multif_perm_a <- array(NA, dim = c(2,1,0)) # create an array with 1 col & 2 rows for null FDis results
fe_perm_a <- array(NA, dim = c(2,2,0)) # create an array with 1 col & 2 rows for FRed indices results

for (i in c(1:n_perm)) { 
  multif_perm_a <- abind(multif_perm_a, multif_perm[[i]])# FDis
  fe_perm_a <- abind(fe_perm_a, fe_perm[[i]])} #FRed 
  
multif_perm_mean <- matrix(NA, ncol = 1, nrow= 2, dimnames = list(rownames(multif_perm_a), colnames(multif_perm_a))) # matrix to save mean values of FDis permutations
multif_perm_sd <- matrix(NA, ncol = 1, nrow= 2, dimnames = list(rownames(multif_perm_a), colnames(multif_perm_a))) # matrix to save sd values of FDis permutations 
for (i in 1:2) {
  for (j in "fdis") {
    multif_perm_mean[i,j] = mean(multif_perm_a[i,j,])
    multif_perm_sd[i,j] = sd(multif_perm_a[i,j,])
  }
} # mean and sd values of FDis permutations

fe_perm_mean <- matrix(NA, ncol = 2, nrow= 2, dimnames = list(rownames(fe_perm_a), colnames(fe_perm_a))) # matrix to save mean values of FRed and FVuln indices permutations
fe_perm_sd <- matrix(NA, ncol = 2, nrow= 2, dimnames = list(rownames(fe_perm_a), colnames(fe_perm_a))) # matrix to save sd values of FRed and FVuln indices permutations 
for (i in 1:2) {
  for (j in 1:2) {
    fe_perm_mean[i,j] = mean(fe_perm_a[i,j,])
    fe_perm_sd[i,j] = sd(fe_perm_a[i,j,]) }
} # mean and sd values of FE indices permutations  


#### SES ####
ses_multif <- (func_ind_sp$functional_diversity_indices[,c("fdis")] - multif_perm_mean) / multif_perm_sd # SES functional multidimensional alpha diversity indices
ses_multif
pnorm(ses_multif)
ses_fe <- (func_ind_fe$asb_fdfe[,c("fred")] - fe_perm_mean) / fe_perm_sd # SES FE indices FRed
ses_fe
pnorm(ses_fe)


################# Identical code as above minus Migration Trait#################
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
sp_traits_nm$dentition <- as.factor(sp_traits_nm$dentition) 
sp_traits_nm$max_diving_depth <- ordered(sp_traits_nm$max_diving_depth, levels=c("epipelagic", "upper_mesopelagic", "lower_mesopelagic", "bathypelagic")) 
sp_traits_nm$average_group_size <- ordered(sp_traits_nm$average_group_size, levels = c("Solitary", "Small", "Medium_Small", "Medium", "Medium_Large", "Large")) 
sp_traits_nm$prey_choice <- as.factor(sp_traits_nm$prey_choice) #setting up species traits as ordinal and nominal variables


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
2*pnorm(1.333395, mean=0, sd=1, lower.tail=FALSE)
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
#means of continuous traits from PCEEZ


morphometric_sd<-cetacean_ind_traits_pceez%>%
  summarise(sd_length=sd(max_length_m), sd_mass=sd(max_mass_kg), sd_ratio=sd(max_mass_max_length_ratio))
morphometric_sd
#sd's of species traits from the PCEEZ

####Null Distributions for the continuous traits####
n_perm=1000
max_length_perm <- as.list(rep(NA, n_perm))
max_mass_perm<-as.list(rep(NA, n_perm))
max_ratio_perm<-as.list(rep(NA, n_perm))

for(i in seq(n_perm)){
  
  max_length_perm_n<-sample(cetacean_ind_traits$max_length_m, size=25, replace=FALSE, prob=NULL)
  max_length_perm_means_n<-mean( max_length_perm_n)
  max_length_perm[[i]]<-max_length_perm_means_n
  
  
  max_mass_perm_n<-sample(cetacean_ind_traits$max_mass_kg, size=25, replace=FALSE, prob=NULL)
  max_mass_perm_means_n<-mean( max_mass_perm_n)
  max_mass_perm[[i]]<-max_mass_perm_means_n
  
  
  max_ratio_perm_n<-sample(cetacean_ind_traits$max_mass_max_length_ratio, size=25, replace=FALSE, prob=NULL)
  max_ratio_perm_means_n<-mean( max_ratio_perm_n)
  max_ratio_perm[[i]]<-max_ratio_perm_means_n
  }
#For loop calculating the mean of each trait for each random selection of 25 species

####Means and SD's of continuous trait null distributions####
view(max_length_perm)
max_length_perm<-unlist(max_length_perm)


mean_max_length_perm<-mean(max_length_perm)
mean_max_length_perm
sd_max_length_perm<-sd(max_length_perm)
sd_max_length_perm

view(max_mass_perm) 
max_mass_perm<-unlist(max_mass_perm)

mean_max_mass_perm<-mean(max_mass_perm)
mean_max_mass_perm
sd_max_mass_perm<-sd(max_mass_perm)
sd_max_mass_perm

view(max_ratio_perm) 
max_ratio_perm<-unlist(max_ratio_perm)

mean_max_ratio_perm<-mean(max_ratio_perm)
mean_max_ratio_perm
sd_max_ratio_perm<-sd(max_ratio_perm)
sd_max_ratio_perm
#Unlist functions convert the list into a usable format for analyses
#calculated means and standard deviations for null distribution of each variable


####SES and p-values for continous traits in PCEEZ####
ses_max_length_pceez<-(morphometric_means$mean_length - mean_max_length_perm)/sd_max_length_perm
ses_max_length_pceez
2*pnorm(ses_max_length_pceez, mean=0, sd=1, lower.tail = FALSE)

ses_max_mass_pceez<-(morphometric_means$mean_mass - mean_max_mass_perm)/sd_max_mass_perm
ses_max_mass_pceez
2*pnorm(ses_max_mass_pceez, mean=0, sd=1, lower.tail = FALSE)

ses_max_ratio_pceez<-(morphometric_means$mean_ratio - mean_max_ratio_perm)/sd_max_ratio_perm
ses_max_ratio_pceez
2*pnorm(ses_max_ratio_pceez, mean=0, sd=1, lower.tail = FALSE)


continuous_trait_plots<-data.frame(max_length_perm, max_mass_perm, max_ratio_perm)
continuous_trait_plots
#Setting up plot creation in Figures_and_tables script

##Left off here trying to Organize. Return here when you go to do the individual trait box and bar plots


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




























####Setups for CWMs Randomization tests####

sp_comm_CWM<-as.data.frame(sp_comm_biog_nm[c(1),])
sp_comm_CWM<-sp_comm_CWM%>%
  rename("PCEEZ"= `sp_comm_biog_nm[c(1), ]`)
sp_comm_CWM<-as.matrix(sp_comm_CWM)
sp_comm_CWM<-t(sp_comm_CWM)
#Creating a species presence absence matrix with 25 species present out of 88.

sp_traits_CWM_length<-as.data.frame(sp_traits_nm[,c(1)])
sp_traits_CWM_length<-sp_traits_CWM_length%>%
  rename("Maximum Length"=`sp_traits_nm[, c(1)]`)
sp_traits_CWM_length<-as.matrix(sp_traits_CWM_length)
rownames(sp_traits_CWM_length)<-c(cetacean_dataset$species_name)
#species trait matrix for length trait

sp_traits_CWM_mass<-as.data.frame(sp_traits_nm[,c(2)])
sp_traits_CWM_mass<-sp_traits_CWM_mass%>%
  rename("Maximum Mass"=`sp_traits_nm[, c(2)]`)
sp_traits_CWM_mass<-as.matrix(sp_traits_CWM_mass)
rownames(sp_traits_CWM_mass)<-c(cetacean_dataset$species_name)
#species trait matrix for mass trait

sp_traits_CWM_ratio<-as.data.frame(sp_traits_nm[,c(3)])
sp_traits_CWM_ratio<-sp_traits_CWM_ratio%>%
  rename("Maximum Mass/Length Ratio"=`sp_traits_nm[, c(3)]`)
sp_traits_CWM_ratio<-as.matrix(sp_traits_CWM_ratio)
rownames(sp_traits_CWM_ratio)<-c(cetacean_dataset$species_name)
#species trait matrix for ratio trait

sp_traits_CWM_dentition<-as.data.frame(sp_traits_nm[,c(4)])
sp_traits_CWM_dentition<-sp_traits_CWM_dentition%>%
  rename("Dentition"=`sp_traits_nm[, c(4)]`)
sp_traits_CWM_dentition<-as.matrix(sp_traits_CWM_dentition)
rownames(sp_traits_CWM_dentition)<-c(cetacean_dataset$species_name)
#species trait matrix for dentition

sp_traits_CWM_diving<-as.data.frame(sp_traits_nm[,c(5)])
sp_traits_CWM_diving<-sp_traits_CWM_diving%>%
  rename("Maximum Diving Depth"=`sp_traits_nm[, c(5)]`)
sp_traits_CWM_diving<-as.matrix(sp_traits_CWM_diving)
rownames(sp_traits_CWM_diving)<-c(cetacean_dataset$species_name)
#species trait matrix for diving

sp_traits_CWM_group_size<-as.data.frame(sp_traits_nm[,c(6)])
sp_traits_CWM_group_size<-sp_traits_CWM_group_size%>%
  rename("Average Group Size"=`sp_traits_nm[, c(6)]`)
sp_traits_CWM_group_size<-as.matrix(sp_traits_CWM_group_size)
rownames(sp_traits_CWM_group_size)<-c(cetacean_dataset$species_name)
#species trait matrix for group size

sp_traits_CWM_prey<-as.data.frame(sp_traits_nm[,c(7)])
sp_traits_CWM_prey<-sp_traits_CWM_prey%>%
  rename("Prey Choice"=`sp_traits_nm[, c(7)]`)
sp_traits_CWM_prey<-as.matrix(sp_traits_CWM_prey)
rownames(sp_traits_CWM_prey)<-c(cetacean_dataset$species_name)
#species trait matrix for prey choice

####Null distribution set-up for each trait####
n_perm = 1000 # number of permutations
sp_traits_CWM_length_perm<- as.list(rep(NA, n_perm)) # list to save permutations results of null distribution
sp_traits_CWM_mass_perm<- as.list(rep(NA, n_perm)) # list to save permutations results of null distribution
sp_traits_CWM_ratio_perm<- as.list(rep(NA, n_perm)) # list to save permutations results of null distribution
sp_traits_CWM_dentition_perm<- as.list(rep(NA, n_perm)) # list to save permutations results of null distribution
sp_traits_CWM_diving_perm<- as.list(rep(NA, n_perm)) # list to save permutations results of null distribution
sp_traits_CWM_group_size_perm<- as.list(rep(NA, n_perm)) # list to save permutations results of null distribution
sp_traits_CWM_prey_perm<- as.list(rep(NA, n_perm)) # list to save permutations results of null distribution


for(i in seq(n_perm)){
sp_comm_CWM_n = randomizeMatrix(sp_comm_CWM, null.model = "richness") # Randomizations maintaining same species richness as in the PCEEZ

sp_traits_CWM_length_n<-functcomp(sp_traits_CWM_length, sp_comm_CWM_n, CWM.type = c("all"))
sp_traits_CWM_length_perm[[i]] <- sp_traits_CWM_length_n 

sp_traits_CWM_mass_n<-functcomp(sp_traits_CWM_mass, sp_comm_CWM_n, CWM.type = c("all"))
sp_traits_CWM_mass_perm[[i]] <- sp_traits_CWM_mass_n

sp_traits_CWM_ratio_n<-functcomp(sp_traits_CWM_ratio, sp_comm_CWM_n, CWM.type = c("all"))
sp_traits_CWM_ratio_perm[[i]] <- sp_traits_CWM_ratio_n

sp_traits_CWM_dentition_n<-functcomp(sp_traits_CWM_dentition, sp_comm_CWM_n, CWM.type = c("all"))
sp_traits_CWM_dentition_perm[[i]] <- sp_traits_CWM_dentition_n

sp_traits_CWM_diving_n<-functcomp(sp_traits_CWM_diving, sp_comm_CWM_n, CWM.type = c("all"))
sp_traits_CWM_diving_perm[[i]] <- sp_traits_CWM_diving_n

sp_traits_CWM_group_size_n<-functcomp(sp_traits_CWM_group_size, sp_comm_CWM_n, CWM.type = c("all"))
sp_traits_CWM_group_size_perm[[i]] <- sp_traits_CWM_group_size_n

sp_traits_CWM_prey_n<-functcomp(sp_traits_CWM_prey, sp_comm_CWM_n, CWM.type = c("all"))
sp_traits_CWM_prey_perm[[i]] <- sp_traits_CWM_prey_n
}

#### Wrap Null Distribution Results for CWMs####
sp_traits_CWM_length_perm_a <- array(NA, dim = c(1,4,0)) # create an array with 4 cols & 1 row for null length CWM results
sp_traits_CWM_mass_perm_a <- array(NA, dim = c(1,4,0))
sp_traits_CWM_ratio_perm_a <- array(NA, dim = c(1,3,0))
sp_traits_CWM_dentition_perm_a <- array(NA, dim = c(1,4,0))
sp_traits_CWM_diving_perm_a <- array(NA, dim = c(1,4,0))
sp_traits_CWM_group_size_perm_a <- array(NA, dim = c(1,6,0))
sp_traits_CWM_prey_perm_a <- array(NA, dim = c(1,5,0))

for (i in c(1:n_perm)) { 
  sp_traits_CWM_length_perm_a <- abind(sp_traits_CWM_length_perm_a, sp_traits_CWM_length_perm[[i]])# length CWMs
  sp_traits_CWM_mass_perm_a <- abind(sp_traits_CWM_mass_perm_a, sp_traits_CWM_mass_perm[[i]])
  sp_traits_CWM_ratio_perm_a <- abind(sp_traits_CWM_ratio_perm_a, sp_traits_CWM_ratio_perm[[i]])
  sp_traits_CWM_dentition_perm_a <- abind(sp_traits_CWM_dentition_perm_a, sp_traits_CWM_dentition_perm[[i]])
  sp_traits_CWM_diving_perm_a <- abind(sp_traits_CWM_diving_perm_a, sp_traits_CWM_diving_perm[[i]])
  sp_traits_CWM_group_size_perm_a <- abind(sp_traits_CWM_group_size_perm_a, sp_traits_CWM_group_size_perm[[i]])
  sp_traits_CWM_prey_perm_a <- abind(sp_traits_CWM_prey_perm_a, sp_traits_CWM_prey_perm[[i]])
}

sp_traits_CWM_length_perm_mean <- matrix(NA, ncol = 4, nrow= 1, dimnames = list(rownames(sp_traits_CWM_length_perm_a), colnames(sp_traits_CWM_length_perm_a))) # matrix to save mean values of length CWM permutations
sp_traits_CWM_length_perm_sd <- matrix(NA, ncol = 4, nrow= 1, dimnames = list(rownames(sp_traits_CWM_length_perm_a), colnames(sp_traits_CWM_length_perm_a))) # matrix to save sd values of length CWM permutations 

sp_traits_CWM_mass_perm_mean <- matrix(NA, ncol = 4, nrow= 1, dimnames = list(rownames(sp_traits_CWM_mass_perm_a), colnames(sp_traits_CWM_mass_perm_a)))
sp_traits_CWM_mass_perm_sd <- matrix(NA, ncol = 4, nrow= 1, dimnames = list(rownames(sp_traits_CWM_mass_perm_a), colnames(sp_traits_CWM_mass_perm_a)))

sp_traits_CWM_ratio_perm_mean <- matrix(NA, ncol = 3, nrow= 1, dimnames = list(rownames(sp_traits_CWM_ratio_perm_a), colnames(sp_traits_CWM_ratio_perm_a)))
sp_traits_CWM_ratio_perm_sd <- matrix(NA, ncol = 3, nrow= 1, dimnames = list(rownames(sp_traits_CWM_ratio_perm_a), colnames(sp_traits_CWM_ratio_perm_a)))

sp_traits_CWM_dentition_perm_mean <- matrix(NA, ncol = 4, nrow= 1, dimnames = list(rownames(sp_traits_CWM_dentition_perm_a), colnames(sp_traits_CWM_dentition_perm_a)))
sp_traits_CWM_dentition_perm_sd <- matrix(NA, ncol = 4, nrow= 1, dimnames = list(rownames(sp_traits_CWM_dentition_perm_a), colnames(sp_traits_CWM_dentition_perm_a)))

sp_traits_CWM_diving_perm_mean <- matrix(NA, ncol = 4, nrow= 1, dimnames = list(rownames(sp_traits_CWM_diving_perm_a), colnames(sp_traits_CWM_diving_perm_a)))
sp_traits_CWM_diving_perm_sd <- matrix(NA, ncol = 4, nrow= 1, dimnames = list(rownames(sp_traits_CWM_diving_perm_a), colnames(sp_traits_CWM_diving_perm_a)))

sp_traits_CWM_group_size_perm_mean <- matrix(NA, ncol = 6, nrow= 1, dimnames = list(rownames(sp_traits_CWM_group_size_perm_a), colnames(sp_traits_CWM_group_size_perm_a)))
sp_traits_CWM_group_size_perm_sd <- matrix(NA, ncol = 6, nrow= 1, dimnames = list(rownames(sp_traits_CWM_group_size_perm_a), colnames(sp_traits_CWM_group_size_perm_a)))

sp_traits_CWM_prey_perm_mean <- matrix(NA, ncol = 5, nrow= 1, dimnames = list(rownames(sp_traits_CWM_prey_perm_a), colnames(sp_traits_CWM_prey_perm_a)))
sp_traits_CWM_prey_perm_sd <- matrix(NA, ncol = 5, nrow= 1, dimnames = list(rownames(sp_traits_CWM_prey_perm_a), colnames(sp_traits_CWM_prey_perm_a)))


for (i in 1) {
  for (j in 1:4) {
    sp_traits_CWM_length_perm_mean[i,j] = mean(sp_traits_CWM_length_perm_a[i,j,])
    sp_traits_CWM_length_perm_sd[i,j] = sd(sp_traits_CWM_length_perm_a[i,j,])
    
    sp_traits_CWM_mass_perm_mean[i,j] = mean(sp_traits_CWM_mass_perm_a[i,j,])
    sp_traits_CWM_mass_perm_sd[i,j] = sd(sp_traits_CWM_mass_perm_a[i,j,])
    
    sp_traits_CWM_dentition_perm_mean[i,j] = mean(sp_traits_CWM_dentition_perm_a[i,j,])
    sp_traits_CWM_dentition_perm_sd[i,j] = sd(sp_traits_CWM_dentition_perm_a[i,j,])
    
    sp_traits_CWM_diving_perm_mean[i,j] = mean(sp_traits_CWM_diving_perm_a[i,j,])
    sp_traits_CWM_diving_perm_sd[i,j] = sd(sp_traits_CWM_diving_perm_a[i,j,])
  }}

for (i in 1) {
  for (j in 1:3) {
    sp_traits_CWM_ratio_perm_mean[i,j] = mean(sp_traits_CWM_ratio_perm_a[i,j,])
    sp_traits_CWM_ratio_perm_sd[i,j] = sd(sp_traits_CWM_ratio_perm_a[i,j,])
  }}

for (i in 1) {
  for (j in 1:5) {
    sp_traits_CWM_prey_perm_mean[i,j] = mean(sp_traits_CWM_prey_perm_a[i,j,])
    sp_traits_CWM_prey_perm_sd[i,j] = sd(sp_traits_CWM_prey_perm_a[i,j,])
  }}

for (i in 1) {
  for (j in 1:6) {
    sp_traits_CWM_group_size_perm_mean[i,j] = mean(sp_traits_CWM_group_size_perm_a[i,j,])
    sp_traits_CWM_group_size_perm_sd[i,j] = sd(sp_traits_CWM_group_size_perm_a[i,j,])
  }}   


####PCEEZ CWMs####
sp_traits_CWM_length_PCEEZ<-functcomp(sp_traits_CWM_length, sp_comm_CWM, CWM.type = c("all"))
sp_traits_CWM_mass_PCEEZ<-functcomp(sp_traits_CWM_mass, sp_comm_CWM, CWM.type = c("all"))
sp_traits_CWM_ratio_PCEEZ<-functcomp(sp_traits_CWM_ratio, sp_comm_CWM, CWM.type = c("all"))
sp_traits_CWM_dentition_PCEEZ<-functcomp(sp_traits_CWM_dentition, sp_comm_CWM, CWM.type = c("all"))
sp_traits_CWM_diving_PCEEZ<-functcomp(sp_traits_CWM_diving, sp_comm_CWM, CWM.type = c("all"))
sp_traits_CWM_group_size_PCEEZ<-functcomp(sp_traits_CWM_group_size, sp_comm_CWM, CWM.type = c("all"))
sp_traits_CWM_prey_PCEEZ<-functcomp(sp_traits_CWM_prey, sp_comm_CWM, CWM.type = c("all"))


####SES for CWMs####
ses_CWM_length <- ( sp_traits_CWM_length_PCEEZ- sp_traits_CWM_length_perm_mean) / sp_traits_CWM_length_perm_sd # SES functional multidimensional alpha diversity indices
ses_CWM_mass <- ( sp_traits_CWM_mass_PCEEZ- sp_traits_CWM_mass_perm_mean) / sp_traits_CWM_mass_perm_sd # SES functional multidimensional alpha diversity indices
ses_CWM_ratio <- ( sp_traits_CWM_ratio_PCEEZ- sp_traits_CWM_ratio_perm_mean) / sp_traits_CWM_ratio_perm_sd # SES functional multidimensional alpha diversity indices
ses_CWM_dentition <- ( sp_traits_CWM_dentition_PCEEZ- sp_traits_CWM_dentition_perm_mean) / sp_traits_CWM_dentition_perm_sd # SES functional multidimensional alpha diversity indices
ses_CWM_diving <- ( sp_traits_CWM_diving_PCEEZ- sp_traits_CWM_diving_perm_mean) / sp_traits_CWM_diving_perm_sd # SES functional multidimensional alpha diversity indices
ses_CWM_group_size <- ( sp_traits_CWM_group_size_PCEEZ- sp_traits_CWM_group_size_perm_mean) / sp_traits_CWM_group_size_perm_sd # SES functional multidimensional alpha diversity indices
ses_CWM_prey <- ( sp_traits_CWM_prey_PCEEZ- sp_traits_CWM_prey_perm_mean) / sp_traits_CWM_prey_perm_sd # SES functional multidimensional alpha diversity indices

####Set up for qualitative trait plots#### 
dentition_CWM_null<-c(sp_traits_CWM_dentition_perm)
dentition_CWM_null<-as.data.frame(dentition_CWM_null)
dentition_CWM_null<-t(dentition_CWM_null)
dentition_CWM_null<-as.data.frame(dentition_CWM_null)
filter_dentition<-rep(c(1,2,3,4), 1000)
dentition_CWM_null$filter<-filter_dentition

dentition_FCT<-dentition_CWM_null%>%
  filter(filter_dentition==1)
dentition_LB<-dentition_CWM_null%>%
  filter(filter_dentition==2)
dentition_NFCT<-dentition_CWM_null%>%
  filter(filter_dentition==3)
dentition_SB<-dentition_CWM_null%>%
  filter(filter_dentition==4)

dentition_null_plots<-data.frame(dentition_FCT$PCEEZ, dentition_NFCT$PCEEZ, dentition_LB$PCEEZ, dentition_SB$PCEEZ)
dentition_null_plots



diving_CWM_null<-c(sp_traits_CWM_diving_perm)
diving_CWM_null<-as.data.frame(diving_CWM_null)
diving_CWM_null<-t(diving_CWM_null)
diving_CWM_null<-as.data.frame(diving_CWM_null)
filter_diving<-rep(c(1,2,3,4), 1000)
diving_CWM_null$filter<-filter_diving

diving_epipelagic<-diving_CWM_null%>%
  filter(filter_diving==2)
diving_upper_mesopelagic<-diving_CWM_null%>%
  filter(filter_diving==4)
diving_lower_mesopelagic<-diving_CWM_null%>%
  filter(filter_diving==3)
diving_bathypelagic<-diving_CWM_null%>%
  filter(filter_diving==1)

diving_null_plots<-data.frame(diving_epipelagic$PCEEZ, diving_upper_mesopelagic$PCEEZ, diving_lower_mesopelagic$PCEEZ, diving_bathypelagic$PCEEZ)
diving_null_plots



group_size_CWM_null<-c(sp_traits_CWM_group_size_perm)
group_size_CWM_null<-as.data.frame(group_size_CWM_null)
group_size_CWM_null<-t(group_size_CWM_null)
group_size_CWM_null<-as.data.frame(group_size_CWM_null)
filter_group_size<-rep(c(1,2,3,4,5,6), 1000)
group_size_CWM_null$filter<-filter_group_size

group_size_solitary<-group_size_CWM_null%>%
  filter(filter_group_size==6)
group_size_small<-group_size_CWM_null%>%
  filter(filter_group_size==5)
group_size_small_medium<-group_size_CWM_null%>%
  filter(filter_group_size==4)
group_size_medium<-group_size_CWM_null%>%
  filter(filter_group_size==2)
group_size_medium_large<-group_size_CWM_null%>%
  filter(filter_group_size==3)
group_size_large<-group_size_CWM_null%>%
  filter(filter_group_size==1)

group_size_null_plots<-data.frame(group_size_solitary$PCEEZ, group_size_small$PCEEZ, group_size_small_medium$PCEEZ, group_size_medium$PCEEZ, group_size_medium_large$PCEEZ, group_size_large$PCEEZ)
group_size_null_plots



prey_CWM_null<-c(sp_traits_CWM_prey_perm)
prey_CWM_null<-as.data.frame(prey_CWM_null)
prey_CWM_null<-t(prey_CWM_null)
prey_CWM_null<-as.data.frame(prey_CWM_null)
filter_prey<-rep(c(1,2,3,4,5), 1000)
prey_CWM_null$filter<-filter_prey

prey_cephalopods<-prey_CWM_null%>%
  filter(filter_prey==1)
prey_fish<-prey_CWM_null%>%
  filter(filter_prey==2)
prey_high_vertebrates<-prey_CWM_null%>%
  filter(filter_prey==3)
prey_non_cephalopod_invertebrates<-prey_CWM_null%>%
  filter(filter_prey==4)
prey_zooplankton<-prey_CWM_null%>%
  filter(filter_prey==5)

prey_null_plots<-data.frame(prey_cephalopods$PCEEZ, prey_fish$PCEEZ, prey_high_vertebrates$PCEEZ, prey_non_cephalopod_invertebrates$PCEEZ, prey_zooplankton$PCEEZ)
prey_null_plots

CWM_dentition_means<-as.data.frame(sp_traits_CWM_dentition_perm_mean)
CWM_dentition_means
CWM_dentition_sds<-as.data.frame(sp_traits_CWM_dentition_perm_sd)
CWM_dentition_sds
CWM_dentition_PCEEZ<-as.data.frame(sp_traits_CWM_dentition_PCEEZ)
CWM_dentition_PCEEZ

####Null distribution plots for dentition####
ggplot(dentition_null_plots, aes(x=dentition_FCT.PCEEZ))+
  geom_histogram(bins=12, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=CWM_dentition_means$Dentition_FCT), colour="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(CWM_dentition_means$Dentition_FCT+(1.96*CWM_dentition_sds$Dentition_FCT))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(CWM_dentition_means$Dentition_FCT-(1.96*CWM_dentition_sds$Dentition_FCT))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=CWM_dentition_PCEEZ$Dentition_FCT), colour="black", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x="Community Weighted Mean", y="Frequency")

ggplot(dentition_null_plots, aes(x=dentition_LB.PCEEZ))+
  geom_histogram(bins=12, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=CWM_dentition_means$Dentition_LB), colour="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(CWM_dentition_means$Dentition_LB+(1.96*CWM_dentition_sds$Dentition_LB))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(CWM_dentition_means$Dentition_LB-(1.96*CWM_dentition_sds$Dentition_LB))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=CWM_dentition_PCEEZ$Dentition_LB), colour="black", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x="Community Weighted Mean", y="Frequency")

ggplot(dentition_null_plots, aes(x=dentition_FCT.PCEEZ))+
  geom_histogram(bins=12, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=CWM_dentition_means$Dentition_FCT), colour="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(CWM_dentition_means$Dentition_FCT+(1.96*CWM_dentition_sds$Dentition_FCT))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(CWM_dentition_means$Dentition_FCT-(1.96*CWM_dentition_sds$Dentition_FCT))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=CWM_dentition_PCEEZ$Dentition_FCT), colour="black", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x="Community Weighted Mean", y="Frequency")

ggplot(dentition_null_plots, aes(x=dentition_SB.PCEEZ))+
  geom_histogram(bins=12, colour="#000000", fill="#999999")+
  theme_bw(base_size = 15)+
  geom_vline(aes(xintercept=CWM_dentition_means$Dentition_SB), colour="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(CWM_dentition_means$Dentition_SB+(1.96*CWM_dentition_sds$Dentition_SB))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=(CWM_dentition_means$Dentition_SB-(1.96*CWM_dentition_sds$Dentition_SB))), colour="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=CWM_dentition_PCEEZ$Dentition_SB), colour="black", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "right") +
  theme(axis.title.y = element_text(vjust = 3)) +
  theme(axis.title.x = element_text(vjust = -1))+
  labs(x="Community Weighted Mean", y="Frequency")


citation()
citation(package="mFD")
citation(package="funrar")
