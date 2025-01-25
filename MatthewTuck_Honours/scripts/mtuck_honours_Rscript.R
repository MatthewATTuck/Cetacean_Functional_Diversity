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
#                                                                              #
################################################################################

#### CSV loading and setting up data frames####

sp_comm_nm <- as.data.frame(cetacean_dataset[,c(19:20)], row.names=cetacean_dataset$species_name)# species occurrence data for PCEEZ without migration
sp_traits_nm <- as.data.frame(cetacean_dataset[,c(3:6, 8:10)], row.names = cetacean_dataset$species_name) #species level trait data without migration
trait_info_nm <- data.frame(trait_name = colnames(sp_traits_nm[,c(1:7)]), trait_type = c("O","O","O","N","O","O","N")) #trait information data frame without migration


sp_traits_nm$max_length_m <- ordered(sp_traits$max_length_m, levels=c("small", "intermediate", "large", "very_large")) 
sp_traits_nm$max_mass_kg <- ordered(sp_traits$max_mass_kg, levels=c("small","intermediate","large","very_large")) 
sp_traits_nm$max_mass_max_length_ratio <- ordered(sp_traits$max_mass_max_length_ratio, levels=c("low", "medium", "high")) 
sp_traits_nm$dentition <- as.factor(sp_traits$dentition) 
sp_traits_nm$max_diving_depth <- ordered(sp_traits$max_diving_depth, levels=c("epipelagic", "upper_mesopelagic", "lower_mesopelagic", "bathypelagic")) 
sp_traits_nm$average_group_size <- ordered(sp_traits$average_group_size, levels = c("Solitary", "Small", "Medium_Small", "Medium", "Medium_Large", "Large")) 
sp_traits_nm$prey_choice <- as.factor(sp_traits$prey_choice) #setting up species traits as ordinal and nominal variables


#### Establishing Gower's Functional Distance between species ####
func_dist_sp_nm <- funct.dist(sp_traits_nm, 
                           tr_cat = trait_info_nm, 
                           metric = "gower", 
                           ordinal_var = "metric",
                           weight_type = "equal") # Functional dissimilarity (Gower) distance between species

#### PCoA construction ####
func_qual_nm <- quality.fspaces(func_dist_sp_nm, 
                             deviation_weighting = c("absolute", "squared"),
                             maxdim_pcoa = 9,
                             fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

apply(func_qual_nm$quality_fspaces, 2, which.min) # best number of dimesions calculated with MAD index (3D) (Maire et al. 2015)
sp_faxes_nm <- func_qual_nm$details_fspaces$sp_pc_coord # species coordinates in the functional space

####Rough PCoA Plots####
ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(2)], colour=cetacean_dataset$Family))+
  geom_point() +
  labs(x="PCoA 1", y="PCoA 2")+
  theme_light()

ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(1)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(3)], colour =cetacean_dataset$Family))+
  geom_point() +
  labs(x="PCoA 1", y="PCoA 3")+
  theme_light()

ggplot(func_qual_nm$details_fspaces$sp_pc_coord, aes(x=func_qual_nm$details_fspaces$sp_pc_coord[,c(2)], y= func_qual_nm$details_fspaces$sp_pc_coord[, c(3)], colour=cetacean_dataset$Family))+
  geom_point() +
  labs(x="PCoA 2", y="PCoA 3")+
  theme_light()

#### Adjustments to species location matrix for later use ####

sp_comm_biog_nm = as.matrix(sp_comm_nm) # create matrix with species classified based on occurrence in the PCEEZ
sp_comm_biog_nm <- decostand(sp_comm_biog_nm, "pa") # presence/abscense transformation (1/0)
sp_comm_biog_nm<-(t(sp_comm_biog_nm)) #transpose matrix for use in the alpha.fd.multidim function

####Calculation of PCEEZ Functional Dispersion (and other indices)####
func_ind_sp_nm <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_nm[,c(1:3)],
                                 asb_sp_w = sp_comm_biog_nm,
                                 scaling = TRUE,
                                 details_returned = T,
                                 verbose = T) # mutidimensinal functional diversity indices in a 3D functional space



func_ind_fdis_nm <- func_ind_sp_nm$functional_diversity_indices$fdis # functional dispersion without migration trait 
func_ind_fdis_nm

#### Unique Trait Combinations ####
func_ent_nm <- sp.to.fe(sp_traits_nm, trait_info_nm, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe_nm <- alpha.fd.fe(sp_comm_biog_nm, sp_to_fe = func_ent_nm) # compute UTC derived indices
func_ind_fe_nm$asb_fdfe # UTC Functional indexes

####Funrar####
sp_traits_unique_nm<-as.matrix(func_dist_sp_nm) #setting the functional dissimilarity data as a matrix
unique_sp_nm<-funrar(pres_matrix = sp_comm_biog_nm, dist_matrix = sp_traits_unique_nm, rel_abund = FALSE) #calculating the number of unique trait assemblages in the data set
uniqueness_nm<-uniqueness(sp_comm_biog_nm, sp_traits_unique_nm) #calculating "uniqueness" for each species in the data set
uniqueness_nm
colSums(uniqueness_nm==0) ## 45 species identified which had uniqueness=0, which means they shared trait assemblage with at least one other species 
#this means that 43 species are unique (~49% of the total species richness)

#### NULL Distribution Setup ####
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
  hist(sapply(multif_perm_nm, "[[", 1))
  multif_perm_nm
  
  func_ind_fe_n_nm <- alpha.fd.fe(sp_comm_n_nm, sp_to_fe = func_ent_nm) # Null FE indices
  fe_perm_nm[[i]] <- func_ind_fe_n_nm$asb_fdfe[,c("fred", "fvuln")]} # store null FRed indices result
fe_perm_nm
fe_perm_nm[,1]
hist(sapply(fe_perm_nm, "[[", 1))
hist(sapply(fe_perm_nm, "[[", 3))
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

#### SES and p-values ####
ses_multif_nm <- (func_ind_sp_nm$functional_diversity_indices[,c("fdis")] - multif_perm_mean_nm) / multif_perm_sd_nm # SES functional multidimensional alpha diversity indices
head(ses_multif_nm)
pnorm(ses_multif_nm)

ses_fe_nm <- (func_ind_fe_nm$asb_fdfe[,c("fred")] - fe_perm_mean_nm) / fe_perm_sd_nm # SES FE indices FRed
ses_fe_nm
pnorm(ses_fe_nm)


####Individual trait analyses####
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

cetacean_ind_traits_pceez <- cetacean_ind_traits %>%
  filter(ind_pceez=="Yes")

morphometric_means<-cetacean_ind_traits_pceez %>%
  summarize(mean_length=mean(max_length_m), mean_mass=mean(max_mass_kg), mean_ratio=mean(max_mass_max_length_ratio))
morphometric_means

morphometric_sd<-cetacean_ind_traits_pceez%>%
  summarise(sd_length=sd(max_length_m), sd_mass=sd(max_mass_kg), sd_ratio=sd(max_mass_max_length_ratio))
morphometric_sd

n_perm=1000
max_length_perm <- as.list(rep(NA, n_perm))
max_mass_perm<-as.list(rep(NA, n_perm))
max_ratio_perm<-as.list(rep(NA, n_perm))

for(i in seq(n_perm)){
  
  max_length_perm_n<-sample(cetacean_ind_traits$max_length_m, size=25, replace=FALSE, prob=NULL)
  max_length_perm_means_n<-mean( max_length_perm_n)
  max_length_perm[[i]]<-max_length_perm_means_n}

max_length_perm 
max_length_perm<-as.vector(t(max_length_perm))

mean_max_length_perm<-mean(max_length_perm)
mean_max_length_perm
sd_max_length_perm<-sd(max_length_perm)
sd_max_length_perm

ses_max_length_pceez<-(morphometric_means$mean_length - mean_max_length_perm)/sd_max_length_perm
ses_max_length_pceez
pnorm(ses_max_length_pceez)
#End of Code



####CWM####
library(FD)
cetacean_CWMs<-data.frame(functcomp(sp_traits, sp_comm_biog, CWM.type = c("all")))
cetacean_CWMs<-data.frame(t(cetacean_CWMs))
cetacean_CWMs<-rownames_to_column(cetacean_CWMs, "trait_designation")
CWMs<-c(cetacean_CWMs$PCEEZ, cetacean_CWMs$Global)
CWMs_plots<-data.frame(CWMs)
CWMs_plots$trait_types<-c("max_length", "max_length", "max_length", "max_length", "max_mass", "max_mass", "max_mass", "max_mass", "max_mass_length_ratio", "max_mass_length_ratio", "max_mass_length_ratio", "dentition", "dentition", "dentition", "dentition", "migration", "migration", "migration", "max_diving_depth", "max_diving_depth", "max_diving_depth", "max_diving_depth", "average_group_size", "average_group_size", "average_group_size", "average_group_size", "average_group_size", "average_group_size", "prey_choice", "prey_choice", "prey_choice", "prey_choice", "prey_choice", "max_length", "max_length", "max_length", "max_length", "max_mass", "max_mass", "max_mass", "max_mass", "max_mass_length_ratio", "max_mass_length_ratio", "max_mass_length_ratio", "dentition", "dentition", "dentition", "dentition", "migration", "migration", "migration", "max_diving_depth", "max_diving_depth", "max_diving_depth", "max_diving_depth", "average_group_size", "average_group_size", "average_group_size", "average_group_size", "average_group_size", "average_group_size", "prey_choice", "prey_choice", "prey_choice", "prey_choice", "prey_choice")
CWMs_plots$location<-c("PCEEZ", "PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ", "PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ","PCEEZ", "Global", "Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global","Global")
CWMs_plots$trait_labels<-c("A", "B", "C", "D", "A", "B", "C", "D", "A", "B", "C", "A", "B", "C", "D", "A", "B", "C", "A", "B", "C", "D", "A", "B", "C", "D", "E", "F", "A", "B", "C", "D", "E", "A", "B", "C", "D", "A", "B", "C", "D", "A", "B", "C", "A", "B", "C", "D", "A", "B", "C", "A", "B", "C", "D", "A", "B", "C", "D", "E", "F", "A", "B", "C", "D", "E")
CWMs_plots$traitdesignations<-c(cetacean_CWMs$trait_designation,cetacean_CWMs$trait_designation) 
CWMs_plots


ggplot(CWMs_plots, aes(fill=trait_labels, x=location, y=CWMs))+
  geom_bar(position="stack", stat = "identity")+
  facet_wrap(~trait_types, ncol=2, nrow=4, scale="free_x")
  theme_light()


