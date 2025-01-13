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


#### Unique Trait Combinations ####
func_ent <- sp.to.fe(sp_traits, trait_info, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe <- alpha.fd.fe(sp_comm_biog, sp_to_fe = func_ent) # compute UTC derived indices
View(func_ind_fe$asb_fdfe) # UTC Functional indexes


library(funrar)

sp_traits_unique<-as.matrix(func_dist_sp) #setting the functional dissimilarity data as a matrix
unique_sp<-funrar(pres_matrix = sp_comm_biog, dist_matrix = sp_traits_unique, rel_abund = FALSE) #calculating the number of unique trait assemblages in the data set
uniqueness<-uniqueness(sp_comm_biog, sp_traits_unique) #calculating "uniqueness" for each species in the data set
colSums(uniqueness==0) ## 44 species identified which had uniqueness=0, which means they shared trait assemblage with at least one other species 
#this means that 44 species are unique (50% of the total species richness)

#### NULL Distribution Setup ####
n_perm = 1000 # number of permutations
multif_perm <- as.list(rep(NA, n_perm)) # list to save permutations results of multidimensional functional alpha-diversity indices (FDis)
sp_perm<-sp_comm_biog[c(1),]

for(i in seq(n_perm)){
  sp_comm_n = randomizeMatrix(sp_comm_biog, null.model = "richness") # Randomizations maintaining species richness in the PCEEZ
  
  multif_n <- alpha.fd.multidim(sp_faxes_coord = sp_faxes[,c(1:4)],
                                asb_sp_w = sp_comm_n,
                                ind_vect= "fdis",
                                scaling = T,
                                details_returned = T,
                                verbose = F) #  Null Distribution for FDis
  
  multif_perm[[i]] <- multif_n$functional_diversity_indices["fdis"] }# store null FDis results 

#### Wrap Null Distribution Results ####
multif_perm_a <- array(NA, dim = c(2,1,0)) # create an array with 1 col & 2 rows for null FDis results

for (i in c(1:n_perm)) { 
  multif_perm_a <- abind(multif_perm_a, multif_perm[[i]])} # FDis
  
multif_perm_mean <- matrix(NA, ncol = 1, nrow= 2, dimnames = list(rownames(multif_perm_a), colnames(multif_perm_a))) # matrix to save mean values of FDis permutations
multif_perm_sd <- matrix(NA, ncol = 1, nrow= 2, dimnames = list(rownames(multif_perm_a), colnames(multif_perm_a))) # matrix to save sd values of FDis permutations 
for (i in 1:2) {
  for (j in "fdis") {
    multif_perm_mean[i,j] = mean(multif_perm_a[i,j,])
    multif_perm_sd[i,j] = sd(multif_perm_a[i,j,])
  }
} # mean and sd values of FDis permutations

#### SES ####
ses_multif <- (func_ind_sp$functional_diversity_indices[,c("fdis")] - multif_perm_mean) / multif_perm_sd # SES functional multidimensional alpha diversity indices


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


#### Unique Trait Combinations ####
func_ent_nm <- sp.to.fe(sp_traits_nm, trait_info_nm, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe_nm <- alpha.fd.fe(sp_comm_biog_nm, sp_to_fe = func_ent_nm) # compute UTC derived indices
View(func_ind_fe_nm$asb_fdfe) # UTC Functional indexes


sp_traits_unique_nm<-as.matrix(func_dist_sp_nm) #setting the functional dissimilarity data as a matrix
unique_sp_nm<-funrar(pres_matrix = sp_comm_biog_nm, dist_matrix = sp_traits_unique_nm, rel_abund = FALSE) #calculating the number of unique trait assemblages in the data set
uniqueness_nm<-uniqueness(sp_comm_biog_nm, sp_traits_unique_nm) #calculating "uniqueness" for each species in the data set
colSums(uniqueness_nm==0) ## 45 species identified which had uniqueness=0, which means they shared trait assemblage with at least one other species 
#this means that 43 species are unique (~49% of the total species richness)

#### NULL Distribution Setup ####
n_perm = 1000 # number of permutations
multif_perm_nm <- as.list(rep(NA, n_perm)) # list to save permutations results of multidimensional functional alpha-diversity indices (FDis)
sp_perm_nm<-sp_comm_biog[c(1),]
for(i in seq(n_perm)){
  
  sp_comm_n_nm = randomizeMatrix(sp_comm_biog_nm, null.model = "richness") # Randomizations maintaining species richness in the PCEEZ
 
   multif_n_nm <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_nm[,c(1:3)],
                                asb_sp_w = sp_comm_n_nm,
                                ind_vect= "fdis",
                                scaling = T,
                                details_returned = T,
                                verbose = F) #  Null Distribution for FDis
  
  multif_perm_nm[[i]] <- multif_n_nm$functional_diversity_indices["fdis"] }# store null FDis results 

#### Wrap Null Distribution Results ####
multif_perm_a_nm <- array(NA, dim = c(2,1,0)) # create an array with 1 col & 2 rows for null FDis results

for (i in c(1:n_perm)) { 
  multif_perm_a_nm <- abind(multif_perm_a_nm, multif_perm_nm[[i]])} # FDis

multif_perm_mean_nm <- matrix(NA, ncol = 1, nrow= 2, dimnames = list(rownames(multif_perm_a_nm), colnames(multif_perm_a_nm))) # matrix to save mean values of FDis permutations
multif_perm_sd_nm <- matrix(NA, ncol = 1, nrow= 2, dimnames = list(rownames(multif_perm_a_nm), colnames(multif_perm_a_nm))) # matrix to save sd values of FDis permutations 
for (i in 1:2) {
  for (j in "fdis") {
    multif_perm_mean_nm[i,j] = mean(multif_perm_a_nm[i,j,])
    multif_perm_sd_nm[i,j] = sd(multif_perm_a_nm[i,j,])
  }
} # mean and sd values of FDis permutations

#### SES ####
ses_multif_nm <- (func_ind_sp_nm$functional_diversity_indices[,c("fdis")] - multif_perm_mean_nm) / multif_perm_sd_nm # SES functional multidimensional alpha diversity indices

#End of Code