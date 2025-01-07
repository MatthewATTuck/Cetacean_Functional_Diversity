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

cetacean_dataset<- read.csv("mtuck_honours_Rdocument.csv") # Import csv dataset

sp_comm <- as.data.frame(cetacean_dataset[,c(19:20)], row.names=cetacean_dataset$species_name)# species occurrence data for PCEEZ
sp_traits <- as.data.frame(cetacean_dataset[,c(3:10)], row.names = cetacean_dataset$species_name)
trait_info <- data.frame(trait_name = colnames(sp_traits[,c(1:8)]), trait_type = c("O","O","O","N","N","O","O","N"))


sp_traits$max_length_m <- ordered(sp_traits$max_length_m, levels=c("A", "B")) 
sp_traits$max_mass_kg <- ordered(sp_traits$max_mass_kg, levels=c("A","B","C","D")) 
sp_traits$max_mass_max_length_ratio <- ordered(sp_traits$max_mass_max_length_ratio, levels=c("A", "B", "C")) 
sp_traits$dentition <- as.factor(sp_traits$dentition) 
sp_traits$migratory_behaviour <- as.factor(sp_traits$migratory_behaviour) 
sp_traits$max_diving_depth <- ordered(sp_traits$max_diving_depth, levels=c("A", "B", "C", "D")) 
sp_traits$group_size <- ordered(sp_traits$group_size, levels = c("A", "B", "C", "D", "E", "F")) 
sp_traits$prey_choice <- as.factor(sp_traits$prey_choice) 

#### Establishing Gower's Functional Distance between species ####
func_dist_sp <- funct.dist(sp_traits, 
                           tr_cat = trait_info, 
                           metric = "gower", 
                           ordinal_var = "metric",
                           weight_type = "equal") # Functional dissimilarity (Gower) distance between species

#### PCoA construction ####
func_qual <- quality.fspaces(func_dist_sp, 
                             deviation_weighting = c("absolute", "squared"),
                             maxdim_pcoa = 9,
                             fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

apply(func_qual$quality_fspaces, 2, which.min) # best number of dimesions calculated with MAD index (4D) (Maire et al. 2015)
sp_faxes <- func_qual$details_fspaces$sp_pc_coord # species coordinates in the functional space

sp_comm_biog = as.matrix(sp_comm) # create new table with species classified in bioprovinces
sp_comm_biog <- decostand(sp_comm_biog, "pa") # presence/abscense transformation (1/0)
sp_comm_biog<-(t(sp_comm_biog))

####Calculation of PCEEZ Functional Dispersion (and other indices)####
func_ind_sp <- alpha.fd.multidim(sp_faxes_coord = sp_faxes[,c(1:4)],
                                 asb_sp_w = sp_comm_biog,
                                 scaling = TRUE,
                                 details_returned = T,
                                 verbose = T) # mutidimensinal functional diversity indices in a 4D functional space of vent regions



func_ind_fdis <- func_ind_sp$functional_diversity_indices$fdis # functional dispersion 


#### Unique Trait Combinations ####
func_ent <- sp.to.fe(sp_traits, trait_info, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe <- alpha.fd.fe(sp_comm_biog, sp_to_fe = func_ent) # compute UTC derived indices
View(func_ind_fe$asb_fdfe) # UTC Functional indexes


install.packages("funrar")
library(funrar)

sp_traits_unique<-as.matrix(func_dist_sp)
unique_sp<-funrar(pres_matrix = sp_comm_biog, dist_matrix = sp_traits_unique, rel_abund = FALSE)
uniqueness<-uniqueness(sp_comm_biog, sp_traits_unique)
colSums(uniqueness==0) ## 44 Unique species Identified

#### NULL Distribution Setup ####
n_perm = 1000 # number of permutations
multif_perm <- as.list(rep(NA, n_perm)) # list to save permutations results of multidimensional functional alpha-diversity indices (FDis)
sp_perm<-sp_comm_biog[c(1),]
for(i in seq(n_perm)){
  
  sp_comm_n = randomizeMatrix(sp_comm_biog, null.model = "richness") # Randomizations maintaining species occurrence frequency (global) and region species richness  
  
  multif_n <- alpha.fd.multidim(sp_faxes_coord = sp_faxes[,c(1:4)],
                                asb_sp_w = sp_comm_n,
                                ind_vect= "fdis",
                                scaling = T,
                                details_returned = T,
                                verbose = F) #  Null FDis
  
  multif_perm[[i]] <- multif_n$functional_diversity_indices["fdis"] }# store null FDis results 

#### Wrap Null Distribution Results ####
multif_perm_a <- array(NA, dim = c(2,1,0)) # create an array with 1 col & 11 rows for null FDis results

for (i in c(1:n_perm)) { 
  multif_perm_a <- abind(multif_perm_a, multif_perm[[i]])} # FRic & FDis
  
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
#### CSV loading and setting up data frames####

sp_comm <- as.data.frame(cetacean_dataset[,c(19:20)], row.names=cetacean_dataset$species_name)# species occurrence data for PCEEZ
sp_traits <- as.data.frame(cetacean_dataset[,c(3:6, 8:10)], row.names = cetacean_dataset$species_name)
trait_info <- data.frame(trait_name = colnames(sp_traits[,c(1:7)]), trait_type = c("O","O","O","N","O","O","N"))


sp_traits$max_length_m <- ordered(sp_traits$max_length_m, levels=c("A", "B")) 
sp_traits$max_mass_kg <- ordered(sp_traits$max_mass_kg, levels=c("A","B","C","D")) 
sp_traits$max_mass_max_length_ratio <- ordered(sp_traits$max_mass_max_length_ratio, levels=c("A", "B", "C")) 
sp_traits$dentition <- as.factor(sp_traits$dentition) 
sp_traits$max_diving_depth <- ordered(sp_traits$max_diving_depth, levels=c("A", "B", "C", "D")) 
sp_traits$group_size <- ordered(sp_traits$group_size, levels = c("A", "B", "C", "D", "E", "F")) 
sp_traits$prey_choice <- as.factor(sp_traits$prey_choice) 

#### Establishing Gower's Functional Distance between species ####
func_dist_sp <- funct.dist(sp_traits, 
                           tr_cat = trait_info, 
                           metric = "gower", 
                           ordinal_var = "metric",
                           weight_type = "equal") # Functional dissimilarity (Gower) distance between species

#### PCoA construction ####
func_qual <- quality.fspaces(func_dist_sp, 
                             deviation_weighting = c("absolute", "squared"),
                             maxdim_pcoa = 9,
                             fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

apply(func_qual$quality_fspaces, 2, which.min) # best number of dimesions calculated with MAD index (4D) (Maire et al. 2015)
sp_faxes <- func_qual$details_fspaces$sp_pc_coord # species coordinates in the functional space

sp_comm_biog = as.matrix(sp_comm) # create new table with species classified in bioprovinces
sp_comm_biog <- decostand(sp_comm_biog, "pa") # presence/abscense transformation (1/0)
sp_comm_biog<-(t(sp_comm_biog))

####Calculation of PCEEZ Functional Dispersion (and other indices)####
func_ind_sp <- alpha.fd.multidim(sp_faxes_coord = sp_faxes[,c(1:3)],
                                 asb_sp_w = sp_comm_biog,
                                 scaling = TRUE,
                                 details_returned = T,
                                 verbose = T) # mutidimensinal functional diversity indices in a 4D functional space of vent regions



func_ind_fdis <- func_ind_sp$functional_diversity_indices$fdis # functional dispersion 


#### Unique Trait Combinations ####
func_ent <- sp.to.fe(sp_traits, trait_info, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe <- alpha.fd.fe(sp_comm_biog, sp_to_fe = func_ent) # compute UTC derived indices
View(func_ind_fe$asb_fdfe) # UTC Functional indexes


sp_traits_unique<-as.matrix(func_dist_sp)
unique_sp<-funrar(pres_matrix = sp_comm_biog, dist_matrix = sp_traits_unique, rel_abund = FALSE)
uniqueness<-uniqueness(sp_comm_biog, sp_traits_unique)
colSums(uniqueness==0) ## 44 Unique species Identified

#### NULL Distribution Setup ####
n_perm = 1000 # number of permutations
multif_perm <- as.list(rep(NA, n_perm)) # list to save permutations results of multidimensional functional alpha-diversity indices (FDis)
sp_perm<-sp_comm_biog[c(1),]
for(i in seq(n_perm)){
  
  sp_comm_n = randomizeMatrix(sp_comm_biog, null.model = "richness") # Randomizations maintaining species occurrence frequency (global) and region species richness  
  
  multif_n <- alpha.fd.multidim(sp_faxes_coord = sp_faxes[,c(1:3)],
                                asb_sp_w = sp_comm_n,
                                ind_vect= "fdis",
                                scaling = T,
                                details_returned = T,
                                verbose = F) #  Null FDis
  
  multif_perm[[i]] <- multif_n$functional_diversity_indices["fdis"] }# store null FDis results 

#### Wrap Null Distribution Results ####
multif_perm_a <- array(NA, dim = c(2,1,0)) # create an array with 1 col & 11 rows for null FDis results

for (i in c(1:n_perm)) { 
  multif_perm_a <- abind(multif_perm_a, multif_perm[[i]])} # FRic & FDis

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

#End of Code