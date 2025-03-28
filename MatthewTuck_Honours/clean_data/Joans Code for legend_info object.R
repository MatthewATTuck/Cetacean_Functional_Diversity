# Trait distribution along the Global Functional Space
######################################################

fe_tr <- as.data.frame(func_ent$fe_tr) # Functional traits per FE
colnames(fe_tr) <- c("Mob", "Size", "Hab.comp.", "Ch.Obl.", "Zone", "Feed")
fe_color <- data.frame(Mob = rep(NA, 167), Size = rep(NA, 167), Hab.comp.= rep(NA, 167), Ch.Obl. = rep(NA, 167), Zone= rep(NA, 167), Feed = rep(NA, 167)) # dataframe for colors
rownames(fe_color) <- rownames(fe_faxes)
legend_info <- as.list(NA) # list for legend colors

color_gr <- sequential_hcl(4, palette = "ag_GrnYl") # color palette for gradients
ord_gr <- c(1:4) # ordinal numbers of mobility
legend_info[["Mob"]] <- matrix(color_gr, nrow = 1, ncol = 4, dimnames = list(1, c(ord_gr)))

for (i in 1:4) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$Mob == ord_gr[i])),])), "Mob"] <- color_gr[i]
}

ord_gr <- c(1,10,100,1000) # ordinal numbers of size
legend_info[["Size"]] <- matrix(color_gr, nrow = 1, ncol = 4, dimnames = list(1, c(ord_gr)))

for (i in 1:4) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$Size == ord_gr[i])),])), "Size"] <- color_gr[i]
}

color_gr <- sequential_hcl(6, palette = "ag_GrnYl")
ord_gr <- c("Does not add", "Burrow forming", "Mat forming", "Bed forming", "Open bush forming", "Dense bush forming") # ordinal numbers of Environment
legend_info[["Hab.comp."]] <- matrix(color_gr, nrow = 1, ncol = 6, dimnames = list(1, c(ord_gr)))

for (i in 1:6) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$Hab.comp. == ord_gr[i])),])), "Hab.comp."] <- color_gr[i]
}

color_gr <- sequential_hcl(3, palette = "ag_GrnYl")
ord_gr <- c("No", "Other CBE", "Vent") # ordinal numbers of Environment
legend_info[["Ch.Obl."]] <- matrix(color_gr, nrow = 1, ncol = 3, dimnames = list(1, c(ord_gr)))

for (i in 1:3) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$Ch.Obl. == ord_gr[i])),])), "Ch.Obl."] <- color_gr[i]
}

color_gr <- sequential_hcl(3, palette = "ag_GrnYl")
ord_gr <- c("Low", "Medium", "High") # ordinal numbers of Environment
legend_info[["Zone"]] <- matrix(color_gr, nrow = 1, ncol = 3, dimnames = list(1, c(ord_gr)))

for (i in 1:3) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$Zone == ord_gr[i])),])), "Zone"] <- color_gr[i]
}

color_gr <- sequential_hcl(4, palette = "ag_GrnYl")
ord_gr <- c("Carnivore", "Detritivore", "Bacterivore", "Symbiotic") # ordinal numbers of Environment
legend_info[["Feed"]] <- matrix(color_gr, nrow = 1, ncol = 4, dimnames = list(1, c(ord_gr)))

for (i in 1:4) {
  fe_color[c(rownames(fe_tr[c(which(fe_tr$Feed == ord_gr[i])),])), "Feed"] <- color_gr[i]
}