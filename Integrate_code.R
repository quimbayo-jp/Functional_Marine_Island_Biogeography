## --------------------------------------------------------------------------------------------------------- ##
##  Project Name:  A trait-based approach to marine island biogeography                                      ##
##                                                                                                           ##
##  Objective:     The Theory of Island Biogeography (TIB) is widely recognized as pivotal for explaining    ## 
##                 species diversity patterns in insular system. However, how this theory can predict the    ## 
##                 functional diversity patterns in oceanic islands, remains a challenge. Here,              ##
##                 we evaluated the influence of past and current island features related to isolation,      ##
##                 area, and geological age in shaping functional diversity of reef fishes across oceanic    ##
##                 islands worldwide.                                                                        ##
##                                                                                                           ##
##  Authors:                                                                                                 ##
##  Debora S. Ferrari (Universidade Federal de Santa Catarina - Brazil)                                      ##
##  Sergio R. Floeter (Universidade Federal de Santa Catarina - Brazil)                                      ##
##  Fabien Leprieur   (Universite de Montpellier - France)                                                   ##
##  Juan P. Quimbayo  (Universidade de São Paulo and MarineGeo Smithsonian Environmental Research Center)    ##
##
##  Date:          2021-11-19                                                                                ##
##                                                                                                           ##
##  Notes:         1. This file is intended to provide a guide to the basic                                  ##
##                    project workflow, attempting to 'integrate_analysis' the steps                         ##
##                    necessary to conduct the analyses & visual outputs.
##
##  Script created by: Juan P. Quimbayo                                                                                                        ##
## --------------------------------------------------------------------------------------------------------- ##


# clean up
rm(list=ls())

source ("R_Code/Functions.R")
source ("R_Code/Packages.R")
source ("R_Code/Maps.R")

# Loading all the database -----------------------------------------

factors_islands <- read.delim ("Datasets/Factors_islands.txt", dec = ",", na.strings=c("", "NA"), sep="\t")
traits_species  <- read.delim ("Datasets/Traits_species.txt", dec = ",", na.strings=c("", "NA"), sep="\t")
local_checklits <- read.delim ("Datasets/Islands_checklists.txt")

# Preparing database to estimated functional index -----------------
local_checklits <- droplevels(local_checklits[!duplicated(local_checklits$Genus_species), ])
traits_species  <- droplevels(traits_species[traits_species$Genus_species %in% local_checklits$Genus_species, ])
traits_species  <- na.omit(traits_species)
local_checklits <- droplevels(local_checklits[local_checklits$Genus_species %in% traits_species$Genus_species, ])

traits_species$Size_class <- NA
traits_species$Size_class <- ifelse (traits_species$Size_cm<=7, "C1", traits_species$Size_class)
traits_species$Size_class <- ifelse (traits_species$Size_cm>7  & traits_species$Size_cm<=15, "C2", traits_species$Size_class)
traits_species$Size_class <- ifelse (traits_species$Size_cm>15 & traits_species$Size_cm<=30, "C3", traits_species$Size_class)
traits_species$Size_class <- ifelse (traits_species$Size_cm>30 & traits_species$Size_cm<=50, "C4", traits_species$Size_class)
traits_species$Size_class <- ifelse (traits_species$Size_cm>50 & traits_species$Size_cm<=80, "C5", traits_species$Size_class)
traits_species$Size_class <- ifelse (traits_species$Size_cm>80, "C6", traits_species$Size_class)
unique (traits_species$Size_class)

dim (local_checklits);dim (traits_species)  

traits  <- traits_species [,c("Genus_species","Size_class","Home_range",
                              "Diet_level","Activity","Level_water",
                              "Schooling")]

rownames(traits) <- traits$Genus_species
traits  <- droplevels (traits [!duplicated(traits$Genus_species), ])
traits  <- traits [,-1]  
traits  <- na.omit(traits)

mobility   <- sapply(traits$Home_range,  function(x) {if (x=="Sed"){1} else if (x=="Mob"){2} else if (x=="VMob"){3}})
activity   <- sapply(traits$Activity,    function(x) {if (x=="Day"){1} else if (x=="Both"){2} else if (x=="Night"){3}})
schooling  <- sapply(traits$Schooling,   function(x) {if (x=="Sol"){1} else if (x=="Pair"){2} else if (x=="SmallG"){3} else if (x=="MedG"){4} else if (x=="LargeG"){5}})
level      <- sapply(traits$Level_water, function(x) {if (x=="Bottom"){1} else if (x=="Low"){2} else if (x=="High"){3}})
sclass     <- sapply(traits$Size_class,  function(x) {if (x=="C1"){1} else if (x=="C2"){2} else if (x=="C3"){3} else if (x=="C4"){4} else if (x=="C5"){5} else if (x=="C6"){6}})

mobility   <- ordered (mobility)
activity   <- ordered (activity)
schooling  <- ordered (schooling)
level      <- ordered (level)
sclass     <- ordered (sclass)

traits     <- data.frame (Genus_species=rownames(traits), Home_range=mobility, Activity=activity,
                       Size_group=schooling, Level_water=level,
                       Size_class=sclass, Diet=traits$Diet_level)

traits     <- traits [order(traits$Genus_species, decreasing = F), ]
rownames(traits) <- traits$Genus_species
traits    <- traits[,-1]

species <- local_checklits [,c(1,3:77)]
species <- droplevels (species [species$Genus_species %in% rownames(traits), ])
species <- species [order(species$Genus_species, decreasing = FALSE), ]
rownames(species) <- species$Genus_species

species <- (species [,-1])
species <- data.matrix(species)
species <- t (species)

rm (mobility, activity, schooling, level)
dim (traits); dim (species)

traits  <- data.matrix(traits)

# Estimating the Gower dissimilarity --------------------------------
gower_matrix <- daisy (traits, metric=c("gower")) 
pcoa <- dudi.pco (quasieuclid(gower_matrix), scannf=F, nf=6)

# Estimation of the quality of functional space ---------------------
source ("R_Code/quality_funct_space_fromdist.R")
quality <- quality_funct_space_fromdist (gower_matrix,  nbdim=6,   plot="Output/Figures/Fig_S1") 
quality$meanSD # the minimal value corresponds to the best space to use, here 6 axes 
Inertia_6axis <- (sum(pcoa$eig[1:6])) /(sum(pcoa$eig)) # percentage of inertia explained by the 6 first axes = 87%
Inertia_4axis <- (sum(pcoa$eig[1:4])) /(sum(pcoa$eig)) # percentage of inertia explained by the 6 first axes = 80%

# Estimation of the functional index -------------------------------
source ("R_Code/multidimFD.R")
FD_Islands <- multidimFD (as.matrix(pcoa$li[,1:4]), species)
FD_Islands <- as.data.frame(FD_Islands)


FD_Islands$Sitename  <- rownames (FD_Islands)
rownames(FD_Islands) <- NULL
FD_Islands <- FD_Islands[,c("Sitename","Nb_sp","FRic","FDiv","FEve")]
# Estimation the functional redundancy and vulnerability ----------
source ("R_Code/species_to_FE.R")
source ("R_Code/FE_metrics.R")

traits <- traits_species [,c("Genus_species","Home_range",
                             "Diet_level","Activity","Level_water",
                             "Schooling", "Size_class")]

rownames(traits) <- traits$Genus_species 
traits           <- traits [,-1]

species <- local_checklits [,c(1,3:77)]
species <- droplevels (species[species$Genus_species %in% rownames(traits), ])        
rownames (species) <- species$Genus_species
species <- species [,-1]
species <- data.matrix(species)
species <- t (species)

F.Entities  <- species_to_FE (traits)
F.Red_F.Vul <- FE_metrics (F.Entities, species)
F.Red_F.Vul <- data.frame (F.Red_F.Vul)
F.Red_F.Vul$Sitename  <- rownames (F.Red_F.Vul)
rownames(F.Red_F.Vul) <- NULL

FD_Islands <- left_join(FD_Islands, F.Red_F.Vul)

# Rescale predictors and testing correlation among them ----------
FDIslands_Factors <- left_join(FD_Islands, factors_islands) 
FDIslands_Factors <- na.omit(FDIslands_Factors)
write.table (FD_Islands, "Output/Tables/FD_Islands.csv", sep = ";", dec = ",", row.names = F)
save.image (file="FIndex_Islands.RData")

# Rescaling the single terms
FDIslands_Factors$Age_Ma_Scale       <- rescale_variables(FDIslands_Factors$Age_Ma)
FDIslands_Factors$Current_Area_Scale <- rescale_variables(log(FDIslands_Factors$Shallow_Reef_area_200mt+1))
FDIslands_Factors$Past_Area_Scale    <- rescale_variables(log(FDIslands_Factors$Area_sub_past_200mt)+1)
FDIslands_Factors$Isolation_Scale    <- rescale_variables(log(FDIslands_Factors$Isol_habitat)+1)

# Including polynomial effect 
FDIslands_Factors$Age_Ma_Scale_2       <- poly (FDIslands_Factors$Age_Ma_Scale, 2)[,2]
FDIslands_Factors$Current_Area_Scale_2 <- poly (FDIslands_Factors$Current_Area_Scale, 2)[,2]
FDIslands_Factors$Past_Area_Scale_2    <- poly (FDIslands_Factors$Past_Area_Scale, 2)[,2]
FDIslands_Factors$Isolation_Scale_2    <- poly (FDIslands_Factors$Isolation_Scale, 2)[,2]

# Exploring correlation between predictors 
predictors <- FDIslands_Factors[,c("Age_Ma_Scale","Current_Area_Scale",
                                   "Past_Area_Scale","Isolation_Scale")]

colnames (predictors) <- c("Age", "C.Area", "P.Area", "Isol")

#pdf (file = "Output/Figures/Fig_S2.pdf")
pairs (predictors, pch=16, diag.panel=panel.hist, lower.panel=panel.cor, panel=panel.smooth)
#dev.off()

rm (predictors)

# Maps of distribution  of functional indices ----------
coord  <-  FDIslands_Factors
rownames (coord) <- coord$Sitename
coord <- coord[,-1]
prj.coord <- project(cbind(coord$Longitude, coord$Latitude), proj = robin_crs)
coord <- cbind (prj.coord, coord)
names (coord)[1:2] <- c("X.prj", "Y.prj")

map_FRic  <- map_relative_richness (coord, FDIslands_Factors$FRic, "Functional richness", "FRic")
map_FEve  <- map_relative_richness (coord, FDIslands_Factors$FEve, "Functional evenness", "FEve")
map_FDiv  <- map_relative_richness (coord, FDIslands_Factors$FDiv, "Functional divergence", "FDiv")
map_FORed <- map_relative_richness (coord, FDIslands_Factors$F_OverRedundancy, "Functional overredundancy", "F_OverRedundancy")
map_FVuln <- map_relative_richness (coord, FDIslands_Factors$F_Vulnerability, "Functional vulnerability", "F_Vulnerability")

png("Output/Figures/Fig_1.png", pointsize=10, width=7400, height=3300, res=300)
ggarrange(map_FRic, map_FEve, map_FDiv, map_FORed, map_FVuln, align="hv", ncol=2, nrow=3,
          font.label = list(size=13, color="black", face="bold", family = "sans"))
dev.off()

# Correlation between species richness and functional indices ----------------------------

head (FDIslands_Factors)
Fig_2a <- Lin_Correlations (FDIslands_Factors, FDIslands_Factors$Nb_sp, 
                  FDIslands_Factors$FRic, FDIslands_Factors$Realm,
                  method_curve = "lm", formula_curve = y ~ splines::bs(x, 2), 
                  name.legend = "Realm",
                  break.scale = c("Atlantic", "East_Pacific", "Indian", "Pacific"), 
                  c("Atlantic", "Eastern Pacific", "Indian", "Pacific"),
                  color.pch = c("darkgreen","#00AFBB","#E7B800", "#FC4E07"),
                  ylabel = "FRichness", xlabel = "Species richness", 
                  axisname1=element_text(size=13, angle=90, family = "sans"),
                  axisname2 = element_blank(),
                  pos.leng="none", 
                  position.text.y="top")
Fig_2b <- Lin_Correlations (FDIslands_Factors, FDIslands_Factors$Nb_sp, 
                            FDIslands_Factors$FEve, FDIslands_Factors$Realm,
                            method_curve = "lm", formula_curve = y ~ splines::bs(x, 2), 
                            name.legend = "Realm",
                            break.scale = c("Atlantic", "East_Pacific", "Indian", "Pacific"), 
                            c("Atlantic", "Eastern Pacific", "Indian", "Pacific"),
                            color.pch = c("darkgreen","#00AFBB","#E7B800", "#FC4E07"),
                            ylabel = "FEvenness", xlabel = "Species richness", 
                            axisname1=element_text(size=13, angle=90, family = "sans"),
                            axisname2 = element_blank(),
                            pos.leng="none", 
                            position.text.y="top")
Fig_2c <- Lin_Correlations (FDIslands_Factors, FDIslands_Factors$Nb_sp, 
                            FDIslands_Factors$FDiv, FDIslands_Factors$Realm,
                            method_curve = "lm", formula_curve = y ~ splines::bs(x, 2), 
                            name.legend = "Realm",
                            break.scale = c("Atlantic", "East_Pacific", "Indian", "Pacific"), 
                            c("Atlantic", "Eastern Pacific", "Indian", "Pacific"),
                            color.pch = c("darkgreen","#00AFBB","#E7B800", "#FC4E07"),
                            ylabel = "FDivergence", xlabel = "Species richness", 
                            axisname1=element_text(size=13, angle=90, family = "sans"),
                            axisname2 = element_blank(),
                            pos.leng="none", 
                            position.text.y="top")
Fig_2d <- Lin_Correlations (FDIslands_Factors, FDIslands_Factors$Nb_sp, 
                            FDIslands_Factors$F_OverRedundancy, FDIslands_Factors$Realm,
                            method_curve = "lm", formula_curve = y ~ splines::bs(x, 2), 
                            name.legend = "Realm",
                            break.scale = c("Atlantic", "East_Pacific", "Indian", "Pacific"), 
                            c("Atlantic", "Eastern Pacific", "Indian", "Pacific"),
                            color.pch = c("darkgreen","#00AFBB","#E7B800", "#FC4E07"),
                            ylabel = "FOver-redundancy", xlabel = "Species richness", 
                            axisname1=element_text(size=13, angle=90, family = "sans"),
                            axisname2 = element_blank(),
                            pos.leng="none", 
                            position.text.y="top")
Fig_2e <- Lin_Correlations (FDIslands_Factors, FDIslands_Factors$Nb_sp, 
                            FDIslands_Factors$F_Vulnerability, FDIslands_Factors$Realm,
                            method_curve = "lm", formula_curve = y ~ splines::bs(x, 2), 
                            name.legend = "Realm",
                            break.scale = c("Atlantic", "East_Pacific", "Indian", "Pacific"), 
                            c("Atlantic", "Eastern Pacific", "Indian", "Pacific"),
                            color.pch = c("darkgreen","#00AFBB","#E7B800", "#FC4E07"),
                            ylabel = "FVulnerability", xlabel = "Species richness", 
                            axisname1=element_text(size=13, angle=90, family = "sans"),
                            axisname2 = element_blank(),
                            pos.leng="none", 
                            position.text.y="top")

pdf("Output/Figures/Fig_2.pdf", height = 5, width = 25, pointsize=30)
ggarrange(Fig_2a, Fig_2b, Fig_2c, Fig_2d, Fig_2e, align="hv", 
          ncol=5, nrow=1, legend = "none",
          font.label = list(size=13, color="black", 
                            face="bold", family = "sans"))
dev.off()


# Models ------------------------------------------------------------

# Exploring distribution variables 
distribution(FDIslands_Factors$FRic, "FRichness")
distribution(FDIslands_Factors$FEve, "FEveness")
distribution(FDIslands_Factors$FDiv, "FDivergence")
distribution(FDIslands_Factors$F_OverRedundancy, "FOverRedundancy")
distribution(FDIslands_Factors$F_Vulnerability, "FVulnerability")

par (mfrow=c(1,2))
plot (FDIslands_Factors$FRic); plot (asin(sqrt(FDIslands_Factors$FRic)))
plot (FDIslands_Factors$FDiv); plot (asin(sqrt(FDIslands_Factors$FDiv)))
plot (FDIslands_Factors$FEve); plot (asin(sqrt(FDIslands_Factors$FEve)))
plot (FDIslands_Factors$F_OverRedundancy); plot (asin(sqrt(FDIslands_Factors$F_OverRedundancy)))
plot (FDIslands_Factors$F_Vulnerability); plot (asin(sqrt(FDIslands_Factors$F_Vulnerability)))


# Build models - Functional richness  
vif (lm(FRic ~ Age_Ma_Scale + Age_Ma_Scale_2 +
          Current_Area_Scale + Current_Area_Scale_2 +
          Past_Area_Scale + Past_Area_Scale_2 +
          Isolation_Scale + Isolation_Scale_2,
        data=FDIslands_Factors))

# Single models

MFric_Age <- betareg(FRic ~ Age_Ma_Scale + Age_Ma_Scale_2, data=FDIslands_Factors)
summary (MFric_Age)

MFric_Area <- betareg(FRic ~ Current_Area_Scale + Current_Area_Scale_2, data=FDIslands_Factors)
summary (MFric_Area)

MFric_Isolation <- betareg(FRic ~ Isolation_Scale + Isolation_Scale_2, data=FDIslands_Factors)
summary (MFric_Isolation)

# Univariate models
# Current Area
MFric <- betareg(FRic ~ Current_Area_Scale,  data=FDIslands_Factors)
summary (MFric)
#Past Area
MFric_past <- betareg(FRic ~ Past_Area_Scale,  data=FDIslands_Factors)
summary (MFric_past)
AIC (MFric, MFric_past)

# Global model
MFric_CArea <- betareg(FRic ~ Age_Ma_Scale + 
                   Current_Area_Scale +
                   Isolation_Scale,
                 data=FDIslands_Factors)
summary (MFric_CArea)

MFric_PArea <- betareg(FRic ~ Age_Ma_Scale + 
                         Past_Area_Scale +
                         Isolation_Scale,
                       data=FDIslands_Factors)
summary (MFric_PArea)

plot_model(MFric_CArea)
Fig1A <- figure_models_current(MFric_CArea, "Functional richness",
                               Values = c(-0.5,0.5),
                               Axis.Inf = element_text(size=13, angle=0, family = "sans"))
Fig1B <- figure_models_past(MFric_PArea,Values = c(-0.5,0.3),
                            Axis.Inf = element_text(size=13, angle=0, family = "sans"))


# Build models - Functional evenness
# Single models
FEve_Age <- betareg(FEve ~ Age_Ma_Scale + Age_Ma_Scale_2,
                     data=FDIslands_Factors)
summary (FEve_Age)

FEve_Area <- betareg(FEve ~ Current_Area_Scale + Current_Area_Scale_2,
                      data=FDIslands_Factors)
summary (FEve_Area)

FEve_Isolation <- betareg(FEve ~ Isolation_Scale + Isolation_Scale_2,
                           data=FDIslands_Factors)
summary (FEve_Isolation)

# Univariate models
# Current Area
FEve <- betareg(FEve ~ Current_Area_Scale,  data=FDIslands_Factors)
summary (FEve)

#Past Area
FEve_past <- betareg(FEve ~ Past_Area_Scale,  data=FDIslands_Factors)
summary (FEve_past)
AIC (FEve, FEve_past)

MFEve_CArea <- betareg(FEve ~ Age_Ma_Scale + 
                         Current_Area_Scale +  
                          Isolation_Scale,
                   data=FDIslands_Factors)
summary (MFEve_CArea)

MFEve_PArea <- betareg(FEve ~ Age_Ma_Scale + 
                         Past_Area_Scale +
                         Isolation_Scale,
                       data=FDIslands_Factors)
summary (MFEve_PArea)
plot_model(MFEve_CArea)

Fig1C <- figure_models_current(MFEve_CArea, "Functional evenness",
                               Values = c(-0.5,0.5),
                               Axis.Inf = element_blank())
Fig1D <- figure_models_past(MFEve_PArea,
                               Values = c(-0.5,0.5),
                               Axis.Inf = element_blank())

# Build models - Functional divergence
# Single models
FDiv_Age <- betareg(FDiv ~ Age_Ma_Scale + Age_Ma_Scale_2,
                    data=FDIslands_Factors)
summary (FDiv_Age)

FDiv_Area <- betareg(FDiv ~ Current_Area_Scale + Current_Area_Scale_2,
                     data=FDIslands_Factors)
summary (FDiv_Area)

FDiv_Isolation <- betareg(FDiv ~ Isolation_Scale + Isolation_Scale_2,
                          data=FDIslands_Factors)
summary (FDiv_Isolation)

MFDiv_CArea <- betareg(FDiv ~ Age_Ma_Scale + 
                   Current_Area_Scale + 
                   Isolation_Scale,
                 data=FDIslands_Factors)
summary (MFDiv_CArea)

MFDiv_PArea <- betareg(FDiv ~ Age_Ma_Scale + 
                         Past_Area_Scale +
                         Isolation_Scale,
                       data=FDIslands_Factors)
summary (MFDiv_PArea)

Fig1E <- figure_models_current(MFDiv_CArea, "Functional divergence", 
                       Values=c(-0.1,0.1),
                       Axis.Inf = element_blank())
Fig1F <- figure_models_past(MFDiv_PArea,
                               Values=c(-0.1,0.1),
                               Axis.Inf = element_blank())

# Build models - Functional Over redundancy
F_OverRedundancy_Age <- betareg(F_OverRedundancy ~ Age_Ma_Scale +
                                  Age_Ma_Scale_2,
                    data=FDIslands_Factors)
summary (F_OverRedundancy_Age)

F_OverRedundancy_Area <- betareg(F_OverRedundancy ~ 
                                   Current_Area_Scale +
                                   Current_Area_Scale_2,
                     data=FDIslands_Factors)
summary (F_OverRedundancy_Area)

F_OverRedundancy_Isolation <- betareg(F_OverRedundancy ~ Isolation_Scale + 
                                        Isolation_Scale_2,
                          data=FDIslands_Factors)
summary (F_OverRedundancy_Isolation)


MFOverRed_CArea <- betareg(F_OverRedundancy ~ Age_Ma_Scale + 
                       Current_Area_Scale +
                   Isolation_Scale,
                 data=FDIslands_Factors)
summary (MFOverRed_CArea)

MFOverRed_PArea <- betareg(F_OverRedundancy ~ Age_Ma_Scale + 
                             Past_Area_Scale +
                             Isolation_Scale,
                           data=FDIslands_Factors)

Fig1G <- figure_models_current(MFOverRed_CArea, "Functional Overredundancy",
                            Values = c(-0.25,0.25),
                            Axis.Inf = element_blank())
Fig1H <- figure_models_past(MFOverRed_PArea,
                               Values = c(-0.25,0.25),
                               Axis.Inf = element_blank())

# Build models - Functional Vulnerability
F_Vulnerability_Age <- betareg(F_Vulnerability ~ Age_Ma_Scale + Age_Ma_Scale_2,
                                data=FDIslands_Factors)
summary (F_Vulnerability_Age)

F_Vulnerability_Area <- betareg(F_Vulnerability ~ Current_Area_Scale 
                                + Current_Area_Scale_2,
                                 data=FDIslands_Factors)
summary (F_Vulnerability_Area)

F_Vulnerability_Isolation <- betareg(F_Vulnerability ~ Isolation_Scale + 
                                       Isolation_Scale_2,
                                      data=FDIslands_Factors)
summary (F_Vulnerability_Isolation)

MFVul_CArea <- betareg(F_Vulnerability ~ Age_Ma_Scale + 
                   Current_Area_Scale +
                   Isolation_Scale,
                     data=FDIslands_Factors)
summary (MFVul_CArea)

MFVul_PArea <- betareg(F_Vulnerability ~ Age_Ma_Scale + 
                         Past_Area_Scale +
                         Isolation_Scale,
                       data=FDIslands_Factors)


Fig1I <- figure_models_current (MFVul_CArea, "Functional vulnerability", 
                                Values = c(-0.25,0.25),
                                Axis.Inf = element_blank())
Fig1J <- figure_models_past (MFVul_PArea,
                                Values = c(-0.25,0.30),
                                Axis.Inf = element_blank())

#write.table (FDIslands_Factors, "Output/Tables/FD_Islands_&_Factors.csv", sep = ";", dec = ",", row.names = F)

# Plotting all figures together -----------------------------------------
pdf("Output/Figures/Fig_3new.pdf", height = 10, width = 25, pointsize=30)
ggarrange(Fig1A, Fig1C, Fig1E, Fig1G, Fig1I,
          Fig1B, Fig1D, Fig1F, Fig1H, Fig1J,
          align="hv", ncol=5, 
          nrow=2, legend = "none",
          font.label = list(size=13, color="black", face="bold", family = "sans"))
dev.off()


# Moran Index ---------------------------------------------

geo_coords <- cbind (FDIslands_Factors$Longitude, FDIslands_Factors$Latitude)
dists      <- as.matrix (dist (cbind (FDIslands_Factors$Longitude, FDIslands_Factors$Latitude)))
dists.inv  <- 1/dists
diag(dists.inv) <- 0 
dists.inv[is.infinite(dists.inv)] <- 0

residual_FRic <- resid (MFric, type="pearson")
Moran.I (residual_FRic, dists.inv, scaled = T) 

residual_FEve <- resid (MFEve, type="pearson")
Moran.I (residual_FEve, dists.inv, scaled = T) 

residual_FDiv <- resid (MFDiv, type="pearson")
Moran.I (residual_FDiv, dists.inv, scaled = T) 

residual_FOverRed <- resid (MFOverRed, type="pearson")
Moran.I (residual_FOverRed, dists.inv, scaled = T) 

residual_FVul <- resid (MFVul, type="pearson")
Moran.I (residual_FVul, dists.inv, scaled = T)



correlogram_residual_FRic <- correlog (coords = geo_coords, z=residual_FRic)
plot (correlogram_residual_FRic)
correlogram_residual_FEve <- correlog (coords = geo_coords, z=residual_FEve)
plot(correlogram_residual_FEve)
correlogram_residual_FDiv <- correlog (coords = geo_coords, z=residual_FDiv)
plot(correlogram_residual_FDiv)
correlogram_residual_FOverRed <- correlog (coords = geo_coords, z=residual_FOverRed)
plot(correlogram_residual_FOverRed)
correlogram_residual_FVul <- correlog (coords = geo_coords, z=residual_FVul)
plot(correlogram_residual_FVul)




# Supplementary material ---------------
### Table S1
Table_S1 <- FDIslands_Factors [, c("Realm","Sitename","Sitecode","Latitude",
                                   "Longitude","Shallow_Reef_area_200mt",
                                   "Area_sub_past_200mt", "Isol_habitat", "Age_Ma")]
colnames (Table_S1) <- c("Realm","Island","Code", "Lat", "Long",
                         "C.area", "P.area","Isol","Age (Ma)")
Table_S1$Realm <- gsub ("_", " ", Table_S1$Realm)
Table_S1$Island <- gsub ("_", " ", Table_S1$Island)
Table_S1$Lat <- round (Table_S1$Lat, 2)
Table_S1$Long <- round (Table_S1$Long, 2)
Table_S1$C.area <- round (Table_S1$C.area, 2)
Table_S1$P.area <- round (Table_S1$P.area, 2)
Table_S1$Isol <- round (Table_S1$Isol, 2)
head (Table_S1)
#write.table(Table_S1, "Output/Tables/Island_features.csv", sep = ";", dec = ",", row.names = F)

# Table S2 
MFric_2 <- betareg(FRic ~ Age_Ma_Scale + Age_Ma_Scale_2 +
                   Past_Area_Scale + Past_Area_Scale_2 +
                   Current_Area_Scale + Current_Area_Scale_2 +
                   Isolation_Scale + Isolation_Scale_2,
                 data=FDIslands_Factors)
summary (MFric_2)
MFEve_2 <- betareg(FEve ~ Age_Ma_Scale + Age_Ma_Scale_2 +
                     Past_Area_Scale + Past_Area_Scale_2 +
                     Current_Area_Scale + Current_Area_Scale_2 +
                     Isolation_Scale + Isolation_Scale_2,
                   data=FDIslands_Factors)
summary (MFEve_2)

MFDiv_2 <- betareg(FDiv ~ Age_Ma_Scale + Age_Ma_Scale_2 +
                     Past_Area_Scale + Past_Area_Scale_2 +
                     Current_Area_Scale + Current_Area_Scale_2 +
                     Isolation_Scale + Isolation_Scale_2,
                   data=FDIslands_Factors)
summary (MFDiv_2)


MFOverRed_2 <- betareg(F_OverRedundancy ~ Age_Ma_Scale + Age_Ma_Scale_2 +
                     Past_Area_Scale + Past_Area_Scale_2 +
                     Current_Area_Scale + Current_Area_Scale_2 +
                     Isolation_Scale + Isolation_Scale_2,
                   data=FDIslands_Factors)
summary (MFOverRed_2)

MFVul_2 <- betareg(F_Vulnerability ~ Age_Ma_Scale + Age_Ma_Scale_2 +
                         Past_Area_Scale + Past_Area_Scale_2 +
                         Current_Area_Scale + Current_Area_Scale_2 +
                         Isolation_Scale + Isolation_Scale_2,
                       data=FDIslands_Factors)
summary (MFVul_2)
