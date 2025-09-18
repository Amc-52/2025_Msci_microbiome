#' # #  @title BacArena Analysis
#' @description Workflow used for creating and running BacArena simulations
#' @author A. J McCarry <https://github.com/Amc-52>


## Installing/loading packages

pkgs <- c("ggplot2", "ggrepel", "RColorBrewer", "devtools", "BacArena", "R.matlab", "parallel", "car", "dplyr")
for (pkg in pkgs) {
  if (pkg %in% rownames(installed.packages())) {
    eval(bquote(library(.(pkg))))
  }
}

## Set working directory to where models are stored
setwd("C:/Users/Ada Mc/Desktop/MSci_documents/JJS_scriptsandresources")

# directory of Genus models
genus_dir <- "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/JJS_scriptsandresources/Genus_Panmodel/Genus_Panmodel/"


# Load gut media (from gapseq docs)
gut <- read.csv("C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/JJS_scriptsandresources/mapped_gutIDs.txt", sep = ",", header = TRUE)
gutnosucr<- subset(gut, gut$name != "Sucrose")
# Load all available pan models for reference (everything in the model dirs except that .Rhistory file.)
modelnamesG <- list.files(path = genus_dir, pattern = "\\.mat$", full.names = TRUE)

# Load lists of bacteria per timepoint
# simulation using top 10 important genera in aGVHD grade from the RF model
top10_grade <- c(
  "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/JJS_scriptsandresources/Genus_Panmodel/Genus_Panmodel/panStreptococcus.mat",
  "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/JJS_scriptsandresources/Genus_Panmodel/Genus_Panmodel/panClostridium.mat",
  "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/JJS_scriptsandresources/Genus_Panmodel/Genus_Panmodel/panHungatella.mat",
  "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/JJS_scriptsandresources/Genus_Panmodel/Genus_Panmodel/panBacteroides.mat",
  "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/JJS_scriptsandresources/Genus_Panmodel/Genus_Panmodel/panEnterococcus.mat",
  "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/JJS_scriptsandresources/Genus_Panmodel/Genus_Panmodel/panParabacteroides.mat",
  "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/JJS_scriptsandresources/Genus_Panmodel/Genus_Panmodel/panLacticaseibacillus.mat",
  "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/JJS_scriptsandresources/Genus_Panmodel/Genus_Panmodel/panLatilactobacillus.mat",
  "C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/JJS_scriptsandresources/Genus_Panmodel/Genus_Panmodel/panPseudomonas.mat"
)
# FUNCTIONS
#' @title prepareModel
#' @description reads in .mat files in a loop and makes them into model objects, then into 
#' Bac objects (bacteria) to be added to the simulation
#' @param  model_path FILEPATH STRING | the filepath to the model
#' @return Bac objects suitable for adding to an arena
prepareModel <- function(model_path) {
  model <- readMATmod(model_path)
  Bac(model, setAllExInf = TRUE, cellweight_sd = 0)
}

#' @title bacSim1
#' @description builds and runs a BacArena simulation 
#' @param  modelfps FILEPATH STRING | filepaths to models (prepareModel is called within this function)
#' @param media CSV | csv file of substances and their fluxes from which the media envrionment of the simulation will be derived. Make sure the namespace of this file and the models match!
#' @param out_prefix STRING |prefix for the output files so that you can tell them apart if you're running multiple in a session
#' @return simulation object and variable substances of the simulation csv. Saved to users computer
bacSim1 <- function(modelfps, media, out_prefix) {
  #  generate panmodels
  panmodels <- lapply(modelfps, prepareModel)
  
  for (model in panmodels) {
    model@type <- model@model@mod_desc
  }
  
  #  Set up arena, add models one at a time so I can control the starting rel abundance of each  (cant simulate predation)
  arena <- Arena(n = 40, m = 40)
  for (i in 1:length(panmodels)) {
    arena <- addOrg(arena, panmodels[[i]], amount = 2)
    arena <- addEssentialMed(arena, panmodels[[i]])
    # removing sucrose!! (they seemed to like that a lot)
    arena<- changeSub(arena, 0, "EX_sucr(e)", unit = "mmol/cell")
  }
  
  
  # Add media to arena
  arena <- addSubs(arena, mediac = media$compounds, smax = media$maxFlux, unit = "mM")
  
  #  Run simulation
  sim <- simEnv(arena, time = 4, with_shadow = TRUE, sec_obj = "mtf") # Adjust time as needed, 4= 4 hours i think?
  
  #  Save simulation and variable substances table to working directory
  saveRDS(sim, file = paste0(out_prefix, "_Sim.rds"))
  df <- getVarSubs(sim)
  write.csv(df, file = paste0(out_prefix, "_Sim_varsubs.rds"))
  
  #  return simulation object for plotting
  return(sim)
}

# do 3 replicates and take the averages!
# CALL FUNCTIONS ----
simG_nosucR3 <- bacSim1(
  modelfps = top10_grade,
  media = gutnosucr, # with/without sucrose
  out_prefix = "simG_nosucR3"
)
#
# pull out cross-feeding csv!
crossfeed_tns <- findFeeding3(simG_grade_no_sucr3, time = 4, mets = names(getVarSubs(simG_grade_no_sucr3)), useNames = TRUE)

crossfeed_tns <- crossfeed_tns[[1]] # to get as dataframe

# making names more informative
crossfeed_t2$prod <- as.factor(crossfeed_t2$prod)
crossfeed_t2$prod <- recode_factor(crossfeed_t2$prod,
                                   model = "Streptococcus", model_1 = "Clostridium", model_2 = "Hungatella", model_3 = "Bacteroides", model_4 = "Enterococcus",
                                   model_5 = "Parabacteroides", model_6 = "Lacticaseibacillus", model_7 = "Latilactobacillus", model_8 = "Pseudomonas"
)

crossfeed_t2$cons <- as.factor(crossfeed_t2$cons)
crossfeed_t2$cons <- recode_factor(crossfeed_t2$cons,
                                   model = "Streptococcus", model_1 = "Clostridium", model_2 = "Hungatella", model_3 = "Bacteroides", model_4 = "Enterococcus",
                                   model_5 = "Parabacteroides", model_6 = "Lacticaseibacillus", model_7 = "Latilactobacillus", model_8 = "Pseudomonas"
)

# other plots that can show interesting things
plotSpecActivity(simG_grade_r2, rm_unused = T, var_nr = 10, useNames = TRUE)
plotGrowthCurve(simG_grade_r2, ret_data = T)[[2]]
plotSubCurve(simG_grade_r2)[[1]]
plotTotFlux(simG_grade_r2, legendpos = "none", num = 20) # plots the time course of reactions with high variation in activity for an Eval object.
plotAbundance(simG_grade_r2, use_biomass = T)
plotFluxVar(simlist = c(simG_grade_r2), metsel = top25)

##plot sub usage with sucrose ----
plotSubUsage(simlist = list(simG_grade_no_sucr, simG_grade_r3, simG_grade_no_sucr2 )) # to get the names of the top substances
SubUsage <- plotSubUsage(simlist = list(simG_grade_no_sucr, simG_grade_r3, simG_grade_no_sucr2), ret_data = T) # to get the names of the top substances
SubUsage$spec <- as.factor(SubUsage$spec)
SubUsage$spec <- recode_factor(SubUsage$spec,
                               model = "Streptococcus", model_1 = "Clostridium", model_2 = "Hungatella", model_3 = "Bacteroides", model_4 = "Enterococcus",
                               model_5 = "Parabacteroides", model_6 = "Lacticaseibacillus", model_7 = "Latilactobacillus", model_8 = "Pseudomonas"
)

# so this is a fun one. You need to load the dplyr package to run recode_factor(), and 
# then to run recode(), unload dplyr and load the package car.
SubUsage$sub <- recode(SubUsage$sub, "'EX_co2(e)'='C02'; 'EX_ac(e)'= 'Acetate'; 'EX_h2(e)'= 'H2'; 'EX_lac_D(e)'= 'D-lactate'; 'EX_lac_L(e)'= 'L-lactate'; 'EX_for(e)'= 'Formate'; 'EX_h2(e)'= 'Hydrogen'; 'EX_h2o(e)'= 'H20'; 'EX_na1(e)'= 'Sodium'; 'EX_sucr(e)'= 'Sucrose'")

# same as in rel abundance plots! for consistency
custom_colours <- c("Clostridium" = "cyan", "Streptococcus" = "green", "Enterococcus" = "red", "Hungatella" = "purple4", "Bacteroides" = "blue", "Parabacteroides" = "yellow", "Lacticaseibacillus" = "darkblue", "Latilactobacillus" = "pink", "Pseudomonas" ="purple")
g1 <- ggplot(SubUsage, aes(sub, mflux)) +
  geom_boxplot(outlier.colour = "black", varwidth = F, aes(colour = spec)) +
  labs(title = "Usage of highly variable substances", colour = "Genus") +
  xlab("Exchange metabolites") +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.1, linetype = "dashed") +
  ylim(-80, 150) +
  scale_color_manual(values = custom_colours) +
  gg_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x=element_blank())

g1

## plot sub usage no sucrose
plotSubUsage(simlist = list(simG_nosucR2, simG_nosucR3, simG_grade_no_sucr3)) # to get the names of the top substances
SubUsage1 <- plotSubUsage(simlist = list(simG_nosucR2, simG_nosucR3, simG_grade_no_sucr3), ret_data = T) # to get the names of the top substances
SubUsage1$spec <- as.factor(SubUsage1$spec)
SubUsage1$spec <- recode_factor(SubUsage1$spec,
                                model = "Streptococcus", model_1 = "Clostridium", model_2 = "Hungatella", model_3 = "Bacteroides", model_4 = "Enterococcus",
                                model_5 = "Parabacteroides", model_6 = "Lacticaseibacillus", model_7 = "Latilactobacillus", model_8 = "Pseudomonas"
)

SubUsage1$sub <- recode(SubUsage$sub, "'EX_co2(e)'='C02'; 'EX_ac(e)'= 'Acetate'; 'EX_h2(e)'= 'H2'; 'EX_lac_D(e)'= 'D-lactate'; 'EX_lac_L(e)'= 'L-lactate'; 'EX_fru(e)'= 'Fructose'; 'EX_h2(e)'= 'Hydrogen'; 'EX_h2o(e)'= 'H20'; 'EX_na1(e)'= 'Sodium'; 'EX_malt(e)'= 'Maltose'")


g <- ggplot(SubUsage1, aes(sub, mflux)) +
  geom_boxplot(outlier.colour = "black", varwidth = F, aes(colour = spec)) +
  labs(colour = "Genus") +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.1, linetype = "dashed") +
  ylim(-80, 150) +
  scale_color_manual(values = custom_colours) +
  gg_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x=element_blank(),
    axis.title.x=element_blank())

g
g<- g + theme(axis.title.x=element_blank())
# combine plots:
wrap_plots(g1,g, ncol = 1) +
  plot_layout(guides = "collect")

# trying to tie reaction mflux to bacteria- doesn't really work
exchanges <- data.frame(simG_grade_r2@exchanges)
exchanges$species <- as.factor(exchanges$species)
exchanges$species <- recode_factor(exchanges$species,
                                   model = "Streptococcus", model_1 = "Clostridium", model_2 = "Hungatella", model_3 = "Bacteroides", model_4 = "Enterococcus",
                                   model_5 = "Parabacteroides", model_6 = "Lacticaseibacillus", model_7 = "Latilactobacillus", model_8 = "Pseudomonas"
)
exchangesNoSpecies <- exchanges[, -c(1)]
top25_ex <- names(sort(colSums(exchangesNoSpecies), decreasing = TRUE))[1:50]
BiocManager::install("reprtree")
save.image()

# making a heatmap for the reactions ----

# getting reactivity data
reactivity <- plotReaActivity(simG_grade_r2, spec_list = c("model", "model_1", "model_2", "model_3", "model_4", "model_5", "model_6", "model_7", "model_8"), reactions = simG_grade_r2@models[[1]]@react_id)

# making the part of the plot with the data into a df
reactivity <- data.frame(reactivity[[2]][["data"]])
# subsetting data (cut down like 90% of the reactions)
reactivity <- subset(reactivity, reactivity$mflux > 0.5)
# fixing species names
reactivity$spec <- as.factor(reactivity$spec)
reactivity$spec <- recode_factor(reactivity$spec,
                                 model = "Streptococcus", model_1 = "Clostridium", model_2 = "Hungatella", model_3 = "Bacteroides", model_4 = "Enterococcus",
                                 model_5 = "Parabacteroides", model_6 = "Lacticaseibacillus", model_7 = "Latilactobacillus", model_8 = "Pseudomonas"
)

# sorting the data (not necessary, just for my own interest)
reactivity$mflux <- sort(reactivity$mflux, decreasing = T)

# pivoting the df wider to coerce it into the correct format for pheatmap
reacWide <- reactivity %>%
  pivot_wider(
    names_from = spec,
    values_from = mflux
  )
# replacing NAs with 0
reacWide <- reacWide %>% replace(is.na(.), 0)
# removing uneccessary columns (time and replicate)
reacWide <- reacWide[, -c(2, 3)]

# R seems to think some of the reactions are non-unique
# all credit and thank-yous to stack overflow user Thomas Mailand for this function! https://stackoverflow.com/questions/49381670/making-non-unique-row-names-unique
uniqify_names <- function(names_vector) {
  names <- unique(names_vector)
  count_table <- rep(0, length(names))
  names(count_table) <- names # works because R has weird symbol lookup
  update_name <- function(name) {
    new_name <- paste0(name, ".", count_table[name])
    count_table[name] <<- count_table[name] + 1
    new_name
  }
  vapply(names_vector, update_name, FUN.VALUE = "character")
}

reacWide_subset <- subset(reacWide, rowSums(reacWide[, c(2:8)]) > 100)

Rowreacs <- uniqify_names(reacWide_subset$rea)
reacWide_subset <- reacWide_subset[, -c(1)]

rownames(reacWide_subset) <- Rowreacs

# normalising the data with z scores a la gene data. z score is a measure of how many
# standard deviations a value is away from the mean of that variable, obvs Strep is
# moving the goalposts a bit but you can see from this plot which reactions are used more/less often by certain bacteria than others
# function from https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
cal_z_score <- function(x) {
  (x - mean(x)) / sd(x)
}

reacWide_subset_norm <- t(apply(reacWide_subset, 1, cal_z_score))

pheatmap(reacWide_subset_norm,
         cluster_rows = T, clustering_method = "ward.D", clustering_distance_rows = "canberra",
         color = viridis(100),
         fontsize_row = 7, main = "simulated maxflux of bacterial reactions"
)

Strep_reactions$mflux <- sort(Strep_reactions$mflux, decreasing = F)
plotReaActivity(simG_grade_r2, spec_list = c("model"), reactions = top20_pos)

top25 <- reactivity$rea[1:25]
