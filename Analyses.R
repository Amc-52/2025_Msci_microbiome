#' # #  @title Analysis
#' @description Collection of analyses performed on the phyloseq/tse objects produced with dadaPipe(), and files produced by MATLAB
#' @author A. J McCarry <https://github.com/Amc-52>

# setup----

# Install packages
# installing and loading bacarena ----
# Installing bacarena for the first time can be a little tricky! What you have to do to get BacArena  to work is download
# the latest version of libSBML and follow the install instructions in the readme. Then install sybil, then
# install BacArena. Otherwise installs will fail because none of these packages are on CRAN.


# after installing manually from sourceforge:
library(libSBML) 

install.packages("remotes")
remotes::install_github("SysBioChalmers/sybil")

library(sybil)
library(devtools)
install_github("euba/bacarena")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#remotes::install_github("jeffreyevans/rfUtilities") # doesnt work, will have to load individual functions

install.packages(c("ggplot2", "devtools", "doParallel", "vegan", "dplyr", "sqldf", "extrafontdb", "Rttf2pt1", "extrafont"))
BiocManager::install(c("dada2", "phyloseq", "phyloseqGraphTest", "Biostrings", "decontam", "MBEC", "microbiome", "caret",
                       "scater", "randomForest", "optRF","BacArena", "R.matlab", "parallel", "mia", "microbiome", "DESeq2"
))

# Load packages
pkgs <- c(
  "devtools", "DECIPHER", "dplyr", "ggplot2", "ggrepel", "ggtext", "dada2", "phyloseq", "Biostrings",
  "doParallel", "vegan", "MBECS", "mia", "microbiome", "RColorBrewer", "sqldf", "patchwork", "caret",
  "randomForest", "optRF", "extrafontdb", "Rttf2pt1", "extrafont", "scater", "BacArena", "R.matlab", "parallel",
  "DESeq2", "viridis", "phyloseqGraphTest", "pairwiseAdonis"
)

for (pkg in pkgs) {
  if (pkg %in% rownames(installed.packages())) {
    eval(bquote(library(.(pkg))))  }
}



# register fonts:
font_import(paths = c("C:/Users/Ada Mc/AppData/Local/Microsoft/Windows/Fonts/"), pattern = ".ttf")
# load them for use:
loadfonts()

# general plot theme
gg_theme <- theme(
  text = element_text(family = "Cambria", size = 12, colour = "grey15"),
  plot.background = element_rect(color = "black", fill = "white", linewidth = 1, linetype = "solid"),
  panel.border = element_rect(fill = NA, colour = "grey10", linewidth = 1),
  axis.title = element_text(size = 12, face = "bold"),
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5),
  legend.title = element_text(face = "bold"),
  legend.key = element_blank(),
  legend.background = element_rect(color = "grey10", fill = NA, linewidth = 0.5, linetype = "solid")
)

# Building TSE/ MAE objects ----

# Load phyloseq objects and mgPipe outputs (phyloseq objects are in a .Rdata file)
load("C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci/DADA2_resultsFINAL/phyloseqObjs.Rdata")
load("C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci/FINAL_msci2025/finalf.Rdata")
reacAb <- read.csv("C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/RWork/mgPipeResults/ReactionAbundance.csv", header = TRUE)
subAb <- read.csv("C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/RWork/mgPipeResults/SubsystemAbundance.csv", header = TRUE)
dietUptake <- read.csv("C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/RWork/mgPipeResults/BRATdiet/inputDiet_net_uptake_fluxes.csv", header = TRUE)
dietSecretion <- read.csv("C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci_documents/RWork/mgPipeResults/BRATdiet/inputDiet_net_secretion_fluxes.csv", header = TRUE)
metadata <- read.csv("C:/Users/Ada Mc/OneDrive - University of Glasgow/MSci/FINAL_msci2025/inputData/I21_metadata.csv", header = T)

# process the mgPipe outputs (had to give samples non-numerical IDs for mgPipe to run; now i have to restore the OG names)
# strip 'sp' from column names
dflist <- list(reacAb, subAb, dietSecretion, dietUptake)
dflist

spStrip <- function(dflist) {
  for (df in dflist) {
    colnames(df) <- gsub("sp", "", colnames(df))
  }
  return(dflist)
}
spStrip(dflist)

colnames(reacAb) <- gsub("sp", "", colnames(reacAb))
colnames(subAb) <- gsub("sp", "", colnames(subAb))

colnames(dietSecretion) <- gsub("sp", "", colnames(dietSecretion))
colnames(dietUptake) <- gsub("sp", "", colnames(dietUptake))

#rownames(modelStats) <- modelStats[, 1]
#rownames(modelStats) <- gsub("sp", "", rownames(modelStats))


colnames(reacAb) <- gsub("X", "Reaction", colnames(reacAb)) # changing to be more informative
colnames(subAb) <- gsub("X", "Reaction", colnames(subAb)) # changing to be more informative


# get/wrangle metadata and other tables for phyloseq TSE

# metadata: sample IDs as row names
rownames(metadata) <- metadata[, 1]
samples <- data.frame(phy21DS@sam_data)
counts <- t(phy21DS@otu_table) # feature IDs (ASVs) as rows, sample IDs as col headers

counts <- as.matrix(counts)
counts <- na.omit(counts)

taxa <- as.data.frame(phy21DS@tax_table) # feature IDs (ASVs) as rows, taxonomic identities in columns

# need to remove taxa that couldn't be mapped so that the mgpipe assay nrow and ncol match the phyloseq ones!
# remove non-modelled taxa from tse
# damn, over half don't have models? I guess I also cut those whose name spaces didn't match and where there were species-level assignments

### 'AGORAtaxa': list of modelled taxa (LONG) ----
AGORAtaxa <- c(
  "Abiotrophia",
  "Achromobacter",
  "Acinetobacter",
  "Actinobacillus",
  "Actinomyces",
  "Adlercreutzia",
  "Afipia",
  "Agathobaculum",
  "Agrobacterium",
  "Akkermansia",
  "Alistipes",
  "Alloprevotella",
  "Amedibacillus",
  "Anaerobutyricum",
  "Anaerococcus",
  "Anaerofustis",
  "Anaerostipes",
  "Anaerotruncus",
  "Bacillus",
  "Bacteroides",
  "Barnesiella",
  "Bifidobacterium",
  "Bilophila",
  "Blautia",
  "Bradyrhizobium",
  "Brevibacterium",
  "Brevundimonas",
  "Brucella",
  "Butyricimonas",
  "Butyrivibrio",
  "Campylobacter",
  "Capnocytophaga",
  "Catenibacterium",
  "Cellulosilyticum",
  "Citrobacter",
  "Clostridioides",
  "Clostridium",
  "Collinsella",
  "Comamonas",
  "Coprobacillus",
  "Coprococcus",
  "Corynebacterium",
  "Cutibacterium",
  "Delftia",
  "Desulfovibrio",
  "Dialister",
  "Dorea",
  "Dysgonomonas",
  "Eggerthella",
  "Enterobacter",
  "Enterocloster",
  "Enterococcus",
  "Escherichia",
  "Eubacterium",
  "Faecalibacterium",
  "Faecalitalea",
  "Finegoldia",
  "Flavonifractor",
  "Fusobacterium",
  "Gemella",
  "Gordonibacter",
  "Granulicatella",
  "Haemophilus",
  "Helicobacter",
  "Holdemanella",
  "Holdemania",
  "Hungatella",
  "Intestinibacter",
  "Kingella",
  "Klebsiella",
  "Kocuria",
  "Lachnoclostridium",
  "Lachnospira",
  "Lacticaseibacillus",
  "Lactobacillus",
  "Lactococcus",
  "Lancefieldella",
  "Latilactobacillus",
  "Lautropia",
  "Leptotrichia",
  "Leuconostoc",
  "Levilactobacillus",
  "Ligilactobacillus",
  "Limosilactobacillus",
  "Listeria",
  "Marvinbryantia",
  "Mediterraneibacter",
  "Megamonas",
  "Megasphaera",
  "Methanobrevibacter",
  "Methylobacterium",
  "Methylorubrum",
  "Methyloversatilis",
  "Microbacterium",
  "Micrococcus",
  "Morganella",
  "Neisseria",
  "Odoribacter",
  "Olsenella",
  "Oribacterium",
  "Oxalobacter",
  "Paenibacillus",
  "Parabacteroides",
  "Paracoccus",
  "Paraprevotella",
  "Parasutterella",
  "Pediococcus",
  "Peptoniphilus",
  "Phascolarctobacterium",
  "Porphyromonas",
  "Prevotella",
  "Proteus",
  "Pseudoflavonifractor",
  "Pseudomonas",
  "Rhodococcus",
  "Roseburia",
  "Roseomonas",
  "Rothia",
  "Ruminiclostridium",
  "Ruminococcus",
  "Salmonella",
  "Scardovia",
  "Schaalia",
  "Schlesneria",
  "Selenomonas",
  "Solobacterium",
  "Sphingobacterium",
  "Staphylococcus",
  "Stenotrophomonas",
  "Streptococcus",
  "Subdoligranulum",
  "Succinivibrio",
  "Sutterella",
  "Terrisporobacter",
  "Trabulsiella",
  "Tyzzerella",
  "Veillonella",
  "Weissella",
  "Yersinia"
)

### end of list ----
# subset taxa by genera for which AGORA2 models are available
taxa <- subset(taxa, taxa$Genus %in% AGORAtaxa)
# do the same for counts
counts<- subset(counts, rownames(counts) %in% rownames(taxa))

# omit these ASVs which still remain (for some reason)
omit<- c( "ASV27",  "ASV64",  "ASV94",  "ASV103", "ASV146", "ASV147", "ASV162", "ASV164", "ASV184", "ASV218", "ASV334",
          "ASV419", "ASV437", "ASV562", "ASV580", "ASV592", "ASV604", "ASV969")
taxa <- subset(taxa, !rownames(taxa) %in% omit)

# check that the row and column names match one-to-one between counts, samples, and taxa:
print(setdiff(union(rownames(samples), colnames(counts)), intersect(rownames(samples), colnames(counts))))
print(setdiff(union(rownames(taxa), rownames(counts)), intersect(rownames(taxa), rownames(counts))))

# 1 sample ID is missing from counts, namely 7. This sample may have been removed due to low abundance?
# remove 7 from samples
samples <- subset(samples, Sample_ID != "7")
# add aGVHD yes no to metadata
samples$aGVHD <- as.factor(samples$aGVHD)

# Match rows and columns (gets the rownames in the same order in taxa and samples!)
counts <- counts[rownames(taxa), rownames(samples)]

# check again that the row and column names match one-to-one between counts, samples, and taxa:
all(rownames(samples) == colnames(counts))
# TRUE

# wrangle mgPipe files so I can add them as Alt Experiments to the TSE

rownames(dietUptake) <- dietUptake[, 1]
dietUptake <- dietUptake[, -1]

rownames(dietSecretion) <- dietSecretion[, 1]
dietSecretion <- dietSecretion[, -1]

rownames(reacAb) <- reacAb[, 1]
reacAb <- reacAb[, -1]

rownames(subAb) <- subAb[, 1]
subAb <- subAb[, -1]

# match order of sample IDs
reacAb <- reacAb[rownames(samples)]
subAb <- subAb[rownames(samples)]
dietSecretion <- dietSecretion[rownames(samples)]
dietUptake <- dietUptake[rownames(samples)]

# convert into SummarisedExperiment objects
reacAbSE <- SummarizedExperiment(reacAb)
subAbSE <- SummarizedExperiment(subAb)
dietSecretionSE <- SummarizedExperiment(dietSecretion)
dietUptakeSE <- SummarizedExperiment(dietUptake)

# Subset phyloseq object into 1 object per patient

# Subset mgPipe outputs into 1 per patient also

# merge into TSEs
ALL_tse <- TreeSummarizedExperiment(
  assays = SimpleList(counts = counts),
  colData = DataFrame(samples),
  rowData = DataFrame(taxa)
)
# ADD ALT EXPERIMENTS (need to be SE objects with same col names as the assay)
altExp(ALL_tse, "reactionAbundances") <- reacAbSE
altExp(ALL_tse, "subsystemAbundances") <- subAbSE
altExp(ALL_tse, "dietUptake") <- dietUptakeSE
altExp(ALL_tse, "dietSecretion") <- dietSecretionSE


# Alpha and beta diversity ----
ALL_tse <- transformAssay(ALL_tse, assay.type = "counts", method = "relabundance")
ALL_tse <- transformAssay(ALL_tse, assay.type = "counts", method = "rclr")

ALL_tse <- addAlpha(
  ALL_tse,
  assay.type = "relabundance",
  detection = 10
)
# beta diversity
ALL_tse <- addMDS(
  ALL_tse,
  FUN = getDissimilarity,
  method = "bray",
  assay.type = "relabundance",
  name = "MDS_bray"
)


p <- plotReducedDim(ALL_tse, "MDS_bray", colour_by = "max_aGVHD_grade", shape_by = "aGVHD", size_by = "shannon_diversity")
# Calculate explained variance


# removing the hashtagged lined will give you a PCoA plot
p <- plotReducedDim(ALL_tse_no24, "NMDS_aitchison", colour_by = "PatientID", size_by = "max_aGVHD_grade", shape_by = "aGVHD")
#e <- attr(reducedDim(ALL_tse, "MDS_bray"), "eig")
#rel_eig <- e / sum(e[e > 0])
# Add explained variance for each axis
p <- p + labs(
  # x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
  #y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = ""),
  title = "NMDS_aitchisoN of aGVHD incidence") + 
  gg_theme + 
  scale_colour_manual(values = palette2) +
  scale_size_binned()
p
wrap_plots(p,p1) +  plot_layout(guides = "collect")

# PERMANOVA TEST TO SEE IF MAX_AGVHD HAS A SIGNIFICANT EFFECT ON ASV ABUNDANCE ----
# since its followup on the MDS plot, i'm looking for significant effect in the timepoints
# i modelled in the MDS. Also so the prediction isn't 'diluted' by the sheer number of ASVs,
# i'm taking the first 250 most abundant overall as a hopefully representative crossection.

phy21DS_pma<- subset_samples(phy21DS, phy21DS@sam_data$timepoint %in% c("d0","w1","w2","w3","m1"))
phy21DS_pma_otu<-phy21DS_pma@otu_table
phy21DS_pma_otu<-data.frame(phy21DS_pma_otu)
phy21DS_pma_otu <- na.omit(phy21DS_pma_otu)

top98 <- names(sort(colSums(phy21DS_pma_otu), decreasing = TRUE))[1:98] # after 98 tis all 0S
phy21DS_pma <- prune_taxa(top98, phy21DS)
phy21DS_pma<- subset_samples(phy21DS_pma, phy21DS_pma@sam_data$timepoint %in% c("d0","w1","w2","w3","m1"))
phy21DS_pma<- subset_samples(phy21DS_pma, phy21DS_pma@sam_data$PatientID != "P24")

phy21DS_pma_sam<- phy21DS_pma@sam_data
phy21DS_pma_sam<- data.frame(phy21DS_pma_sam)
phy21DS_pma_otu<-phy21DS_pma@otu_table

# combining aGVHD incidence and grade variables

adonis2(phy21DS_pma_otu~max_aGVHD_grade,data=phy21DS_pma_sam, permutations=9999, method="bray")

# post-hoc validation of PERMANOVA
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

pairwise.adonis(phy21DS_pma_otu, phy21DS_pma_sam$max_aGVHD_grade, sim.method = "bray",
                p.adjust.m = "bonferroni")

# more beta diversity methods to compare
ALL_tse <- addNMDS(
  ALL_tse,
  FUN = getDissimilarity,
  method = "euclidean",
  assay.type = "rclr",
  name = "NMDS_aitchison"
)


ALL_tse <- addNMDS(
  ALL_tse,
  FUN = getDissimilarity,
  method = "bray",
  assay.type = "relabundance",
  name = "NMDS_bray"
)


# plotting multiple beta diversity methods
plots <- lapply(
  c("MDS_bray", "NMDS_bray", "MDS_aitchison","NMDS_aitchison"),
  plotReducedDim,
  object = ALL_tse,
  colour_by = "PatientID"
)

wrap_plots(plots) +
  plot_layout(guides = "collect")

# MDS bray seems to partition out samples/patients the best

## alpha diversity and effective number of species ----
plotColData(
  ALL_tse,
  "shannon_diversity",
  "timepoint",
  colour_by = "aGVHD"
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "timepoint", y = "Shannon")

# plotting richness- altering an internal phyloseq function to plot a boxplot
# instead of a scatter graph
# geom_violin(scale= "width",alpha=0.9) + geom_boxplot(width= 0.2, linewidth=0.2) extra bit for a violin plot
plot_richnessBox <- function(physeq, x = "samples", color = NULL, shape = NULL, title = NULL,
                             scales = "free_y", nrow = 1, shsi = NULL, measures = NULL, sortby = NULL) {
  # Calculate the relevant alpha-diversity measures
  erDF <- estimate_richness(physeq, split = TRUE, measures = measures)
  # Measures may have been renamed in `erDF`. Replace it with the name from erDF
  measures <- colnames(erDF)
  # Define "measure" variables and s.e. labels, for melting.
  ses <- colnames(erDF)[grep("^se\\.", colnames(erDF))]
  # Remove any S.E. from `measures`
  measures <- measures[!measures %in% ses]
  # Make the plotting data.frame.
  # This coerces to data.frame, required for reliable output from reshape2::melt()
  if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    # Include the sample data, if it is there.
    DF <- data.frame(erDF, sample_data(physeq))
  } else {
    # If no sample data, leave it out.
    DF <- data.frame(erDF)
  }
  if (!"samples" %in% colnames(DF)) {
    # If there is no "samples" variable in DF, add it
    DF$samples <- sample_names(physeq)
  }
  # sample_names used to be default, and should also work.
  # #backwardcompatibility
  if (!is.null(x)) {
    if (x %in% c("sample", "samples", "sample_names", "sample.names")) {
      x <- "samples"
    }
  } else {
    # If x was NULL for some reason, set it to "samples"
    x <- "samples"
  }
  # melt to display different alpha-measures separately
  mdf <- reshape2::melt(DF, measure.vars = measures)
  # Initialize the se column. Helpful even if not used.
  mdf$se <- NA_integer_
  if (length(ses) > 0) {
    ## Merge s.e. into one "se" column
    # Define conversion vector, `selabs`
    selabs <- ses
    # Trim the "se." from the names
    names(selabs) <- substr(selabs, 4, 100)
    # Make first letter of selabs' names uppercase
    substr(names(selabs), 1, 1) <- toupper(substr(names(selabs), 1, 1))
    # use selabs conversion vector to process `mdf`
    mdf$wse <- sapply(as.character(mdf$variable), function(i, selabs) {
      selabs[i]
    }, selabs)
    for (i in 1:nrow(mdf)) {
      if (!is.na(mdf[i, "wse"])) {
        mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
      }
    }
    # prune the redundant columns
    mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
  }
  ## Interpret measures
  # If not provided (default), keep all
  if (!is.null(measures)) {
    if (any(measures %in% as.character(mdf$variable))) {
      # If any measures were in mdf, then subset to just those.
      mdf <- mdf[as.character(mdf$variable) %in% measures, ]
    } else {
      # Else, print warning about bad option choice for measures, keeping all.
      warning("Argument to `measures` not supported. All alpha-diversity measures (should be) included in plot.")
    }
  }
  if (!is.null(shsi)) {
    # Deprecated:
    # If shsi is anything but NULL, print a warning about its being deprecated
    warning("shsi no longer supported option in plot_richness. Please use `measures` instead")
  }
  # Address `sortby` argument
  if (!is.null(sortby)) {
    if (!all(sortby %in% levels(mdf$variable))) {
      warning("`sortby` argument not among `measures`. Ignored.")
    }
    if (!is.discrete(mdf[, x])) {
      warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
    }
    if (all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[, x])) {
      # Replace x-factor with same factor that has levels re-ordered according to `sortby`
      wh.sortby <- which(mdf$variable %in% sortby)
      mdf[, x] <- factor(mdf[, x],
                         levels = names(sort(tapply(
                           X = mdf[wh.sortby, "value"],
                           INDEX = mdf[wh.sortby, x],
                           mean,
                           na.rm = TRUE, simplify = TRUE
                         )))
      )
    }
  }
  # Define variable mapping
  richness_map <- aes_string(x = x, y = "value", colour = color, shape = shape)
  # Make the ggplot.
  p <- ggplot(mdf, richness_map) +
    geom_boxplot(varwidth = F) 
  # Add error bars if mdf$se is not all NA
  if (any(!is.na(mdf[, "se"]))) {
    p <- p + geom_errorbar(aes(ymax = value + se, ymin = value - se), width = 0.1)
  }
  # Rotate horizontal axis labels, and adjust
  p <- p + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))
  # Add y-label
  p <- p + ylab("Alpha Diversity Measure")
  # Facet wrap using user-options
  p <- p + facet_wrap(~variable, nrow = nrow, scales = scales)
  # Optionally add a title to the plot
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


phy21DS@sam_data$aGVHD<- as.factor(phy21DS@sam_data$aGVHD)
phy21DS@sam_data$aGVHD<- recode_factor(phy21DS@sam_data$aGVHD, "0" = "no", "1"= "yes")
phy21DS@sam_data$aGVHD <- with(phy21DS@sam_data, ifelse(max_aGVHD_grade > 0, "yes", "no"))
phy21DS@sam_data$timepoint_long <- factor(phy21DS@sam_data$timepoint)
phy21DS@sam_data$timepoint_long <- recode_factor(phy21DS@sam_data$timepoint_long,
                                                 pre= "pre-exam", conditioning= "conditioning", d0= "day_of_HSCT", w1= "week_1", w2= "week_2",
                                                 w3= "week_3", m1= "month_1", m3= "month_3", m6= "month_6", m12= "month_12")

# plot richness
p <- plot_richnessBox(phy21DS, x = "timepoint_long", color = "aGVHD", measures = c("Shannon", "Observed", "InvSimpson"))
p <- p + scale_colour_manual(values = c("no" = "#108888", "yes" = "darkred")) + gg_theme + labs(x="Timepoint")
p
# + labs(title = "Alpha Diversity over time by aGVHD incidence")

# Network Analysis ----
# another way of visualising how similar patients are to each other
# based on ASV abundance, and how well patient samples cluster!
set.seed(4) # there's randomness in how the plot looks, so set seed for consistency
phy21DS@sam_data$max_aGVHD_grade<- as.factor(phy21DS@sam_data$max_aGVHD_grade)

# aGVHD cohort were developing the condition
phy21DS_no24<- subset_samples(phy21DS, phy21DS@sam_data$"timepoint" %in% c("d0","w1","w2","w3","m1"))
phy21DS_no24<- subset_samples(phy21DS_no24, phy21DS@sam_data$"PatientID" != "P24")

plot_net(phy21DS, distance = "bray", color = "PatientID", shape = "aGVHD",  title = "Network plot") +
  gg_theme + scale_color_manual(values = palette2) + labs(title = "Network Plot")

gt <- graph_perm_test(phy21DS, "PatientID",
                      distance = "jaccard",
                      type = "mst", nperm = 1000
)
plot_test_network(gt) + gg_theme + scale_color_manual(values = palette2)
gt
# Output from graph_perm_test

#   Observed test statistic: 43 pure edges
# 97 total edges in the graph
# Permutation p-value: 0.000999000999000999

gt <- graph_perm_test(phy21DS, "aGVHD",
                      grouping = "patientID",
                      distance = "jaccard", type = "mst", nperm = 1000
)


# Relative abundance ----
# Plotting relative abundances for each timepoint in 2 batches: no aGVHD and yes aGVHD. Tried to make function, does not work, idk why

# this is the 'plot_bar' function from phyloseq. I edited it so it would make a box plot instead
plot_box <- function(physeq, x = "Sample", y = "Abundance", fill = NULL,
                     title = NULL, facet_grid = NULL) {
  # Start by melting the data in the "standard" way using psmelt.
  mdf <- psmelt(physeq)
  
  # Build the plot data structure
  p <- ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  
  # Add the bar geometric object. Creates a basic graphic. Basis for the rest.
  # Test weather additional
  p <- p + geom_boxplot(outlier.colour = "black")
  
  # By default, rotate the x-axis labels (they might be long)
  p <- p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  
  # Add faceting, if given
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  
  # Optionally add a title to the plot
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

# PLOTTING top taxa in each patient with aGVHD at the time they developed it
phy_deSamples <- subset_samples(phy21DS, phy21DS@sam_data$aGVHD == "no")
# or (can't just put aGVHD= 1 bc have to remove P24)
phy_deSamples <- subset_samples(phy21DS, phy21DS@sam_data$PatientID == "P11")
#,"P03","P09","P12","P17","P20"
#phy_deSamples <- subset_samples(phy21DS, phy21DS@sam_data$Sample_ID == "211")

phy_deSamples <- subset_samples(phy_deSamples, phy_deSamples@sam_data$timepoint == "w1")

# transform compositional
phy_deSamples <- transform_sample_counts(phy_deSamples, function(x) x * 100 / sum(x))
phy_deSamples <- transform_sample_counts(phy_deSamples,function(x) (x - mean(x)) / sd(x))

#phy_deSamples@sam_data$timepoint <- factor(phy_deSamples@sam_data$timepoint, levels = c("pre", "conditioning", "d0", "w1", "w2", "w3", "w4", "w5", "m1", "m3", "m6", "m12"))
REL_MERGEDotu <- as.data.frame(phy_deSamples@otu_table)
REL_MERGEDotuNA <- na.omit(REL_MERGEDotu)
top25 <- names(sort(colSums(REL_MERGEDotuNA), decreasing = TRUE))[1:20]
ps.top25 <- prune_taxa(top25, phy_deSamples)
zero<-plot_bar(ps.top25, x = "timepoint", fill = "Genus") +
  labs(title = "Top Genera in patient 11, week 1 (non-aGVHD)", x = "Genus", y = "Abundance (compositional)") +
  scale_fill_genus() +
  gg_theme 
zero
# subtitle = "T.O.O at +13d; sampled at +21d"
# theme(legend.position = "none")
a22
a1 # patient 2
a11 # p2 alt
a2 # patient 3
a22 #p3 alt
a3 # P9
a4 # p12
a44 # patient 12 alt

a6 # p20
# pre_a0, con_a0, d0_a0, w1_a0, w2_a0,w3_a0, m1_a0, m3_a0, m6_a0, m12_a0
# pre_a1, con_a1, d0_a1, w1_a1, w2_a1,w3_a1, m1_a1, m3_a1, m6_a1, m12_a1





# overall changes to relative abundance
# compositional transformation
phy21DS <- transform_sample_counts(phy21DS, function(x) x * 100 / sum(x))

phy21DS@sam_data$timepoint <- factor(phy21DS@sam_data$timepoint, levels = c("pre", "conditioning", "d0", "w1", "w2", "w3", "w4", "w5", "m1", "m3", "m6", "m12"))
REL_MERGEDotu <- as.data.frame(phy21DS@otu_table)
REL_MERGEDotuNA <- na.omit(REL_MERGEDotu)
top25_overall <- names(sort(colSums(REL_MERGEDotuNA), decreasing = TRUE))[1:50]
ps.top25_overall <- prune_taxa(top25, phy21DS)

plot_vio(ps.top25_overall, x = "timepoint", fill = "Genus") + labs(title = "m1", x = "Genus", y = "Abundance (compositional)") +
  scale_fill_genus() +
  gg_theme

# affixing genus name to colours
phy21DS_df_top$Genus <- as.factor(phy21DS_df_top$Genus)
phy21DS_df_genus <- unique(phy21DS_df_genus)
## kind of a horrible function but it works ----
scale_fill_genus <- function(phylo,...) {
  phylo_df<- psmelt(phylo)
  phylo_df$Genus<- as.factor(phylo_df$Genus)
  phyloGenus<- unique(phylo_df$Genus)
  
  ggplot2:::manual_scale(
    "fill",
    values = setNames(c("blue",	"pink",	"yellow",	"green",	"red",	"cyan",
                        "purple",	"#432c8f",	"#56cd92",	"#fca420",	"#291d4f",	"#bddf66",
                        "#6844d3",	"#d8c237",	"#cf3b43",	"#613876",	"#dcd489",	"#64a79f",
                        "#c7e535",	"#8d3193",	"#67e0d9",	"#601f1c",	   "#99e48e",    	"#d33d78",
                        "#832d50",	"#72e750",	"#e35324",	"#998a99",	"#3c6c29",	"#bcdbbc",
                        "#47b645",	"#c378de",	"#81bce0",	"#6a7dda",   	"#bb6aa3",	"#411c28",
                        "#dc42ab",	"#77a036",	"#d39435",	"#546b9d",	"#948f47",	"#e0885e",
                        "#273141",	"#d0ac9c",	"#383d1f",	"#cea9da",	"#546b9d",	"#578361",
                        "#d57884",	"#3b646b",	"#856338",	"#8c92a8",	"#806166", "#dacedd"
                        
                        
    ), phyloGenus[1:46]),
    ...
  )
}

# Agglomerate by genus and subset by prevalence
tse <- subsetByPrevalent(ALL_tse, rank = "Genus", prevalence = 5/100)

# Transform count assay to relative abundances
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")

tse@colData$aGVHD_status<- factor(tse@colData$aGVHD_status, levels= c("Control", "aGVHD"))

tsep<- convertToPhyloseq(tse)
tsep@sam_data$aGVHD_status<- ifelse(tsep@sam_data$aGVHD == "1", "aGVHD", "no aGVHD")
tsep <- transform_sample_counts(tsep, function(x) x * 100 / sum(x))
tsep_filt <- transform_sample_counts(tsep_filt,function(x) (x - mean(x)) / sd(x)) # z-scores?

#taxa<- data.frame(tsep@tax_table)
# making nicer legend labels
tsep@sam_data$timepoint_long<- factor(tsep@sam_data$timepoint)
tsep@sam_data$timepoint_long <- recode_factor(tsep@sam_data$timepoint_long,
                                              pre= "pre-exam", conditioning= "conditioning", d0= "day_of_HSCT", w1= "week_1", w2= "week_2",
                                              w3= "week_3", m1= "month_1", m3= "month_3", m6= "month_6", m12= "month_12")
tsep_filt<- subset_samples(tsep, tsep@sam_data$PatientID == "P04")

plot_bar(tsep, x="timepoint_long", fill = "Genus") +
  geom_bar(aes(fill = Genus), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n", title = "All patients~aGVHD") +
  facet_wrap(~aGVHD_status, scales = "free") +
  scale_fill_genus(tsep) +
  gg_theme +
  # ylim(0,800) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x=element_blank())

# , subtitle= "TOO: +10d/w1 | CS: +9d/w1"


# mgPipe prep ----

# Preparing data for mgPipe

# making the merged OTU table into a dataframe
species_counts <- as.data.frame(phy21DS@otu_table)
taxonomies <- as.data.frame(phy21DS@tax_table)

# transposing it so that sample names are columns and ASVs are rows
species_counts_t <- as.data.frame(t(species_counts))

# using the packages microbiome and microbiomeutilities to sort data for mgPipe

mergedotu_tab <- microbiome::abundances(phy21DS)
tax_tab <- phyloseq::tax_table(phy21DS)
tax_tab[1:5, 1:5]

merged_family <- tax_glom(phy21DS, "Family", NArm = TRUE)
taxa_names(merged_family) <- tax_table(merged_family)[, "Family"]
merged_genus <- tax_glom(phy21DS, "Genus", NArm = TRUE)
taxa_names(merged_genus) <- tax_table(merged_genus)[, "Genus"]

genus_taxa <- as.data.frame(merged_genus@tax_table)
family_taxa <- as.data.frame(merged_family@tax_table)
merged_family@otu_table[, colSums(merged_family@otu_table != 0) > 0] # removing columns with no values
# normalising abundances:
norm_merged_genus <- transform(merged_genus@otu_table, "compositional")
norm_merged_family <- transform(merged_family, "compositional")

# taxonmic abundance table:
genus_counts_t <- as.data.frame(t(norm_merged_genus))
family_counts_t <- as.data.frame(t(norm_merged_family@sotu_table))
family_counts_t <- family_counts_t[, colSums(family_counts_t != 0) > 0] # removing columns with no values
genus_counts_t <- genus_counts_t[, colSums(genus_counts_t != 0) > 0] # removing columns with no values


# export as .csv
write.csv(genus_counts_t, file = "genus_Abs12_6.csv", row.names = TRUE)

write.csv(family_counts_t, file = "fam_abundances.csv", row.names = TRUE)
#############################
# AFTER MGPIPE #
#############################

# Random Forest models of diet reactions, subsystems and ASVs ----

# subsystems
# setting up the data
subAb_t <- data.frame(t(subAb)) # make reaction abundance into a dataframe
subAb_t <- na.omit(subAb_t) # remove NAs
subAb_t$aGVHD <- samples[rownames(subAb_t), "aGVHD"]

# set seed
set.seed(612)

# run model
RF_subsys_aGVHD <- randomForest(x = subAb_t[, 1:(ncol(subAb_t) - 1)], y = subAb_t[, ncol(subAb_t)], ntree = 10000, importance = TRUE, proximities = TRUE)
RF_subsys_aGVHD 


# plotting error/ntrees
plot(RF_subsys_aGVHD)

# plotting top 10 important subsystems
par(mfrow = c(1, 2))
RF_state_classify_imp <- as.data.frame(RF_subsys_aGVHD$importance)
RF_state_classify_imp$features <- rownames(RF_state_classify_imp)
RF_state_classify_imp$`%IncMSE` <- sort(-RF_state_classify_imp$`%IncMSE`)
MSE_top10 <- names(RF_state_classify_imp$features)[1:10]
RF_state_classify_imp_top10 <- RF_state_classify_imp[1:10, ]

ggplot(data = RF_state_classify_imp_top10, aes(x = features, y = `%IncMSE`, fill = features)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Spectral") +
  labs(y = "% Increase in Mean Squared Error", x = "Subsystems", title = "RF Important Subsystems in aGVHD Incidence") +
  gg_theme +
  theme(axis.text.x = element_blank())
barplot(RF_state_classify_imp[1:10, "%IncMSE"], names.arg = RF_state_classify_imp[1:10, "features"], ylab = "% Increase in Mean Squared Error (Variable Importance)", main = "Regression RF")


# RF for ASVS

tse_Phylo <- convertToPhyloseq(ALL_tse)

# need to get data into the right format
otu_table_scaled_state <- data.frame(t(counts))
otu_table_scaled_maxaGVHD <- data.frame(t(counts))
samples$aGVHD <- with(samples, ifelse(max_aGVHD_grade == 0, 0, 1))
otu_table_scaled_state$aGVHD <- samples[rownames(otu_table_scaled_state), "aGVHD"]
otu_table_scaled_maxaGVHD$max_aGVHD <- samples[rownames(otu_table_scaled_maxaGVHD), "max_aGVHD_grade"]

# set seed to make  RF reproducible
set.seed(414)

RF_aGVHD_classify <- randomForest(x = otu_table_scaled_state[, 1:(ncol(otu_table_scaled_state) - 1)], y = otu_table_scaled_state[, ncol(otu_table_scaled_state)], ntree = 501, importance = TRUE, proximities = TRUE)

# looking at the models

RF_aGVHD_classify
# Type of random forest: regression
# Number of trees: 501
# No. of variables tried at each split: 485
#
# Mean of squared residuals: 1.451227
# % Var explained: 40.14



# plotting important features
par(mfrow = c(1, 2))
RF_state_classify_imp <- as.data.frame(RF_aGVHD_classify$importance)
RF_state_classify_imp$features <- rownames(RF_state_classify_imp)
RF_state_classify_imp$`%IncMSE` <- sort(-RF_state_classify_imp$`%IncMSE`)
MSE_top10 <- names(RF_state_classify_imp$features)[1:10]
RF_state_classify_imp_top10 <- RF_state_classify_imp[1:10, ]
rownames(RF_state_classify_imp_top10)

ggplot(data = RF_state_classify_imp_top10, aes(x = features, y = `%IncMSE`, fill = features)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Spectral") +
  labs(y = "% Increase in Mean Squared Error", x = "ASVs", title = "RF Important ASVs in aGVHD Incidence") +
  gg_theme +
  theme(axis.text.x = element_blank())

# RF models of intake/secreted diet flux exchange reactions
# need to get data into the right format
otu_table_scaled_state <- data.frame(t(dietUptake))
otu_table_scaled_state <- na.omit(otu_table_scaled_state) # remove NAs

otu_table_scaled_maxaGVHD <- data.frame(t(dietUptake))
otu_table_scaled_maxaGVHD <- na.omit(otu_table_scaled_maxaGVHD) # remove NAs

samples$aGVHD <- with(samples, ifelse(max_aGVHD_grade == 0, 0, 1))
otu_table_scaled_state$aGVHD <- samples[rownames(otu_table_scaled_state), "aGVHD"]
otu_table_scaled_maxaGVHD$max_aGVHD <- samples[rownames(otu_table_scaled_maxaGVHD), "max_aGVHD_grade"]
RF_aGVHD_dUp
set.seed(212)
RF_aGVHD_dUp <- randomForest(x = otu_table_scaled_state[, 1:(ncol(otu_table_scaled_state) - 1)], y = otu_table_scaled_state[, ncol(otu_table_scaled_state)], ntree = 1000, importance = TRUE, proximities = TRUE)
RF_max_aGVHD_dUp <- randomForest(x = otu_table_scaled_maxaGVHD[, 1:(ncol(otu_table_scaled_maxaGVHD) - 1)], y = otu_table_scaled_maxaGVHD[, ncol(otu_table_scaled_maxaGVHD)], ntree = 1000, importance = TRUE, proximities = TRUE)

RF_aGVHD_dUp
# Type of random forest: regression
# Number of trees: 28001
# No. of variables tried at each split: 1120
# 
# Mean of squared residuals: 1.682143
# % Var explained: 31.31

RF_max_aGVHD_dUp

par(mfrow = c(1, 2))
RF_state_classify_imp <- as.data.frame(RF_max_aGVHD_dUp$importance)
RF_state_classify_imp$features <- rownames(RF_state_classify_imp)
RF_state_classify_imp$`%IncMSE` <- sort(RF_state_classify_imp$`%IncMSE`)
MSE_top10 <- names(RF_state_classify_imp$features)[1:10]
RF_state_classify_imp_top10 <- RF_state_classify_imp[1:10, ]
rownames(RF_state_classify_imp_top10)

ggplot(data = RF_state_classify_imp_top10, aes(x = features, y = `%IncMSE`, fill = features)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Spectral") +
  labs(y = "% Increase in Mean Squared Error", x = "EX_reactions", title = "RF Important dietUptake reactions in aGVHD grade") +
  gg_theme +
  theme(axis.text.x = element_blank())

otu_table_scaled_state <- data.frame(t(dietSecretion))
otu_table_scaled_state <- na.omit(otu_table_scaled_state) # remove NAs

otu_table_scaled_maxaGVHD <- data.frame(t(dietSecretion))
otu_table_scaled_maxaGVHD <- na.omit(otu_table_scaled_maxaGVHD) # remove NAs

samples$aGVHD <- with(samples, ifelse(max_aGVHD_grade == 0, 0, 1))
otu_table_scaled_state$aGVHD <- samples[rownames(otu_table_scaled_state), "aGVHD"]
otu_table_scaled_maxaGVHD$max_aGVHD <- samples[rownames(otu_table_scaled_maxaGVHD), "max_aGVHD_grade"]

RF_aGVHD_dSec <- randomForest(x = otu_table_scaled_state[, 1:(ncol(otu_table_scaled_state) - 1)], y = otu_table_scaled_state[, ncol(otu_table_scaled_state)], ntree = 1000, importance = TRUE, proximities = TRUE)
RF_max_aGVHD_dSec <- randomForest(x = otu_table_scaled_maxaGVHD[, 1:(ncol(otu_table_scaled_maxaGVHD) - 1)], y = otu_table_scaled_maxaGVHD[, ncol(otu_table_scaled_maxaGVHD)], ntree = 1000, importance = TRUE, proximities = TRUE)

# pretty much the same...

par(mfrow = c(1, 2))
RF_state_classify_imp <- as.data.frame(RF_aGVHD_dSec$importance)
RF_state_classify_imp$features <- rownames(RF_state_classify_imp)
RF_state_classify_imp$`%IncMSE` <- sort(-RF_state_classify_imp$`%IncMSE`)
MSE_top10 <- names(RF_state_classify_imp$features)[1:10]
RF_state_classify_imp_top10 <- RF_state_classify_imp[1:10, ]
rownames(RF_state_classify_imp_top10)

ggplot(data = RF_state_classify_imp_top10, aes(x = features, y = `%IncMSE`, fill = features)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Spectral") +
  labs(y = "% Increase in Mean Squared Error", x = "EX_reactions", title = "RF Important dietSecretion reactions in aGVHD grade") +
  gg_theme +
  theme(axis.text.x = element_blank())
# uptake and secretion are reactions are the same



