### MSstat analysis for whole proteome and phosphoproteome of mistranslating yeast strains
### This script analyze the DIANN ouputs using MSstats to determine protein and phosphosite abundance as well as test differential abundance in each strain relative to the control strain
### Input is DIANN Whole Proteome and Phosphoproteome output (the report.parquet files)
### Outputs are CSV files containing protein/phosphosite level abundances in each yeast strain and protein/phosphosite differential abundance relative to the control strain

# Load packages needed for analysis
library(MSstats)
library(MSstatsPTM)
library(tidyverse)
library(arrow)

# Set working directory to where DIANN files are located
setwd("")

## Read in DIANN v2.1 proteome report, filter and prepare whole proteome data for MSstats
WholeProteome.DIANNInput <- read_parquet("WholeProteome_DIANNreport.parquet")

# Filter out decoys, filter Q-values and filter quantification quality
WholeProteome.DIANNInput.Filtered <- WholeProteome.DIANNInput %>%
  filter(Decoy == 0) %>%
  filter(Q.Value <= 0.01,
         Global.Q.Value <= 0.01) %>%
  filter(Quantity.Quality > 0.5)

rm(WholeProteome.DIANNInput)

# Set up data in format compatible with MSstats
MSstats.Columns <- WholeProteome.DIANNInput.Filtered %>%
  separate(Run, into = c("RawFile", "Experiment", "Condition", "BioReplicate", "Measure", "Measure2"), remove = F) %>%
  select(ProteinName = Genes,
         PeptideSequence = Modified.Sequence,
         PrecursorCharge = Precursor.Charge,
         Intensity = Precursor.Quantity,
         Run,
         BioReplicate,
         Condition) %>%
  filter(Intensity != 0) %>%
  mutate(FragmentIon = NA,
         ProductCharge = 0,
         IsotopeLabelType = "L")

# Read in DIANN v2.1 report phosphoproteome data, filter and prepare phosphoproteome data for MSstatsPTM
PhosphoProteome.DIANNInput <- read_parquet("Phosphoproteome_DIANNreport.parquet")

# Filter out decoys, filter on Q.value, global Q.value, PTM localization and PTM Q.value, filter on quantity quality and only keep peptides that are phosphorylation
PhosphoProteome.DIANNInput.Filtered <- PhosphoProteome.DIANNInput %>%
  filter(Decoy == 0) %>%
  filter(Q.Value <= 0.01,
         Global.Q.Value <= 0.01,
         PTM.Site.Confidence >= 0.75,
         Global.Peptidoform.Q.Value <= 0.01) %>%
  filter(Quantity.Quality > 0.5) %>%
  filter(str_detect(Modified.Sequence, "(UniMod:21)"))

rm(PhosphoProteome.DIANNInput)

# Set up data in format compatible with MSstats
MSstats.Columns.Phospho <- PhosphoProteome.DIANNInput.Filtered %>%
  separate(Run, into = c("RawFile", "Experiment", "Condition", "BioReplicate", "Measure", "Measure2"), remove = F) %>%
  select(ProteinName = Genes,
         PeptideSequence = Modified.Sequence,
         PrecursorCharge = Precursor.Charge,
         Intensity = Precursor.Quantity,
         Run,
         BioReplicate,
         Condition,
         Protein.Sites) %>%
  filter(Intensity != 0) %>%
  mutate(FragmentIon = NA,
         ProductCharge = 0,
         IsotopeLabelType = "L")

# Extract the exact residue number where the phosphorylation occurs (to roll up to phosphosite level)
MSstats.Phospho <- MSstats.Columns.Phospho %>%
  mutate(Sites = str_extract(Protein.Sites, "(?<=:)[^\\]]+"),
         Sites.List = str_split(Sites, ","),
         STY.Sites = lapply(Sites.List, function(x) x[str_starts(x, "S|T|Y")]),
         STY.Sites = sapply(STY.Sites, function(x) paste(x, collapse = ","))) %>%
  mutate(ProteinName = paste0(ProteinName, "_", STY.Sites)) %>%
  select(-Sites,
         -Sites.List,
         -STY.Sites)

# Make individual entries for each phosphosite (if one peptide has two phosphorylation events, that will be represented on two lines in the data frame)
MSstats.Phospho.Split <- MSstats.Phospho %>%
  separate(ProteinName, into = c("Protein", "Sites"), sep = "_", remove = FALSE) %>%
  separate_rows(Sites, sep = ",") %>%
  mutate(ProteinName = paste0(Protein, "_", Sites)) %>%
  select(-Protein, -Sites)

# Merge proteome and phosphoproteome into a single list for MSstats
MSstatsPTM.Input <- list(PTM = MSstats.Phospho.Split, PROTEIN = MSstats.Columns)

# Summarize the data using MSstats
MSstatsPTM.Summary <- dataSummarizationPTM(MSstatsPTM.Input,
                                           logTrans = 2,
                                           normalization = "equalizeMedians",
                                           normalization.PTM = "equalizeMedians",
                                           featureSubset = "all",
                                           featureSubset.PTM = "all",
                                           min_feature_count = 2,
                                           min_feature_count.PTM = 1,
                                           summaryMethod = "TMP",
                                           censoredInt = "NA",
                                           MBimput = T,
                                           MBimpute.PTM = T,
                                           verbose = T,
                                           use_log_file = F,
                                           append = F)

# Write protein level and phosphoproteome level output to csv files
ProteinLevel <- MSstatsPTM.Summary$PROTEIN$ProteinLevelData
PhosphoLevel <- MSstatsPTM.Summary$PTM$ProteinLevelData

write.csv(ProteinLevel, "../FINAL_Data/ProteinLevel_MSstatsSummary.csv")
write.csv(PhosphoLevel, "../FINAL_Data/PhosphoLevel_MSstatsSummary.csv")

# Differential protein and phosphosite abundance analysis with MSstats
comparison <- matrix(c(
  -1,  1,  0,  0,  # ProAla vs Control
  -1,  0,  1,  0,  # ProSer vs Control
  -1,  0,  0,  1   # ArgSer vs Control
), nrow = 3, byrow = TRUE)

colnames(comparison) <- c("Ctrl", "ProAla", "ProSer", "ArgSer")
rownames(comparison) <- c("ProAla_vs_Ctrl", "ProSer_vs_Ctrl", "ArgSer_vs_Ctrl")

Proteome.Model <- groupComparisonPTM(MSstatsPTM.Summary, data.type = "LabelFree", contrast.matrix = comparison)

# Write differential abundance protein and phosphosite output to csv file
Proteme.DA <- Proteome.Model$PROTEIN.Model
Phospho.DA <- Proteome.Model$ADJUSTED.Model

write.csv(Proteme.DA, "../FINAL_Data/ProteinDifferentialAbundance_MSstats.csv")
write.csv(Phospho.DA, "../FINAL_Data/ProteinAdjustedPhosphositeDifferentialAbundance_MSstats.csv")