### DDA analysis of mistranslating tRNA Strains digested with LysC
### This code analyzes how much mistranslation is occuring in each strain and creates plots that make up Figure 1 and Figure S1

## Packages needed to complete analysis and create plots
library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(cowplot)
library(ggsci)
library(GGally)
library(ggpubr)
library(scales)
library(rmotifx)
library(ggseqlogo)
library(dagLogo)
library(data.table)
library(ggridges)
library(Biostrings)

## My R Theme - set theme for all plots
theme_matt <- function(){
  theme_classic() %+replace%
    theme(panel.background = element_rect(fill = "transparent", color = NA),
          axis.text = element_text(color = "black", size = rel(1.2)),
          axis.title.x = element_text(margin = margin(t = 5), size = rel(1.5)),
          axis.title.y = element_text(angle = 90, vjust = 1, margin = margin(r = 8), size = rel(1.5)),
          axis.line = element_line(linewidth = rel(0)),
          axis.ticks = element_line(linewidth = rel(1.3)),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.4))
}

## Set working directory for where the MSFragger searched csv files are located
setwd("")

## Read in MSFragger searched combined modified peptide tsv file
WholeProteome.DDA <- read.csv("WholeProteome_LysC_DDA_MSFragger_CombinedPeptide.tsv", sep = "\t")

# Convert data to tidy format and filter
WholeProteome.DDA.Tidy <- WholeProteome.DDA %>%
  dplyr::select(Peptide.Sequence, Modified.Sequence, Start, End, Gene, Mapped.Genes, ArgSer_1.Intensity:ProSer_6.Intensity, Protein.ID) %>%
  gather(value = "Intensity", key = "Sample", ArgSer_1.Intensity:ProSer_6.Intensity) %>%
  mutate(ProSer = ifelse(str_detect(Modified.Sequence, "P\\[-10.0207\\]"), T, F),
         ProAla = ifelse(str_detect(Modified.Sequence, "P\\[-26.0156\\]"), T, F),
         ArgSer = ifelse(str_detect(Modified.Sequence, "R\\[-69.0691\\]"), T, F)) %>%
  mutate(ProSer.Count = str_count(Modified.Sequence, "P\\[-10.0207\\]"),
         ProAla.Count = str_count(Modified.Sequence, "P\\[-26.0156\\]"),
         ArgSer.Count = str_count(Modified.Sequence, "R\\[-69.0691\\]")) %>%
  separate(Sample, into = c("Sample", "Replicate"))

# Filter out peptides with two or more types of mistranslation and for peptides with no signal intensity
WholeProteome.DDA.Tidy.Filter <- WholeProteome.DDA.Tidy %>%
  filter(((ProSer + ProAla + ArgSer) == 1) | ((ProSer + ProAla + ArgSer) == 0),
         Intensity > 0) %>%
  mutate(NumMuts = ProAla.Count + ProSer.Count + ArgSer.Count)

# Plot number of mutations per peptide for proline to alanine substitution
NumMuts.ProAla <- WholeProteome.DDA.Tidy.Filter %>%
  group_by(Sample, Replicate, ProAla.Count) %>%
  filter(ProAla.Count > 0)

NumMuts.ProAla$Sample <- factor(NumMuts.ProAla$Sample, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))

ggplot(NumMuts.ProAla, aes(x = Sample, fill = Sample, z = Replicate)) + 
  geom_bar(color = "black", position = "dodge") + 
  facet_wrap(~ProAla.Count) +
  ylab("Unique Peptides") + 
  xlab("") +
  scale_fill_brewer(palette = "Set1") +
  theme_matt() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1300))

# Plot number of mutations per peptide for proline to serine substitution
NumMuts.ProSer <- WholeProteome.DDA.Tidy.Filter %>%
  group_by(Sample, Replicate, ProSer.Count) %>%
  filter(ProSer.Count > 0)

NumMuts.ProSer$Sample <- factor(NumMuts.ProSer$Sample, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))

ggplot(NumMuts.ProSer, aes(x = Sample, fill = Sample, z = Replicate)) + 
  geom_bar(color = "black", position = "dodge") + 
  facet_wrap(~ProSer.Count) +
  ylab("Unique Peptides") + 
  xlab("") +
  scale_fill_brewer(palette = "Set1") +
  theme_matt()+
  scale_y_continuous(expand = c(0,0), limits = c(0,1300))

# Plot number of mutations per peptide for proline to serine substitution
NumMuts.ArgSer <- WholeProteome.DDA.Tidy.Filter %>%
  group_by(Sample, Replicate, ArgSer.Count) %>%
  filter(ArgSer.Count > 0)

NumMuts.ArgSer$Sample <- factor(NumMuts.ArgSer$Sample, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))

ggplot(NumMuts.ArgSer, aes(x = Sample, fill = Sample, z = Replicate)) + 
  geom_bar(color = "black", position = "dodge") + 
  facet_wrap(~ArgSer.Count) +
  ylab("Unique Peptides") + 
  xlab("") +
  scale_fill_brewer(palette = "Set1") +
  theme_matt() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1300))

# Filter to only keep peptides with none or one mistranslation event and create a stripped sequence with no mistranslation
Peptides.OneMist <- WholeProteome.DDA.Tidy.Filter %>%
  filter(NumMuts == 0 | NumMuts == 1) %>%
  mutate(Stripped.Sequence.ProAla = str_replace_all(Modified.Sequence, "P\\[-10.0207\\]", "P"),
         Stripped.Sequence.ProSer = str_replace_all(Stripped.Sequence.ProAla, "P\\[-26.0156\\]", "P"),
         Stripped.Sequence = str_replace_all(Stripped.Sequence.ProSer, "R\\[-69.0691\\]", "R")) %>%
  mutate(PepType = ifelse(NumMuts == 0, "WT", "Mutant")) %>%
  dplyr::select(Sample, Replicate, Modified.Sequence, Stripped.Sequence, NumMuts, PepType, Intensity, Gene)

#Separate mistranslated peptides from WT peptides
WT.Peptides <- Peptides.OneMist %>%
  ungroup() %>%
  filter(PepType == "WT")

Mistranslated.Peptides <- Peptides.OneMist %>%
  ungroup() %>%
  filter(PepType == "Mutant")

# Join to filter out any mutant peptide with no WT counterpart
Filtered.MistranslatedPeptides <- left_join(WT.Peptides, Mistranslated.Peptides, by = c("Sample", "Replicate", "Stripped.Sequence"), suffix = c(".wt", ".mut")) %>%
  mutate(MistranslationType = ifelse(str_detect(Modified.Sequence.mut, "P\\[-10.0207\\]"), "ProSer", 
                                     ifelse(str_detect(Modified.Sequence.mut, "P\\[-26.0156\\]"), "ProAla",
                                            ifelse(str_detect(Modified.Sequence.mut, "R\\[-69.0691\\]"), "ArgSer", NA))))

# Calculate the percent mistranslation based on peptide counts
PtoA.PercentMistranslation <- Filtered.MistranslatedPeptides  %>%
  filter(str_detect(Stripped.Sequence, "P"),
         (MistranslationType == "ProAla" | is.na(MistranslationType))) %>%
  dplyr::select(Replicate, Sample, PepType.wt, PepType.mut) %>%
  group_by(Replicate, Sample) %>%
  summarise(countWT = n(), countMutant = sum(!is.na(PepType.mut))) %>%
  mutate(PercentMistranslation = (countMutant/countWT)*100,
         MistranslationEvent = "P to A")

PtoA.PercentMistranslation$Sample <- factor(PtoA.PercentMistranslation$Sample, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))

PtoA <- ggplot(PtoA.PercentMistranslation, aes(x = Sample, y = PercentMistranslation, fill = Sample, color = Sample)) +
  geom_bar(stat = "summary", alpha = 0.35) +
  geom_jitter(width = 0.25, size = 3) +
  ylab("% Proline to Alanine") +
  xlab("") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_matt() +
  scale_y_continuous(lim = c(0,10), expand = c(0,0)) +
  theme(legend.position = "none")

PtoS.PercentMistranslation <- Filtered.MistranslatedPeptides  %>%
  filter(str_detect(Stripped.Sequence, "P"),
         (MistranslationType == "ProSer" | is.na(MistranslationType))) %>%
  dplyr::select(Replicate, Sample, PepType.wt, PepType.mut) %>%
  group_by(Replicate, Sample) %>%
  summarise(countWT = n(), countMutant = sum(!is.na(PepType.mut))) %>%
  mutate(PercentMistranslation = (countMutant/countWT)*100,
         MistranslationEvent = "P to S")

PtoS.PercentMistranslation$Sample <- factor(PtoS.PercentMistranslation$Sample, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))

PtoS <- ggplot(PtoS.PercentMistranslation, aes(x = Sample, y = PercentMistranslation, fill = Sample, color = Sample)) +
  geom_bar(stat = "summary", alpha = 0.35) +
  geom_jitter(width = 0.25, size = 3) +
  ylab("% Proline to Serine") +
  xlab("") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_matt() +
  scale_y_continuous(lim = c(0,10), expand = c(0,0)) +
  theme(legend.position = "none")

RtoS.PercentMistranslation <- Filtered.MistranslatedPeptides  %>%
  filter(str_detect(Stripped.Sequence, "R"),
         (MistranslationType == "ArgSer" | is.na(MistranslationType))) %>%
  dplyr::select(Replicate, Sample, PepType.wt, PepType.mut) %>%
  group_by(Replicate, Sample) %>%
  summarise(countWT = n(), countMutant = sum(!is.na(PepType.mut))) %>%
  mutate(PercentMistranslation = (countMutant/countWT)*100,
         MistranslationEvent = "R to S")

RtoS.PercentMistranslation$Sample <- factor(RtoS.PercentMistranslation$Sample, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))

RtoS <- ggplot(RtoS.PercentMistranslation, aes(x = Sample, y = PercentMistranslation, fill = Sample, color = Sample)) +
  geom_bar(stat = "summary", alpha = 0.35) +
  geom_jitter(width = 0.25, size = 3) +
  ylab("% Arginine to Serine") +
  xlab("") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_matt() +
  scale_y_continuous(lim = c(0,10), expand = c(0,0)) +
  theme(legend.position = "none")

plot_grid(PtoA, PtoS, RtoS, nrow = 1)

PtoA.AverageMistranslation <- PtoA.PercentMistranslation %>%
  group_by(Sample) %>%
  summarise(MeanMistranslation = mean(PercentMistranslation))

PtoS.AverageMistranslation <- PtoS.PercentMistranslation %>%
  group_by(Sample) %>%
  summarise(MeanMistranslation = mean(PercentMistranslation))

RtoS.AverageMistranslation <- RtoS.PercentMistranslation %>%
  group_by(Sample) %>%
  summarise(MeanMistranslation = mean(PercentMistranslation))

# Summarize mistranslation frequency based on counts
PtoA.SummaryCounts <- PtoA.PercentMistranslation %>%
  ungroup() %>%
  group_by(Sample) %>%
  summarise(MeanMis = mean(PercentMistranslation),
            SDMis = sd(PercentMistranslation)) %>%
  mutate(MistranslationEvent = "P to A")

PtoA.SummaryCounts

PtoA.Ttest <- PtoA.PercentMistranslation %>%
  filter(Sample == "Ctrl" |
           Sample == "ProAla")

write.csv(PtoA.Ttest, "PtoA.Mistranslation.csv")

t.test(PercentMistranslation ~ Sample, data = PtoA.Ttest, var.equal = TRUE) ## p-value < 2.2e-16

PtoS.SummaryCounts <- PtoS.PercentMistranslation %>%
  ungroup() %>%
  group_by(Sample) %>%
  summarise(MeanMis = mean(PercentMistranslation),
            SDMis = sd(PercentMistranslation)) %>%
  mutate(MistranslationEvent = "P to S")

PtoS.SummaryCounts

PtoS.Ttest <- PtoS.PercentMistranslation %>%
  filter(Sample == "Ctrl" |
           Sample == "ProSer")

write.csv(PtoS.Ttest, "PtoS.Mistranslation.csv")

t.test(PercentMistranslation ~ Sample, data = PtoS.Ttest, var.equal = TRUE) ## p-value < 2.2e-16

RtoS.SummaryCounts <- RtoS.PercentMistranslation %>%
  ungroup() %>%
  group_by(Sample) %>%
  summarise(MeanMis = mean(PercentMistranslation),
            SDMis = sd(PercentMistranslation)) %>%
  mutate(MistranslationEvent = "R to S")

RtoS.SummaryCounts

RtoS.Ttest <- RtoS.PercentMistranslation %>%
  filter(Sample == "Ctrl" |
           Sample == "ArgSer")

write.csv(RtoS.Ttest, "RtoS.Mistranslation.csv")

t.test(PercentMistranslation ~ Sample, data = RtoS.Ttest, var.equal = TRUE) ## p-value = 4.922e-13

## Plot with all the bars together (Figure 1B)
Together.PercentMistranslation <- PtoA.PercentMistranslation %>%
  rbind(PtoS.PercentMistranslation, RtoS.PercentMistranslation) %>%
  filter((MistranslationEvent == "P to A" & (Sample == "ProAla" | Sample == "Ctrl")) | (MistranslationEvent == "P to S" & (Sample == "ProSer" | Sample == "Ctrl")) | (MistranslationEvent == "R to S" & (Sample == "ArgSer" | Sample == "Ctrl")))

Together.Summary <- Together.PercentMistranslation %>%
  ungroup() %>%
  group_by(Sample, MistranslationEvent) %>%
  summarise(Mean = mean(PercentMistranslation))

ggplot(Together.PercentMistranslation, aes(x = MistranslationEvent, y = PercentMistranslation, fill = Sample)) +
  geom_bar(stat = "summary", position = "dodge", alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 2.5, shape = 21) +
  labs(y = "% Mistranslated",
       x = "Mistranslation Event",
       fill = "tRNA") +
  scale_fill_manual(values = c("#6496bf", "#58c87e", "#f8c754", "#a165a7"), labels = c("EV", expression(tRNA[G3:U70]^{"Pro"}), expression(tRNA[UGG,G26A]^{"Ser"}), expression(tRNA[UCU,G26A]^{"Ser"}))) +
  theme_matt() +
  scale_y_continuous(lim = c(0,9), expand = c(0,0)) +
  scale_x_discrete(labels = c("Pro\u2192Ala", "Pro\u2192Ser", "Arg\u2192Ser"))

ggsave("20240812_MistranslationFrequencyCount_MSFragger.pdf", device = cairo_pdf, width = 5, height = 4, units = "in")

## Sum signal intensity for mistranslated peptides and normalize it to total signal of all identified peptides (for Figure 1E)
# Pro->Ala
Mistranslation.SignalIntensity.ProAla <- WholeProteome.DDA.Tidy.Filter %>%
  mutate(Mistranslated = ifelse(ProAla.Count > 0, "Yes", "No")) %>%
  group_by(Sample, Replicate, Mistranslated) %>%
  summarise(SumInt = sum(Intensity))

NotMistranslatedSignal.ProAla <- Mistranslation.SignalIntensity.ProAla %>%
  filter(Mistranslated == "No")

MistranslatedSignal.ProAla <- Mistranslation.SignalIntensity.ProAla %>%
  filter(Mistranslated == "Yes") %>%
  left_join(NotMistranslatedSignal.ProAla, by = c("Sample", "Replicate")) %>%
  mutate(PercentSignalMistranslated = (SumInt.x/(SumInt.y+SumInt.x))*100)

Summary.ProAla <- MistranslatedSignal.ProAla %>%
  ungroup() %>%
  group_by(Sample) %>%
  summarise(Mean = mean(PercentSignalMistranslated), SD = sd(PercentSignalMistranslated))

# Pro->Ser
Mistranslation.SignalIntensity.ProSer <- WholeProteome.DDA.Tidy.Filter %>%
  mutate(Mistranslated = ifelse(ProSer.Count > 0, "Yes", "No")) %>%
  group_by(Sample, Replicate, Mistranslated) %>%
  summarise(SumInt = sum(Intensity))

NotMistranslatedSignal.ProSer <- Mistranslation.SignalIntensity.ProSer %>%
  filter(Mistranslated == "No")

MistranslatedSignal.ProSer <- Mistranslation.SignalIntensity.ProSer %>%
  filter(Mistranslated == "Yes") %>%
  left_join(NotMistranslatedSignal.ProSer, by = c("Sample", "Replicate")) %>%
  mutate(PercentSignalMistranslated = (SumInt.x/(SumInt.y+SumInt.x))*100)

Summary.ProSer <- MistranslatedSignal.ProSer %>%
  ungroup() %>%
  group_by(Sample) %>%
  summarise(Mean = mean(PercentSignalMistranslated), SD = sd(PercentSignalMistranslated))

# Arg->Ser
Mistranslation.SignalIntensity.ArgSer <- WholeProteome.DDA.Tidy.Filter %>%
  mutate(Mistranslated = ifelse(ArgSer.Count > 0, "Yes", "No")) %>%
  group_by(Sample, Replicate, Mistranslated) %>%
  summarise(SumInt = sum(Intensity))

NotMistranslatedSignal.ArgSer <- Mistranslation.SignalIntensity.ArgSer %>%
  filter(Mistranslated == "No")

MistranslatedSignal.ArgSer <- Mistranslation.SignalIntensity.ArgSer %>%
  filter(Mistranslated == "Yes") %>%
  left_join(NotMistranslatedSignal.ArgSer, by = c("Sample", "Replicate")) %>%
  mutate(PercentSignalMistranslated = (SumInt.x/(SumInt.y+SumInt.x))*100)

Summary.ArgSer <- MistranslatedSignal.ArgSer %>%
  ungroup() %>%
  group_by(Sample) %>%
  summarise(Mean = mean(PercentSignalMistranslated), SD = sd(PercentSignalMistranslated))

# Merge the mistranslation signal calculations for each strain together
MistranslatedSignalTogether.ProAla <- MistranslatedSignal.ProAla %>%
  filter(Sample == "ProAla") %>%
  select(Sample, Replicate, PercentSignalMistranslated)
  
MistranslatedSignalTogether.ProSer <- MistranslatedSignal.ProSer %>%
  filter(Sample == "ProSer") %>%
  select(Sample, Replicate, PercentSignalMistranslated)

MistranslatedSignalTogether.ArgSer <- MistranslatedSignal.ArgSer %>%
  filter(Sample == "ArgSer") %>%
  select(Sample, Replicate, PercentSignalMistranslated)

MistranslatedSignalTogether <- rbind(MistranslatedSignalTogether.ArgSer,
                                     MistranslatedSignalTogether.ProAla,
                                     MistranslatedSignalTogether.ProSer)


# Import relative growth data and plot against mistranslation frequency to create Figure 1E
growth <- read.csv("../GrowthtRNAs/NormalizedGrowth.csv")

growth.mistranslation <- growth %>% 
  filter(tRNA != "Y69") %>%
  dplyr::select(Replicate, tRNA, NormalizedGrowth = Normalized) %>%
  mutate(Replicate = case_when(Replicate == "A" ~ "1",
                               Replicate == "B" ~ "2",
                               Replicate == "C" ~ "3",
                               Replicate == "D" ~ "4",
                               Replicate == "E" ~ "5",
                               Replicate == "F" ~ "6")) %>%
  left_join(MistranslatedSignalTogether %>% 
              dplyr::select(Replicate, tRNA = Sample, PercentSignalMistranslated) %>%
              filter(tRNA != "Ctrl")) %>%
  mutate(tRNA = fct_relevel(tRNA, c("ProAla", "ProSer", "ArgSer")))

ggplot(growth.mistranslation, aes(x = PercentSignalMistranslated, y = NormalizedGrowth, color = tRNA)) +
  geom_point(size = 3.5) +
  theme_matt() +
  scale_y_continuous(limits = c(0.8, 1)) +
  scale_x_continuous(limits = c(0,1))  +
  scale_color_manual(values = c("#58c87e", "#f8c754", "#a165a7"), labels = c("Pro\u2192Ala", "Pro\u2192Ser", "Arg\u2192Ser")) +
  ylab("Relative Growth") +
  xlab("% Proteome Mistranslated")

ggsave("PercentMistranslationVSRelativeGrowth.pdf", device = cairo_pdf, width = 4, height = 3)

## Determine the number of mistranslated peptides observed for each codon (For Figure S1)
# Filter for peptides containing only one mistranslation event (because of localization issues)

# Read in yeast FASTA file
YeastCDSFASTA <- readDNAStringSet("202411_SGD_orf_coding.fasta")
headers <- names(YeastCDSFASTA)
gene_names <- str_extract(headers, "(?<= )[A-Z0-9]+")
names(YeastCDSFASTA) <- gene_names

# Function to extract codon identities for proline residues
extract_codon_for_proline <- function(row, YeastCDSFASTA) {
  
  # Extract relevant data from the row
  gene_id <- row$Gene
  proline_cds_start <- row$CDSStart
  
  # Check if gene exists in the FASTA sequences
  if (!(gene_id %in% names(YeastCDSFASTA))) {
    return(NA)  # Return NA if gene not found
  }
  
  # Extract the gene sequence
  gene_sequence <- as.character(YeastCDSFASTA[gene_id])
  
  # Calculate the CDS position of the proline
  proline_cds_end <- proline_cds_start + 2
  
  # Ensure positions are within bounds
  if (proline_cds_end > nchar(gene_sequence)) {
    return(NA)
  }
  
  # Extract the codon
  codon <- substr(gene_sequence, proline_cds_start, proline_cds_end)
  return(codon)
}

ProAla.UniqueCodon <- WholeProteome.DDA.Tidy.Filter %>%
  mutate(NumPro = str_count(Peptide.Sequence, "P")) %>%
  filter(NumMuts == 0 | NumMuts == 1) %>%
  filter(Sample == "Ctrl" | Sample == "ProAla") %>%
  filter(NumPro == 1) %>%
  filter(Mapped.Genes == "") %>%
  filter(ProSer.Count == 0) %>%
  filter(ArgSer.Count == 0) %>%
  select(Peptide.Sequence, Start, End, Gene) %>%
  unique() %>%
  mutate(ProlinePosition = str_locate(Peptide.Sequence, "P")[,"start"],
         ProteinPosition = (Start + ProlinePosition - 1),
         CDSStart = (ProteinPosition*3)-2) %>%
  rowwise() %>%
  mutate(codon = extract_codon_for_proline(cur_data(), YeastCDSFASTA))


ProAla.Codon <- WholeProteome.DDA.Tidy.Filter %>%
  mutate(NumPro = str_count(Peptide.Sequence, "P")) %>%
  filter(NumMuts == 0 | NumMuts == 1) %>%
  filter(Sample == "Ctrl" | Sample == "ProAla") %>%
  filter(NumPro == 1) %>%
  filter(Mapped.Genes == "") %>%
  filter(ProSer.Count == 0) %>%
  filter(ArgSer.Count == 0) %>%
  left_join(ProAla.UniqueCodon %>% select(-ProlinePosition, -ProteinPosition, -CDSStart), by = c("Peptide.Sequence", "Start", "End", "Gene")) %>% 
  na.omit() %>%
  mutate(PepType = ifelse(NumMuts == 0, "WT", "Mutant"),
         Stripped.Sequence.ProAla = str_replace_all(Modified.Sequence, "P\\[-26.0156\\]", "P"))

#Separate mistranslated peptides from WT peptides
WT.Peptides.ProAla <- ProAla.Codon %>%
  ungroup() %>%
  filter(PepType == "WT")

Mistranslated.Peptides.ProAla <- ProAla.Codon %>%
  ungroup() %>%
  filter(PepType == "Mutant")

# Join to filter out any mutant peptide with no WT counterpart
Filtered.MistranslatedPeptides <- left_join(WT.Peptides.ProAla, Mistranslated.Peptides.ProAla, by = c("Sample", "Replicate", "Stripped.Sequence.ProAla", "codon"), suffix = c(".wt", ".mut"))

# Calculate the percent mistranslation based on peptide counts
PtoA.PercentMistranslation <- Filtered.MistranslatedPeptides  %>%
  filter(codon == "CCA" |
           codon == "CCT" |
           codon == "CCC" |
           codon == "CCG") %>%
  filter(str_detect(Stripped.Sequence.ProAla, "P")) %>%
  dplyr::select(Replicate, Sample, PepType.wt, PepType.mut, codon) %>%
  group_by(Replicate, Sample, codon) %>%
  summarise(countWT = n(), countMutant = sum(!is.na(PepType.mut))) %>%
  mutate(PercentMistranslation = (countMutant/countWT)*100,
         MistranslationEvent = "P to A")

PtoA.PercentMistranslation$Sample <- factor(PtoA.PercentMistranslation$Sample, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))

ggplot(PtoA.PercentMistranslation, aes(x = codon, y = PercentMistranslation, fill = Sample)) +
  geom_bar(stat = "summary", alpha = 0.35, position = position_dodge(width = 0.9)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 2) +
  ylab("% Proline to Alanine") +
  xlab("") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_matt() +
  scale_y_continuous(lim = c(0,8), expand = c(0,0)) +
  xlab("Codon") +
  scale_fill_manual(values = c("#6496bf", "#58c87e"), labels = c("EV", expression(tRNA[G3:U70]^{"Pro"})))

ggsave("CodonUsage_ProAla.pdf", device = cairo_pdf, width = 5, height = 4)

# Pro->Ser
ProSer.UniqueCodon <- WholeProteome.DDA.Tidy.Filter %>%
  mutate(NumPro = str_count(Peptide.Sequence, "P")) %>%
  filter(NumMuts == 0 | NumMuts == 1) %>%
  filter(Sample == "Ctrl" | Sample == "ProSer") %>%
  filter(NumPro == 1) %>%
  filter(Mapped.Genes == "") %>%
  filter(ProAla.Count == 0) %>%
  filter(ArgSer.Count == 0) %>%
  select(Peptide.Sequence, Start, End, Gene) %>%
  unique() %>%
  mutate(ProlinePosition = str_locate(Peptide.Sequence, "P")[,"start"],
         ProteinPosition = (Start + ProlinePosition - 1),
         CDSStart = (ProteinPosition*3)-2) %>%
  rowwise() %>%
  mutate(codon = extract_codon_for_proline(cur_data(), YeastCDSFASTA))

ProSer.Codon <- WholeProteome.DDA.Tidy.Filter %>%
  mutate(NumPro = str_count(Peptide.Sequence, "P")) %>%
  filter(NumMuts == 0 | NumMuts == 1) %>%
  filter(Sample == "Ctrl" | Sample == "ProSer") %>%
  filter(NumPro == 1) %>%
  filter(Mapped.Genes == "") %>%
  filter(ProAla.Count == 0) %>%
  filter(ArgSer.Count == 0) %>%
  left_join(ProSer.UniqueCodon %>% select(-ProlinePosition, -ProteinPosition, -CDSStart), by = c("Peptide.Sequence", "Start", "End", "Gene")) %>% 
  na.omit() %>%
  mutate(PepType = ifelse(NumMuts == 0, "WT", "Mutant"),
         Stripped.Sequence.ProSer = str_replace_all(Modified.Sequence, "P\\[-10.0207\\]", "P"))

#Separate mistranslated peptides from WT peptides
WT.Peptides.ProSer <- ProSer.Codon %>%
  ungroup() %>%
  filter(PepType == "WT")

Mistranslated.Peptides.ProSer <- ProSer.Codon %>%
  ungroup() %>%
  filter(PepType == "Mutant")

# Join to filter out any mutant peptide with no WT counterpart
Filtered.MistranslatedPeptides.ProSer <- left_join(WT.Peptides.ProSer, Mistranslated.Peptides.ProSer, by = c("Sample", "Replicate", "Stripped.Sequence.ProSer", "codon"), suffix = c(".wt", ".mut"))

# Calculate the percent mistranslation based on peptide counts
PtoS.PercentMistranslation <- Filtered.MistranslatedPeptides.ProSer  %>%
  filter(codon == "CCA" |
           codon == "CCT" |
           codon == "CCC" |
           codon == "CCG") %>%
  filter(str_detect(Stripped.Sequence.ProSer, "P")) %>%
  dplyr::select(Replicate, Sample, PepType.wt, PepType.mut, codon) %>%
  group_by(Replicate, Sample, codon) %>%
  summarise(countWT = n(), countMutant = sum(!is.na(PepType.mut))) %>%
  mutate(PercentMistranslation = (countMutant/countWT)*100,
         MistranslationEvent = "P to S")

PtoS.PercentMistranslation$Sample <- factor(PtoS.PercentMistranslation$Sample, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))

ggplot(PtoS.PercentMistranslation, aes(x = codon, y = PercentMistranslation, fill = Sample)) +
  geom_bar(stat = "summary", alpha = 0.35, position = position_dodge(width = 0.9)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 2) +
  ylab("% Proline to Serine") +
  xlab("") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_matt() +
  scale_y_continuous(lim = c(0,8), expand = c(0,0)) +
  xlab("Codon") +
  scale_fill_manual(values = c("#6496bf", "#f8c754"), labels = c("EV", expression(tRNA[UGG,G26A]^{"Ser"})))

ggsave("CodonUsage_ProSer.pdf", device = cairo_pdf, width = 5, height = 4)

t.test(PercentMistranslation ~ Sample, data = PtoS.PercentMistranslation %>% filter(codon == "CCA"), var.equal = TRUE) # < 2.2e-16
t.test(PercentMistranslation ~ Sample, data = PtoS.PercentMistranslation %>% filter(codon == "CCC"), var.equal = TRUE) # < 7e-6
t.test(PercentMistranslation ~ Sample, data = PtoS.PercentMistranslation %>% filter(codon == "CCG"), var.equal = TRUE) # < 3.691e-10
t.test(PercentMistranslation ~ Sample, data = PtoS.PercentMistranslation %>% filter(codon == "CCT"), var.equal = TRUE) # < 4.539e-14

# Arg->Ser

# Function to extract codon identities for arginine residues
extract_codon_for_arginine <- function(row, YeastCDSFASTA) {
  
  # Extract relevant data from the row
  gene_id <- row$Gene
  arginine_cds_start <- row$CDSStart
  
  # Check if gene exists in the FASTA sequences
  if (!(gene_id %in% names(YeastCDSFASTA))) {
    return(NA)  # Return NA if gene not found
  }
  
  # Extract the gene sequence
  gene_sequence <- as.character(YeastCDSFASTA[gene_id])
  
  # Calculate the CDS position of the proline
  arginine_cds_end <- arginine_cds_start + 2
  
  # Ensure positions are within bounds
  if (arginine_cds_end > nchar(gene_sequence)) {
    return(NA)
  }
  
  # Extract the codon
  codon <- substr(gene_sequence, arginine_cds_start, arginine_cds_end)
  return(codon)
}

ArgSer.UniqueCodon <- WholeProteome.DDA.Tidy.Filter %>%
  mutate(NumArg = str_count(Peptide.Sequence, "R")) %>%
  filter(NumMuts == 0 | NumMuts == 1) %>%
  filter(Sample == "Ctrl" | Sample == "ArgSer") %>%
  filter(NumArg == 1) %>%
  filter(Mapped.Genes == "") %>%
  filter(ProAla.Count == 0) %>%
  filter(ProSer.Count == 0) %>%
  select(Peptide.Sequence, Start, End, Gene) %>%
  unique() %>%
  mutate(ArgininePosition = str_locate(Peptide.Sequence, "R")[,"start"],
         ProteinPosition = (Start + ArgininePosition - 1),
         CDSStart = (ProteinPosition*3)-2) %>%
  rowwise() %>%
  mutate(codon = extract_codon_for_arginine(cur_data(), YeastCDSFASTA))

ArgSer.Codon <- WholeProteome.DDA.Tidy.Filter %>%
  mutate(NumArg = str_count(Peptide.Sequence, "R")) %>%
  filter(NumMuts == 0 | NumMuts == 1) %>%
  filter(Sample == "Ctrl" | Sample == "ArgSer") %>%
  filter(NumArg == 1) %>%
  filter(Mapped.Genes == "") %>%
  filter(ProAla.Count == 0) %>%
  filter(ProSer.Count == 0) %>%
  left_join(ArgSer.UniqueCodon %>% select(-ArgininePosition, -ProteinPosition, -CDSStart), by = c("Peptide.Sequence", "Start", "End", "Gene")) %>% 
  na.omit() %>%
  mutate(PepType = ifelse(NumMuts == 0, "WT", "Mutant"),
         Stripped.Sequence.ArgSer = str_replace_all(Modified.Sequence, "R\\[-69.0691\\]", "R"))

#Separate mistranslated peptides from WT peptides
WT.Peptides.ArgSer <- ArgSer.Codon %>%
  ungroup() %>%
  filter(PepType == "WT")

Mistranslated.Peptides.ArgSer <- ArgSer.Codon %>%
  ungroup() %>%
  filter(PepType == "Mutant")

# Join to filter out any mutant peptide with no WT counterpart
Filtered.MistranslatedPeptides.ArgSer <- left_join(WT.Peptides.ArgSer, Mistranslated.Peptides.ArgSer, by = c("Sample", "Replicate", "Stripped.Sequence.ArgSer", "codon"), suffix = c(".wt", ".mut"))

# Calculate the percent mistranslation based on peptide counts
RtoS.PercentMistranslation <- Filtered.MistranslatedPeptides.ArgSer  %>%
  filter(codon == "CGT" |
           codon == "CGG" |
           codon == "CGC" |
           codon == "CGA" |
           codon == "AGA" |
           codon == "AGG") %>%
  filter(str_detect(Stripped.Sequence.ArgSer, "R")) %>%
  dplyr::select(Replicate, Sample, PepType.wt, PepType.mut, codon) %>%
  group_by(Replicate, Sample, codon) %>%
  summarise(countWT = n(), countMutant = sum(!is.na(PepType.mut))) %>%
  mutate(PercentMistranslation = (countMutant/countWT)*100,
         MistranslationEvent = "R to S")

RtoS.PercentMistranslation$Sample <- factor(RtoS.PercentMistranslation$Sample, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))

ggplot(RtoS.PercentMistranslation, aes(x = codon, y = PercentMistranslation, fill = Sample, color = Sample)) +
  geom_bar(stat = "summary", alpha = 0.35, position = position_dodge(width = 0.9)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 2) +
  ylab("% Arginine to Serine") +
  xlab("") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_matt() +
  scale_y_continuous(lim = c(0,4), expand = c(0,0)) +
  xlab("Codon") +
  scale_fill_manual(values = c("#6496bf", "#a165a7"), labels = c("EV", expression(tRNA[UCU,G26A]^{"Ser"})))

ggsave("CodonUsage_ArgSer.pdf", device = cairo_pdf, width = 6, height = 4)

t.test(PercentMistranslation ~ Sample, data = RtoS.PercentMistranslation %>% filter(codon == "AGA"), var.equal = TRUE) # 8.005e-12
t.test(PercentMistranslation ~ Sample, data = RtoS.PercentMistranslation %>% filter(codon == "AGG"), var.equal = TRUE) # 0.3951
t.test(PercentMistranslation ~ Sample, data = RtoS.PercentMistranslation %>% filter(codon == "CGA"), var.equal = TRUE) # 0
t.test(PercentMistranslation ~ Sample, data = RtoS.PercentMistranslation %>% filter(codon == "CGC"), var.equal = TRUE) # 0.3409
t.test(PercentMistranslation ~ Sample, data = RtoS.PercentMistranslation %>% filter(codon == "CGG"), var.equal = TRUE) # 0
t.test(PercentMistranslation ~ Sample, data = RtoS.PercentMistranslation %>% filter(codon == "CGT"), var.equal = TRUE) # 0.9897
