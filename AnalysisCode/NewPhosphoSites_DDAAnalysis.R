### DDA Analysis of new potential phosphosites created due to mis-incorporation of serine at either proline or arginine residues
### Input file is the NewPhosphosites_TrypsinLysC_DDA_MSFragger_CombinedModifiedPeptide.tsv output file from MSFragger search of DDA phosphoproteome measurements
## Also requires "Yeast_SGD_all_protSeqs.txt" for the phosphomotif enrichment analysis

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
library(foreach)

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
NewPhospho.DDA <- read.csv("combined_modified_peptide.tsv", sep = "\t")

# Convert data to tidy format and filter
NewPhospho.DDA.Tidy <- NewPhospho.DDA %>%
  dplyr::select(Peptide.Sequence, Modified.Sequence, Start, End, Gene, Mapped.Genes, ArgSer_1.Spectral.Count:ProSer_6.Spectral.Count, Protein.ID) %>%
  gather(value = "Spectral.Count", key = "Sample", ArgSer_1.Spectral.Count:ProSer_6.Spectral.Count) %>%
  mutate(ProPhosphoSer = ifelse(str_detect(Modified.Sequence, "P\\[69.9456\\]"), T, F),
         ArgPhosphoSer = ifelse(str_detect(Modified.Sequence, "R\\[10.8972\\]"), T, F)) %>%
  mutate(ProPhosphoSer.Count = str_count(Modified.Sequence, "P\\[69.9456\\]"),
         ArgPhosphoSer.Count = str_count(Modified.Sequence, "R\\[10.8972\\]")) %>%
  separate(Sample, into = c("Sample", "Replicate"))

# Filter out peptides with two or more types of mistranslation and for peptides with no signal intensity
NewPhospho.DDA.Tidy.Filter <- NewPhospho.DDA.Tidy %>%
  filter((ProPhosphoSer + ArgPhosphoSer) == 1,
         Spectral.Count > 0) %>%
  mutate(NumMuts = ProPhosphoSer.Count + ArgPhosphoSer.Count)

# New phosphosite must be seen in all six replicates
NewPhospho.DDA.Tidy.Filter.Counts <- NewPhospho.DDA.Tidy.Filter %>%
  group_by(Modified.Sequence, Sample, ProPhosphoSer.Count, ArgPhosphoSer.Count) %>%
  summarise(Count = n())

NewPhospho.DDA.Tidy.FilterCounts <- NewPhospho.DDA.Tidy.Filter %>%
  left_join(NewPhospho.DDA.Tidy.Filter.Counts) %>%
  filter(Count >= 6)

## Plot number of modifications per peptide for proline to phosphoserine (Figure 7A)
NumMuts.ProPhosphoSer <- NewPhospho.DDA.Tidy.FilterCounts %>%
  filter(ProPhosphoSer.Count > 0) %>%
  ungroup() %>%
  dplyr::select(Modified.Sequence, Protein.ID, Sample) %>%
  unique()

NumMuts.ProPhosphoSer$Sample <- factor(NumMuts.ProPhosphoSer$Sample, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))
condition_colors <- c("Ctrl" = "#6496bf", "ProAla" = "#58c87e", "ProSer" = "#f8c754", "ArgSer" = "#a165a7")
condition_labels <- c("Ctrl" = "EV", "ProAla" = "P\u2192A", "ProSer" = "P\u2192S", "ArgSer" = "R\u2192S")

ggplot(NumMuts.ProPhosphoSer, aes(x = Sample, fill = Sample)) + 
  geom_bar(color = "black", position = "dodge") + 
  ylab("Unique Peptides") + 
  xlab("") +
  theme_matt() +
  scale_y_continuous(expand = c(0,0), limits = c(0,250))  +
  scale_fill_manual(values = condition_colors) +
  scale_x_discrete(labels = condition_labels) +
  theme(legend.position = "none") +
  ggtitle("Proline to Phosphoserine")

ggsave("../FINAL_Figures/Figure7/UniquePeptidesProlinePhosphoSerine.pdf", device = cairo_pdf, width = 3, height = 4)

## Plot number of modifications per peptide for arginine to phosphoserine (Figure 7A)
NumMuts.ArgPhosphoSer <- NewPhospho.DDA.Tidy.FilterCounts %>%
  filter(ArgPhosphoSer.Count > 0) %>%
  ungroup() %>%
  dplyr::select(Modified.Sequence, Protein.ID, Sample) %>%
  unique() %>%
  group_by(Sample) %>%
  summarise(Count = n()) %>%
  bind_rows(tibble(Sample = "ProSer", Count = 0))

NumMuts.ArgPhosphoSer$Sample <- factor(NumMuts.ArgPhosphoSer$Sample, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))
condition_colors <- c("Ctrl" = "#6496bf", "ProAla" = "#58c87e", "ProSer" = "#f8c754", "ArgSer" = "#a165a7")
condition_labels <- c("Ctrl" = "EV", "ProAla" = "P\u2192A", "ProSer" = "P\u2192S", "ArgSer" = "R\u2192S")

ggplot(NumMuts.ArgPhosphoSer, aes(x = Sample, y = Count, fill = Sample)) + 
  geom_bar(color = "black", stat = "identity") + 
  ylab("Unique Peptides") + 
  xlab("") +
  theme_matt() +
  scale_y_continuous(expand = c(0,0), limits = c(0,40)) +
  scale_fill_manual(values = condition_colors) +
  scale_x_discrete(labels = condition_labels) +
  theme(legend.position = "none") +
  ggtitle("Arginine to Phosphoserine")

ggsave("../FINAL_Figures/Figure7/UniquePeptidesArgininePhosphoSerine.pdf", device = cairo_pdf, width = 3, height = 4)

# Site summary
Site.Summary <- NewPhospho.DDA.Tidy.FilterCounts %>%
  ungroup() %>%
  select(Modified.Sequence, Protein.ID, Start, End, Gene, Sample, ProPhosphoSer, ArgPhosphoSer) %>%
  unique()

Summary.ProPhosphoSer <- Site.Summary %>%
  ungroup() %>%
  filter(ProPhosphoSer == T) %>%
  group_by(Sample) %>%
  summarise(Counts = n())

Summary.ArgPhosphoSer <- Site.Summary %>%
  ungroup() %>%
  filter(ArgPhosphoSer == T) %>%
  group_by(Sample) %>%
  summarise(Counts = n())

ArgSer.Peptides <- Site.Summary %>%
  filter(ArgPhosphoSer == T) %>%
  filter(Sample == "ArgSer") %>%
  select(Modified.Sequence) %>%
  unique() %>%
  mutate(SPMotif = ifelse(str_detect(Modified.Sequence, "R\\[10.8972\\]P"), T, F))

## Determine if there is a motif around observed phosphorylated mistranslation site (Figure 7D)
# Determine where mistranslation is in the protein, take 15 resiudes up and downstream as foreground
sequence <- read_tsv("..//annotations/Yeast_SGD_all_protSeqs.txt")

MistranslatedPositions <- NewPhospho.DDA.Tidy.Filter %>%
  filter(NumMuts == 1) %>%
  mutate(Modified.Sequence.NoMetOx = str_replace_all(Modified.Sequence, "M\\[15.9949\\]", "M"),
         Modified.Sequence.NonAc = str_replace_all(Modified.Sequence.NoMetOx, "n\\[42.0106\\]", ""),
         Modified.Sequence.NoCys = str_replace_all(Modified.Sequence.NonAc, "C\\[57.0215\\]", "C"),
         Modified.Sequence.NoPhosS = str_replace_all(Modified.Sequence.NoCys, "S\\[79.9663\\]", "S"),
         Modified.Sequence.NoPhosT = str_replace_all(Modified.Sequence.NoPhosS, "T\\[79.9663\\]", "T"),
         Modified.Sequence.NoPhosY = str_replace_all(Modified.Sequence.NoPhosT, "Y\\[79.9663\\]", "Y")) %>%
  mutate(mod.position.peptide = str_locate(Modified.Sequence.NoPhosY, "\\[")[,1] - 1,
         modified.aa = str_sub(Modified.Sequence.NoPhosY, mod.position.peptide, mod.position.peptide),
         mod.position.protein = Start - 1 + mod.position.peptide) %>%
  left_join(sequence %>% dplyr::rename(Gene = GeneName)) %>%
  drop_na(Sequence) %>%
  mutate(FrontExtra = "__________________________________________________",
         EndExtra = "__________________________________________________") %>%
  unite("Sequence", c(FrontExtra, Sequence, EndExtra), sep = "") %>%
  mutate(Seq.Window = str_sub(string = Sequence, start = (mod.position.protein-7+50), end = (mod.position.protein+7+50))) %>%
  mutate(length = str_length(Seq.Window)) %>%
  filter(length == 15)  %>%
  distinct(Gene, Protein.ID, Sample, modified.aa, mod.position.protein, Reference, Seq.Window)

# Create a background from all phosphorylated peptides detected with 15 residues up and downstream
# For proline mistranslation
Proline.Positions <- NewPhospho.DDA.Tidy.Filter %>%
  mutate(Modified.Sequence.NoMetOx = str_replace_all(Modified.Sequence, "M\\[15.9949\\]", "M"),
         Modified.Sequence.NonAc = str_replace_all(Modified.Sequence.NoMetOx, "n\\[42.0106\\]", ""),
         Modified.Sequence.NoCys = str_replace_all(Modified.Sequence.NonAc, "C\\[57.0215\\]", "C"),
         Modified.Sequence.NoProSer = str_replace_all(Modified.Sequence.NoCys, "P\\[69.9456\\]", "P"),
         Modified.Sequence.NoArgSer = str_replace_all(Modified.Sequence.NoProSer, "R\\[10.8972\\]", "R"),
         Modified.Sequence.NoPhosS = str_replace_all(Modified.Sequence.NoArgSer, "S\\[79.9663\\]", "S"),
         Modified.Sequence.NoPhosT = str_replace_all(Modified.Sequence.NoPhosS, "T\\[79.9663\\]", "T"),
         Modified.Sequence.NoPhosY = str_replace_all(Modified.Sequence.NoPhosT, "Y\\[79.9663\\]", "Y")) %>%
  filter(str_detect(Modified.Sequence.NoPhosY, "P")) %>%
  rowwise() %>%
  mutate(Pro.Positions = list(str_locate_all(Modified.Sequence.NoPhosY, "P")[[1]][, 1])) %>%
  unnest(Pro.Positions)%>%
  mutate(modified.aa = str_sub(Modified.Sequence.NoPhosY, Pro.Positions, Pro.Positions),
         mod.position.protein = Start - 1 + Pro.Positions) %>%
  mutate(mod.position.peptide = str_locate(Modified.Sequence.NoPhosY, "P")[,1],
         modified.aa = str_sub(Modified.Sequence.NoPhosY, mod.position.peptide, mod.position.peptide),
         mod.position.protein = Start - 1 + mod.position.peptide) %>%
  left_join(sequence %>% dplyr::rename(Gene = GeneName)) %>%
  drop_na(Sequence) %>%
  mutate(FrontExtra = "__________________________________________________",
         EndExtra = "__________________________________________________") %>%
  unite("Sequence", c(FrontExtra, Sequence, EndExtra), sep = "") %>%
  mutate(Seq.Window = str_sub(string = Sequence, start = (mod.position.protein-7+50), end = (mod.position.protein+7+50))) %>%
  mutate(length = str_length(Seq.Window)) %>%
  filter(length == 15) %>%
  distinct(Gene, Protein.ID, Sample, modified.aa, mod.position.protein, Reference, Seq.Window)

# For arginine mistranslation
Arginine.Positions <- NewPhospho.DDA.Tidy.Filter %>%
  mutate(Modified.Sequence.NoMetOx = str_replace_all(Modified.Sequence, "M\\[15.9949\\]", "M"),
         Modified.Sequence.NonAc = str_replace_all(Modified.Sequence.NoMetOx, "n\\[42.0106\\]", ""),
         Modified.Sequence.NoCys = str_replace_all(Modified.Sequence.NonAc, "C\\[57.0215\\]", "C"),
         Modified.Sequence.NoProSer = str_replace_all(Modified.Sequence.NoCys, "P\\[69.9456\\]", "P"),
         Modified.Sequence.NoArgSer = str_replace_all(Modified.Sequence.NoProSer, "R\\[10.8972\\]", "R"),
         Modified.Sequence.NoPhosS = str_replace_all(Modified.Sequence.NoArgSer, "S\\[79.9663\\]", "S"),
         Modified.Sequence.NoPhosT = str_replace_all(Modified.Sequence.NoPhosS, "T\\[79.9663\\]", "T"),
         Modified.Sequence.NoPhosY = str_replace_all(Modified.Sequence.NoPhosT, "Y\\[79.9663\\]", "Y")) %>%
  filter(str_detect(Modified.Sequence.NoPhosY, "R")) %>%
  rowwise() %>%
  mutate(Arg.Positions = list(str_locate_all(Modified.Sequence.NoPhosY, "R")[[1]][, 1])) %>%
  unnest(Arg.Positions)%>%
  mutate(modified.aa = str_sub(Modified.Sequence.NoPhosY, Arg.Positions, Arg.Positions),
         mod.position.protein = Start - 1 + Arg.Positions) %>%
  left_join(sequence %>% dplyr::rename(Gene = GeneName)) %>%
  drop_na(Sequence) %>%
  mutate(FrontExtra = "__________________________________________________",
         EndExtra = "__________________________________________________") %>%
  unite("Sequence", c(FrontExtra, Sequence, EndExtra), sep = "") %>%
  mutate(Seq.Window = str_sub(string = Sequence, start = (mod.position.protein-7+50), end = (mod.position.protein+7+50))) %>%
  mutate(length = str_length(Seq.Window)) %>%
  filter(length == 15) %>%
  distinct(Gene, Protein.ID, Sample, modified.aa, mod.position.protein, Reference, Seq.Window)

# Compare background and foreground sequences around mistranslated site for proline to alanine
ProSer.foreground <- MistranslatedPositions %>%
  filter(Sample == "ProSer",
         modified.aa == "P")

ProSer.background <- Proline.Positions %>%
  filter(Sample == "ProSer",
         modified.aa == "P")

ProSer.foreground.seqs <- ProSer.foreground$Seq.Window
ProSer.background.seqs <- ProSer.background$Seq.Window

ProSer.motif <- motifx(ProSer.foreground.seqs, ProSer.background.seqs, central.res = "P", min.seqs = 50, pval.cutoff = 1e-6)

# Compare background and foreground sequences around mistranslated site for arginine to serine
ArgSer.foreground <- MistranslatedPositions %>%
  filter(Sample == "ArgSer",
         modified.aa == "R")

ArgSer.background <- Arginine.Positions %>%
  filter(Sample == "ArgSer",
         modified.aa == "R")

ArgSer.foreground.seqs <- ArgSer.foreground$Seq.Window
ArgSer.background.seqs <- ArgSer.background$Seq.Window

ArgSer.motif <- motifx(ArgSer.foreground.seqs, ArgSer.background.seqs, central.res = "R", min.seqs = 50, pval.cutoff = 1e-6)

# Sequence logo for proline to serine sites
ProSer.df <- ProSer.foreground %>% 
  dplyr::select(PTM_window = Seq.Window) %>% 
  mutate(test = 1, regulated = "up") %>% 
  full_join(ProSer.background %>% dplyr::select(PTM_window = Seq.Window) %>% 
              mutate(test =1, regulated = "ns"))

ProSer.res_exp <- as.data.table(ProSer.df)

ProSer.seq_nchar <- c("desired" = 15L, "current" = ProSer.res_exp[, nchar(PTM_window[1])])
ProSer.flank_nchar <- (ProSer.seq_nchar-1)/2
ProSer.seq_sub <- c((1+(ProSer.seq_nchar[2]-ProSer.seq_nchar[1])/2), (ProSer.seq_nchar[2]-(ProSer.seq_nchar[2]-ProSer.seq_nchar[1])/2))
ProSer.res_exp[, PTM_window := str_sub(PTM_window, ProSer.seq_sub[1], ProSer.seq_sub[2])]

## Create sequence logos for each type of mistranslation
# Function for plotting dagLogo with ggseqlogo
dagLogo_ggseqlogo <-
  function(seq_reg, seq_ns, background = c("manual", "ztest", "fisher")){
    #check variables
    if(!is.character(seq_reg)) stop("'seq_reg' needs to be a character vector")
    if(!is.character(seq_ns) & background=="manual") stop("'seq_ns' needs to be a character vector")
    logo_type <- match.arg(background)
    
    #manually create S4 dagPeptides foreground object
    data_fg <- data.frame(anchor = str_sub(seq_reg, ProSer.flank_nchar[1]+1, ProSer.flank_nchar[1]+1),
                          upstream = str_sub(seq_reg, 1, ProSer.flank_nchar[1]),
                          downstream = str_sub(seq_reg, ProSer.flank_nchar[1]+2, ProSer.seq_nchar[1]))
    peptides_fg <- matrix(data = unlist(strsplit(seq_reg, "")), ncol = ProSer.seq_nchar[1], byrow = T)
    test_fg <- new("dagPeptides", data = data_fg, peptides = peptides_fg, upstreamOffset = 7, downstreamOffset = 7, type = "fasta")
    
    if(background %in% c("ztest", "fisher")){
      
      #prepare yeast proteome
      proteome <- prepareProteome(fasta = "..//..//0_data/cerevisiae_orf_trans_all.fasta")
      proteome@species <- "Saccharomyces cerevisiae"
      
      #prepare background object
      if(background == "fisher"){
        test_bg <- buildBackgroundModel(test_fg, background = "wholeProteome",
                                        proteome = proteome, testType = "fisher")
      }
      if(background == "ztest"){
        test_bg <- buildBackgroundModel(test_fg, background = "wholeProteome", 
                                        proteome = proteome, testType = "ztest")
      }
    }
    
    if(background == "manual"){
      #manually build background object
      peptides_bg <- matrix(data = unlist(strsplit(seq_ns, "")), ncol = ProSer.seq_nchar[1], byrow = T)
      test_bg <-
        new("dagBackground", background = list(peptides_bg),
            numSubsamples = 1L, testType = "fisher")
    }
    
    #logo enrichment test
    test_logo <- testDAU(test_fg, dagBackground = test_bg)
    # dagHeatmap(test_logo)
    # dagLogo(test_logo)
    return(test_logo)
  }

ProSer.logo_list <- foreach(i = ProSer.res_exp[, unique(test)], .combine = c) %do% {
  #extract regulated sequences based on test regulatory column
  seq_reg <- unique(ProSer.res_exp[test==1 & regulated=="up", PTM_window])
  seq_ns <- unique(ProSer.res_exp[test==1 & regulated!="up", PTM_window])
  
  #create and return plot
  temp <- dagLogo_ggseqlogo(seq_reg, seq_ns, background = "manual")
  p <- list(temp@difference * (temp@pvalue <= 0.05) * 100)
  names(p) <- i
  return(p)
}

ProSer.logo_list <- ProSer.logo_list[order(as.numeric(names(ProSer.logo_list)))]

ggplot() + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = (ProSer.seq_nchar[1]+1)/2, linetype = "dashed") +
  geom_logo(ProSer.logo_list, method='custom', seq_type='aa') + 
  ylab("Enrichment Score") +
  ylim(-20,20)+
  scale_x_continuous(breaks = seq(1, ProSer.seq_nchar[1], 1),
                     labels = seq(-((ProSer.seq_nchar[1]-1)/2), ((ProSer.seq_nchar[1]-1)/2), 1)) +
  theme_matt() +
  theme(panel.grid.major = element_blank(), legend.position = "bottom")

ggsave("../FINAL_Figures/Figure6/ProSerMistranslation_SequenceLogo.pdf", units = "in", width = 6, height = 4)

# Sequence logo for arginine to serine sites
ArgSer.df <- ArgSer.foreground %>% 
  dplyr::select(PTM_window = Seq.Window) %>% 
  mutate(test =1, regulated = "up") %>% 
  full_join(ArgSer.background %>% dplyr::select(PTM_window = Seq.Window) %>% 
              mutate(test =1, regulated = "ns"))

ArgSer.res_exp <- as.data.table(ArgSer.df)

ArgSer.seq_nchar <- c("desired" = 15L, "current" = ArgSer.res_exp[, nchar(PTM_window[1])])
ArgSer.flank_nchar <- (ArgSer.seq_nchar-1)/2
ArgSer.seq_sub <- c((1+(ArgSer.seq_nchar[2]-ArgSer.seq_nchar[1])/2), (ArgSer.seq_nchar[2]-(ArgSer.seq_nchar[2]-ArgSer.seq_nchar[1])/2))
ArgSer.res_exp[, PTM_window := str_sub(PTM_window, ArgSer.seq_sub[1], ArgSer.seq_sub[2])]

## Create sequence logos for each type of mistranslation
# Function for plotting dagLogo with ggseqlogo
dagLogo_ggseqlogo <-
  function(seq_reg, seq_ns, background = c("manual", "ztest", "fisher")){
    #check variables
    if(!is.character(seq_reg)) stop("'seq_reg' needs to be a character vector")
    if(!is.character(seq_ns) & background=="manual") stop("'seq_ns' needs to be a character vector")
    logo_type <- match.arg(background)
    
    #manually create S4 dagPeptides foreground object
    data_fg <- data.frame(anchor = str_sub(seq_reg, ArgSer.flank_nchar[1]+1, ArgSer.flank_nchar[1]+1),
                          upstream = str_sub(seq_reg, 1, ArgSer.flank_nchar[1]),
                          downstream = str_sub(seq_reg, ArgSer.flank_nchar[1]+2, ArgSer.seq_nchar[1]))
    peptides_fg <- matrix(data = unlist(strsplit(seq_reg, "")), ncol = ArgSer.seq_nchar[1], byrow = T)
    test_fg <- new("dagPeptides", data = data_fg, peptides = peptides_fg, upstreamOffset = 7, downstreamOffset = 7, type = "fasta")
    
    if(background %in% c("ztest", "fisher")){
      
      #prepare yeast proteome
      proteome <- prepareProteome(fasta = "..//..//0_data/cerevisiae_orf_trans_all.fasta")
      proteome@species <- "Saccharomyces cerevisiae"
      
      #prepare background object
      if(background == "fisher"){
        test_bg <- buildBackgroundModel(test_fg, background = "wholeProteome",
                                        proteome = proteome, testType = "fisher")
      }
      if(background == "ztest"){
        test_bg <- buildBackgroundModel(test_fg, background = "wholeProteome", 
                                        proteome = proteome, testType = "ztest")
      }
    }
    
    if(background == "manual"){
      #manually build background object
      peptides_bg <- matrix(data = unlist(strsplit(seq_ns, "")), ncol = ArgSer.seq_nchar[1], byrow = T)
      test_bg <-
        new("dagBackground", background = list(peptides_bg),
            numSubsamples = 1L, testType = "fisher")
    }
    
    #logo enrichment test
    test_logo <- testDAU(test_fg, dagBackground = test_bg)
    # dagHeatmap(test_logo)
    # dagLogo(test_logo)
    return(test_logo)
  }


ArgSer.logo_list <- foreach(i = ArgSer.res_exp[, unique(test)], .combine = c) %do% {
  #extract regulated sequences based on test regulatory column
  seq_reg <- unique(ArgSer.res_exp[test==i & regulated=="up", PTM_window])
  seq_ns <- unique(ArgSer.res_exp[test==i & regulated!="up", PTM_window])
  
  #create and return plot
  temp <- dagLogo_ggseqlogo(seq_reg, seq_ns, background = "manual")
  p <- list(temp@difference * (temp@pvalue <= 0.05) * 100)
  names(p) <- i
  return(p)
}

ArgSer.logo_list <- ArgSer.logo_list[order(as.numeric(names(ArgSer.logo_list)))]

ggplot() + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = (ArgSer.seq_nchar[1]+1)/2, linetype = "dashed") +
  geom_logo(ArgSer.logo_list, method='custom', seq_type='aa') + 
  ylab("Enrichment Score") +
  ylim(-20,20)+
  scale_x_continuous(breaks = seq(1, ArgSer.seq_nchar[1], 1),
                     labels = seq(-((ArgSer.seq_nchar[1]-1)/2), ((ArgSer.seq_nchar[1]-1)/2), 1)) +
  theme_matt() +
  theme(panel.grid.major = element_blank(), legend.position = "bottom")

ggsave("ArgSerMistranslation_SequenceLogo.pdf", units = "in", width = 6, height = 4)