## Data analaysis and figures for whole proteome and phosphoproteome analysis of mistranslating yeast strains
## Input for this analysis is MSstats protein and phosphosite abundance summary data and MSstats differential abundance for protein and phosphosite
## Also requires yeast fasta from uniprot (Yeast_uniprot.csv) and currated yeast kinase-substrate data (kinase_phos_substrate_yeast_2.csv) for the kinase-substrate enrichment analysis

# Load theme
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

# Load in required packages
library(tidyverse)
library(cowplot)
library(pheatmap)
library(ggpubr)
library(viridis)
library(ggforce)
library(Biostrings)

# Set working directory to MSstats output data
setwd("")

# Read in MSstats data
ProteinLevel <- read.csv("ProteinLevel_MSstatsSummary.csv")
PhosphoLevel <- read.csv("PhosphoLevel_MSstatsSummary.csv")
Proteome.DA <- read.csv("ProteinDifferentialAbundance_MSstats.csv")
Phospho.DA <- read.csv("ProteinAdjustedPhosphositeDifferentialAbundance_MSstats.csv")

# Analysis of the number of unique proteins and create plots (Figure 2B)
ProteinCount.Summary <- ProteinLevel %>%
  dplyr::select(GROUP, SUBJECT, Protein) %>%
  unique() %>%
  group_by(GROUP, SUBJECT) %>%
  summarise(ProteinCount = n())

ProteinCount.Summary$GROUP <- factor(ProteinCount.Summary$GROUP, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))

ggplot(ProteinCount.Summary, aes(x = GROUP, y = ProteinCount, fill = GROUP, group = SUBJECT)) +
  geom_col(color = "black", position = position_dodge(width = 0.9), width = 0.9) +
  ylab("Number of Proteins") + 
  xlab("") +
  scale_fill_manual(values = c("#6496bf", "#58c87e", "#f8c754", "#a165a7")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 5000)) +
  theme_matt() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("EV", "Pro\u2192Ala", "Pro\u2192Ser", "Arg\u2192Ser"))

ggsave("../FINAL_Figures/Figure2/Proteome_UniqueProteinsCounts.pdf", device = cairo_pdf, width = 5, height = 3.5, units = "in")

# Determine the average number of proteins identified per sample
# There are 4661 +/- 14.1 unique proteins identified on average across the samples
Summary.ProteinCount.Summary <- ProteinCount.Summary %>% 
  ungroup() %>%
  summarise(MeanCount = mean(ProteinCount),
            SDCount = sd(ProteinCount))

# Analysis of the number of unique phosphosites and create plots (Figure 2B)
PhosphoCount.Summary <- PhosphoLevel %>%
  dplyr::select(GROUP, SUBJECT, Protein) %>%
  unique() %>%
  group_by(GROUP, SUBJECT) %>%
  summarise(PhosphoCount = n())

PhosphoCount.Summary$GROUP <- factor(PhosphoCount.Summary$GROUP, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer"))

ggplot(PhosphoCount.Summary, aes(x = GROUP, y = PhosphoCount, fill = GROUP, group = SUBJECT)) +
  geom_col(color = "black", position = position_dodge(width = 0.9), width = 0.9) +
  ylab("Number of Phosphosites") + 
  xlab("") +
  scale_fill_manual(values = c("#6496bf", "#58c87e", "#f8c754", "#a165a7")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 11000)) +
  theme_matt() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("EV", "Pro\u2192Ala", "Pro\u2192Ser", "Arg\u2192Ser"))

ggsave("../FINAL_Figures/Figure2/Phospho_UniquePhosphoCounts.pdf", device = cairo_pdf, width = 5, height = 3.5, units = "in")

# Determine the average number of proteins identified per sample
# There are 10290 +/- 280 unique proteins identified on average across the samples
Summary.PhosphoCount.Summary <- PhosphoCount.Summary %>%
  ungroup() %>%
  summarise(MeanCount = mean(PhosphoCount),
            SDCount = sd(PhosphoCount))

# Determine CVs of protein level quantifications and make plot (Figure 2C)
CVSummarizedProtein <- ProteinLevel %>%
  distinct(Protein, GROUP, LogIntensities) %>%
  drop_na(LogIntensities) %>%
  mutate(RawAbundance = 2^LogIntensities) %>%
  group_by(Protein) %>%
  mutate(cv_combined = (stats::sd(RawAbundance)/mean(RawAbundance)) *100)%>%
  group_by(GROUP, Protein) %>%
  mutate(cv = (stats::sd(RawAbundance)/mean(RawAbundance)) * 100) %>%
  distinct(GROUP, Protein, .data$cv_combined, .data$cv) %>%
  drop_na() %>%
  group_by(GROUP) %>%
  mutate(median_cv = stats::median(.data$cv))%>%
  ungroup() %>%
  mutate(median_cv_combined = stats::median(.data$cv_combined)) %>%
  dplyr::select(-Protein, -c("cv_combined", "cv")) %>%
  distinct()

CVResultProtein <- ProteinLevel %>%
  distinct(Protein, GROUP, LogIntensities) %>%
  drop_na(LogIntensities) %>%
  mutate(RawAbundance = 2^LogIntensities) %>%
  group_by(Protein) %>%
  mutate(cv_combined = (stats::sd(RawAbundance)/mean(RawAbundance)) *100)%>%
  group_by(GROUP, Protein) %>%
  mutate(cv = (stats::sd(RawAbundance)/mean(RawAbundance)) * 100) %>%
  ungroup() %>%
  distinct(GROUP, Protein, .data$cv_combined, .data$cv) %>%
  drop_na() %>%
  pivot_longer(cols = starts_with("cv"), names_to = "type", values_to = "values") %>%
  mutate(type = ifelse(.data$type == "cv", paste0(GROUP), "combined"))%>%
  mutate(type = forcats::fct_relevel(as.factor(.data$type), "combined")) %>%
  dplyr::select(-GROUP) %>%
  group_by(.data$type) %>%
  mutate(median = stats::median(.data$values)) %>%
  distinct()

CVResultProtein$type <- factor(CVResultProtein$type, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer", "combined"))

ggplot(CVResultProtein, aes(x = type, y = values, fill = type)) +
  geom_violin(na.rm = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", na.rm = TRUE, alpha = 0.6, outlier.shape = NA) +
  labs(
    x = "",
    y = "CV (%)",
    fill = "Condition"
  ) +
  scale_fill_manual(values = c("#6496bf", "#58c87e", "#f8c754", "#a165a7", "darkgrey")) +
  theme_matt() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("EV", "P\u2192A", "P\u2192S", "R\u2192S", "All")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))

ggsave("../FINAL_Figures/Figure2/Protein_CVs.pdf", device = cairo_pdf, width = 4, height = 3.5, units = "in")

# Determine CVs of the Phosphosite quantification and make plots (Figure 2C)
CVSummarizedPhospho <- PhosphoLevel %>%
  distinct(Protein, GROUP, LogIntensities) %>%
  drop_na(LogIntensities) %>%
  mutate(RawAbundance = 2^LogIntensities) %>%
  group_by(Protein) %>%
  mutate(cv_combined = (stats::sd(RawAbundance)/mean(RawAbundance)) *100)%>%
  group_by(GROUP, Protein) %>%
  mutate(cv = (stats::sd(RawAbundance)/mean(RawAbundance)) * 100) %>%
  distinct(GROUP, Protein, .data$cv_combined, .data$cv) %>%
  drop_na() %>%
  group_by(GROUP) %>%
  mutate(median_cv = stats::median(.data$cv))%>%
  ungroup() %>%
  mutate(median_cv_combined = stats::median(.data$cv_combined)) %>%
  dplyr::select(-Protein, -c("cv_combined", "cv")) %>%
  distinct()

CVResultPhospho <- PhosphoLevel%>%
  distinct(Protein, GROUP, LogIntensities) %>%
  drop_na(LogIntensities) %>%
  mutate(RawAbundance = 2^LogIntensities) %>%
  group_by(Protein) %>%
  mutate(cv_combined = (stats::sd(RawAbundance)/mean(RawAbundance)) *100)%>%
  group_by(GROUP, Protein) %>%
  mutate(cv = (stats::sd(RawAbundance)/mean(RawAbundance)) * 100) %>%
  ungroup() %>%
  distinct(GROUP, Protein, .data$cv_combined, .data$cv) %>%
  drop_na() %>%
  pivot_longer(cols = starts_with("cv"), names_to = "type", values_to = "values") %>%
  mutate(type = ifelse(.data$type == "cv", paste0(GROUP), "combined"))%>%
  mutate(type = forcats::fct_relevel(as.factor(.data$type), "combined")) %>%
  dplyr::select(-GROUP) %>%
  group_by(.data$type) %>%
  mutate(median = stats::median(.data$values)) %>%
  distinct()

CVResultPhospho$type <- factor(CVResultPhospho$type, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer", "combined"))

ggplot(CVResultPhospho, aes(x = type, y = values, fill = type)) +
  geom_violin(na.rm = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", na.rm = TRUE, alpha = 0.6, outlier.shape = NA) +
  labs(
    x = "",
    y = "CV (%)",
    fill = "Condition"
  ) +
  scale_fill_manual(values = c("#6496bf", "#58c87e", "#f8c754", "#a165a7", "darkgrey")) +
  theme_matt() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("EV", "P\u2192A", "P\u2192S", "R\u2192S", "All")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,200), breaks = c(0,50,100,150,200))

ggsave("../FINAL_Figures/Figure2/Phospho_CVs.pdf", device = cairo_pdf, width = 4, height = 3.5, units = "in")

# PCA of protein level data and make plots (Figure 2D)
pca_input <- ProteinLevel %>%
  distinct(originalRUN, Protein, LogIntensities) %>%
  group_by(originalRUN, Protein) %>%
  summarise(intensity = sum(LogIntensities)) %>%
  pivot_wider(names_from = originalRUN, values_from = intensity) %>%
  tidyr::drop_na() %>%
  dplyr::select(-Protein) %>%
  t(.)

annotation <- ProteinLevel %>%
  distinct(originalRUN, GROUP)

pca <- prcomp(pca_input, center = TRUE)

pca_df <- as.data.frame(pca$x) %>%
  mutate(originalRUN := factor(row.names(.))) %>%
  left_join(annotation, by = rlang::as_name("originalRUN"))

pca_sdev_df <- as.data.frame(pca$sdev)

pca_sdev_df <- pca_sdev_df %>%
  mutate(
    percent_variance = (pca$sdev^2 / sum(pca$sdev^2) * 100),
    dimension = row.names(.)
  ) %>%
  mutate(dimension = factor(.data$dimension,
                            levels = unique(stringr::str_sort(.data$dimension, numeric = TRUE))
  ))

ggplot(pca_df, aes(x = PC1, y = PC2, col = factor(GROUP, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer")), fill = factor(GROUP, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer")))) +
  geom_point(size = 5) +
  labs(x = paste0("PC1: ", round(pca_sdev_df$percent_variance[1],1), "%"),
       y = paste0("PC2: ", round(pca_sdev_df$percent_variance[2],1), "%"),
       color = ""
  ) +
  scale_color_manual(values= c("#6496bf", "#58c87e", "#f8c754", "#a165a7"), labels = c("EV", "Pro\u2192Ala", "Pro\u2192Ser", "Arg\u2192Ser")) +
  scale_fill_manual(values= c("#6496bf", "#58c87e", "#f8c754", "#a165a7"), labels = c("EV", "Pro\u2192Ala", "Pro\u2192Ser", "Arg\u2192Ser")) +
  theme_matt() +
  theme(legend.position = "none") +
  stat_ellipse(type = "norm", level = 0.68, geom = "polygon", alpha = 0.2)

ggsave("../FINAL_Figures/Figure2/Proteome_PCA.pdf", device = cairo_pdf, width = 4, height = 4, units = "in")

# PCA of phosphosite data and make plots (Figure 2D)
pca_input <- PhosphoLevel %>%
  distinct(originalRUN, Protein, LogIntensities) %>%
  group_by(originalRUN, Protein) %>%
  summarise(intensity = sum(LogIntensities)) %>%
  pivot_wider(names_from = originalRUN, values_from = intensity) %>%
  tidyr::drop_na() %>%
  dplyr::select(-Protein) %>%
  t(.)

annotation <-  PhosphoLevel %>%
  distinct(originalRUN, GROUP)

pca <- prcomp(pca_input, center = TRUE, scale = T)

pca_df <- as.data.frame(pca$x) %>%
  mutate(originalRUN := factor(row.names(.))) %>%
  left_join(annotation, by = rlang::as_name("originalRUN"))

pca_sdev_df <- as.data.frame(pca$sdev)

pca_sdev_df <- pca_sdev_df %>%
  mutate(
    percent_variance = (pca$sdev^2 / sum(pca$sdev^2) * 100),
    dimension = row.names(.)
  ) %>%
  mutate(dimension = factor(.data$dimension,
                            levels = unique(stringr::str_sort(.data$dimension, numeric = TRUE))
  ))

ggplot(pca_df, aes(x = PC1, y = PC2, col = factor(GROUP, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer")), fill = factor(GROUP, levels = c("Ctrl", "ProAla", "ProSer", "ArgSer")))) +
  geom_point(size = 5) +
  labs(x = paste0("PC1: ", round(pca_sdev_df$percent_variance[1],1), "%"),
       y = paste0("PC2: ", round(pca_sdev_df$percent_variance[2],1), "%"),
       color = ""
  ) +
  scale_color_manual(values= c("#6496bf", "#58c87e", "#f8c754", "#a165a7"), labels = c("EV", "Pro\u2192Ala", "Pro\u2192Ser", "Arg\u2192Ser")) +
  scale_fill_manual(values= c("#6496bf", "#58c87e", "#f8c754", "#a165a7"), labels = c("EV", "Pro\u2192Ala", "Pro\u2192Ser", "Arg\u2192Ser")) +
  theme_matt() +
  theme(legend.position = "none") +
  stat_ellipse(type = "norm", level = 0.68, geom = "polygon", alpha = 0.2)

ggsave("../FINAL_Figures/Figure2/Phospho_PCA.pdf", device = cairo_pdf, width = 4, height = 4, units = "in")


# Create volcano plots from the differential abundance protein analysis (Figure 3A)
ProAlaVolcano <- ggplot(Proteome.DA %>% filter(Label == "ProAla_vs_Ctrl"), aes (x = log2FC, y = -log10(adj.pvalue), text = Protein, label = Protein)) +
  geom_point(color = "lightgrey") +
  geom_point(data = Proteome.DA %>% filter(Label == "ProAla_vs_Ctrl") %>% filter(log2FC > 0.6 | log2FC < -0.6) %>% filter(adj.pvalue < 0.05), color = "#58c87e", size = 3) + 
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = -0.6, linetype = "dashed", color = "darkgrey") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  xlab("") +
  ylab("") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  theme_matt() +
  ggtitle("Pro\u2192Ala") + 
  theme(plot.title = element_text(size = 12, face = "bold"))

ProSerVolcano <- ggplot(Proteome.DA %>% filter(Label == "ProSer_vs_Ctrl"), aes (x = log2FC, y = -log10(adj.pvalue), text = Protein, label = Protein)) +
  geom_point(color = "lightgrey") +
  geom_point(data = Proteome.DA %>% filter(Label == "ProSer_vs_Ctrl") %>% filter(log2FC > 0.6 | log2FC < -0.6) %>% filter(adj.pvalue < 0.05), color = "#f8c754", size = 3) + 
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = -0.6, linetype = "dashed", color = "darkgrey") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  xlab("") +
  ylab("") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  theme_matt() +
  ggtitle("Pro\u2192Ser") + 
  theme(plot.title = element_text(size = 12, face = "bold"))

ArgSerVolcano <- ggplot(Proteome.DA %>% filter(Label == "ArgSer_vs_Ctrl"), aes (x = log2FC, y = -log10(adj.pvalue), text = Protein, label = Protein)) +
  geom_point(color = "grey") +
  geom_point(data = Proteome.DA %>% filter(Label == "ArgSer_vs_Ctrl") %>% filter(log2FC > 0.6 | log2FC < -0.6) %>% filter(adj.pvalue < 0.05), color = "#a165a7", size = 3) + 
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = -0.6, linetype = "dashed", color = "darkgrey") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  xlab("") +
  ylab("") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  theme_matt() +
  ggtitle("Arg\u2192Ser") + 
  theme(plot.title = element_text(size = 12, face = "bold")) 

plot_grid(ProAlaVolcano, ProSerVolcano, ArgSerVolcano, nrow = 1, scale = 0.999)
ggsave("../FINAL_Figures/Figure3/WholeProteomeVolcano.pdf", device = cairo_pdf, width = 10, height = 3.5, units = "in")

# There are 142 upregulated proteins in ProSer
ProSer.Up <- Proteome.DA %>%
  filter(Label == "ProSer_vs_Ctrl") %>%
  filter(adj.pvalue <= 0.05) %>%
  filter(log2FC > 0.6)

# There are 134 downregulated proteins in ProSer
ProSer.Dn <- Proteome.DA %>%
  filter(Label == "ProSer_vs_Ctrl") %>%
  filter(adj.pvalue <= 0.05) %>%
  filter(log2FC < -0.6)

# There are 32 upregulated proteins in ArgSer
ArgSer.Up <- Proteome.DA %>%
  filter(Label == "ArgSer_vs_Ctrl") %>%
  filter(adj.pvalue <= 0.05) %>%
  filter(log2FC > 0.6)

# There are 26 downregulated proteins in ArgSer
ArgSer.Dn <- Proteome.DA %>%
  filter(Label == "ArgSer_vs_Ctrl") %>%
  filter(adj.pvalue <= 0.05) %>%
  filter(log2FC < -0.6)

# There are 34 upregulated proteins in ProAla
ProAla.Up <- Proteome.DA %>%
  filter(Label == "ProAla_vs_Ctrl") %>%
  filter(adj.pvalue <= 0.05) %>%
  filter(log2FC > 0.6)

# There are 31 downregulated proteins in ProAla
ProAla.Dn <- Proteome.DA %>%
  filter(Label == "ProAla_vs_Ctrl") %>%
  filter(adj.pvalue <= 0.05) %>%
  filter(log2FC < -0.6)

# How many differential abundant proteins across all mistranslating strains?
AffectedProteins <- Proteome.DA %>% filter(adj.pvalue <= 0.05) %>% filter(log2FC > 0.6 | log2FC < -0.6)
length(unique(AffectedProteins$Protein)) # 316 proteins

## Make individual plots of different groups of protein quality control proteins (Figure 3C)
# Heat shock and chaperone proteins
# Hsp70 proteins
Hsp70.proteins <- c("SSA1", "SSA2", "SSA3", "SSA4", "SSB1", "SSB2", "SSE1", "SSE2", "SSZ1", "KAR2", "LHS1", "SSC1", "SSQ1", "ECM10")

Hsp70.proteins.FC <- data.frame(Proteome.DA) %>%
  filter(Protein %in% Hsp70.proteins) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("SSA1", "SSA2", "SSA3", "SSA4", "SSB1", "SSB2", "SSE1", "SSE2", "SSZ1", "KAR2", "LHS1", "SSC1", "SSQ1")))))

ggplot(Hsp70.proteins.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.01, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)")) +
  ggtitle("Hsp70s")

ggsave("../FINAL_Figures/Figure3/Hsp70Abundance.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 20)

# Hsp40s
Hsp40.proteins <- c("YDJ1", "XDJ1", "APJ1", "SIS1", "DJP1", "ZUO1", "SWA2", "JJJ1", "JJJ2", "JJJ3", "CAJ1", "CWC23", "MDJ1", "MDJ2", "PAM18", "JAC1", "JID1", "SCJ1", "HLJ1", "JEM1", "SEC63", "ERJ5")

Hsp40.proteins.FC <- data.frame(Proteome.DA) %>%
  filter(Protein %in% Hsp40.proteins) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("YDJ1", "XDJ1", "APJ1", "SIS1", "DJP1", "ZUO1", "SWA2", "JJJ1", "CAJ1", "CWC23", "MDJ1", "PAM18", "JAC1", "SCJ1", "HLJ1", "JEM1", "SEC63", "ERJ5")))))

ggplot(Hsp40.proteins.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.01, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)")) +
  ggtitle("Hsp40s")

ggsave("../FINAL_Figures/Figure3/Hsp40Abundance.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 20)

# Hsp90s
Hsp90.proteins <- c("HSP82", "HSC82", "STI1", "CDC37", "CNS1", "SBA1", "CPR6", "CPR7", "HCH1", "AHA1")

Hsp90.proteins.FC <- data.frame(Proteome.DA) %>%
  filter(Protein %in% Hsp90.proteins) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("HSP82", "HSC82", "STI1", "CDC37", "CNS1", "SBA1", "CPR6", "CPR7", "HCH1", "AHA1")))))

ggplot(Hsp90.proteins.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.01, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)")) +
  ggtitle("Hsp90s")

ggsave("../FINAL_Figures/Figure3/Hsp90Abundance.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 20)

# Other heatshock
Other.proteins <- c("HSP104", "HSP42", "HSP26")

Other.proteins.FC <- data.frame(Proteome.DA) %>%
  filter(Protein %in% Other.proteins) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("HSP104", "HSP42", "HSP26")))))

ggplot(Other.proteins.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.01, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)")) +
  ggtitle("Other heat shock proteins")

ggsave("../FINAL_Figures/Figure3/OtherHeatShockAbundance.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 20)

# Make plot to show abundance change of proteoasome subunits across all three mistranslating strains (Figure S2A)
Proteasome.proteins <- c("PRE1", "PRE2", "PRE3", "PRE4", "PRE5", "PRE6", "PRE7", "PRE8", "PRE9", "PRE10", "PUP1", "PUP2", "PUP3", "RPN1", "RPN2", "RPN3", "RPN5", "RPN6", "RPN7", "RPN8", "RPN9", "RPN10", "RPN11", "RPN12", "RPN13", "RPT1", "RPT2", "RPT3", "RPT4", "RPT5", "RPT6", "SCL1", "SEM1")

Proteasome.FC <- data.frame(Proteome.DA) %>%
  filter(Protein %in% Proteasome.proteins) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein))

ggplot(Proteasome.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.01, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)")) +
  ggtitle("Other heat shock proteins")

ggplot(Proteasome.FC, aes(x = Mistranslation, y = log2FC, fill = Mistranslation)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.5) +
  theme_matt() +
  ylab(expression(Log[2] * "(Fold-change)")) +
  scale_fill_manual(values= c("#58c87e", "#f8c754", "#a165a7")) +
  coord_cartesian(ylim = c(-1,1)) +
  geom_hline(yintercept = 0, linewidth = 1, color = "red", linetype = "dashed") +
  theme(legend.position = "none") +
  ggtitle("Proteasome Subunits")

ggsave("../FINAL_Figures/FigureS2/ProteasomeAbundance.pdf", device = cairo_pdf, units = "in", width = 3, height = 4)

# Make plot to show abundance change of autophagy proteins across all three mistranslating strains (Figure S2B)
Autophagy.proteins <- c("ATG1", "ATG13", "ATG17", "ATG29", "ATG31", "ATG9", "VPS15", "VPS30", "VPS34", "ATG14", "ATG38", "ATG2", "ATG18", "ATG21", "ATG5", "ATG10", "ATG12", "ATG16", "ATG7", "ATG3", "ATG4", "ATG8")

Autophagy.FC <- data.frame(Proteome.DA) %>%
  filter(Protein %in% Autophagy.proteins) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("ATG1", "ATG13", "ATG17", "ATG29", "ATG31", "ATG9", "VPS15", "VPS30", "VPS34", "ATG38", "ATG2", "ATG18", "ATG21", "ATG5", "ATG12", "ATG16", "ATG7", "ATG3", "ATG4", "ATG8")))))

ggplot(Autophagy.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.01, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)")) +
  ggtitle("Autophagy proteins")

ggsave("../FINAL_Figures/FigureS2/AutophagyProteinAbundance.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 20)


## Heatmap of differentially regulate proteins from all mistranslating strains (Figure 3B)
# Take only the columns needed and split into two data frames - one for each tRNA
ProSer.DA <- Proteome.DA %>%
  filter(Label == "ProSer_vs_Ctrl") %>%
  dplyr::select(Protein, log2FC, adj.pvalue) %>%
  mutate(FoldChange = ifelse(adj.pvalue > 0.05, 0, log2FC)) %>% 
  filter(is.finite(FoldChange))

ArgSer.DA <- Proteome.DA %>%
  filter(Label == "ArgSer_vs_Ctrl") %>%
  dplyr::select(Protein, log2FC, adj.pvalue) %>%
  mutate(FoldChange = ifelse(adj.pvalue > 0.05, 0, log2FC)) %>% 
  filter(is.finite(FoldChange))

ProAla.DA <-Proteome.DA %>%
  filter(Label == "ProAla_vs_Ctrl") %>%
  dplyr::select(Protein, log2FC, adj.pvalue) %>%
  mutate(FoldChange = ifelse(adj.pvalue > 0.05, 0, log2FC)) %>% 
  filter(is.finite(FoldChange))

# Merge the data frames
DiffAbund <- ProSer.DA %>%
  left_join(ArgSer.DA, by = c("Protein"), suffix = c(".ProSer", ".ArgSer")) %>%
  left_join(ProAla.DA, by = c("Protein"))%>%
  dplyr::rename(log2FC.ProAla = log2FC,
                adj.pvalue.ProAla = adj.pvalue,
                FoldChange.ProAla = FoldChange) %>%
  na.omit() %>%
  mutate(Significant = ifelse(adj.pvalue.ProSer < 0.05 | adj.pvalue.ArgSer < 0.05 | adj.pvalue.ProAla < 0.05, T, F))

# Heatmap of differentially abundant proteins
Heatmap <- DiffAbund %>%
  dplyr::select(Protein, FoldChange.ProSer, FoldChange.ArgSer, FoldChange.ProAla)

colnames(Heatmap) <- c("Gene", "ProSer", "ArgSer", "ProAla")
rownames(Heatmap) <- Heatmap$Gene
Heatmap <- Heatmap[-1]

CleanHeatmap <- t(Heatmap[which(rowSums(Heatmap) != 0), ])
breaks <- seq(-1, 1, by = 0.02)

Heatmap.Figure <- pheatmap(CleanHeatmap, 
                           show_colnames = F, 
                           col = colorRampPalette(c("blue3", "grey93", "red3"))(100),
                           annotation_legend = T, 
                           clustering_distance_cols = "euclidean", 
                           clustering_method = "ward.D2", 
                           cutree_cols = 4, 
                           annotation_col = annot,
                           breaks = breaks)

svg(filename = "../FINAL_Figures/Figure3/WholeProteomeHeatMap.svg", width = 10, height = 5)
Heatmap.Figure
dev.off()

heatmap.clust <- cbind(t(CleanHeatmap), cluster = cutree(Heatmap.Figure$tree_col, k = 4))
annot <- as.data.frame(cutree(Heatmap.Figure$tree_col, k = 4))
annot$`cutree(Heatmap.Figure$tree_col, k = 4)` <- factor(annot$`cutree(Heatmap.Figure$tree_col, k = 4)`)

write.csv(heatmap.clust, "ProteomeHeatmapClusters.csv")

write.csv(Proteome.DA, "AllProteinsInProteomeDataset.csv")


## PHOSPHO ANALYSIS ##

# Make volcano plots from MSstats analysis
ProAlaVolcano <- ggplot(Phospho.DA %>% filter(Label == "ProAla_vs_Ctrl"), aes (x = log2FC, y = -log10(adj.pvalue), text = Protein, label = Protein)) +
  geom_point(color = "lightgrey") +
  geom_point(data = Phospho.DA %>% filter(Label == "ProAla_vs_Ctrl") %>% filter(log2FC > 0.6 | log2FC < -0.6) %>% filter(adj.pvalue < 0.05), color = "#58c87e", size = 3) + 
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = -0.6, linetype = "dashed", color = "darkgrey") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  xlab("") +
  ylab("") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  theme_matt() +
  ggtitle("Pro\u2192Ala") + 
  theme(plot.title = element_text(size = 12, face = "bold"))

ProSerVolcano <- ggplot(Phospho.DA %>% filter(Label == "ProSer_vs_Ctrl"), aes (x = log2FC, y = -log10(adj.pvalue), text = Protein, label = Protein)) +
  geom_point(color = "lightgrey") +
  geom_point(data = Phospho.DA %>% filter(Label == "ProSer_vs_Ctrl") %>% filter(log2FC > 0.6 | log2FC < -0.6) %>% filter(adj.pvalue < 0.05), color = "#f8c754", size = 3) + 
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = -0.6, linetype = "dashed", color = "darkgrey") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  xlab("") +
  ylab("") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  theme_matt() +
  ggtitle("Pro\u2192Ser") + 
  theme(plot.title = element_text(size = 12, face = "bold"))

ArgSerVolcano <- ggplot(Phospho.DA %>% filter(Label == "ArgSer_vs_Ctrl"), aes (x = log2FC, y = -log10(adj.pvalue), text = Protein, label = Protein)) +
  geom_point(color = "grey") +
  geom_point(data = Phospho.DA %>% filter(Label == "ArgSer_vs_Ctrl") %>% filter(log2FC > 0.6 | log2FC < -0.6) %>% filter(adj.pvalue < 0.05), color = "#a165a7", size = 3) + 
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = -0.6, linetype = "dashed", color = "darkgrey") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  xlab("") +
  ylab("") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  theme_matt() +
  ggtitle("Arg\u2192Ser") + 
  theme(plot.title = element_text(size = 12, face = "bold")) 

plot_grid(ProAlaVolcano, ProSerVolcano, ArgSerVolcano, nrow = 1, scale = 0.999)
ggsave("../FigureX/PhosphoProteomeVolcano.pdf", device = cairo_pdf, width = 10, height = 3.5, units = "in")

# There are 405 upregulated phosphosites in ProSer
ProSer.PhosphoUp <- Phospho.DA %>%
  filter(Label == "ProSer_vs_Ctrl") %>%
  filter(adj.pvalue <= 0.05) %>%
  filter(log2FC > 0.6)

# There are 253 downregulated phosphosites in ProSer
ProSer.PhosphoDn <- Phospho.DA %>%
  filter(Label == "ProSer_vs_Ctrl") %>%
  filter(adj.pvalue <= 0.05) %>%
  filter(log2FC < -0.6)

# There are 38 upregulated phosphosites in ArgSer
ArgSer.PhosphoUp <- Phospho.DA %>%
  filter(Label == "ArgSer_vs_Ctrl") %>%
  filter(adj.pvalue <= 0.05) %>%
  filter(log2FC > 0.6)

# There are 15 downregulated phosphosites in ArgSer
ArgSer.PhosphoDn <- Phospho.DA %>%
  filter(Label == "ArgSer_vs_Ctrl") %>%
  filter(adj.pvalue <= 0.05) %>%
  filter(log2FC < -0.6)

# There are 51 upregulated phosphosites in ProAla
ProAla.PhosphoUp <- Phospho.DA %>%
  filter(Label == "ProAla_vs_Ctrl") %>%
  filter(adj.pvalue <= 0.05) %>%
  filter(log2FC > 0.6)

# There are 41 downregulated phosphosites in ProAla
ProAla.PhosphoDn <- Phospho.DA %>%
  filter(Label == "ProAla_vs_Ctrl") %>%
  filter(adj.pvalue <= 0.05) %>%
  filter(log2FC < -0.6)

# How many proteins have at least one differential phosphosites across all mistranslating strains?
AffectedPhosphoProteins <- Phospho.DA %>% filter(adj.pvalue <= 0.05) %>% filter(log2FC > 0.6 | log2FC < -0.6)
length(unique(AffectedPhosphoProteins$GlobalProtein)) # 456 proteins


## Heatmap of differentially regulate phosphosites from all mistranslating strains (Figure 4A)
# Take only the columns needed and split into two data frames - one for each tRNA
ProSer.PhosphoDA <- Phospho.DA %>%
  filter(Label == "ProSer_vs_Ctrl") %>%
  dplyr::select(Protein, log2FC, adj.pvalue) %>%
  mutate(FoldChange = ifelse(adj.pvalue > 0.05, 0, log2FC)) %>% 
  filter(is.finite(FoldChange))

ArgSer.PhosphoDA <- Phospho.DA %>%
  filter(Label == "ArgSer_vs_Ctrl") %>%
  dplyr::select(Protein, log2FC, adj.pvalue) %>%
  mutate(FoldChange = ifelse(adj.pvalue > 0.05, 0, log2FC)) %>% 
  filter(is.finite(FoldChange))

ProAla.PhosphoDA <-Phospho.DA %>%
  filter(Label == "ProAla_vs_Ctrl") %>%
  dplyr::select(Protein, log2FC, adj.pvalue) %>%
  mutate(FoldChange = ifelse(adj.pvalue > 0.05, 0, log2FC)) %>% 
  filter(is.finite(FoldChange))

# Merge the data frames
DiffAbundPhospho <- ProSer.PhosphoDA %>%
  left_join(ArgSer.PhosphoDA, by = c("Protein"), suffix = c(".ProSer", ".ArgSer")) %>%
  left_join(ProAla.PhosphoDA, by = c("Protein"))%>%
  dplyr::rename(log2FC.ProAla = log2FC,
                adj.pvalue.ProAla = adj.pvalue,
                FoldChange.ProAla = FoldChange) %>%
  na.omit() %>%
  mutate(Significant = ifelse(adj.pvalue.ProSer < 0.05 | adj.pvalue.ArgSer < 0.05 | adj.pvalue.ProAla < 0.05, T, F))

# Heatmap of differentially abundant phosphosites
HeatmapPhospho <- DiffAbundPhospho %>%
  dplyr::select(Protein, FoldChange.ProSer, FoldChange.ArgSer, FoldChange.ProAla)

colnames(HeatmapPhospho) <- c("Gene", "ProSer", "ArgSer", "ProAla")
rownames(HeatmapPhospho) <- HeatmapPhospho$Gene
HeatmapPhospho <- HeatmapPhospho[-1]

CleanHeatmapPhospho <- t(HeatmapPhospho[which(rowSums(HeatmapPhospho) != 0), ])
breaks <- seq(-2, 2, by = 0.04)

Heatmap.Figure.Phospho <- pheatmap(CleanHeatmapPhospho, 
                           show_colnames = F, 
                           col = colorRampPalette(c("blue3", "grey93", "red3"))(100),
                           annotation_legend = T, 
                           clustering_distance_cols = "euclidean", 
                           clustering_method = "ward.D2", 
                           cutree_cols = 3, 
                           #annotation_col = annot.phospho,
                           breaks = breaks)

svg(filename = "../FINAL_Figures/Figure4/PhosphoProteomeHeatMap.svg", width = 10, height = 5)
Heatmap.Figure.Phospho
dev.off()

heatmap.clust.phospho <- cbind(t(CleanHeatmapPhospho), cluster = cutree(Heatmap.Figure.Phospho$tree_col, k = 3))
annot.phospho <- as.data.frame(cutree(Heatmap.Figure.Phospho$tree_col, k = 3))
annot.phospho$`cutree(Heatmap.Figure.Phospho$tree_col, k = 3)` <- factor(annot.phospho$`cutree(Heatmap.Figure.Phospho$tree_col, k = 3)`)

write.csv(heatmap.clust.phospho, "../FINAL_data/PhosphoproteomeHeatmapClusters.csv")

write.csv(Phospho.DA, "../FINAL_data/AllPhosphoProteinsInProteomeDataset.csv")

# Kinase substrate enrichment analysis and plot (Figure 4C)
# Read in Kinase-PhosphoSite file
KS.data <- Phospho.DA %>%
  separate(Protein, into = c("Protein", "ResidueSite")) %>%
  dplyr::select(test_label = Label,
                gene = Protein,
                log2_fc = log2FC,
                adj_pval = adj.pvalue,
                p_position = ResidueSite) %>%
  unite("ref", c(gene, p_position), remove = F) %>%
  mutate(significant = ifelse(adj_pval < 0.05, "Yes", "No"),
         regulated = ifelse(significant == "No", "ns",
                            ifelse(log2_fc > 0, "up", "dn"))) %>%
  dplyr::select(-adj_pval, -significant)

# Fisher Test Functions
## Function to perform fisher test and get p.val
exact_pval_fun <- function(data){
  dataframe <- data.frame(data) %>% 
    column_to_rownames(var = "selection")
  matrix <- as.matrix(dataframe) 
  exact_test_result <- fisher.test(matrix, alternative = "two.sided")
  return(exact_test_result$p.value)
}

## Function to perform fisher test and get enrichment
enrich_fun<-function(data){
  dataframe<-data.frame(data) %>% 
    mutate(frac = pos/neg) %>% 
    dplyr::select(-c(pos,neg)) %>% 
    pivot_wider(names_from = selection, values_from = frac) %>% 
    summarise(enrichment = foreground/background) 
  return(dataframe$enrichment[1])
}

uniprot <- read_csv("input_data/Yeast_uniprot.csv") %>% 
  dplyr::select(gene = Gene, reference = Reference) %>% 
  separate(gene, "gene")

ks <- read_csv("input_data/kinase_phos_substrate_yeast_2.csv") %>% 
  separate(ref, c("reference","B","C")) %>% 
  unite("A", c(B,C), sep = "") %>% 
  left_join(uniprot) %>% 
  unite("ref", c(gene, A)) %>% 
  dplyr::select(reference = ref, kinase.name = kinase_activity)

clust.annot <- KS.data %>% 
  filter(regulated != "ns") %>% 
  unite("cluster", c(test_label, regulated)) %>% 
  dplyr::select(reference = ref, cluster) %>% 
  distinct() %>% 
  left_join(ks)

# Plot of sites for each kinase (Supplemental Figure and Figure 4B)
Site.Kinase.FC <- KS.data %>%
  dplyr::rename(reference = ref) %>%
  left_join(ks, by = c("reference")) %>%
  na.omit() %>%
  mutate(tRNA = ifelse(test_label == "ProAla_vs_Ctrl", "Pro\u2192Ala",
                       ifelse(test_label == "ProSer_vs_Ctrl", "Pro\u2192Ser",
                              ifelse(test_label == "ArgSer_vs_Ctrl", "Arg\u2192Ser", NA)))) %>%
  mutate(tRNA = fct_relevel(tRNA, "Pro\u2192Ala", "Pro\u2192Ser", "Arg\u2192Ser"))

ggplot(Site.Kinase.FC, aes(x = kinase.name, y = log2_fc, color = tRNA)) +
  geom_boxplot() +
  theme_matt() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  facet_wrap(~tRNA, ncol = 1) +
  ylab("log2(Fold Change)") +
  xlab("") +
  scale_color_manual(values= c("Pro\u2192Ala" = "#58c87e", "Pro\u2192Ser" = "#f8c754", "Arg\u2192Ser" = "#a165a7"))

ggsave("../FINAL_Figures/Figure4/AllKinaseSites.pdf", device = cairo_pdf, width = 17, height = 12)

# Plot boxplots for ATG1, CDC28 and Tor1 specifically
SpecificKinases <- Site.Kinase.FC %>%
  filter(kinase.name == "ATG1" |
           kinase.name == "CDC28" |
           kinase.name == "TORC1") %>%
  mutate(tRNA = ifelse(test_label == "ProAla_vs_Ctrl", "Pro\u2192Ala",
                       ifelse(test_label == "ProSer_vs_Ctrl", "Pro\u2192Ser",
                              ifelse(test_label == "ArgSer_vs_Ctrl", "Arg\u2192Ser", NA)))) %>%
  mutate(tRNA = fct_relevel(tRNA, "Pro\u2192Ala", "Pro\u2192Ser", "Arg\u2192Ser"))


ggplot(SpecificKinases, aes(x = kinase.name, y = log2_fc, color = tRNA)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_violin(position = position_dodge(width = 0.8), trim = FALSE) +
  geom_boxplot(width = 0.15, position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(data = SpecificKinases %>% filter(regulated != "ns"), aes(color = tRNA), position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), alpha = 0.75) +
  theme_matt() +
  ylab("Log2(Fold Change)") +
  xlab("") +
  scale_color_manual(values= c("Pro\u2192Ala" = "#58c87e", "Pro\u2192Ser" = "#f8c754", "Arg\u2192Ser" = "#a165a7"))

ggsave("../FINAL_Figures/Figure4/SpecificKinasesBoxplot.pdf", device = cairo_pdf, width = 6, height = 5)

kinase.df <- ks %>%
  ungroup %>% 
  drop_na()

kinase.bg <- kinase.df %>% 
  mutate(all = n_distinct(reference)) %>% 
  drop_na() %>% 
  group_by(kinase.name) %>% 
  summarise("pos" = n_distinct(reference),
            "neg" = all - pos) %>% 
  distinct() %>% 
  mutate(selection = "background") %>% 
  ungroup() %>% 
  filter(pos > 4) #remove go with less then 15 hits in background

kinase.df <- kinase.bg %>% 
  dplyr::select(kinase.name) %>% 
  left_join(kinase.df) %>% 
  drop_na()

kinase.fg <- clust.annot %>% 
  dplyr::select(reference, cluster) %>% 
  mutate(selection = "foreground") %>% 
  dplyr::select(reference, selection, cluster) %>% 
  right_join(kinase.df, relationship = "many-to-many") %>% 
  group_by(selection, cluster) %>% 
  mutate(term_all = n_distinct(reference)) %>% 
  group_by(kinase.name, selection, cluster) %>% 
  summarise("pos" = n_distinct(reference),
            "neg" = term_all - pos)  %>% 
  ungroup() %>% 
  drop_na(cluster) %>% 
  distinct()


#### multiple comparison #####
kinase.nest <- kinase.bg %>% 
  left_join(kinase.fg %>% distinct(kinase.name, cluster)) %>% 
  distinct() %>% 
  full_join(kinase.fg) %>% 
  drop_na(cluster) %>% 
  drop_na(kinase.name) %>% 
  distinct() %>% 
  group_by(kinase.name, cluster) %>% 
  nest() 
             
fish.kinase.df <- kinase.nest %>% 
  mutate(p.val = data %>% map_dbl(exact_pval_fun))%>% #calc p val
  mutate(enrich = data %>% map_dbl(enrich_fun)) %>%   #calc enrichment
  dplyr::select(-c(data)) %>% 
  mutate(p.val.adj = p.adjust(p.val, method = "fdr")) %>% 
  ungroup() %>% 
  left_join(kinase.fg %>% dplyr::select(cluster, kinase.name, fg = pos))   #foregournd

ks.enrich <- fish.kinase.df %>% 
  filter(p.val.adj < 0.05) %>% 
  filter(fg > 4,
         enrich > 1)

ks.enrich2 <- ks.enrich %>% 
  separate(cluster, c("A", "B", "C", "D"), sep = "_") %>% 
  mutate(enrich = log2(enrich)) %>% 
  dplyr::rename(strain = A,
         reg = D) %>%
  distinct()

ggplot(ks.enrich2) +
  geom_point(aes(y= strain, x = kinase.name, size = -log10(p.val.adj), fill = enrich), shape = 21, color="black", stroke = 0.5) +
  scale_fill_viridis(limits = c(1, 4), breaks = c(0, 1, 2, 3, 4), option = "plasma")+
  scale_size(name = waiver(),
             limits = c(1,15), 
             range = c(1, 10),
             trans = "identity",
             guide = "legend") +
  ggforce::facet_col(facets = "reg", scales = "free_y", space = "free") +
  theme_bw() +
  theme(legend.position="top",  plot.margin = margin(50, 50, 50, 50)) +
  theme(axis.text.x = element_text(size = 10,angle = 90), axis.text.y = element_text(size = 10), legend.text =element_text(size = 10) )

ggsave("../FINAL_Figures/Figure4/DownRegulatedKinase_plot.pdf", device = cairo_pdf, width = 4, height = 4, units = "in")

# Heat map for CDC5 protein level and phospho activation site, Bfa1 phospho, Nud1 phospho, Cdc14 phospho, Net1 phospho (Figure 4D)
CDC5Protein.FC <- data.frame(Proteome.DA) %>%
  filter(Protein == "CDC5") %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("CDC5")))))

ggplot(CDC5Protein.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.05, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)"))

ggsave("../FINAL_Figures/Figure4/Cdc5_Protein.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 10)

CDC5Phospho.FC <- data.frame(Phospho.DA) %>%
  filter(Protein == "CDC5_T70") %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("CDC5_T70")))))

ggplot(CDC5Phospho.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.01, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)"))

ggsave("../FINAL_Figures/Figure4/Cdc5_T70.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 10)

CDC14Phospho.FC <- data.frame(Phospho.DA) %>%
  filter(Protein == "CDC14_S408" |
           Protein == "CDC14_S429") %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("CDC14_S408", "CDC14_S429")))))

ggplot(CDC14Phospho.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.05, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)"))

ggsave("../FINAL_Figures/Figure4/Cdc14_Phosphorylation.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 10)

Net1Sites <- c("NET1_T676", "NET1_T357", "NET1_T288", "NET1_T1017", "NET1_S830", "NET1_S803", "NET1_S744", "NET1_S679", "NET1_S385", "NET1_S166", "NET1_S1032")

NET1Phospho.FC <- data.frame(Phospho.DA) %>%
  filter(Protein %in% Net1Sites) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("NET1_S166", "NET1_T288", "NET1_T357", "NET1_S385", "NET1_T676", "NET1_S679", "NET1_S744", "NET1_S803", "NET1_S830", "NET1_T1017", "NET1_S1032")))))

ggplot(NET1Phospho.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.05, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)"))

ggsave("../FINAL_Figures/Figure4/NET1_Phosphorylation.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 10)

Bfa1Sites <- c("BFA1_S274", "BFA1_T340")

BFA1Phospho.FC <- data.frame(Phospho.DA) %>%
  filter(Protein %in% Bfa1Sites) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("BFA1_S274", "BFA1_T340")))))

ggplot(BFA1Phospho.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.05, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)"))

ggsave("../FINAL_Figures/Figure4/Bfa1_Phosphorylation.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 10)

Nud1Sites <- c("NUD1_S63", "NUD1_S240", "NUD1_S460")

NUD1Phospho.FC <- data.frame(Phospho.DA) %>%
  filter(Protein %in% Nud1Sites) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein) ,
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("NUD1_S63", "NUD1_S240", "NUD1_S460")))))

ggplot(NUD1Phospho.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.05, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)"))

ggsave("../FINAL_Figures/Figure4/Nud1_Phosphorylation.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 10)

# Plot of Clb2 levels (Figure 4E)
Clb2 <- ProteinLevel %>% 
  filter(Protein == "CLB2") %>%
  mutate(MistranslationShort = case_when(GROUP == "ProAla" ~ "P\u2192A",
                                         GROUP == "ProSer" ~ "P\u2192S",
                                         GROUP == "ArgSer" ~ "R\u2192S",
                                         GROUP == "Ctrl" ~ "EV",))

condition_colors <- c("EV" = "#6496bf", "P\u2192A" = "#58c87e", "P\u2192S" = "#f8c754", "R\u2192S" = "#a165a7")
condition_colors_alpha <- alpha(condition_colors, 0.25) 

ggplot(Clb2, aes(x = MistranslationShort, y = LogIntensities, color = MistranslationShort, fill = MistranslationShort)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  theme_matt() +
  theme(legend.position = "none") +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors_alpha) +
  scale_y_continuous(limits = c(17, 20)) +
  labs(y = expression(log[2]("Protein Abundance")),
       x = "") +
  ggtitle("Clb2p Levels")

ggsave("../FINAL_Figures/Figure4/Clb2LevelsBoxplot.pdf", device = cairo_pdf, width = 3.5, height = 4)

# Plot of Sic1 levels (Figure 4E)
Sic1 <- ProteinLevel %>% 
  filter(Protein == "SIC1") %>%
  mutate(MistranslationShort = case_when(GROUP == "ProAla" ~ "P\u2192A",
                                         GROUP == "ProSer" ~ "P\u2192S",
                                         GROUP == "ArgSer" ~ "R\u2192S",
                                         GROUP == "Ctrl" ~ "EV",))

condition_colors <- c("EV" = "#6496bf", "P\u2192A" = "#58c87e", "P\u2192S" = "#f8c754", "R\u2192S" = "#a165a7")
condition_colors_alpha <- alpha(condition_colors, 0.25) 

ggplot(Sic1, aes(x = MistranslationShort, y = LogIntensities, color = MistranslationShort, fill = MistranslationShort)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  theme_matt() +
  theme(legend.position = "none") +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors_alpha) +
  scale_y_continuous(limits = c(18.5, 20.5)) +
  labs(y = expression(log[2]("Protein Abundance")),
       x = "") +
  ggtitle("Sic1p Levels")

ggsave("../FINAL_Figures/Figure4/Sic1LevelsBoxplot.pdf", device = cairo_pdf, width = 3.5, height = 4)
  

# Lollipop plot of Eft1 (Figure 5A)
# Plot phospho-site differential abundance for different proteins as lollipop plots
# Read in yeast FASTA file
YeastFASTA <- readAAStringSet("2023-02-23-decoys-reviewed-contam-UP000002311.fas")

# Extract the accession numbers and rename AAStringSet
accession_numbers <- sapply(YeastFASTA@ranges@NAMES, function(x) str_extract(as.character(x), "(?<=\\|)[^|]+(?=\\|)"))
names(YeastFASTA) <- accession_numbers

# Merge with earlier data set to annotate each protein with its full length
YeastProteinLength <- data.frame(YeastFASTA@ranges@width, YeastFASTA@ranges@NAMES) %>%
  dplyr::select(ProteinLength = YeastFASTA.ranges.width, 
                Protein.Group = YeastFASTA.ranges.NAMES)

sgd <- read_tsv("20250727_UniprotSGD.txt") %>%
  dplyr::select(Protein.Group = Entry,
                Gene = `Gene Names`,
                Length) %>%
  separate_rows(Gene, sep = " ")

ProteinLength <- PhosphoLevel %>%
  dplyr::select(Protein)%>%
  unique() %>%
  separate(Protein, into = c("Gene", "Site"), sep = "_") %>%
  filter(!str_detect(Gene, "cRAP")) %>%
  dplyr::select(Gene) %>%
  unique() %>%
  separate_rows(Gene, sep = ";") %>%
  full_join(sgd) %>%
  na.omit()

Phospho.DA.SiteLength <- Phospho.DA %>%
  separate(Protein, into = c("Gene", "Site"), remove = F) %>%
  left_join(ProteinLength, by = c("Gene")) %>%
  na.omit() %>%
  mutate(Site = gsub("[^0-9]", "", Site)) %>%
  mutate(Site = as.numeric(Site))

Eft1PhosphoSites <- Phospho.DA.SiteLength %>% filter(Gene == "EFT1")

condition_colors <- c("ProAla_vs_Ctrl" = "#58c87e", "ProSer_vs_Ctrl" = "#f8c754", "ArgSer_vs_Ctrl" = "#a165a7")
condition_labels <- c("ProAla_vs_Ctrl" = "Pro\u2192Ala", "ProSer_vs_Ctrl" = "Pro\u2192Ser", "ArgSer_vs_Ctrl" = "Arg\u2192Ser")

ggplot(Eft1PhosphoSites, aes(x = Site, y = log2FC, color = Label)) +
  geom_segment(data = Eft1PhosphoSites %>% filter(adj.pvalue < 0.05), aes(xend = Site, yend = 0), color = "black") +
  geom_segment(x = 0, y = 0, xend = max(Eft1PhosphoSites$Length), yend = 0, linewidth = 2, color = "black") +
  geom_point(data = Eft1PhosphoSites %>% filter(adj.pvalue < 0.05), size = 4) +
  theme_matt() +
  scale_x_continuous(limits = c(0, Eft1PhosphoSites$Length)) +
  ylab("Log2(Fold-Change)") +
  scale_color_manual(values= condition_colors, labels = condition_labels) +
  geom_point(data = Eft1PhosphoSites, x = Eft1PhosphoSites$Site, y = 0, color = "red", size = 3) +
  ggtitle("Eft1 PhosphoSites")

ggsave("../FINAL_Figures/Figure5/Eft1PhosphoSites.pdf", device = cairo_pdf, width = 5, height = 3)

# Phosphosites on proteins involved in translation other than EFT1 (Figure 5C)
# Ded1 regulated phosphosites
Ded1Sites <- c("DED1_T282", "DED1_S576", "DED1_S572", "DED1_S55", "DED1_S543", "DED1_S369", "DED1_S263")

DED1Phospho.FC <- data.frame(Phospho.DA) %>%
  filter(Protein %in% Ded1Sites) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein) ,
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("DED1_S55", "DED1_S263" , "DED1_T282", "DED1_S369", "DED1_S543", "DED1_S572", "DED1_S576")))))

ggplot(DED1Phospho.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.05, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)"))

ggsave("../FINAL_Figures/Figure5/Ded1_Phosphorylation.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 10)

# Sup45 regulated phosphosites
Sup45Sites <- c("SUP45_S67")

SUP45Phospho.FC <- data.frame(Phospho.DA) %>%
  filter(Protein %in% Sup45Sites) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein), ,
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("SUP45_S67")))))

ggplot(SUP45Phospho.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.05, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)"))

ggsave("../FINAL_Figures/Figure5/Sup45_Phosphorylation.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 10)

# Tef1 regulated phosphosites
Tef1Sites <- c("TEF1_S237", "TEF1_T240")

Tef1Phospho.FC <- data.frame(Phospho.DA) %>%
  filter(Protein %in% Tef1Sites) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("TEF1_S237", "TEF1_T240")))))

ggplot(Tef1Phospho.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.05, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)"))

ggsave("../FINAL_Figures/Figure5/Tef1_Phosphorylation.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 10)

# Sui2 regulated phosphosites
Sui2Sites <- c("SUI2_S52", "SUI2_S91")

Sui2Phospho.FC <- data.frame(Phospho.DA) %>%
  filter(Protein %in% Sui2Sites) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("SUI2_S52", "SUI2_S91")))))

ggplot(Sui2Phospho.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.05, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)"))

ggsave("../FINAL_Figures/Figure5/Sui2_Phosphorylation.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 10)

# Tif4631 regulated phosphosites
Tif4631Sites <- c("TIF4631_T276", "TIF4631_S916", "TIF4631_S920", "TIF4631_S502")

Tif4631Phospho.FC <- data.frame(Phospho.DA) %>%
  filter(Protein %in% Tif4631Sites) %>%
  mutate(Mistranslation = case_when(Label == "ProAla_vs_Ctrl" ~ "P\u2192A",
                                    Label == "ProSer_vs_Ctrl" ~ "P\u2192S",
                                    Label == "ArgSer_vs_Ctrl" ~ "R\u2192S"),
         Proteins = str_to_title(Protein),
         Proteins = fct_relevel(Proteins, str_to_title(rev(c("TIF4631_T276", "TIF4631_S502", "TIF4631_S916", "TIF4631_S920")))))

ggplot(Tif4631Phospho.FC, aes(y = Proteins, x = Mistranslation, fill = log2FC)) +
  geom_tile(colour = "black", size = 0.75)+
  scale_fill_gradient2(low = "blue3", mid = "gray93", high = "red3", limits = c(-2, 2), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))  +
  geom_point(aes(shape=ifelse(adj.pvalue < 0.05, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none")  +
  theme_matt()  +
  ylab("") +
  xlab("") +
  coord_fixed(ratio = 0.6) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.ontop = T) +
  labs(fill = expression(Log[2] * "(Fold-change)"))

ggsave("../FINAL_Figures/Figure5/Tif4631_Phosphorylation.pdf", device = cairo_pdf, units = "in", width = 4.5, height = 10)