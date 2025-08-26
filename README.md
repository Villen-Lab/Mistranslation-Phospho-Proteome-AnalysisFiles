# Analysis of yeast proteome and phosphoproteome in response to mistranslation

This repository contains outputs and data analysis scripts to accompany analysis of the proteome and phosphoproteome response to mistranslating tRNAs. The raw mass spectrometry proteomics data have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifiers XXX. The folders in this repository contain the following:

**LogFiles:** Log files produced when analyzing the raw files with either FragPipe or DIA-NN containing the parameters used.

**DIALibraries:** Proteome and phosphoproteome libraries used to analyze DIA raw files created with FragPipe. The proteome library was made from all proteome DIA and DDA runs along with DIA-GPF runs from a pooled sample. The phosphoproteome library was created with all the phosphoproteome DIA and DDA files supplemented with deep yeast phosphoproteome raw files from Leutert et al. (The regulatory landscape of the yeast phosphoproteome. 2023. Nat Struct Mol Biol 30:1761-1773). 

**SearchOutputs:** Output files of either FragPipe (for DDA analysis of mistranslation frequency and new phosphosites) or DIA-NN (for DIA analysis of proteome and phosphoproteome samples) that are further analyzed in the R scripts containing in the AnalysisCode folder.

**AnalysisCode:** R scripts used to analyze the DIA-NN and FragPipe outputs and create the figures presented in the manuscript. Also contains additional input files (e.g. yeast kinase-substrate data) required for some of the analyses.


