# Mertzetal2024
## This respository contains R code associated with "Gut microbiome provides amino acids to host to compensate for dietary deficiencies in deer mice (Peromyscus)"

**16s_processing.R:** Initial processing of bacterial 16s amplicon sequences using DADA2 pipeline. 16S sequence data are deposited in the National Center for Biotechnology Information with the primary accession code PRJNA1108251. Inputs are fastq.gz files from each sequencing run. Output is phyloseq object used in downstream analyses.

**16S_DecontamRarefy.Rmd** This code details how to decontaminate and rarefy 16S microbiome and food data. Input is un-rarefied phyloseq object, output is cleaned and rarefied microbiome phyloseq object.

**16S_DiversityComposition.Rmd** Script for alpha and beta diversity analyses for microbiome and diet data. Includes code for 16S alpha and beta diversity figures, and relative abundance figures.

**qPCR_BacterialBiomass.Rmd** Script used to perform statistical and visual analyses of qPCR bacterial biomass data. Code for supplemental qPCR figure is included.

**MixSIAR_d13C_Models.Rmd** Script for running compound-specific stable isotope mixing models of amino acid δ13C values from mice muscle (mixture/consumer) and dietary or microbial protein sources (casein for synthetic diet; cornmeal and cricket powder for semi-natural diet) to estimate microbial contributions of amino acids to host tissue.

**MetagenomicsProcessingPipeline** Code detailing pipeline used for metagenomic data processing and analysis of cecum samples sequenced at 7 million reads. Inputs are raw fastq files, output include assembled MAGs and their resulting proteomes.

**MAGS_SankeyDataVisualization** This script details how to make a phyloseq object from metagenomic processing pipeline output and includes code on how to visualize data using Sankey diagrams.
