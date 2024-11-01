# Mertzetal2024
## This respository contains R code associated with "Gut microbiome provides amino acids to host to compensate for dietary deficiencies in deer mice (Peromyscus)"

**16s_processing.R:** Initial processing of bacterial 16s amplicon sequences using DADA2 pipeline. 16S sequence data are deposited in the National Center for Biotechnology Information with the primary accession code PRJNA1108251. Inputs are fastq.gz files from each sequencing run. Output is phyloseq object used in downstream analyses.

**16S_DecontamRarefy.Rmd** This code details how to decontaminate and rarefy 16S microbiome and food data. Input is un-rarefied phyloseq object, output is cleaned and rarefied microbiome phyloseq object.

**16S_DiversityComposition.Rmd** Script for alpha and beta diversity analyses for microbiome and diet data. Includes code for 16S alpha and beta diversity figures, and relative abundance figures.

**qPCR_BacterialBiomass.Rmd** Script used to perform statistical and visual analyses of qPCR bacterial biomass data. Code for supplemental qPCR figure is included.
