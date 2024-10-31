#Initial processing of 16S rRNA amplicon sequences
#starts with fastq files, and then produces single phyloseq object for the project
#ran on R studio 

#code based on:
#Dada2: https://benjjneb.github.io/dada2/tutorial.html
#combining sequencing runs: https://benjjneb.github.io/dada2/bigdata.html
#Phyloseq:https://joey711.github.io/phyloseq/


#Step 1- create seqtab object from each run (danny2021,fred2022,Bonnie2022)
#Step2- merge all of the seqtabobjects 
#Step3- proceed to chimera step

#load packages
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(seqinr)

#####Make seqtab for Danny_run2021########
path <- "/Users/connermertz/Documents/Sequence Analysis/E1_E2_All_Samples/GUTS_Danny_run2021/GUTS_DannyFastQ"

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
forward <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
reverse <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(forward), "_"), `[`, 1)

#plotQualityProfile(forward[1:6])
#plotQualityProfile(reverse[1:6])

# Place filtered files in filtered/ subdirectory
filteredForward <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filteredReverse <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filteredForward) <- sample.names
names(filteredReverse) <- sample.names

#messed w/parameters, best are (260,200) and keep at (2,2), and removed primers using trimLeft function.  
out <- filterAndTrim(forward, filteredForward, reverse, filteredReverse, truncLen=c(260,200), 
                     trimLeft= c(19, 20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
out

errF <- learnErrors(filteredForward, multithread=TRUE)
errR <- learnErrors(filteredReverse, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaForward <- dada(filteredForward, err=errF, multithread=TRUE)
dadaReverse <- dada(filteredReverse, err=errR, multithread=TRUE)
dadaForward[[1]]

mergers <- mergePairs(dadaForward, filteredForward, dadaReverse, filteredReverse, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

saveRDS(seqtab, "/Users/connermertz/Documents/Sequence Analysis/E1_E2_All_Samples/GUTS_Danny_run2021/seqtab.rds")





#####Make seqtab for Fred_run2022#######
path <- "/Users/connermertz/Documents/Sequence Analysis/E1_E2_All_Samples/GUTS_Fred_run2022/GUTS_FredFastQ"

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
forward <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
reverse <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(forward), "_"), `[`, 1)

#plotQualityProfile(forward[1:6])
#plotQualityProfile(reverse[1:6])

# Place filtered files in filtered/ subdirectory
filteredForward <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filteredReverse <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filteredForward) <- sample.names
names(filteredReverse) <- sample.names

#messed w/parameters, best are (260,200) and keep at (2,2), and removed primers using trimLeft function.  
out <- filterAndTrim(forward, filteredForward, reverse, filteredReverse, truncLen=c(260,200), 
                     trimLeft= c(19, 20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
out

errF <- learnErrors(filteredForward, multithread=TRUE)
errR <- learnErrors(filteredReverse, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaForward <- dada(filteredForward, err=errF, multithread=TRUE)
dadaReverse <- dada(filteredReverse, err=errR, multithread=TRUE)
dadaForward[[1]]

mergers <- mergePairs(dadaForward, filteredForward, dadaReverse, filteredReverse, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

saveRDS(seqtab, "/Users/connermertz/Documents/Sequence Analysis/E1_E2_All_Samples/GUTS_Fred_run2022/seqtab.rds")




####make seqtab for Bonnie_run2022###
path <- "/Users/connermertz/Documents/Sequence Analysis/E1_E2_All_Samples/GUTS_Bonnie_run2022/GUTS_BonnieFastQ"

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
forward <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
reverse <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(forward), "_"), `[`, 1)

#plotQualityProfile(forward[1:6])
#plotQualityProfile(reverse[1:6])

# Place filtered files in filtered/ subdirectory
filteredForward <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filteredReverse <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filteredForward) <- sample.names
names(filteredReverse) <- sample.names

#messed w/parameters, best are (260,200) and keep at (2,2), and removed primers using trimLeft function.  
out <- filterAndTrim(forward, filteredForward, reverse, filteredReverse, truncLen=c(260,200), 
                     trimLeft= c(19, 20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
out

errF <- learnErrors(filteredForward, multithread=TRUE)
errR <- learnErrors(filteredReverse, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaForward <- dada(filteredForward, err=errF, multithread=TRUE)
dadaReverse <- dada(filteredReverse, err=errR, multithread=TRUE)
dadaForward[[1]]

mergers <- mergePairs(dadaForward, filteredForward, dadaReverse, filteredReverse, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

saveRDS(seqtab, "/Users/connermertz/Documents/Sequence Analysis/E1_E2_All_Samples/GUTS_Bonnie_run2022/seqtab.rds")


####make seqtab for Cindy_run2023###
#My laptop ran out of room so I moved all of my fastq's to Vlab dropbox on 6/27/23. I am going to try to upload fastq's from dropbox instead. 
path <- "/Users/connermertz/Dropbox (vlab)/Conner Mertz/Fastqs/E1_E2_16S/GUTS_Cindy_run2023/GUTS_CindyFastQ"


list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
forward <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
reverse <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(forward), "_"), `[`, 1)

#plotQualityProfile(forward[1:6])
#plotQualityProfile(reverse[1:6])

# Place filtered files in filtered/ subdirectory
filteredForward <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filteredReverse <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filteredForward) <- sample.names
names(filteredReverse) <- sample.names

#messed w/parameters, best are (260,200) and keep at (2,2), and removed primers using trimLeft function.  
out <- filterAndTrim(forward, filteredForward, reverse, filteredReverse, truncLen=c(260,200), 
                     trimLeft= c(19, 20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
out

errF <- learnErrors(filteredForward, multithread=TRUE)
errR <- learnErrors(filteredReverse, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaForward <- dada(filteredForward, err=errF, multithread=TRUE)
dadaReverse <- dada(filteredReverse, err=errR, multithread=TRUE)
dadaForward[[1]]

mergers <- mergePairs(dadaForward, filteredForward, dadaReverse, filteredReverse, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

saveRDS(seqtab, "/Users/connermertz/Dropbox (vlab)/Conner Mertz/Fastqs/E1_E2_16S/GUTS_Cindy_run2023/seqtab.rds")






#6.27.23 stopped here :) continue tomorrow
#Then combine the seqtab objects like this

st1 <- readRDS("/Users/connermertz/Dropbox (vlab)/Conner Mertz/Fastqs/E1_E2_16S/GUTS_Danny_run2021/seqtab.rds")
st2 <- readRDS("/Users/connermertz/Dropbox (vlab)/Conner Mertz/Fastqs/E1_E2_16S/GUTS_Fred_run2022/seqtab.rds")
st3 <- readRDS("/Users/connermertz/Dropbox (vlab)/Conner Mertz/Fastqs/E1_E2_16S/GUTS_Bonnie_run2022/seqtab.rds")
st4 <- readRDS("/Users/connermertz/Dropbox (vlab)/Conner Mertz/Fastqs/E1_E2_16S/GUTS_Cindy_run2023/seqtab.rds")


##repeat for each dataset sequenced separately
st.all <- mergeSequenceTables(st1, st2, st3,st4)



#remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(st.all)

#I don't need these stats b/c I already generated them for each individual run. When I try I get an error b/c it id calling the "out"object from DADA 2 workflow. Since I have many "out" objects from mult runs, I am going to skip this and just assign taxonomy. That is what Tinas code does. 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaForward, getN), sapply(dadaReverse, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
#head only shows top 6 samples
#track

#assign taxom=nomy using Silva 138.1 prokaryotic SSU taxonomic training data formateed for DADA 2.
#These training fasta files are derived from the Silva Project's version 138.1 release and formatted for use with DADA2. These files are #intended for use in classifying prokaryotic 16S sequencing data and are not appropriate for classifying eukaryotic ASVs. 
fn <- "/Users/connermertz/Dropbox (vlab)/Conner Mertz/SILVAdb/silva_nr99_v138.1_train_set.fa.gz"
file.exists(fn)
taxa_Dada <- assignTaxonomy(seqtab.nochim, "/Users/connermertz/Dropbox (vlab)/Conner Mertz/SILVAdb/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa_Dada # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)




#How to save DADA2 tables as CSVs
#create excel sheet with naming conventions-> saved in Sequencing Analysis

Meta <- read.csv("/Users/connermertz/Dropbox (vlab)/Conner Mertz/UNM biology/Chapter 1 E1 vs E2/Methods/16S rRNA gene sequencing/Sequencing_workup/E1_E2_Meta_ALL.csv")
#Meta <- Meta[c(1:6), c(1:6)]
row.names(Meta) <- Meta$Sample_ID

#here are the csv's used to make phyloseq objects

ps <- phyloseq(otu_table(seqtab.nochim,taxa_are_rows=FALSE), 
               sample_data(Meta), 
               tax_table(taxa_Dada))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#Save ps object as rds
saveRDS(ps, "/Users/connermertz/Dropbox (vlab)/Conner Mertz/UNM biology/Chapter 1 E1 vs E2/Methods/16S rRNA gene sequencing/Sequencing_workup/DADA2_output/raw_ps.rds")




#this was done so that I can now have data frames of the two things I make in dada2 and DON'T HAVE TO RUN IT AGAIN.
otu_export <- as.data.frame(otu_table(ps))
tax_export <- as.data.frame(tax_table(ps))
richness_export <- estimate_richness(ps, split=TRUE, measures=NULL)



#write the csvs for counts, taxa, and richness, ## these out after you run it once
counts <- write.csv(otu_export, "/Users/connermertz/Dropbox (vlab)/Conner Mertz/UNM biology/Chapter 1 E1 vs E2/Methods/16S rRNA gene sequencing/Sequencing_workup/DADA2_output/ALL_counts.csv")
taxa <- write.csv(tax_export, "/Users/connermertz/Dropbox (vlab)/Conner Mertz/UNM biology/Chapter 1 E1 vs E2/Methods/16S rRNA gene sequencing/Sequencing_workup/DADA2_output/ALL_taxa.csv")
richness <- write.csv(richness_export, "/Users/connermertz/Dropbox (vlab)/Conner Mertz/UNM biology/Chapter 1 E1 vs E2/Methods/16S rRNA gene sequencing/Sequencing_workup/DADA2_output/ALL_richness.csv")


ASV <- read.csv("/Users/connermertz/Dropbox (vlab)/Conner Mertz/UNM biology/Chapter 1 E1 vs E2/Methods/16S rRNA gene sequencing/Sequencing_workup/DADA2_output/ALL_counts.csv", row.names = "Group")
ASV <- t(ASV)
ASV <- data.matrix(ASV)

taxa <- read.csv("/Users/connermertz/Dropbox (vlab)/Conner Mertz/UNM biology/Chapter 1 E1 vs E2/Methods/16S rRNA gene sequencing/Sequencing_workup/DADA2_output/ALL_taxa.csv", row.names = "ASV")
taxa<- as.matrix(taxa)

diversity <- read.csv("/Users/connermertz/Dropbox (vlab)/Conner Mertz/UNM biology/Chapter 1 E1 vs E2/Methods/16S rRNA gene sequencing/Sequencing_workup/DADA2_output/ALL_richness.csv", row.names = "Sample_ID")
diversity<- as.matrix(diversity)

