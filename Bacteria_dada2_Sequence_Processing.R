####dada2####
library(dada2)
packageVersion("dada2") #confirm that DADA2 version is 1.6 or later

#Follow the dada2 pipeline 
#https://benjjneb.github.io/dada2/tutorial_1_8.html

path <- "/F:/Sonoran Desert Microbiome/Downloads/idemp-master/demultiplexed/"

#Set working directory to the unzipped mock dataset folder, here located in downloads
#Here this is done by defining the variable path
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#Read in the files in forward and reverse 
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 6)
sample.names <- sapply(strsplit(basename(sample.names), "\\."), `[`, 1)

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

#Visualize the quality profiles of the forward reads
#forward reads are generally in good, so trim forward reads at 240 bp
plotQualityProfile(fnFs[1:2])

#Visualize the quality profiles of the reverse reads
#reverse reads tend to be poorer, so trim reverse reads at 160 bp
plotQualityProfile(fnRs[1:2])


#Assign the filenames for the filtered fastq.gz files
filt_path <- file.path(path, 'filtered')
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#Trim first 10 base pairs (trimLeft = 10) because the first 10 tend to be poor quality
#Trim reads at 160 reads for both fwd and rev based on QC graphs
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=10, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
View(out)


#Randomly choosing 8 samples to look at and to learn errors ("1:8")
#Predict and calculate error rates from forward and reverse reads
#Computationally intensive step for all samples, so use 25% of the samples
errF <- learnErrors(filtFs[1:8], multithread=TRUE)
errR <- learnErrors(filtRs[1:8], multithread=TRUE)

#Plot the error rates
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

#DEREPLICATION - this step matches exact sequences from all samples to each other and records their abundance
#This step finds and determines the exact sequence variants
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Check the number of unique sequences from each sample 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#Describe data from the forward reads of the first sample
dadaFs[[1]]

#Merge the trimmed paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[1])

#Create sequence ('ASV') read table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#Determine the percentage of total sequence reads that are NOT chimeras 
#Percentage of reads that are chimeras should not be too high
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline to look how many reads went through the pipeline
#trimmed, filtered, chimeria removal, etc.
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
View(track)

#ASSIGN TAXONOMY - dada gives the option of using RDP, greengeens, or Silva as the database for assigning taxonomy
#IMPORTANT - go to this site https://benjjneb.github.io/dada2/training.html and download the training set
#and place it in the same working directory 'path' before starting the command. For now just use the silva set.
#NOTE: Provide the full/absolute path to the Tax file otherwise you will get an error.
taxa <- assignTaxonomy(seqtab.nochim, "/F:/Sonoran Desert Microbiome/Downloads/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

#SPECIES ASSIGNMENT -available for silva and RDP
#Be sure to download the appropriate files and add to the same directory as 'path'
taxa <- addSpecies(taxa, "/F:/Sonoran Desert Microbiome/Downloads/silva_species_assignment_v132.fa.gz")

#view the assigned taxonomic table
#NOTE:This is the end of the DADA2 pipeline 
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
View(taxa.print)

#Export ASV and taxonomy tables
setwd("/F:/Sonoran Desert Microbiome/Results")
write.table(taxa, "All_Seqs_silva.txt", sep="\t", quote=F, row.names=T, col.names=NA)
write.table(seqtab.nochim, "All_Seqs_asv_dada.txt", sep="\t", quote=F, row.names=T, col.names=NA)
saveRDS(seqtab.nochim, "All_Seqs_asv_dada.rds")
write.table(track, "All_Seqs_Reads.txt", sep="\t", quote=F, row.names=T, col.names=NA)


