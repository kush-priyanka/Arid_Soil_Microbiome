####dada2####
#Follow the dada2 pipeline for ITS:
#https://benjjneb.github.io/dada2/ITS_workflow.html

install.packages("dada")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.10")
library(dada2)
packageVersion("dada2") #confirm that DADA2 version is 1.8 or later (1.10)
library(ShortRead)
packageVersion("ShortRead") #1.40.0
library(Biostrings)
packageVersion("Biostrings") #2.50.1

path <- "/F:/Sonoran Desert Microbiome/Downloads/idemp-master/demultiplexed_aridityITS/"

#Set working directory to the unzipped mock dataset folder, here located in downloads
#Here this is done by defining the variable path
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#Identify primers
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  #forward primer sequence
REV <- "GCTGCGTTCTTCATCGATGC"  #reverse primer sequence

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

FWD.orients
REV.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#Count the number of times the primers appear in the forward and reverse read, 
#while considering all possible primer orientations
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[2]]))

#Remove primers (cutadapt tool)
#CHANGE the cutadapt path on your machine
#Run shell commands from R (version 1.18)
cutadapt <- "/usr/local/bin/cutadapt" 
system2(cutadapt, args = "--version") 

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#Count the presence of primers in the first cutadapt-ed samples
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[2]]))

#Primers should not be detected in the cutadapted reads now.

#Forward and reverse fastq filenames must have the format:
cutFs <- sort(list.files(path.cut, pattern="_R1_001.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern="_R2_001.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have same format:
sample.names <- sapply(strsplit(basename(cutFs), "_"), `[`, 6)
sample.names <- sapply(strsplit(basename(sample.names), "\\.f"), `[`, 1)
head(sample.names)

#Inspect read quality profiles
plotQualityProfile(cutFs[1:4])
plotQualityProfile(cutRs[1:4])

#Filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

head(out)
View(out)

#Learn the error rates
errF <- learnErrors(filtFs[11:25], multithread = TRUE)
errR <- learnErrors(filtRs[11:25], multithread = TRUE)

#Plot error rates
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

#Dereplicate identical reads
derepFs <- derepFastq(filtFs[-which(out[,2]==0)], verbose = TRUE)
derepRs <- derepFastq(filtRs[-which(out[,2]==0)], verbose = TRUE)

#Name the derep-class objects by the sample names
names(derepFs) <- sample.names[-which(out[,2]==0)]
names(derepRs) <- sample.names[-which(out[,2]==0)]

#Sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#Construct sequence (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

#Determine the percentage of total sequence reads that are NOT chimeras 
#Percentage of reads that are chimeras should not be too high
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
#trimmed, filtered, chimeria removal, etc.
getN <- function(x) sum(getUniques(x))
track <- cbind(out[-which(out[,2]==0),], sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names[-which(out[,2]==0)]
head(track)

#Assign taxonomy using UNITE ITS database
#IMPORTANT - go to this site https://benjjneb.github.io/dada2/training.html and download the training set
#Place it in the same working directory 'path' before starting the command. 
taxa <- assignTaxonomy(seqtab.nochim, "/F:/Sonoran Desert Microbiome/Downloads/sh_general_release_dynamic_s_01.12.2017.fasta", multithread=TRUE, tryRC=T)

#Removing sequence rownames for display only
taxa.print <- taxa  
rownames(taxa.print) <- NULL
head(taxa.print)
View(taxa.print)

#Export ASV and taxonomy tables to path for upload and analysis
setwd("/F:/Sonoran Desert Microbiome/ITS/")
write.table(taxa, "All_Seqs_unite.txt", sep="\t", quote=F, row.names=T, col.names=NA)
write.table(seqtab.nochim, "All_Seqs_asv_dada.txt", sep="\t", quote=F, row.names=T, col.names=NA)
saveRDS(seqtab.nochim, "All_Seqs_asv_dada.rds")

write.table(track, "All_Seqs_Reads.txt", sep="\t", quote=F, row.names=T, col.names=NA)
