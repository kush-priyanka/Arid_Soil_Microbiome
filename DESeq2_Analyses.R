install.packages("DESeq2")
#Might have to install DESeq2 using BiocManager 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

#load libraries
library(DESeq2)

#Import bacterial/archaeal or fungal ASV table and metadata file
bac.asv <- read.table("F:/Sonoran Desert Microbiome/otu_table_wTax_dada2_16S.txt", sep = "\t", header = T, row.names=1)
metadata.deseq<-read.table("F:/Sonoran Desert Microbiome/Mapping.File.Richness_16S.txt", sep = "\t", header = T, row.names=1)

#Remove samples not related to the project

#Remove samples not related to the project
bac.asv<-bac.asv[rownames(metadata.deseq),]
dim(bac.asv) #28 18554
all(rownames(metadata.deseq)==rownames(bac.asv)) #TRUE

#Separate the taxonomy information into another variable
bac.tax<-bac.asv
dim(bac.tax)# 18554 41
bac.tax<- bac.tax[, c(35:41)]
dim(bac.tax)# 18554 7

#Remove empty ASVs
bac.asv<-bac.asv[rownames(metadata),]
bac.asv<-subset(bac.asv, select=colSums(bac.asv)!=0) #28 16484

#Match ASV with the taxonomy table
bac.tax<-bac.tax[colnames(bac.asv),]
dim(bac.tax)# 16484 7
bac.tax<-as.data.frame(bac.tax)

bac.asv<-t(bac.asv)
bac.asv<-as.data.frame(bac.asv)


#Metadata for DESeq analysis
metadata.deseq$Microsite <- factor(metadata.deseq$Microsite, levels = c("C", "G"))
metadata.deseq$Geolocation <- factor(metadata.deseq$Geolocation, levels = c("1", "2", "3", "4"))

#Generate Phylum table
bac.phylum <- bac.asv
bac.phylum$Phylum <- bac.tax$Phylum
bac.phylum <- aggregate(.~Phylum, data = bac.phylum, sum)
rownames(bac.phylum) <- bac.phylum$Phylum
bac.phylum <- bac.phylum[,-1]

#Generate Class table
bac.class <- bac.asv
bac.class$Class <- bac.tax$Class
bac.class <- aggregate(.~Class, data = bac.class, sum)
rownames(bac.class) <- bac.class$Class
bac.class <- bac.class[,-1]

#Generate Order table
bac.order <- bac.asv
bac.order$Order <- bac.tax$Order
bac.order <- aggregate(.~Order, data = bac.order, sum)
rownames(bac.order) <- bac.order$Order
bac.order <- bac.order[,-1]

#Generate Family table
bac.family <- bac.asv
bac.family$Family <- bac.tax$Family
bac.family <- aggregate(.~Family, data = bac.family, sum)
rownames(bac.family) <- bac.family$Family
bac.family <- bac.family[,-1]

#Generate Genus table
bac.genus <- bac.asv
bac.genus$Genus <- bac.tax$Genus
bac.genus <- aggregate(.~Genus, data = bac.genus, sum)
rownames(bac.genus) <- bac.genus$Genus
bac.genus <- bac.genus[,-1]

# Generate species table
bac.species <- bac.asv
bac.species$Species <- bac.tax$Species
bac.species <- aggregate(.~Species, data = bac.species, sum)
rownames(bac.species) <- bac.species$Species
bac.species <- bac.species[,-1]

#DESeq analysis at Phylum
bac.phylum.dds <- DESeqDataSetFromMatrix(countData = bac.phylum, 
                                         colData = metadata.deseq,
                                         design = ~ Geolocation + Microsite)

bac.phylum.dds <- DESeq(bac.phylum.dds)

bac.phylum.deseq.res <- results(bac.phylum.dds)
bac.phylum.deseq.res <- results(bac.phylum.dds,contrast=c("Microsite", "C","G"))
head(bac.phylum.deseq.res)

bac.phylum.deseq.res <- as.data.frame(bac.phylum.deseq.res)

bac.phylum.deseq.sig <- subset(bac.phylum.deseq.res, bac.phylum.deseq.res$padj < 0.05)

#DESeq analysis at Class
bac.class.dds <- DESeqDataSetFromMatrix(countData = bac.class, 
                                        colData = metadata.deseq,
                                        design = ~ Geolocation + Microsite)

bac.class.dds <- DESeq(bac.class.dds)

bac.class.deseq.res <- results(bac.class.dds ,contrast=c("Microsite", "C","G")) 

bac.class.deseq.res <- as.data.frame(bac.class.deseq.res)

bac.class.deseq.sig <- subset(bac.class.deseq.res, bac.class.deseq.res$padj < 0.05)


#DESeq analysis at Order
bac.order.dds <- DESeqDataSetFromMatrix(countData = bac.order, 
                                        colData = metadata.deseq,
                                        design = ~ Geolocation + Microsite)

bac.order.dds <- DESeq(bac.order.dds)

bac.order.deseq.res <- results(bac.order.dds,contrast=c("Microsite", "C","G")) 

bac.order.deseq.res <- as.data.frame(bac.order.deseq.res)

bac.order.deseq.sig <- subset(bac.order.deseq.res, bac.order.deseq.res$padj < 0.05)



#DESeq analysis at Family
bac.family.dds <- DESeqDataSetFromMatrix(countData = bac.family, 
                                         colData = metadata.deseq,
                                         design = ~ Geolocation + Microsite)

bac.family.dds <- DESeq(bac.family.dds)

bac.family.deseq.res <- results(bac.family.dds ,contrast=c("Microsite", "C","G")) 

bac.family.deseq.res <- as.data.frame(bac.family.deseq.res)

bac.family.deseq.sig <- subset(bac.family.deseq.res, bac.family.deseq.res$padj < 0.05)


#DESeq analysis at Genus
bac.genus.dds <- DESeqDataSetFromMatrix(countData = bac.genus, 
                                        colData = metadata.deseq,
                                        design = ~ Geolocation + Microsite)

bac.genus.dds <- DESeq(bac.genus.dds)

bac.genus.deseq.res <- results(bac.genus.dds,contrast=c("Microsite", "C","G")) 

bac.genus.deseq.res <- as.data.frame(bac.genus.deseq.res)

bac.genus.deseq.sig <- subset(bac.genus.deseq.res, bac.genus.deseq.res$padj < 0.05)


#DESeq at species
bac.species.dds <- DESeqDataSetFromMatrix(countData = bac.species, 
                                          colData = metadata.deseq,
                                          design = ~ Geolocation + Microsite)

bac.species.dds <- DESeq(bac.species.dds)

bac.species.deseq.res <- results(bac.species.dds,contrast=c("Microsite", "C","G")) 

bac.species.deseq.res <- as.data.frame(bac.species.deseq.res)

bac.species.deseq.sig <- subset(bac.species.deseq.res, bac.species.deseq.res$padj < 0.05)

#Check the results
bac.phylum.deseq.sig
bac.order.deseq.sig
bac.class.deseq.sig
bac.family.deseq.sig
bac.genus.deseq.sig
bac.species.deseq.sig

#Save DESeq2 results 
write.table(bac.phylum.deseq.sig, file="bac.phylum.deseq.sig.txt", sep="\t", quote=F, row.names=T, col.names=NA)
write.table(bac.order.deseq.sig, file="bac.order.deseq.sig.txt", sep="\t", quote=F, row.names=T, col.names=NA)
write.table(bac.class.deseq.sig, file="bac.class.deseq.sig.txt", sep="\t", quote=F, row.names=T, col.names=NA)
write.table(bac.family.deseq.sig, file="bac.family.deseq.sig.txt", sep="\t", quote=F, row.names=T, col.names=NA)
write.table(bac.genus.deseq.sig, file="bac.genus.deseq.sig.txt", sep="\t", quote=F, row.names=T, col.names=NA)
write.table(bac.species.deseq.sig, file="bac.species.deseq.sig.txt", sep="\t", quote=F, row.names=T, col.names=NA)

