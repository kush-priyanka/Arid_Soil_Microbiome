#Install packages
install.packages("vegan")
install.packages("ggplot2")
install.packages("devtools")
install.packages("reshape")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("metagenomeSeq")

#Load libraries
library(vegan) #for multivariate statistics
library(ggplot2) #for visualization and plotting
library(devtools)
library(reshape)
library(metagenomeSeq) #for ASV table normalization

#Set working directory to the folder where your data files are
setwd("/home/pree/Arid_Soil_Microbiome")

#If continuing analyeses from dada2 analyses, then use the following variables
#seqtab.nochim (ASV table)
#taxa (Taxonomic table)
#Else upload the files from your desktop
bac.otu<-read.table('All_Seqs_otu_dada.txt', sep='\t', header=T, row.names=1)
bac.tax.table<-read.table('All_Seqs_silva.txt', sep='\t', header=T, row.names=1)
#Check the dimensions of your file
#Regular check can help with catching errors
dim(bac.otu) 
dim(bac.tax.table)

#Fix the rownames as they dont match the mapping file
#Strip the .fastq.gz extension and keep only name of sample @ position1
bac.rownames <- sapply(strsplit(rownames(bac.aj.otu), ".fastq.gz"), `[`, 1)
rownames(bac.aj.otu)<-bac.rownames


#Read metadata
map<-read.table("/MappingFile.txt", sep="\t", header=T, row.names=1, check.names=F)
dim(map)
#Remove samples not related to the project by matching samples IDs from mapping file
bac.otu<-bac.otu[rownames(map),] 
dim(bac.otu)
#Check if all sample names in ASV table match mapping file
rownames(map)==rownames(bac.otu)

#As you have removed samples,you may have empty OTU columns 
#Remove empty OTUs from taxonomy table
bac.otu<-subset(bac.otu, select=colSums(bac.otu)!=0)
dim(bac.otu)

#Remove any samples from taxonomy table that are not in ASV table 
bac.tax.table<-as.data.frame(bac.tax.table[colnames(bac.otu),])
dim(bac.tax.table)
##Check if all sample names in ASV table match taxonomy sample
rownames(bac.tax.table)==colnames(bac.otu)

#Remove contaminants:chloroplasts, mitochondria and Eukaryota
grep("Chloroplast", bac.tax.table$Order)
grep("Mitochondria", bac.tax.table$Family)
grep("Eukaryota", bac.tax.table$Kingdom)

#Remove chloroplasts from ASV and taxonomy table
bac.otu<-bac.otu[,-grep("Chloroplast", bac.tax.table$Order)]
bac.tax.table<-bac.tax.table[-grep("Chloroplast", bac.tax.table$Order),]
dim(bac.otu)
dim(bac.tax.table)

#Remove Mitochondria from ASV and taxonomy table
bac.otu<-bac.otu[,-grep("Mitochondria", bac.tax.table$Family)]
bac.tax.table<-bac.tax.table[-grep("Mitochondria", bac.tax.table$Family),]
dim(bac.otu)
dim(bac.tax.table) 
#Remove Eukaryota from ASV and taxonomy table
bac.otu<-bac.otu[,-grep("Eukaryota", bac.tax.table$Kingdom)]
bac.tax.table<-bac.tax.table[-grep("Eukaryota", bac.tax.table$Kingdom),]
dim(bac.otu)
dim(bac.tax.table)


#Check your blank DNA extracts and remove potential lab contamination
#Assign a variable to each of your blanks
bac.ext01.control<-bac.otu[which(rownames(bac.otu)=="C9"),]
bac.ext02.control<-bac.otu[which(rownames(bac.otu)=="D9"),]
bac.ext03.control<-bac.otu[which(rownames(bac.otu)=="E9"),]
bac.ext04.control<-bac.otu[which(rownames(bac.otu)=="F9"),]

#Check the number of ASVs in each blank
sum(bac.ext01.control)
sum(bac.ext02.control)
sum(bac.ext03.control)
sum(bac.ext04.control)

#Now check, which ASV is greater than 1% of the total ASVs in that blank
which(bac.ext01.control>0.01*sum(bac.ext01.control))
which(bac.ext02.control>0.01*sum(bac.ext01.control))
which(bac.ext03.control>0.01*sum(bac.ext01.control))
which(bac.ext04.control>0.01*sum(bac.ext01.control))

#Now check, taxa information for the ASV that is greater than 1% of the total ASVs
bac.tax.table[names(which(bac.ext01.control>0.01*sum(bac.ext01.control))),]
bac.tax.table[names(which(bac.ext02.control>0.01*sum(bac.ext02.control))),]
bac.tax.table[names(which(bac.ext03.control>0.01*sum(bac.ext03.control))),]
bac.tax.table[names(which(bac.ext04.control>0.01*sum(bac.ext04.control))),]

#Substract the ASVs that are greater than 1% of the total sum from all the samples in ASV table 
bac.otu<-bac.otu[,-which(colnames(bac.otu)%in%names(which(bac.ext01.control>0.01*sum(bac.ext01.control))))]
bac.otu<-bac.otu[,-which(colnames(bac.otu)%in%names(which(bac.ext02.control>0.01*sum(bac.ext02.control))))]
bac.otu<-bac.otu[,-which(colnames(bac.otu)%in%names(which(bac.ext03.control>0.01*sum(bac.ext03.control))))]
bac.otu<-bac.otu[,-which(colnames(bac.otu)%in%names(which(bac.ext04.control>0.01*sum(bac.ext04.control))))]

#Substract the ASVs that are greater than 1% of the total sum from all the samples in taxonomy table
bac.tax.table<-bac.tax.table[-which(rownames(bac.tax.table)%in%names(which(bac.ext01.control>0.01*sum(bac.ext01.control)))),]
bac.tax.table<-bac.tax.table[-which(rownames(bac.tax.table)%in%names(which(bac.ext02.control>0.01*sum(bac.ext02.control)))),]
bac.tax.table<-bac.tax.table[-which(rownames(bac.tax.table)%in%names(which(bac.ext03.control>0.01*sum(bac.ext03.control)))),]
bac.tax.table<-bac.tax.table[-which(rownames(bac.tax.table)%in%names(which(bac.ext04.control>0.01*sum(bac.ext04.control)))),]

#Completely remove the blanks
bac.aj.otu<-bac.aj.otu[-which(rownames(bac.aj.otu)=="C9"),]
bac.aj.otu<-bac.aj.otu[-which(rownames(bac.aj.otu)=="D9"),]
bac.aj.otu<-bac.aj.otu[-which(rownames(bac.aj.otu)=="E9"),]
bac.aj.otu<-bac.aj.otu[-which(rownames(bac.aj.otu)=="F9"),]

dim(bac.aj.otu)
dim(bac.aj.tax.table) 

#Remove any empty ASVs
bac.otu<-subset(bac.otu, select=colSums(bac.otu)!=0)
bac.tax.table<-bac.tax.table[colnames(bac.otu),]
dim(bac.aj.otu)
dim(bac.aj.tax.table)

#Merge and taxonomic information
bac.otu.tax<-cbind(t(bac.otu), bac.tax.table)
dim(bac.otu.tax)

#Save the file
write.table(bac.otu.tax, file="otu_table_wTax_dada2.txt", sep="\t", quote=F, row.names=T, col.names=NA)

#Check number of sequences per sample 
hist(rowSums(bac.otu))
sort(rowSums(bac.otu))
summary(rowSums(bac.otu)) 



#Normalization using metagenomeSeq
bac.MR <- newMRexperiment(t(bac.otu))
p <- cumNormStat(bac.MR)
bac.MR <- cumNorm(bac.MR, p=p)
bac.norm <- t(MRcounts(bac.MR, norm=T, log=F))

#Check number of sequences per normalized sample
hist(rowSums(bac.norm))
sort(rowSums(bac.norm))
summary(rowSums(bac.norm)) 

dim(bac.norm) 

#Calculate Richness and Shannon
map$Richness<-specnumber(bac.norm) 
map$Shannon<-diversity(bac.norm, index="shannon")

#Overwrite the mapping filing with Richness and Shannon
write.table(map, file="Mapping.File.txt", sep="\t", quote=F, row.names=T, col.names=NA)
mean(map$Richness) 
mean(map$Shannon) 

#Save Richness plot
pdf("bac_richness.pdf", width=8)
ggplot(map, aes(Geolocation, Richness), color=Microsite)+
  geom_jitter(aes(color=Geolocation, shape=Microsite), alpha=0.8, size=4)+
  geom_boxplot(aes(fill=Geolocation), alpha=0.3, outlier.shape = NA)+
  facet_wrap(~Microsite, nrow=1)+
  theme_bw()
dev.off()

#Test differences among treatments (richness)
anova(lm(map$Richness~map$Geolocation))


#NMDS ordination using bray-curtis distance
bac.dist<-vegdist(bac.norm, method="bray")
bac.nmds<-metaMDS(bac.dist, k=2)
bac.nmds$stress 

map$Axis01<-bac.nmds$points[,1]
map$Axis02<-bac.nmds$points[,2]

#Save the NMDS plots
pdf("bac_nmds_dada2.pdf", width=8)
ggplot(map, aes(Axis01, Axis02))+
  geom_point(aes(color=Microsite, shape=Geolocation), alpha=0.8, size=4)+
  theme_bw()+
  theme(legend.position="right")
dev.off()

#Perform PERMANOVA to see the variation in microbial diversity
#using the variable Geolocation 
adonis(bac.dist~map$Geolocation)

