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
setwd("/F:/Sonoran Desert Microbiome/Bacteria")

#If continuing analyses from dada2 analyses, then use the following variables
#seqtab.nochim (ASV table)
#taxa (Taxonomy table)

#Else upload the files from your desktop
bac.asv<-read.table('All_Seqs_asv_dada.txt', sep='\t', header=T, row.names=1)
bac.tax.table<-read.table('All_Seqs_silva.txt', sep='\t', header=T, row.names=1)
#Check the dimensions of your file
#Regular check can help with catching errors
dim(bac.asv) 
dim(bac.tax.table)

#Fix the rownames as they don't match the mapping file
#Strip the .fastq.gz extension and keep only name of sample @ position1
bac.rownames <- sapply(strsplit(rownames(bac.asv), ".fastq.gz"), `[`, 1)
rownames(bac.asv)<-bac.rownames


#Read metadata
map<-read.table("/MappingFile.txt", sep="\t", header=T, row.names=1, check.names=F)
dim(map)
#Remove samples not related to the project by matching samples IDs from mapping file
bac.asv<-bac.asv[rownames(map),] 
dim(bac.asv)
#Check if all sample names in ASV table match mapping file
rownames(map)==rownames(bac.asv)

#As you have removed samples,you may have empty ASV columns 
#Remove empty ASVs from taxonomy table
bac.asv<-subset(bac.asv, select=colSums(bac.asv)!=0)
dim(bac.asv)

#Remove any samples from taxonomy table that are not in ASV table 
bac.tax.table<-as.data.frame(bac.tax.table[colnames(bac.asv),])
dim(bac.tax.table)
##Check if all sample names in ASV table match taxonomy sample
rownames(bac.tax.table)==colnames(bac.asv)

#Remove contaminants:chloroplasts, mitochondria and Eukaryota
grep("Chloroplast", bac.tax.table$Order)
grep("Mitochondria", bac.tax.table$Family)
grep("Eukaryota", bac.tax.table$Kingdom)

#Remove chloroplasts from ASV and taxonomy table
bac.asv<-bac.asv[,-grep("Chloroplast", bac.tax.table$Order)]
bac.tax.table<-bac.tax.table[-grep("Chloroplast", bac.tax.table$Order),]
dim(bac.asv)
dim(bac.tax.table)

#Remove Mitochondria from ASV and taxonomy table
bac.asv<-bac.asv[,-grep("Mitochondria", bac.tax.table$Family)]
bac.tax.table<-bac.tax.table[-grep("Mitochondria", bac.tax.table$Family),]
dim(bac.asv)
dim(bac.tax.table) 

#Remove Eukaryota from ASV and taxonomy table
bac.asv<-bac.asv[,-grep("Eukaryota", bac.tax.table$Kingdom)]
bac.tax.table<-bac.tax.table[-grep("Eukaryota", bac.tax.table$Kingdom),]
dim(bac.asv)
dim(bac.tax.table)


#Check your blank DNA extracts and remove potential lab contamination
#Assign a variable to each of your blanks
bac.ext01.control<-bac.asv[which(rownames(bac.asv)=="C9"),]
bac.ext02.control<-bac.asv[which(rownames(bac.asv)=="D9"),]
bac.ext03.control<-bac.asv[which(rownames(bac.asv)=="E9"),]
bac.ext04.control<-bac.asv[which(rownames(bac.asv)=="F9"),]

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
bac.asv<-bac.asv[,-which(colnames(bac.asv)%in%names(which(bac.ext01.control>0.01*sum(bac.ext01.control))))]
bac.asv<-bac.asv[,-which(colnames(bac.asv)%in%names(which(bac.ext02.control>0.01*sum(bac.ext02.control))))]
bac.asv<-bac.asv[,-which(colnames(bac.asv)%in%names(which(bac.ext03.control>0.01*sum(bac.ext03.control))))]
bac.asv<-bac.asv[,-which(colnames(bac.asv)%in%names(which(bac.ext04.control>0.01*sum(bac.ext04.control))))]


#Substract the ASVs that are greater than 1% of the total sum from all the samples in taxonomy table
bac.tax.table<-bac.tax.table[-which(rownames(bac.tax.table)%in%names(which(bac.ext01.control>0.01*sum(bac.ext01.control)))),]
bac.tax.table<-bac.tax.table[-which(rownames(bac.tax.table)%in%names(which(bac.ext02.control>0.01*sum(bac.ext02.control)))),]
bac.tax.table<-bac.tax.table[-which(rownames(bac.tax.table)%in%names(which(bac.ext03.control>0.01*sum(bac.ext03.control)))),]
bac.tax.table<-bac.tax.table[-which(rownames(bac.tax.table)%in%names(which(bac.ext04.control>0.01*sum(bac.ext04.control)))),]

#Completely remove the blanks
bac.asv<-bac.asv[-which(rownames(bac.asv)=="C9"),]
bac.asv<-bac.asv[-which(rownames(bac.asv)=="D9"),]
bac.asv<-bac.asv[-which(rownames(bac.asv)=="E9"),]
bac.asv<-bac.asv[-which(rownames(bac.asv)=="F9"),]

dim(bac.asv)
dim(bac.tax.table) 

#Remove any empty ASVs
bac.asv<-subset(bac.asv, select=colSums(bac.asv)!=0)
bac.tax.table<-bac.tax.table[colnames(bac.asv),]
dim(bac.asv)
dim(bac.tax.table)


#Merge and taxonomic information
bac.asv.tax<-cbind(t(bac.asv), bac.tax.table)
dim(bac.asv.tax)

#Save the file
write.table(bac.asv.tax, file="asv_table_wTax_dada2.txt", sep="\t", quote=F, row.names=T, col.names=NA)

hist(rowSums(bac.asv))
sort(rowSums(bac.asv))
summary(rowSums(bac.asv)) 


#Normalization using metagenomeSeq
bac.MR <- newMRexperiment(t(bac.asv))
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
pdf("Bac_final_Richness.pdf", width=8)
ggplot(map, aes(Micorsite, Richness), color=Micosite)+
  geom_jitter(aes(color=Micosite, shape=Geolocation), alpha=0.8, size=4)+
  geom_boxplot(aes(fill=Micosite), alpha=0.3, outlier.shape = NA)+
  scale_fill_manual(values=c("C" = "#9AC801", "G" = "#3552C5"))+
  scale_color_manual(values=c("C" = "#9AC801", "G" = "#3552C5"))+
  facet_wrap(~Geolocation, nrow=1)+
  theme_bw()
dev.off()

#Test differences among treatments (richness)
anova(lm(map$Richness~map$Geolocation+map$Geolocation/map$Microsite))


#NMDS ordination using Bray-Curtis distance
bac.dist<-vegdist(bac.norm, method="bray")
bac.nmds<-metaMDS(bac.dist, k=2)
bac.nmds$stress 

map$Axis01<-bac.nmds$points[,1]
map$Axis02<-bac.nmds$points[,2]

#Save the NMDS plots
pdf("Bac_final_NMDS.pdf", width=8)
ggplot(map, aes(Axis01, Axis02))+
  geom_point(aes(color=Micosite, shape=Geolocation), alpha=0.8, size=4)+
  scale_fill_manual(values=c("C" = "#9AC801", "G" = "#3552C5"))+
  scale_color_manual(values=c("C" = "#9AC801", "G" = "#3552C5"))+
  theme_bw()+
  theme(legend.position="right")
dev.off()

#Perform PERMANOVA to test community similarity differences
adonis(bac.dist~map$Geolocation/map$Microsite, strata=map$Geolocation)

#Calculate Beta-dispersivity
bac.dispersion<-betadisper(bac.dist, group=map$Microsite,type = "centroid")
#Check Average distance to centroid values
bac.dispersion
#Test whether the beta-dispersivity across Microsite is different
permutest(bac.dispersion)

#Plot the distance to centroid values
pdf("Bac_betadispersion_microsite_centroid.pdf", width=8)
boxplot(bac.dispersion, xlab="Microsite", ylab="Distance to centroid",
        col=c("Canopy" = "#9AC801", "Gap" = "#3552C5"))
dev.off()