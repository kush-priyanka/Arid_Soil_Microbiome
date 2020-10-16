#load libraries
library(vegan) #for multivariate statistics
library(ggplot2) #for visualization and plotting
library(metagenomeSeq) #for normalization

#set working directory
setwd("/F:/Sonoran Desert Microbiome")
#If continuing analyses from dada2 analyses, then use the following variables
#seqtab.nochim (ASV table)
#taxa (Taxonomic table)

#Else upload the files from your desktop
fungi.asv<-read.table('All_Seqs_asv_dada1.10.txt', sep='\t', header=T, row.names=1)
fungi.tax.table<-read.table('All_Seqs_unite.txt', sep='\t', header=T, row.names=1)
#Check the dimensions of your file
#Regular check can help with catching errors
dim(fungi.asv) 
dim(fungi.tax.table)

#Read metadata
map<-read.table("Desert2017_mapping_feb_ITS.txt", sep="\t", header=T, row.names=1, check.names=F)

#Remove "k__Eukaryota_kgd_Incertae_sedis" "k__Metazoa" "k__Rhizaria"
fungi.asv<-fungi.asv[,grep("k__Fungi", fungi.tax.table$Kingdom)]
fungi.tax.table<-fungi.tax.table[grep("k__Fungi", fungi.tax.table$Kingdom),]

#Remove samples not related to the project
fungi.asv<-fungi.asv[rownames(map),]
rownames(map)==rownames(fungi.asv)

#Remove Barren (B) samples
map<-subset(map, map$Eco_type!="B")

#Remove Geo_location 5 (Tumamoc Hill) samples
map<-subset(map, map$Geo_location!=5)
map$Geo_location<-as.factor(map$Geo_location)

map$Eco_type<-droplevels(map$Eco_type)
map$Geo_location<-droplevels(map$Geo_location)

fungi.asv<-fungi.asv[rownames(map),]
fungi.asv<-subset(fungi.asv, select=colSums(fungi.asv)!=0) #28 2237

#Remove empty OTUs from taxonomy table
fungi.tax.table<-as.data.frame(fungi.tax.table[colnames(fungi.asv),])
rownames(fungi.tax.table)==colnames(fungi.asv)

#Remove potential lab contamination
#Check your blank DNA extracts and remove potential lab contamination
#Assign a variable to each of your blanks
fungi.ext01.control<-fungi.asv[which(rownames(fungi.asv)=="JulieITS.00106"),]
fungi.ext02.control<-fungi.asv[which(rownames(fungi.asv)=="JulieITS.00107"),]
fungi.ext03.control<-fungi.asv[which(rownames(fungi.asv)=="JulieITS.00108"),]
fungi.ext04.control<-fungi.asv[which(rownames(fungi.asv)=="JulieITS.00109"),]

#Check total number of ASVs in each blank
sum(fungi.ext01.control)
sum(fungi.ext02.control)
sum(fungi.ext03.control)
sum(fungi.ext04.control)

#Now check, which ASV is greater than 1% of the total ASVs in that blank
which(fungi.ext01.control>0.01*sum(fungi.ext01.control))
which(fungi.ext02.control>0.01*sum(fungi.ext02.control))
which(fungi.ext03.control>0.01*sum(fungi.ext03.control))
which(fungi.ext04.control>0.01*sum(fungi.ext04.control))

#who are the contaminants?
fungi.tax.table[which(fungi.ext01.control>0.01*sum(fungi.ext01.control)),]
fungi.tax.table[which(fungi.ext02.control>0.01*sum(fungi.ext02.control)),]
fungi.tax.table[which(fungi.ext03.control>0.01*sum(fungi.ext03.control)),]
fungi.tax.table[which(fungi.ext04.control>0.01*sum(fungi.ext04.control)),]

#Substract the ASVs that are greater than 1% of the total sum from all the samples in ASV table 
fungi.asv<-fungi.asv[,-which(fungi.ext01.control>0.01*sum(fungi.ext01.control))]
fungi.asv<-fungi.asv[,-which(fungi.ext03.control>0.01*sum(fungi.ext03.control))]


#Substract the ASVs that are greater than 1% of the total sum from all the samples in taxonomy table
fungi.tax.table<-fungi.tax.table[-which(fungi.ext01.control>0.01*sum(fungi.ext01.control)),]
fungi.tax.table<-fungi.tax.table[-which(fungi.ext03.control>0.01*sum(fungi.ext03.control)),]

#Now we can remove control samples
fungi.asv<-fungi.asv[-which(rownames(fungi.asv)=="JulieITS.00106"),]
fungi.asv<-fungi.asv[-which(rownames(fungi.asv)=="JulieITS.00107"),]
fungi.asv<-fungi.asv[-which(rownames(fungi.asv)=="JulieITS.00108"),]
fungi.asv<-fungi.asv[-which(rownames(fungi.asv)=="JulieITS.00109"),]

#Remove empty ASVs
fungi.asv<-subset(fungi.asv, select=colSums(fungi.asv)!=0) #24 2234
fungi.tax.table<-fungi.tax.table[colnames(fungi.asv),]

#Create table with ASV and taxonomic information
fungi.asv.tax<-cbind(t(fungi.asv), fungi.tax.table)

#write.table(fungi.asv.tax, file="ITS_asv_table_wTax_dada2_1.10.txt", sep="\t", quote=F, row.names=T, col.names=NA)

#check number of sequences per sample
hist(rowSums(fungi.asv))
sort(rowSums(fungi.asv))

summary(rowSums(fungi.asv)) #min 43008 - max 206480

#Normalization (metagenomeSeq)
fungi.MR <- newMRexperiment(t(fungi.asv))
p <- cumNormStat(fungi.MR)
fungi.MR <- cumNorm(fungi.MR, p=p)
fungi.norm <- t(MRcounts(fungi.MR, norm=T, log=F))

#match map with OTU table
map<-map[rownames(fungi.asv),]

#richness and shannon
map$Richness<-specnumber(fungi.norm)
map$Shannon<-diversity(fungi.norm, index="shannon")

#write.table(map, file="ITS.Mapping.Richness.txt", sep="\t", quote=F, row.names=T, col.names=NA)

mean(map$Richness) #217
mean(map$Shannon) #3.64

#Overwrite the mapping filing with Richness and Shannon
write.table(map, file="ITS.Mapping.File.txt", sep="\t", quote=F, row.names=T, col.names=NA)

#Save Richness plot
pdf("Fungi_final_Richness.pdf", width=8)
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

#NMDS ordination
fungi.dist<-vegdist(fungi.norm, method="bray")
fungi.nmds<-metaMDS(fungi.dist, k=2)
fungi.nmds$stress #0.1343737

map$Axis01<-fungi.nmds$points[,1]
map$Axis02<-fungi.nmds$points[,2]

pdf("Fungi_final_NMDS.pdf", width=8)
ggplot(map, aes(Axis01, Axis02))+
  geom_point(aes(color=Microsite, shape=Geolocation), alpha=0.8, size=4)+
  scale_fill_manual(values=c("C" = "#9AC801", "G" = "#3552C5"))+
  scale_color_manual(values=c("C" = "#9AC801", "G" = "#3552C5"))+
  theme_bw()+
  theme(legend.position="right")
dev.off()

#Perform PERMANOVA to test community similarity differences
adonis(fungi.dist~map$Geolocation/map$Microsite, strata=map$Geolocation)

#Calculate beta-dispersivity
fungi.dispersion<-betadisper(fungi.dist, group=map$Microsite,type = "centroid")
#Check Average distance to centroid values
fungi.dispersion
#Test whether the beta-dispersivity across Microsite is different
permutest(fungi.dispersion)

#Plot the distance to centroid values
pdf("Fungi_betadispersion_microsite_centroid.pdf", width=8)
boxplot(fungi.dispersion, xlab="Microsite", ylab="Distance to centroid",
        col=c("Canopy" = "#9AC801", "Gap" = "#3552C5"))
dev.off()

#Create table for FUNGuild
fun.funguild<-fungi.tax.table #after processing
fun.funguild$taxonomy<-paste(fungi.tax.table$Kingdom, fungi.tax.table$Phylum, fungi.tax.table$Class, 
                             fungi.tax.table$Order, fungi.tax.table$Family, fungi.tax.table$Genus,
                             fungi.tax.table$Species, sep=";")

fun.asv.funguild<-cbind(t(fungi.asv), fun.funguild$taxonomy)
colnames(fun.asv.funguild)[ncol(fun.asv.funguild)]<-'taxonomy'

#Save file
write.table(fun.asv.funguild, file="aridity_asv_table_wTax_dada2_funguild.txt", sep="\t", quote=F, row.names=T, col.names=NA)

#Read funguild results
funguild<-t(read.table("aridity_asv_table_wTax_dada2_funguild.guilds.txt", sep="\t", header=T, row.names=1, check.names=F))
fungi.asv<-as.data.frame(funguild[1:24,])
fungi.asv<-apply(fungi.asv, 2, function (x) {as.numeric(as.character(x))})
rownames(fungi.asv)<-rownames(funguild[1:24,])

guilds<-as.data.frame(apply((apply(fungi.asv, 1, function (x) by(x, as.factor(funguild[29,]), sum))), 2, function (x) {x/sum(x)}))

#Save file
write.table(guilds, file="funguild_results.txt", sep="\t", quote=F, row.names=T, col.names=NA)

funguild<-t(read.table("funguild_results_filt.txt", sep="\t", header=T, row.names=1))
funguild<-subset(funguild, select=colSums(funguild)!=0)

#Combine the funguild results with the mapping file
map.funguild<-cbind(map, funguild)

#Save file
write.table(map.funguild, file="funguild_results_mapping.txt", sep="\t", quote=F, row.names=T, col.names=NA)



