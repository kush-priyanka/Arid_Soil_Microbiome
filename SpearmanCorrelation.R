#Calculate Spearman's correlation

#Import data file
dat<-read.table(file.choose(), sep="\t", header=T) 

#Spearman correlations between TOC and different Nitrogen forms
a<-cor.test(dat$TOC,dat$TN,method="spearman")
b<-cor.test(dat$TOC,dat$DN,method="spearman")
c<-cor.test(dat$TOC,dat$Ammonium,method="spearman")
d<-cor.test(dat$TOC,dat$Nitrate,method="spearman")

#Spearman correlations between different Nitrogen forms and gene abundance
e<-cor.test(dat$TN,dat$log_copy_number_16S,method="spearman")
f<-cor.test(dat$TN,dat$log_copy_number_amoA,method="spearman")
g<-cor.test(dat$TN,dat$log_copy_number_ureC,method="spearman")

h<-cor.test(dat$DN,dat$log_copy_number_16S,method="spearman")
i<-cor.test(dat$DN,dat$log_copy_number_amoA,method="spearman")
j<-cor.test(dat$DN,dat$log_copy_number_ureC,method="spearman")

k<-cor.test(dat$Ammonium,dat$log_copy_number_16S,method="spearman")
l<-cor.test(dat$Ammonium,dat$log_copy_number_amoA,method="spearman")
m<-cor.test(dat$Ammonium,dat$log_copy_number_ureC,method="spearman")

n<-cor.test(dat$Nitrate,dat$log_copy_number_16S,method="spearman")
o<-cor.test(dat$Nitrate,dat$log_copy_number_amoA,method="spearman")
p<-cor.test(dat$Nitrate,dat$log_copy_number_ureC,method="spearman")

q<-cor.test(dat$Urease,dat$log_copy_number_ureC,method="spearman")

#Create a variable with all the p-values
p.value<-c(a$p.value,b$p.value,c$p.value,d$p.value,
           e$p.value,f$p.value,g$p.value,h$p.value,i$p.value,
           j$p.value,k$p.value,l$p.value,m$p.value,n$p.value, o$p.value,
           p$p.value, q$p.value)

#False discovery rate correction for p-values
fdr.corr<-p.adjust(p.value,method="BH")
fdr.corr
