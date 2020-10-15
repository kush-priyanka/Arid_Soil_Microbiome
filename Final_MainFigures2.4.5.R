library(vegan) #for multivariate statistics
library(ggplot2)
library(dplyr)
install.packages("ggpubr")
library("ggpubr")

#Import data file
dat<-read.table(file.choose() , sep='\t', header=T, row.names=1, check.names=F)
dim(dat)

dat$Geolocation <- as.factor(dat$Geolocation)

#Figure 2A & C, Bacterial/Archaeal and Fungal Richness
bac_rich<-dat %>%
  group_by(Geolocation, Microsite) %>%
  summarise(
    sd = sd(Bacterial_Richness, na.rm = TRUE),
    len = mean(Bacterial_Richness)
  )
bac_rich

ITS_rich<-dat %>%
  group_by(Geolocation, Microsite) %>%
  summarise(
    sd = sd(ITS_Richness, na.rm = TRUE),
    len = mean(ITS_Richness)
  )
ITS_rich

pdf("16S_Richness.pdf", width=8)
ggplot(bac_rich, aes(Microsite, len)) +
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, color = Microsite),
    position = position_dodge(0.3), width = 0.3,size=0.8)+
  geom_point(aes(color = Microsite), alpha=0.8, size=5, position = position_dodge(0.3)) +
  scale_color_manual(values = c("#9AC801", "#3552C5"))+
  facet_wrap(~Geolocation, nrow=1)+
  theme_bw()
dev.off()

pdf("ITS_Richness.pdf", width=8)
ggplot(ITS_rich, aes(Microsite, len)) +
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, color = Microsite),
    position = position_dodge(0.3), width = 0.3,size=0.8)+
  geom_point(aes(color = Microsite), alpha=0.8, size=5, position = position_dodge(0.3)) +
  scale_color_manual(values = c("#9AC801", "#3552C5"))+
  facet_wrap(~Geolocation, nrow=1)+
  theme_bw()
dev.off()


#Figure 4, FUNGuild Results
#Import table with Funguild Propotion Data
funguild<-read.table(file.choose() , sep='\t', header=T, row.names=1, check.names=F)
dim(funguild)
colnames(funguild)

funguild$Geolocation <- as.factor(funguild$Geolocation)
mycorrhizal<-funguild %>%
  group_by(Geolocation, Microsite) %>%
  summarise(
    sd = sd(Arbuscular_Mycorrhizal, na.rm = TRUE),
    AM=mean(Arbuscular_Mycorrhizal)
  )
mycorrhizal

pathogen<-funguild %>%
  group_by(Geolocation, Microsite) %>%
  summarise(
    sd = sd(Plant_Pathogen, na.rm = TRUE),
  pathogen=mean(Plant_Pathogen)
  )
pathogen

Wood_Sap<-funguild %>%
  group_by(Geolocation, Microsite) %>%
  summarise(
    sd = sd(Wood_Saprotroph, na.rm = TRUE),
    Wood_Sap=mean(Wood_Saprotroph)
  )
Wood_Sap

pdf("Arbuscular_Fungi.pdf",width=8)
ggplot(mycorrhizal, aes(Microsite, AM)) +
  geom_errorbar(
    aes(ymin = AM-sd, ymax = AM+sd, color = Microsite),
    position = position_dodge(0.3), width = 0.3,size=0.8)+
  geom_point(aes(color = Microsite), alpha=0.8, size=5, position = position_dodge(0.3)) +
  scale_color_manual(values = c("#9AC801", "#3552C5"))+
  facet_wrap(~Geolocation, nrow=1)+
  theme_bw()
dev.off()

pdf("Pathogen_Fungi.pdf",width=8)
ggplot(pathogen, aes(Microsite, pathogen)) +
  geom_errorbar(
    aes(ymin = pathogen-sd, ymax = pathogen+sd, color = Microsite),
    position = position_dodge(0.3), width = 0.3,size=0.8)+
  geom_point(aes(color = Microsite), alpha=0.8, size=5, position = position_dodge(0.3)) +
  scale_color_manual(values = c("#9AC801", "#3552C5"))+
  facet_wrap(~Geolocation, nrow=1)+
  theme_bw()
dev.off()


pdf("Wood_Saprotroph_Fungi.pdf",width=8)
ggplot(Wood_Sap, aes(Microsite, Wood_Sap)) +
  geom_errorbar(
    aes(ymin = Wood_Sap-sd, ymax = Wood_Sap+sd, color = Microsite),
    position = position_dodge(0.3), width = 0.3,size=0.8)+
  geom_point(aes(color = Microsite), alpha=0.8, size=5, position = position_dodge(0.3)) +
  scale_color_manual(values = c("#9AC801", "#3552C5"))+
  facet_wrap(~Geolocation, nrow=1)+
  theme_bw()
dev.off()

#Figure 5, qPCR Gene abundance
bac.16S <- dat %>%
  group_by(Geolocation, Microsite) %>%
  summarise(
    sd = sd(log_copy_number_16S, na.rm = TRUE),
    len = mean(log_copy_number_16S)
  )
bac.16S

ureC <- dat %>%
  group_by(Geolocation, Microsite) %>%
  summarise(
    sd = sd(log_copy_number_ureC, na.rm = TRUE),
    len = mean(log_copy_number_ureC)
  )
ureC



pdf("16S_geneabundance.pdf", width=8)
a=ggplot(bac.16S, aes(Microsite, len)) +
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, color = Microsite),
    position = position_dodge(0.3), width = 0.3,size=0.8)+
  geom_point(aes(color = Microsite), alpha=0.8, size=5, position = position_dodge(0.3)) +
  scale_color_manual(values = c("#9AC801", "#3552C5"))+
  facet_wrap(~Geolocation, nrow=1)+
  theme_bw()
dev.off()


pdf("ureC_geneabundance.pdf", width=8)
c=ggplot(ureC, aes(Microsite, len)) +
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, color = Microsite),
    position = position_dodge(0.3), width = 0.3,size=0.8)+
  geom_point(aes(color = Microsite), alpha=0.8, size=5, position = position_dodge(0.3)) +
  scale_color_manual(values = c("#9AC801", "#3552C5"))+
  facet_wrap(~Geolocation, nrow=1)+
  theme_bw()
dev.off()

