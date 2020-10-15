#Install libraries
install.packages("lmerTest")
install.packages("performance")
install.packages("MuMIn")
install.packages("nlme")


# Load required libraries
library(lmerTest)
library(performance)
library(MuMIn)
library(nlme)

#Import the .txt file with soil properties and richness data
dat <- read.table(file.choose(), header = T, row.names = 1)
dat$Geolocation <- as.factor(dat$Geolocation)

# Linear mixed-effects model: Fixed effect: Microsite (canopy/gap)
#Random effect: Geolocation (site)

#Alpha diversity
anova(lmer(Bacterial_Richness ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(Bacterial_Shannon ~ Microsite + (1|Geolocation), data = dat))

anova(lmer(ITS_Richness ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(ITS_Shannon ~ Microsite + (1|Geolocation), data = dat))

anova(lmer(DNA_Biomass ~ Microsite + (1|Geolocation), data = dat))

#Gene abundance copy number
anova(lmer(log_copy_number_16S ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(log_copy_number_ureC ~ Microsite + (1|Geolocation), data = dat))

#Soil properties
anova(lmer(pH ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(EC ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(DOC ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(DIC ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(TOC ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(TC ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(TN ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(DN ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(Ammonium ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(Nitrate ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(BAP ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(Sand ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(Silt ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(Clay ~ Microsite + (1|Geolocation), data = dat))
anova(lmer(Urease ~ Microsite + (1|Geolocation), data = dat))

#Model R2
#Alpha diversity
r2(lmer(Bacterial_Richness ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(Bacterial_Shannon ~ Microsite + (1|Geolocation), data = dat))

r2(lmer(ITS_Richness ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(ITS_Shannon ~ Microsite + (1|Geolocation), data = dat))

r2(lmer(DNA_Biomass ~ Microsite + (1|Geolocation), data = dat))

#Gene abundance copy number
r2(lmer(log_copy_number_16S ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(log_copy_number_ureC ~ Microsite + (1|Geolocation), data = dat))

# Soil properties
r2(lmer(pH ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(EC ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(DOC ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(DIC ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(TOC ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(TC ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(TN ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(DN ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(Ammonium ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(Nitrate ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(BAP ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(Sand ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(Silt ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(Clay ~ Microsite + (1|Geolocation), data = dat))
r2(lmer(Urease~ Microsite + (1|Geolocation), data = dat))

r.squaredGLMM(model.lme)

#Alpha diversity
r.squaredGLMM(lme(Bacterial_Richness ~ Microsite, ~1|Geolocation, data = dat))
r.squaredGLMM(lme(Bacterial_Shannon ~ Microsite, ~1|Geolocation, data = dat))

r.squaredGLMM(lme(ITS_Richness ~ Microsite, ~1|Geolocation, data = dat))
r.squaredGLMM(lme(ITS_Shannon ~ Microsite, ~1|Geolocation, data = dat))

r.squaredGLMM(lme(DNA_Biomass ~ Microsite, ~1|Geolocation, data = dat))

#Gene abundance copy number
r.squaredGLMM(lme(log_copy_number_16S ~ Microsite, ~ 1|Geolocation, data = dat))
r.squaredGLMM(lme(log_copy_number_ureC ~ Microsite , ~1|Geolocation, data = dat))

#Soil properties
r.squaredGLMM(lme(TOC ~ Microsite, ~ 1|Geolocation, data = dat))
r.squaredGLMM(lme(TC ~ Microsite, ~1|Geolocation, data = dat))
r.squaredGLMM(lme(TN ~ Microsite, ~1|Geolocation, data = dat))
r.squaredGLMM(lme(DN ~ Microsite, ~1|Geolocation, data = dat))
r.squaredGLMM(lme(Ammonium ~ Microsite, ~1|Geolocation, data = dat))
r.squaredGLMM(lme(Nitrate ~ Microsite, ~1|Geolocation, data = dat))
r.squaredGLMM(lme(BAP ~ Microsite, ~1|Geolocation, data = dat))