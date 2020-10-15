##Install necessary packages
install.packages("raster")
install.packages("sp")
##load libraries
library(raster)
library(sp)
#Download data files for global Aridity and PET data on your computer from the weblink: 
#https://figshare.com/articles/Global_Aridity_Index_and_Potential_Evapotranspiration_ET0_Climate_Database_v2/7504448/3
#Extract all the files on your computer: ET_annual, ET_monthly, AI_annual
#You need .tif files to extract data

#Now import .tif file containing Global Aridity Index (AI) data
##Make sure to use the correct location of your file
ai_et_0<-raster(x="C:/Users/pree/Documents/R_worldClimData/global-ai_et0/ai_et0/ai_et0.tif")
#Multiply AI values by 0.0001 for actual data
ai_values<-ai_et_0*0.0001
#Plot global AI
plot(ai_values)

#Create a variable with longitude and latitude outline of AZ boundaries 
#latitude and longitude valuesspecific for AZ
az <- as(extent(-114.866667, -109,31.333333, 37), 'SpatialPolygons')
crs(az) <- "+proj=longlat +datum=WGS84 +no_defs"
#Crop AI values for Arizona limits
az.map.ai<-crop(ai_values, az)
#Plot AI map for AZ
plot(az.map.ai)
#write/save the sub-setted data as a GTiff file
writeRaster(az.map.ai, "AZ_Aridity.tiff", "GTiff")

#Import the AZ image from the saved file as a raster 
##if starting fresh else use az.map.ai
az.raster<-raster(x="C:/Users/pree/Documents/Arid_Soil'_Microbiome/AZ_Aridity.tif")
#Plot the az.raster
plot(az.raster)
#Plot the AI map for AZ with titles
plot(az.raster, main="Aridity Index across state of Arizona")

#Adding gps points to the map specific to sampling sites
long= c(-114.19702, -113.92156,-114.20385,-114.30834,-112.074036, -111.651299,-110.911789) 
lat= c(33.37893,33.72123,33.75074,34.05754, 33.448376, 35.198284, 32.253460)
df<-as.data.frame(cbind(long, lat))

##Plot the graph and save it as a pdf 
pdf("AZmap_AI_02.20.2020.pdf", width=8)
plot(az.raster, main="Aridity Index across state of Arizona",axes=FALSE,box=FALSE, legend=T)
#Add symbols with site locations on it
points(df, pch=17, cex=0.8, col="black")
#Add sites name
name2=c("1 Quartzsite", "2", "3", "4 Parker", "Phoenix", "Flagstaff", "Tucson")
text(df,name2,cex=0.5, pos=3.5, offset=0.4)
dev.off()






