#Read in BBS data

routes <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_species_20170712.csv")
weather <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")

#species used in Huang 2017 GCB
huang_species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\huang-2017-bbs-spp.txt", header = TRUE, sep = "\t")

#Remove routes weather runtype = 0
routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("stateroute", "latitude", "longitude")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

#for 1st spp- all speciestotal values by route/lat/lon with weather runtype = 1
spp1 <- subset(counts, aou == huang_species$ID[1], select = c("stateroute", "year", "aou", "speciestotal"))
spp1.routes <- merge(spp1, RT1.routes[,c("stateroute", "latitude","longitude")], by = "stateroute", all.x = TRUE)

#Average routes for 1 spp over five year time windows, calculate centroid
#install.packages("dplyr")
#install.packages("DBI")
#install.packages("lazyeval")
library("dplyr")

years <- seq(1966, 2016, by = 5)
results <- matrix(nrow = 10, ncol = 4)
for (i in 1:10) {
yr <- years[i]  
results[i,1] <- yr
bobwhite <- distinct(filter(spp1.routes, year >= yr & year <= yr + 4, aou == 2890))
bobwhite.mean <- bobwhite %>%
  group_by(stateroute, latitude, longitude) %>%
  summarize(mean_abund = mean(speciestotal, na.rm = TRUE))
results[i, 4] <- mean(bobwhite.mean$mean_abund)
results[i,2] <- sum(bobwhite.mean$latitude*bobwhite.mean$mean_abund, na.rm=TRUE)/sum(bobwhite.mean$mean_abund, na.rm=TRUE)
results[i,3] <- sum(bobwhite.mean$longitude*bobwhite.mean$mean_abund, na.rm=TRUE)/sum(bobwhite.mean$mean_abund, na.rm = TRUE)
}
bobwhite.centroid <- data.frame(year = results[,1], centroid.lat = results[,2], centroid.lon = results[,3], mean_abund = results[,4])

