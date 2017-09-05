#Read in BBS data

routes <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_species_20170712.csv")
weather <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")

#species used in Huang 2017 GCB
huang_species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\huang-2017-bbs-spp.txt", header = TRUE, sep = "\t")

#Remove routes weather runtype = 0
routes$stateroute <- routes$countrynum*100000 + routes$statenum*1000 + routes$route
weather$stateroute <- weather$countrynum*100000 + weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("stateroute", "latitude", "longitude")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$countrynum*100000 + counts$statenum*1000 + counts$route

#for 1st spp- all speciestotal values by route/lat/lon with weather runtype = 1
spp1 <- subset(counts, aou == huang_species$ID[1], select = c("stateroute", "year", "aou", "speciestotal"))
spp1.routes <- merge(spp1, RT1.routes[,c("stateroute", "latitude","longitude")], by = "stateroute", all.x = TRUE)

#next need to average routes over 5 year windows for time range used in Huang paper
