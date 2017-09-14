###### Data cleaning ######

#install.packages("dplyr")
#install.packages("DBI")
#install.packages("lazyeval")
library("dplyr")

setwd("C:/Users/gdicecco/Documents/bbs-centroid/results")

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
routes.short <- subset(RT1.routes, year >= 1969, select = c("stateroute", "year", "latitude", "longitude"))
counts$stateroute <- counts$statenum*1000 + counts$route

#Group routes by 5 year windows and spatial windows, check for routes in every time/spatial window
scale <- 5
routes.short$time.window <- scale*round(routes.short$year/scale) - 1
routes.short$lat.window <- scale*floor(routes.short$latitude/scale) + scale/2
routes.short$lon.window <- scale*floor(routes.short$longitude/scale) + scale/2
routes.short$spatial.window <- paste(routes.short$lat.window, routes.short$lon.window, sep = "")

##Routes that have at least one year per five year window
route.windows <- routes.short %>%
  select(stateroute, time.window)

unique.rtes <- unique(route.windows)

route.coverage.time <- unique.rtes %>%
  group_by(stateroute) %>%
  summarize(total = n())

time.windows <- route.coverage.time %>%
  filter(total == 10)

route.windows.s <- routes.short %>%
  select(stateroute, time.window, spatial.window)

unique.rtes2 <- unique(route.windows.s)

##Spatial grids that have observed routes for at least one year per time window
route.coverage.spatial <- unique.rtes2 %>%
  group_by(spatial.window) %>%
  select(time.window)

unique.rtes3 <- unique(route.coverage.spatial)

windows <- unique.rtes3 %>%
  group_by(spatial.window) %>%
  summarize(total= n()) %>%
  filter(total == 10)

routes.subs <- filter(routes.short, stateroute %in% time.windows$stateroute | spatial.window %in% windows$spatial.window)
write.csv(routes.subs, "cleaned_bbs_routes.csv")

#Average routes for each spp over five year time windows, calculate centroid, record shifted distance, velocity, shift bearing, population change
##Pull relevant species/species totals from counts df
counts.subs <- counts %>%
  filter(aou %in% huang_species$ID) %>%
  filter(stateroute %in% routes.subs$stateroute) %>%
  select(year, aou, speciestotal, stateroute)
write.csv(counts.subs, "cleaned_bbs_counts.csv")
