library(dplyr)
library(geosphere)
library(maps)
library(ggplot2)
library(RColorBrewer)
library(rgdal)
library(rgeos)

####### Data cleaning ######

#Read in BBS data
routes <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_species_20170712.csv")
weather <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")
bcrshp <- readOGR("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\bcr_terrestrial_shape\\BCR_Terrestrial_master.shp") #BCRs
##Typo in bcrshp - has TENNESSE

#species used in Huang 2017 GCB
huang_species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\huang-2017-bbs-species.csv", header = TRUE)

#Remove routes weather runtype = 0
routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
routes.short <- subset(RT1.routes, year >= 1969, select = c("statenum","stateroute", "year", "latitude", "longitude", "bcr"))
counts$stateroute <- counts$statenum*1000 + counts$route

#Plot centroid movement for each analysis
par(mfrow=c(1,1))
longs = c(-125,-60)
lats = c(26,50)

#calcuate distance ratio- determine if centroid moves in consistent direction over time
#input- centroids df with centroids per time window for all spp. (10 lat lon coordinates of centroid over 5 year time windows)
#will subset centroids df for spp of interest by AOU
distance.ratio <- function(df, AOU, randomizations = FALSE, n = 1000) {
  
  x = filter(df, aou == AOU)
  
  years = seq(1969, 2014, by = 5)
  
  if (randomizations) {
    rand.ratios = c()
    for (r in 1:n) {
      sampleyears = sample(years, 10)
      distances <- c()
      for(i in 1:9) {
        distances <- c(distances, distGeo(x[x$time.window == sampleyears[i], c("centroid_lon", "centroid_lat")], 
                                          x[x$time.window == sampleyears[i+1], c("centroid_lon", "centroid_lat")]))
        sumdist <- sum(distances)
      }
      fl <- distGeo(x[1,c("centroid_lon", "centroid_lat")], x[10,c("centroid_lon", "centroid_lat")])
      rand.ratios = c(rand.ratios, fl/sumdist)
    }
  } else {
    rand.ratios = NULL
  }
   
  sampleyears = years
  distances <- c()
  for(i in 1:9) {
    distances <- c(distances, distGeo(x[x$time.window == sampleyears[i], c("centroid_lon", "centroid_lat")], 
                                      x[x$time.window == sampleyears[i+1], c("centroid_lon", "centroid_lat")]))
    sumdist <- sum(distances)
  }
  fl <- distGeo(x[1,c("centroid_lon", "centroid_lat")], x[10,c("centroid_lon", "centroid_lat")])
  obs.ratio = fl/sumdist
  
  return(list(obs = obs.ratio, rand = rand.ratios))
}    
    
####### Centroids by routes ############
#Subset counts df for relevant spp and time period
scale5 <- 5
routes.short$time.window <- scale5*round(routes.short$year/scale5) - 1
counts.short <- counts %>%
  filter(year >= 1969) %>%
  filter(aou %in% huang_species$aou) %>%
  select(year, aou, speciestotal, stateroute)
counts.short.merged <- merge(routes.short, counts.short, by = c("stateroute", "year"))

#mean species abundance for each route over time.windows
route_spp_means <- counts.short.merged %>%
  group_by(stateroute, latitude, longitude, aou, time.window) %>%
  summarize(avg_abund = mean(speciestotal, na.rm = TRUE))

#mean centroid for each species in five year time windows
centroids <- route_spp_means %>%
  group_by(aou, time.window) %>%
  summarize(centroid_lat = sum(latitude*avg_abund, na.rm = TRUE)/sum(avg_abund, na.rm = TRUE), centroid_lon = sum(longitude*avg_abund, na.rm = TRUE)/sum(avg_abund, na.rm = TRUE),
            mean_total_abund = mean(avg_abund))
#setwd("C:/Users/gdicecco/Documents/bbs-centroid/centroids-by-routes/")
#write.csv(centroids, "centroids-by-routes.csv", row.names=F)

#Plot centroid movement
#map(database="world",xlim = longs, ylim = lats) #plot just on states
#map(database = "state", add = TRUE)
plot(bcrshp[bcrshp$WATER == 3,], ylim = lats, xlim = longs, border = "gray73", col = "gray95") #plot centroids on BCR map
mtext("Centroids by route",3,cex=2,line=.5)

#shifted distance, velocity, bearing of shift, population change, shift direction regression
results <- matrix(nrow = 35, ncol = 12)
for(i in 1:35) {
  species <- huang_species$aou[i]
  results[i,1] <- species
  df <- centroids %>%
    filter(aou == species)
  
  results[i,2] <- distGeo(df[1,4:3], df[10,4:3])/1000
  results[i,3] <- results[i,2]/(2016-1969)
  results[i,4] <- bearing(df[1,4:3], df[10,4:3])
  results[i,5] <- log(mean(df$mean_total_abund[9:10])/mean(df$mean_total_abund[1:2]))
  distratios <- distance.ratio(centroids, species, T)
  results[i,10] <- distratios$obs
  results[i,11] <- mean(distratios$rand)
  results[i,12] <- sd(distratios$rand)
  
  mod.test <- lm(df$centroid_lat ~ df$time.window)
  sum <- summary(mod.test)
  results[i,7] <- sum$coefficients[2,4]
  results[i,6] <- sum$coefficients[2,1]
  
  mod.test.2 <- lm(df$centroid_lon ~ df$time.window)
  sum.2 <- summary(mod.test.2)
  results[i,9] <- sum.2$coefficients[2,4]
  results[i,8] <- sum.2$coefficients[2,1]
  
  points(df$centroid_lon[1], df$centroid_lat[1], pch = 16, col = 'red')
  points(df$centroid_lon, df$centroid_lat, type = 'l')
  points(df$centroid_lon[10], df$centroid_lat[10], pch = 17, col = 'blue')
  text(df$centroid_lon[10], df$centroid_lat[10]+1, species, cex = .75)
}
results.df <- data.frame(aou = results[,1], 
                               shiftdist = results[,2], 
                               velocity = results[,3], 
                               bearing = results[,4], 
                               r = results[,5],
                               lat_slope = results[,6], 
                               lat_pval = results[,7], 
                               lon_slope = results[,8], 
                               lon_pval = results[,9],
                         distance_ratio = results[, 10],
                         dratio_rand_mean = results[, 11],
                         dratio_rand_sd = results[,12])

#assign shift directions
direction.lat <- c()
for(i in 1:35){
  slope <- results.df$lat_slope[i]
  pval <- results.df$lat_pval[i]
  if (slope > 0 & pval < 0.05) {
    direction.lat <- c(direction.lat, "north")
  } else if (slope < 0 & pval < 0.05) {
    direction.lat <- c(direction.lat, "south")
  } else 
    direction.lat <- c(direction.lat, "")
}
results.df <- cbind(results.df, direction.lat)

direction.lon <- c()
for(i in 1:35){
  slope <- results.df$lon_slope[i]
  pval <- results.df$lon_pval[i]
  if (slope > 0 & pval < 0.05) {
    direction.lon <- c(direction.lon, "east")
  } else if (slope < 0 & pval < 0.05) {
    direction.lon <- c(direction.lon, "west")
  } else 
    direction.lon <- c(direction.lon, "")
}
results.df <- cbind(results.df,direction.lon)
results.df$shiftdir <- paste(results.df$direction.lat,results.df$direction.lon, sep = "")

popchange <- c()
for(i in 1:35){
  r <- results.df$r[i]
  if (r > 0) {
    popchange <- c(popchange, "increasing")
  } else popchange <- c(popchange, "decreasing")
}
results.df <- cbind(results.df, popchange)
#setwd("C:/Users/gdicecco/Documents/bbs-centroid/centroids-by-routes/")
#write.csv(results.df, "centroids-by-routes-results.csv", row.names=F)
final.df <- data.frame(species = huang_species$species, 
                             aou= huang_species$aou, 
                             shiftdistance = results.df$shiftdist, 
                             velocity = results.df$velocity, 
                             direction = results.df$shiftdir, 
                             abundance = results.df$popchange,
                       distanceratio = results.df$distance_ratio)
#write.csv(final.df, "results-routes-summarized.csv", row.names=F)

####### Centroids by 1 degree grids ######
#Group routes by 5 year windows and spatial windows, check for routes in every time/spatial window
scale1 <- 1
routes.short$lat.window <- scale1*floor(routes.short$latitude/scale1) + scale1/2
routes.short$lon.window <- scale1*floor(routes.short$longitude/scale1) + scale1/2
routes.short$spatial.window <- paste(routes.short$lat.window, routes.short$lon.window, sep = "")

##Spatial grids that have observed routes for at least one year per time window
route.windows <- routes.short %>%
  select(stateroute, time.window, spatial.window) %>%
  unique() %>%
  group_by(spatial.window) %>%
  select(time.window) %>%
  unique() %>%
  group_by(spatial.window) %>%
  summarize(total = n()) %>%
  filter(total == 10)

routes.subs <- filter(routes.short, spatial.window %in% route.windows$spatial.window)

#Pull relevant species/species totals from counts df
counts.subs <- counts %>%
  filter(year >= 1969) %>%
  filter(aou %in% huang_species$aou) %>%
  filter(stateroute %in% routes.subs$stateroute) %>%
  select(year, aou, speciestotal, stateroute)

counts.merged <- merge(routes.subs, counts.subs, by = c("stateroute", "year"))

#mean species abundance for each route over time.windows
spp_abund_means <- counts.merged %>%
  group_by(stateroute, latitude, longitude, aou, time.window) %>%
  summarize(avg_abund = mean(speciestotal, na.rm = TRUE))

#mean centroid for each species in five year time windows
centroids.grids <- spp_abund_means %>%
  group_by(aou, time.window) %>%
  summarize(centroid_lat = sum(latitude*avg_abund, na.rm = TRUE)/sum(avg_abund, na.rm = TRUE), centroid_lon = sum(longitude*avg_abund, na.rm = TRUE)/sum(avg_abund, na.rm = TRUE),
            mean_total_abund = mean(avg_abund))
#setwd("C:/Users/gdicecco/Documents/bbs-centroid/centroids-by-grids/")
#write.csv(centroids.grids, "centroids-by-grids.csv", row.names=F)

#Plot centroid movement
#map(database="world",xlim = longs, ylim = lats) #plot just on states
#map(database = "state", add = TRUE)
plot(bcrshp[bcrshp$WATER == 3,], ylim = lats, xlim = longs, border = "gray73", col = "gray95") #plot centroids on BCR map
mtext("Centroids by 1 deg grid",3,cex=2,line=.5)

#shifted distance, velocity, bearing of shift, population change, shift direction regression
results.grids <- matrix(nrow = 35, ncol = 12)
for(i in 1:35) {
  species <- huang_species$aou[i]
  results.grids[i,1] <- species
  df <- centroids.grids %>%
    filter(aou == species)
  
  results.grids[i,2] <- distGeo(df[1,4:3], df[10,4:3])/1000
  results.grids[i,3] <- results.grids[i,2]/(2016-1969)
  results.grids[i,4] <- bearing(df[1,4:3], df[10,4:3])
  results.grids[i,5] <- log(mean(df$mean_total_abund[9:10])/mean(df$mean_total_abund[1:2]))
  distratios <- distance.ratio(centroids.grids, species, T)
  results.grids[i,10] <- distratios$obs
  results.grids[i,11] <- mean(distratios$rand)
  results.grids[i,12] <- sd(distratios$rand)
  
  
  mod.test <- lm(df$centroid_lat ~ df$time.window)
  sum <- summary(mod.test)
  results.grids[i,7] <- sum$coefficients[2,4]
  results.grids[i,6] <- sum$coefficients[2,1]
  
  mod.test.2 <- lm(df$centroid_lon ~ df$time.window)
  sum.2 <- summary(mod.test.2)
  results.grids[i,9] <- sum.2$coefficients[2,4]
  results.grids[i,8] <- sum.2$coefficients[2,1]
  
  points(df$centroid_lon[1], df$centroid_lat[1], pch = 16, col = 'red')
  points(df$centroid_lon, df$centroid_lat, type = 'l')
  points(df$centroid_lon[10], df$centroid_lat[10], pch = 17, col = 'blue')
  text(df$centroid_lon[10], df$centroid_lat[10]+1, species, cex = .75)
}
results.grids.df <- data.frame(aou = results.grids[,1], 
                         shiftdist = results.grids[,2], 
                         velocity = results.grids[,3], 
                         bearing = results.grids[,4], 
                         r = results.grids[,5],
                         lat_slope = results.grids[,6], 
                         lat_pval = results.grids[,7], 
                         lon_slope = results.grids[,8], 
                         lon_pval = results.grids[,9],
                         distance_ratio = results.grids[,10],
                         dratio_rand_mean = results.grids[,11],
                         dratio_rand_sd = results.grids[,12])

#assign shift directions
direction.lat <- c()
for(i in 1:35){
  slope <- results.grids.df$lat_slope[i]
  pval <- results.grids.df$lat_pval[i]
  if (slope > 0 & pval < 0.05) {
    direction.lat <- c(direction.lat, "north")
  } else if (slope < 0 & pval < 0.05) {
    direction.lat <- c(direction.lat, "south")
  } else 
    direction.lat <- c(direction.lat, "")
}
results.grids.df <- cbind(results.grids.df, direction.lat)

direction.lon <- c()
for(i in 1:35){
  slope <- results.grids.df$lon_slope[i]
  pval <- results.grids.df$lon_pval[i]
  if (slope > 0 & pval < 0.05) {
    direction.lon <- c(direction.lon, "east")
  } else if (slope < 0 & pval < 0.05) {
    direction.lon <- c(direction.lon, "west")
  } else 
    direction.lon <- c(direction.lon, "")
}
results.grids.df <- cbind(results.grids.df,direction.lon)
results.grids.df$shiftdir <- paste(results.grids.df$direction.lat,results.grids.df$direction.lon, sep = "")

popchange <- c()
for(i in 1:35){
  r <- results.grids.df$r[i]
  if (r > 0) {
    popchange <- c(popchange, "increasing")
  } else popchange <- c(popchange, "decreasing")
}
results.grids.df <- cbind(results.grids.df, popchange)
#setwd("C:/Users/gdicecco/Documents/bbs-centroid/centroids-by-grids/")
#write.csv(results.grids.df, "centroids-by-grids-results.csv", row.names=F)
final.grids.df <- data.frame(species = huang_species$species, 
                           aou= huang_species$aou, 
                           shiftdistance = results.grids.df$shiftdist, 
                           velocity = results.grids.df$velocity, 
                           direction = results.grids.df$shiftdir, 
                           abundance = results.grids.df$popchange,
                           distanceratio = results.grids.df$distance_ratio)
#write.csv(final.grids.df, "results-grids-summarized.csv", row.names=F)


####### Centroids by strata ########

#Calculate BCR/state strata centroids
regioncodes <- read.table("C:/Users/gdicecco/Desktop/regioncodes.txt", sep = "\t")
colnames(regioncodes) <- c("countrynum", "statenum", "PROVINCE_S")

polys.df <- data.frame(BCR = c(0), BCRNAME = c(0), PROVINCE_S = c(0), COUNTRY = c(0), REGION = c(0),
                       WATER = c(0), Shape_Leng = c(0), Shape_Area = c(0), x = c(0), y = c(0))
for(state in regioncodes$PROVINCE_S) {
  shape <- subset(bcrshp, PROVINCE_S == state)
  centers <- gCentroid(shape, byid = TRUE)
  df.temp <- data.frame(cbind(shape@data, centers@coords))
  polys.df <- rbind(polys.df, df.temp)
}
polys.df <- polys.df[-1,] #remove zero row
polys.merged <- merge(regioncodes, polys.df, by = "PROVINCE_S", all.x = TRUE) #merge PROVINCE_S with statenum
polys.merged.land <- subset(polys.merged, WATER == 3) #remove centroids for bodies of water

#merge strata centroids with routes
routes.short$statebcr <- routes.short$statenum*1000 + routes.short$bcr
polys.merged.land$statebcr <- polys.merged.land$statenum*1000 + polys.merged.land$BCR
routes.short.centers <- merge(routes.short, polys.merged.land, by = "statebcr")
counts.merged.centers <- merge(routes.short.centers, counts.short, by = c("stateroute","year"))

subs.routes4 <- counts.merged.centers %>%
  group_by(aou, statebcr) %>%
  distinct(stateroute) %>%
  summarize(total = n()) %>%
  filter(total > 4)

counts.merged.subs.strata <- merge(counts.merged.centers, subs.routes4, by = c("aou","statebcr"))

#mean centroid for each species in five year time windows per bcr/state strata
spp_abund_means_bcr <- counts.merged.subs.strata %>%
  group_by(statebcr, x, y, stateroute, bcr, latitude, longitude, aou, time.window) %>%
  summarize(avg_abund = mean(speciestotal, na.rm = TRUE))

centroids.bcr <- spp_abund_means_bcr %>%
  group_by(aou, time.window, statebcr) %>%
  summarize(centroid_lat = sum(y*avg_abund, na.rm = TRUE)/sum(avg_abund, na.rm = TRUE), centroid_lon = sum(x*avg_abund, na.rm = TRUE)/sum(avg_abund, na.rm = TRUE),
            mean_total_abund = mean(avg_abund))

#mean centroid weighted from bcr centroids
centroids.bcr2 <- centroids.bcr %>%
  group_by(aou, time.window) %>%
  summarize(centroid_lat = sum(centroid_lat*mean_total_abund, na.rm = TRUE)/sum(mean_total_abund, na.rm = TRUE), centroid_lon = sum(centroid_lon*mean_total_abund, na.rm = TRUE)/sum(mean_total_abund, na.rm = TRUE),
            mean_total_abund = mean(mean_total_abund))
#setwd("C:/Users/gdicecco/Desktop/git/bbs-centroid/centroids-by-strata/")
#write.csv(centroids.bcr2, "centroids-by-strata.csv", row.names=F)

#Plot centroid movement
#map(database="world",xlim = longs, ylim = lats) #plot just on states
#map(database = "state", add = TRUE)
plot(bcrshp[bcrshp$WATER == 3,], ylim = lats, xlim = longs, border = "gray73", col = "gray95") #plot centroids on BCR map
mtext("Centroids by strata",3,cex=2,line=.5)

#shifted distance, velocity, bearing of shift, population change, shift direction regression
results.bcr <- matrix(nrow = 35, ncol = 12)
for(i in 1:35) {
  species <- huang_species$aou[i]
  results.bcr[i,1] <- species
  df <- centroids.bcr2 %>%
    filter(aou == species)
  
  results.bcr[i,2] <- distGeo(df[1,4:3], df[10,4:3])/1000
  results.bcr[i,3] <- results.bcr[i,2]/(2016-1969)
  results.bcr[i,4] <- bearing(df[1,4:3], df[10,4:3])
  results.bcr[i,5] <- log(mean(df$mean_total_abund[9:10])/mean(df$mean_total_abund[1:2]))
  distratios <- distance.ratio(centroids.bcr2, species, T)
  results.bcr[i,10] <- distratios$obs
  results.bcr[i,11] <- mean(distratios$rand)
  results.bcr[i,12] <- sd(distratios$rand)
  
  mod.test <- lm(df$centroid_lat ~ df$time.window)
  sum <- summary(mod.test)
  results.bcr[i,7] <- sum$coefficients[2,4]
  results.bcr[i,6] <- sum$coefficients[2,1]
  
  mod.test.2 <- lm(df$centroid_lon ~ df$time.window)
  sum.2 <- summary(mod.test.2)
  results.bcr[i,9] <- sum.2$coefficients[2,4]
  results.bcr[i,8] <- sum.2$coefficients[2,1]
  
  points(df$centroid_lon[1], df$centroid_lat[1], pch = 16, col = 'red')
  points(df$centroid_lon, df$centroid_lat, type = 'l')
  points(df$centroid_lon[10], df$centroid_lat[10], pch = 17, col = 'blue')
  text(df$centroid_lon[10], df$centroid_lat[10]+1, species, cex = .75)
  
  
}
results.bcr.df <- data.frame(aou = results.bcr[,1], 
                         shiftdist = results.bcr[,2], 
                         velocity = results.bcr[,3], 
                         bearing = results.bcr[,4], 
                         r = results.bcr[,5],
                         lat_slope = results.bcr[,6], 
                         lat_pval = results.bcr[,7], 
                         lon_slope = results.bcr[,8], 
                         lon_pval = results.bcr[,9],
                         distance_ratio = results.bcr[,10],
                         dratio_rand_mean = results.bcr[,11],
                         dratio_rand_sd = results.bcr[,12])

#assign shift directions
direction.lat <- c()
for(i in 1:35){
  slope <- results.bcr.df$lat_slope[i]
  pval <- results.bcr.df$lat_pval[i]
  if (slope > 0 & pval < 0.05) {
    direction.lat <- c(direction.lat, "north")
  } else if (slope < 0 & pval < 0.05) {
    direction.lat <- c(direction.lat, "south")
  } else 
    direction.lat <- c(direction.lat, "")
}
results.bcr.df <- cbind(results.bcr.df, direction.lat)

direction.lon <- c()
for(i in 1:35){
  slope <- results.bcr.df$lon_slope[i]
  pval <- results.bcr.df$lon_pval[i]
  if (slope > 0 & pval < 0.05) {
    direction.lon <- c(direction.lon, "east")
  } else if (slope < 0 & pval < 0.05) {
    direction.lon <- c(direction.lon, "west")
  } else 
    direction.lon <- c(direction.lon, "")
}
results.bcr.df <- cbind(results.bcr.df,direction.lon)
results.bcr.df$shiftdir <- paste(results.bcr.df$direction.lat,results.bcr.df$direction.lon, sep = "")

popchange <- c()
for(i in 1:35){
  r <- results.bcr.df$r[i]
  if (r > 0) {
    popchange <- c(popchange, "increasing")
  } else popchange <- c(popchange, "decreasing")
}
results.bcr.df <- cbind(results.bcr.df, popchange)

#setwd("C:/Users/gdicecco/Desktop/git/bbs-centroid/centroids-by-strata/")
#write.csv(results.bcr.df, "centroids-by-strata-results.csv", row.names=F)

final.bcr.df <- data.frame(species = huang_species$species, 
                       aou= huang_species$aou, 
                       shiftdistance = results.bcr.df$shiftdist, 
                       velocity = results.bcr.df$velocity, 
                       direction = results.bcr.df$shiftdir, 
                       abundance = results.bcr.df$popchange,
                       distanceratio = results.bcr.df$distance_ratio)
#write.csv(final.bcr.df, "results-strata-summarized.csv", row.names=F)

####### Comparison figures #######
#centroids by route
final.df$analysis <- rep(x = "By route", times = 35)
#centroids by grid
final.grids.df$analysis <- rep(x = "By grid", times = 35)
#centroids by bcr
final.bcr.df$analysis <- rep(x = "By strata", times = 35)
#huang data
huang_species$distanceratio <- rep(NA)
huang_species$analysis <- rep(x = "Huang et al.", times = 35)

sppnumbers <- c(1:15, 17:36)

compiled.df <- rbind(final.df, final.grids.df, final.bcr.df, huang_species)
compiled.df$sppnumber <- rep(x = sppnumbers, times = 4)
compiled.df$aou <- as.factor(compiled.df$aou)
 
blank <- theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#compare shift distance
##boxplot
ggplot(compiled.df, aes(analysis, shiftdistance)) + theme_bw() + geom_boxplot(fill = "gray") + labs(x = "Analysis", y = "Shift Distance") + blank

##multipanel regression
par(bg = "white")
par(new = T)
plot.new()
par(oma = c(0,0,4,0))
split.screen(c(3,3))
par(mar = c(5.1,4.1,4.1,2.1))

screen(1)
par(mar = c(3.1,4.1,3.1,2.1))
plot(compiled.df$shiftdistance[compiled.df$analysis == "By route"], compiled.df$shiftdistance[compiled.df$analysis == "By grid"], xlab = "", ylab = "By grid", type="n", ylim = c(0,750), xlim = c(0,750))
text(compiled.df$shiftdistance[compiled.df$analysis == "By route"], compiled.df$shiftdistance[compiled.df$analysis == "By grid"], label = compiled.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.df$shiftdistance[compiled.df$analysis == "By grid"] ~ compiled.df$shiftdistance[compiled.df$analysis == "By route"]), col = "blue")

screen(4)
par(mar = c(3.1,4.1,2.1,2.1))
plot(compiled.df$shiftdistance[compiled.df$analysis == "By route"], compiled.df$shiftdistance[compiled.df$analysis == "By strata"], xlab = "", ylab = "By strata", type = "n", ylim = c(0,750), xlim = c(0,750))
text(compiled.df$shiftdistance[compiled.df$analysis == "By route"], compiled.df$shiftdistance[compiled.df$analysis == "By strata"], label = compiled.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.df$shiftdistance[compiled.df$analysis == "By strata"] ~ compiled.df$shiftdistance[compiled.df$analysis == "By route"]), col = "blue")

screen(5)
par(mar = c(3.1,4.1,2.1,2.1))
plot(compiled.df$shiftdistance[compiled.df$analysis == "By grid"], compiled.df$shiftdistance[compiled.df$analysis == "By strata"],xlab = "", ylab = "", type = "n", ylim = c(0,750), xlim = c(0,750))
text(compiled.df$shiftdistance[compiled.df$analysis == "By grid"], compiled.df$shiftdistance[compiled.df$analysis == "By strata"], label = compiled.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.df$shiftdistance[compiled.df$analysis == "By strata"] ~ compiled.df$shiftdistance[compiled.df$analysis == "By grid"]), col = "blue")

screen(7)
par(mar = c(5.1,4.1,2.1,2.1))
plot(compiled.df$shiftdistance[compiled.df$analysis == "By route"], compiled.df$shiftdistance[compiled.df$analysis == "Huang et al."],xlab = "By route", ylab = "Huang et al.", type = "n", ylim = c(0,750), xlim = c(0,750))
text(compiled.df$shiftdistance[compiled.df$analysis == "By route"], compiled.df$shiftdistance[compiled.df$analysis == "Huang et al."], label = compiled.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.df$shiftdistance[compiled.df$analysis == "Huang et al."] ~ compiled.df$shiftdistance[compiled.df$analysis == "By route"]), col = "blue")

screen(8)
par(mar = c(5.1,4.1,2.1,2.1))
plot(compiled.df$shiftdistance[compiled.df$analysis == "By grid"], compiled.df$shiftdistance[compiled.df$analysis == "Huang et al."],xlab = "By grid", ylab = "", type = "n", ylim = c(0,750), xlim = c(0,750))
text(compiled.df$shiftdistance[compiled.df$analysis == "By grid"], compiled.df$shiftdistance[compiled.df$analysis == "Huang et al."], label = compiled.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.df$shiftdistance[compiled.df$analysis == "Huang et al."] ~ compiled.df$shiftdistance[compiled.df$analysis == "By grid"]), col = "blue")

screen(9)
par(mar = c(5.1,4.1,2.1,2.1))
plot(compiled.df$shiftdistance[compiled.df$analysis == "By strata"], compiled.df$shiftdistance[compiled.df$analysis == "Huang et al."],xlab = "By strata", ylab = "", type = "n", ylim = c(0,750), xlim = c(0,750))
text(compiled.df$shiftdistance[compiled.df$analysis == "By strata"], compiled.df$shiftdistance[compiled.df$analysis == "Huang et al."], label = compiled.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.df$shiftdistance[compiled.df$analysis == "Huang et al."] ~ compiled.df$shiftdistance[compiled.df$analysis == "By strata"]), col = "blue")
close.screen(all = T)
title(main = "Shift distance (km)", outer = T)


#compare shift direction
results.df$analysis <- rep(x = "By route", times = 35)
#centroids by grid
results.grids.df$analysis <- rep(x = "By grid", times = 35)
#centroids by bcr
results.bcr.df$analysis <- rep(x = "By strata", times = 35)

compiled.results.df <- rbind(results.df, results.grids.df, results.bcr.df)
compiled.results.df$sppnumber <- rep(x = sppnumbers, times = 3)
compiled.results.df$aou <- as.factor(compiled.results.df$aou)

#multipanel shift direction regression plot
par(bg = "white")
par(new = T)
plot.new()
par(oma = c(0,0,2,0))
split.screen(c(2,2))
par(mar = c(5.1,4.1,4.1,2.1))

screen(1)
par(mar = c(3.1,4.1,3.1,2.1))
plot(compiled.results.df$bearing[compiled.results.df$analysis == "By route"], compiled.results.df$bearing[compiled.results.df$analysis == "By grid"], xlab = "", ylab = "By grid", type="n", ylim = c(-180,180), xlim = c(-180,180))
text(compiled.results.df$bearing[compiled.results.df$analysis == "By route"], compiled.results.df$bearing[compiled.results.df$analysis == "By grid"], label = compiled.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.results.df$bearing[compiled.results.df$analysis == "By grid"] ~ compiled.results.df$bearing[compiled.results.df$analysis == "By route"]), col = "blue")

screen(3)
par(mar = c(5.1,4.1,2.1,2.1))
plot(compiled.results.df$bearing[compiled.results.df$analysis == "By route"], compiled.results.df$bearing[compiled.results.df$analysis == "By strata"], xlab = "By route", ylab = "By strata", type = "n", ylim = c(-180,180), xlim = c(-180,180))
text(compiled.results.df$bearing[compiled.results.df$analysis == "By route"], compiled.results.df$bearing[compiled.results.df$analysis == "By strata"], label = compiled.results.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.results.df$bearing[compiled.results.df$analysis == "By strata"] ~ compiled.results.df$bearing[compiled.results.df$analysis == "By route"]), col = "blue")

screen(4)
par(mar = c(5.1,4.1,2.1,2.1))
plot(compiled.results.df$bearing[compiled.results.df$analysis == "By grid"], compiled.results.df$bearing[compiled.results.df$analysis == "By strata"],xlab = "By grid", ylab = "", type = "n",ylim = c(-180,180), xlim = c(-180,180))
text(compiled.results.df$bearing[compiled.results.df$analysis == "By grid"], compiled.results.df$bearing[compiled.results.df$analysis == "By strata"], label = compiled.results.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.results.df$bearing[compiled.results.df$analysis == "By strata"] ~ compiled.results.df$bearing[compiled.results.df$analysis == "By grid"]), col = "blue")

close.screen(all = T)
title(main = "Bearing (degrees)", outer = T)

#shift direction multipanel
scale45 <- 45
compiled.results.df$roundbearing <- scale45*round(compiled.results.df$bearing/scale45)
bearinghuang <- c()
for(x in compiled.df$direction[compiled.df$analysis == "Huang et al."]) {
  if(x == "north") {
    bearinghuang <- c(bearinghuang, 0)
  } else if (x == "northeast") {
    bearinghuang <- c(bearinghuang, 45)
  } else if (x == "east") {
    bearinghuang <- c(bearinghuang,90)
  } else if (x == "southeast") {
    bearinghuang <- c(bearinghuang, 135)
  } else if (x == "southwest") {
    bearinghuang <- c(bearinghuang, -135)
  } else if (x == "south") {
    bearinghuang <- c(bearinghuang, -180)
  } else if (x == "southwest") {
    bearinghuang <- c(bearinghuang, -135)
  } else if (x == "west") {
    bearinghuang <- c(bearinghuang, -90)
  } else if (x == "northwest") {
    bearinghuang <- c(bearinghuang, -45)
  }
}

#plot code
par(bg = "white")
par(new = T)
plot.new()
par(oma = c(0,0,4,0))
split.screen(c(3,3))
par(mar = c(5.1,4.1,4.1,2.1))

screen(1)
par(mar = c(3.1,4.1,3.1,2.1))
plot(compiled.results.df$roundbearing[compiled.results.df$analysis == "By route"], compiled.results.df$roundbearing[compiled.results.df$analysis == "By grid"], xlab = "", ylab = "By grid", type="n", ylim = c(-180,180), xlim = c(-180,180))
text(compiled.results.df$roundbearing[compiled.results.df$analysis == "By route"], compiled.results.df$roundbearing[compiled.results.df$analysis == "By grid"], label = compiled.results.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.results.df$roundbearing[compiled.results.df$analysis == "By grid"] ~ compiled.results.df$roundbearing[compiled.results.df$analysis == "By route"]), col = "blue")

#by route, by strata
screen(4)
par(mar = c(3.1,4.1,2.1,2.1))
plot(compiled.results.df$roundbearing[compiled.results.df$analysis == "By route"], compiled.results.df$roundbearing[compiled.results.df$analysis == "By strata"], xlab = "", ylab = "By strata", type="n", ylim = c(-180,180), xlim = c(-180,180))
text(compiled.results.df$roundbearing[compiled.results.df$analysis == "By route"], compiled.results.df$roundbearing[compiled.results.df$analysis == "By strata"], label = compiled.results.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.results.df$roundbearing[compiled.results.df$analysis == "By strata"] ~ compiled.results.df$roundbearing[compiled.results.df$analysis == "By route"]), col = "blue")

#by grid, by strata
screen(5)
par(mar = c(3.1,4.1,2.1,2.1))
plot(compiled.results.df$roundbearing[compiled.results.df$analysis == "By grid"], compiled.results.df$roundbearing[compiled.results.df$analysis == "By strata"], xlab = "", ylab = "", type="n", ylim = c(-180,180), xlim = c(-180,180))
text(compiled.results.df$roundbearing[compiled.results.df$analysis == "By grid"], compiled.results.df$roundbearing[compiled.results.df$analysis == "By strata"], label = compiled.results.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.results.df$roundbearing[compiled.results.df$analysis == "By strata"] ~ compiled.results.df$roundbearing[compiled.results.df$analysis == "By grid"]), col = "blue")

#by route, huang et al.
screen(7)
par(mar = c(5.1,4.1,2.1,2.1))
plot(compiled.results.df$roundbearing[compiled.results.df$analysis == "By route"], bearinghuang, xlab = "By route", ylab = "Huang et al.", type="n", ylim = c(-180,180), xlim = c(-180,180))
text(compiled.results.df$roundbearing[compiled.results.df$analysis == "By route"], bearinghuang, label = compiled.results.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(bearinghuang ~ compiled.results.df$roundbearing[compiled.results.df$analysis == "By route"]), col = "blue")

#by grid, huang et al.
screen(8)
par(mar = c(5.1,4.1,2.1,2.1))
plot(compiled.results.df$roundbearing[compiled.results.df$analysis == "By grid"], bearinghuang, xlab = "By grid", ylab = "", type="n", ylim = c(-180,180), xlim = c(-180,180))
text(compiled.results.df$roundbearing[compiled.results.df$analysis == "By grid"], bearinghuang, label = compiled.results.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(bearinghuang ~ compiled.results.df$roundbearing[compiled.results.df$analysis == "By grid"]), col = "blue")

#by strata, huang et al.
screen(9)
par(mar = c(5.1,4.1,2.1,2.1))
plot(compiled.results.df$roundbearing[compiled.results.df$analysis == "By strata"], bearinghuang, xlab = "By strata", ylab = "", type="n", ylim = c(-180,180), xlim = c(-180,180))
text(compiled.results.df$roundbearing[compiled.results.df$analysis == "By strata"], bearinghuang, label = compiled.results.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(bearinghuang ~ compiled.results.df$roundbearing[compiled.results.df$analysis == "By route"]), col = "blue")

close.screen(all = T)
title(main = "Shift direction", outer = T)


#compare pop change
ggplot(compiled.df, aes(x = abundance, fill = analysis))+stat_count()+theme_bw()+labs(x = "Population Status", y = "Number of Species") +
  facet_wrap(~analysis, ncol = 2) + scale_fill_brewer(palette = "Set1") + theme(legend.position = "none")

#compare distance ratio
par(bg = "white")
par(new = T)
plot.new()
par(oma = c(0,0,2,0))
split.screen(c(2,2))
par(mar = c(5.1,4.1,4.1,2.1))

screen(1)
par(mar = c(3.1,4.1,3.1,2.1))
plot(compiled.results.df$distance_ratio[compiled.results.df$analysis == "By route"], compiled.results.df$distance_ratio[compiled.results.df$analysis == "By grid"], xlab = "", ylab = "By grid", type="n", ylim = c(0,1), xlim = c(0,1))
text(compiled.results.df$distance_ratio[compiled.results.df$analysis == "By route"], compiled.results.df$distance_ratio[compiled.results.df$analysis == "By grid"], label = compiled.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.results.df$distance_ratio[compiled.results.df$analysis == "By grid"] ~ compiled.results.df$distance_ratio[compiled.results.df$analysis == "By route"]), col = "blue")

screen(3)
par(mar = c(5.1,4.1,2.1,2.1))
plot(compiled.results.df$distance_ratio[compiled.results.df$analysis == "By route"], compiled.results.df$distance_ratio[compiled.results.df$analysis == "By strata"], xlab = "By route", ylab = "By strata", type = "n", ylim = c(0,1), xlim = c(0,1))
text(compiled.results.df$distance_ratio[compiled.results.df$analysis == "By route"], compiled.results.df$distance_ratio[compiled.results.df$analysis == "By strata"], label = compiled.results.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.results.df$distance_ratio[compiled.results.df$analysis == "By strata"] ~ compiled.results.df$distance_ratio[compiled.results.df$analysis == "By route"]), col = "blue")

screen(4)
par(mar = c(5.1,4.1,2.1,2.1))
plot(compiled.results.df$distance_ratio[compiled.results.df$analysis == "By grid"], compiled.results.df$distance_ratio[compiled.results.df$analysis == "By strata"],xlab = "By grid", ylab = "", type = "n",ylim = c(0,1), xlim = c(0,1))
text(compiled.results.df$distance_ratio[compiled.results.df$analysis == "By grid"], compiled.results.df$distance_ratio[compiled.results.df$analysis == "By strata"], label = compiled.results.df$sppnumber[1:35])
abline(0,1,col = "black")
abline(lm(compiled.results.df$distance_ratio[compiled.results.df$analysis == "By strata"] ~ compiled.results.df$distance_ratio[compiled.results.df$analysis == "By grid"]), col = "blue")

close.screen(all = T)
title(main = "Distance ratio", outer = T)
