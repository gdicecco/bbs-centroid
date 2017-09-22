###### Data cleaning ######

library(dplyr)
library(geosphere)
library(maps)

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
RT1.routes <- merge(RT1, routes[ , c("stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
routes.short <- subset(RT1.routes, year >= 1969, select = c("stateroute", "year", "latitude", "longitude", "bcr"))
counts$stateroute <- counts$statenum*1000 + counts$route

#Plot centroid movement for each analysis
par(mfrow=c(1,1))
longs = c(-125,-60)
lats = c(26,50)

####### Centroids by routes ############
#Subset counts df for relevant spp and time period
scale5 <- 5
routes.short$time.window <- scale5*round(routes.short$year/scale5) - 1
counts.short <- counts %>%
  filter(year >= 1969) %>%
  filter(aou %in% huang_species$ID) %>%
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
setwd("C:/Users/gdicecco/Documents/bbs-centroid/centroids-by-routes/")
write.csv(centroids, "centroids-by-routes.csv")

#Plot centroid movement
map(database="world",xlim = longs, ylim = lats)
map(database = "state", add = TRUE)
mtext("Centroids by route",3,cex=2,line=.5)

#shifted distance, velocity, bearing of shift, population change, shift direction regression
results <- matrix(nrow = 35, ncol = 9)
for(i in 1:35) {
  species <- huang_species$ID[i]
  results[i,1] <- species
  df <- centroids %>%
    filter(aou == species)
  
  results[i,2] <- distGeo(df[1,4:3], df[10,4:3])/1000
  results[i,3] <- results[i,2]/(2016-1969)
  results[i,4] <- bearing(df[1,4:3], df[10,4:3])
  results[i,5] <- log(mean(df$mean_total_abund[9:10])/mean(df$mean_total_abund[1:2]))
  
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
                               lon_pval = results[,9])

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
setwd("C:/Users/gdicecco/Documents/bbs-centroid/centroids-by-routes/")
write.csv(results.df, "centroids-by-routes-results.csv")
final.df <- data.frame(species = huang_species$Species, 
                             aou= huang_species$ID, 
                             shiftdistance = results.df$shiftdist, 
                             velocity = results.df$velocity, 
                             direction = results.df$shiftdir, 
                             abundance = results.df$popchange)
write.csv(final.df, "results-routes-summarized.csv")

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
  filter(aou %in% huang_species$ID) %>%
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
setwd("C:/Users/gdicecco/Documents/bbs-centroid/centroids-by-grids/")
write.csv(centroids.grids, "centroids-by-grids.csv")

#Plot centroid movement
map(database="world",xlim = longs, ylim = lats)
map(database = "state", add = TRUE)
mtext("Centroids by 1 deg grid",3,cex=2,line=.5)

#shifted distance, velocity, bearing of shift, population change, shift direction regression
results.grids <- matrix(nrow = 35, ncol = 9)
for(i in 1:35) {
  species <- huang_species$ID[i]
  results.grids[i,1] <- species
  df <- centroids.grids %>%
    filter(aou == species)
  
  results.grids[i,2] <- distGeo(df[1,4:3], df[10,4:3])/1000
  results.grids[i,3] <- results.bcr[i,2]/(2016-1969)
  results.grids[i,4] <- bearing(df[1,4:3], df[10,4:3])
  results.grids[i,5] <- log(mean(df$mean_total_abund[9:10])/mean(df$mean_total_abund[1:2]))
  
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
                         lon_pval = results.grids[,9])

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
setwd("C:/Users/gdicecco/Documents/bbs-centroid/centroids-by-grids/")
write.csv(results.grids.df, "centroids-by-grids-results.csv")
final.grids.df <- data.frame(species = huang_species$Species, 
                           aou= huang_species$ID, 
                           shiftdistance = results.grids.df$shiftdist, 
                           velocity = results.grids.df$velocity, 
                           direction = results.grids.df$shiftdir, 
                           abundance = results.grids.df$popchange)
write.csv(final.grids.df, "results-grids-summarized.csv")


####### Centroids by 1 deg grids & BCR ########

#mean centroid for each species in five year time windows/bcr
spp_abund_means_bcr <- counts.merged %>%
  group_by(stateroute, bcr, latitude, longitude, aou, time.window) %>%
  summarize(avg_abund = mean(speciestotal, na.rm = TRUE))

centroids.bcr <- spp_abund_means_bcr %>%
  group_by(aou, time.window, bcr) %>%
  summarize(centroid_lat = sum(latitude*avg_abund, na.rm = TRUE)/sum(avg_abund, na.rm = TRUE), centroid_lon = sum(longitude*avg_abund, na.rm = TRUE)/sum(avg_abund, na.rm = TRUE),
            mean_total_abund = mean(avg_abund))

#mean centroid weighted from bcr centroids
centroids.bcr2 <- centroids.bcr %>%
  group_by(aou, time.window) %>%
  summarize(centroid_lat = sum(centroid_lat*mean_total_abund, na.rm = TRUE)/sum(mean_total_abund, na.rm = TRUE), centroid_lon = sum(centroid_lon*mean_total_abund, na.rm = TRUE)/sum(mean_total_abund, na.rm = TRUE),
            mean_total_abund = mean(mean_total_abund))
setwd("C:/Users/gdicecco/Documents/bbs-centroid/centroids-by-grids-bcr/")
write.csv(centroids.bcr2, "centroids-by-grids-bcr.csv")

#Plot centroid movement
map(database="world",xlim = longs, ylim = lats)
map(database = "state", add = TRUE)
mtext("Centroids by 1 deg grid and BCR",3,cex=2,line=.5)

#shifted distance, velocity, bearing of shift, population change, shift direction regression
results.bcr <- matrix(nrow = 35, ncol = 9)
for(i in 1:35) {
  species <- huang_species$ID[i]
  results.bcr[i,1] <- species
  df <- centroids.bcr2 %>%
    filter(aou == species)
  
  results.bcr[i,2] <- distGeo(df[1,4:3], df[10,4:3])/1000
  results.bcr[i,3] <- results.bcr[i,2]/(2016-1969)
  results.bcr[i,4] <- bearing(df[1,4:3], df[10,4:3])
  results.bcr[i,5] <- log(mean(df$mean_total_abund[9:10])/mean(df$mean_total_abund[1:2]))
  
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
                         lon_pval = results.bcr[,9])

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

setwd("C:/Users/gdicecco/Documents/bbs-centroid/centroids-by-grids-bcr/")
write.csv(results.bcr.df, "centroids-by-grids-bcr-results.csv")

final.bcr.df <- data.frame(species = huang_species$Species, 
                       aou= huang_species$ID, 
                       shiftdistance = results.bcr.df$shiftdist, 
                       velocity = results.bcr.df$velocity, 
                       direction = results.bcr.df$shiftdir, 
                       abundance = results.bcr.df$popchange)
write.csv(final.bcr.df, "results-bcr-summarized.csv")
