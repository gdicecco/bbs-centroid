# Analysis/plots for Bayesian centroids output

# Libraries
library(dplyr)
library(tidyr)
library(geosphere)
library(maps)
library(ggplot2)
library(RColorBrewer)
library(rgdal)
library(rgeos)

# Read in data
routes <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
weather <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")

bcrshp <- readOGR("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\bcr_terrestrial_shape\\BCR_Terrestrial_master.shp") #BCRs
counts.strata <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\BBS-centroids\\bayesian-model-output\\counts_w_modeloutput.csv", stringsAsFactors = F)
model.grids <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\BBS-centroids\\bayesian-model-output\\model_parameters_grids.csv", stringsAsFactors = F)

huang_species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\BBS-centroids\\huang-2017-bbs-species.csv", header = TRUE)
huang_species <- huang_species %>%
  filter(!aou == 7030, !aou == 7270)

# make counts.grids
routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
routes.short <- subset(RT1.routes, year >= 1969, select = c("statenum","stateroute", "year", "latitude", "longitude", "bcr"))
counts$stateroute <- counts$statenum*1000 + counts$route

#Subset counts by species
counts.short <- counts %>%
  filter(year >= 1969) %>%
  filter(aou %in% huang_species$aou) 
counts.grids <- merge(routes.short, counts.short, by = c("stateroute", "year"))

#Add obs-route ID and dummy variable for first year observers
counts.grids$obsroute <- paste0(counts.grids$obsn, counts.grids$stateroute, sep = "")

counts.grids <- arrange(counts.grids, year)

s <- as.character(paste0(counts.grids$statenum.x, counts.grids$spatial.window, sep = ""))
sall <- sort(unique(s))
is <- match(s, sall)
counts.grids$strata <- is
t <- counts.grids$year - 1968
counts.grids$t <- t

# get fixed effects for strat
model.betas <- model.grids %>% 
  mutate(coef = substr(param, 1, 2)) %>%
  filter(coef == "b1" | coef == "b2")

betas.results <- matrix(0, nrow = 1, ncol = ncol(model.betas) + 1)
colnames(betas.results) <- c(colnames(model.betas), "strata")

for(i in 1:length(unique(model.betas$aou))) {
  model.short <- filter(model.betas, aou == unique(model.betas$aou)[i], coef == "b1", model == "fixedEffects")
  model.short$strata <- as.numeric(rownames(model.short))
  
  model.b2 <- filter(model.betas, aou == unique(model.betas$aou)[i], coef == "b2", model == "fixedEffects")
  model.b2$strata <- as.numeric(rownames(model.b2))
  
  model.short.re <- filter(model.betas, aou == unique(model.betas$aou)[i], coef == "b1", model == "randomEffects")
  model.short.re$strata <- as.numeric(rownames(model.short.re))
  
  model.b2.re <- filter(model.betas, aou == unique(model.betas$aou)[i], coef == "b2", model == "randomEffects")
  model.b2.re$strata <- as.numeric(rownames(model.b2.re))
  
  
  tmp <- rbind(model.short, model.b2, model.short.re, model.b2.re)
  betas.results <- rbind(betas.results, tmp)
}
betas.results <- as.data.frame(betas.results[-1, ])

abunds <- betas.results %>%
  select(aou, model, mean, coef, strata) %>%
  spread(key = coef, value = mean)

counts.grids.ind <- counts.grids %>%
  left_join(abunds, by = c("aou", "strata")) %>%
  mutate(abundind = exp(b1 + b2*t))
counts.grids.re <- filter(counts.grids.ind, model == "randomEffects")

#Plot centroid movement for each analysis
par(mfrow=c(1,1))
longs = c(-125,-60)
lats = c(26,50)

#calcuate distance ratio- determine if centroid moves in consistent direction over time
#input- centroids df with centroids per time window for all spp. (10 lat lon coordinates of centroid over 5 year time windows)
#will subset centroids df for spp of interest by AOU
distance.ratio <- function(df, AOU, randomizations = FALSE, n = 1000) {
  
  x = filter(df, aou == AOU)
  
  years = seq(1, 48)
  
  if (randomizations) {
    rand.ratios = c()
    for (r in 1:n) {
      sampleyears = sample(years, 10)
      distances <- c()
      for(i in 1:48) {
        distances <- c(distances, distGeo(x[x$time.window == sampleyears[i], c("centroid_lon", "centroid_lat")], 
                                          x[x$time.window == sampleyears[i+1], c("centroid_lon", "centroid_lat")]))
        sumdist <- sum(distances)
      }
      fl <- distGeo(x[1,c("centroid_lon", "centroid_lat")], x[48,c("centroid_lon", "centroid_lat")])
      rand.ratios = c(rand.ratios, fl/sumdist)
    }
  } else {
    rand.ratios = NULL
  }
  
  sampleyears = years
  distances <- c()
  for(i in 1:48) {
    distances <- c(distances, distGeo(x[x$time.window == sampleyears[i], c("centroid_lon", "centroid_lat")], 
                                      x[x$time.window == sampleyears[i+1], c("centroid_lon", "centroid_lat")]))
    sumdist <- sum(na.omit(distances))
  }
  fl <- distGeo(x[1,c("centroid_lon", "centroid_lat")], x[48,c("centroid_lon", "centroid_lat")])
  obs.ratio = fl/sumdist
  
  return(list(obs = obs.ratio, rand = rand.ratios))
} 

####### Centroids by 1 degree grids ######
#Group routes by spatial windows
scale1 <- 1
counts.grids.re$lat.window <- scale1*floor(counts.grids.re$latitude/scale1) + scale1/2
counts.grids.re$lon.window <- scale1*floor(counts.grids.re$longitude/scale1) + scale1/2
counts.grids.re$spatial.window <- paste(counts.grids.re$lat.window, counts.grids.re$lon.window, sep = "")

counts.grids.re$ceiling.lat <- ceiling(counts.grids.re$lat.window)
counts.grids.re$floor.lon <- floor(counts.grids.re$lon.window)

#mean centroid for each species in five year time windows
centroids.grids <- counts.grids.re %>%
  group_by(aou, t) %>% ##for each grid cell
  summarize(centroid_lat = sum(ceiling.lat*abundind, na.rm = TRUE)/sum(abundind, na.rm = TRUE),
            centroid_lon = sum(floor.lon*abundind, na.rm = TRUE)/sum(abundind, na.rm = TRUE),
            mean_total_abund = mean(abundind)) 

colnames(centroids.grids)[2] <- "time.window"

setwd("C:/Users/gdicecco/Desktop/git/bbs-centroid/bayesian-output/")
write.csv(centroids.grids, "centroids-by-grids-raw.csv", row.names=F)

#Plot centroid movement
#map(database="world",xlim = longs, ylim = lats) #plot just on states
#map(database = "state", add = TRUE)
plot(bcrshp[bcrshp$WATER == 3,], ylim = lats, xlim = longs, border = "gray73", col = "gray95") #plot centroids on BCR map
mtext("Centroids by 1 deg grid",3,cex=2,line=.5)

#shifted distance, velocity, bearing of shift, population change, shift direction regression
results.grid <- matrix(nrow = length(unique(centroids.grids$aou)), ncol = 12)
for(i in 1:length(unique(centroids.grids$aou))) {
  species <- unique(centroids.grids$aou)[i]
  results.grid[i,1] <- species
  df <- centroids.grids %>%
    filter(aou == species)
  lastrow <- nrow(df)
  
  results.grid[i,2] <- distGeo(df[1,4:3], df[48,4:3])/1000
  results.grid[i,3] <- results.grid[i,2]/(2016-1969)
  results.grid[i,4] <- bearing(df[1,4:3], df[48,4:3])
  results.grid[i,5] <- log(mean(df$mean_total_abund[47:48])/mean(df$mean_total_abund[1:2]))
  distratios <- distance.ratio(centroids.grids, species, F)
  results.grid[i,10] <- distratios$obs
  results.grid[i,11] <- mean(distratios$rand)
  results.grid[i,12] <- sd(distratios$rand)
  
  
  mod.test <- lm(df$centroid_lat ~ df$time.window)
  sum <- summary(mod.test)
  results.grid[i,7] <- sum$coefficients[2,4]
  results.grid[i,6] <- sum$coefficients[2,1]
  
  mod.test.2 <- lm(df$centroid_lon ~ df$time.window)
  sum.2 <- summary(mod.test.2)
  results.grid[i,9] <- sum.2$coefficients[2,4]
  results.grid[i,8] <- sum.2$coefficients[2,1]
  
  points(df$centroid_lon[1], df$centroid_lat[1], pch = 16, col = 'red')
  points(df$centroid_lon, df$centroid_lat, type = 'l', col = "gray45")
  points(df$centroid_lon[48], df$centroid_lat[48], pch = 17, col = 'blue')
  text(df$centroid_lon[48], df$centroid_lat[48]+1, species, cex = .75)
}
results.grids.df <- data.frame(aou = results.grids[,1], 
                               shiftdist = results.grid[,2], 
                               velocity = results.grid[,3], 
                               bearing = results.grid[,4], 
                               r = results.grid[,5],
                               lat_slope = results.grid[,6], 
                               lat_pval = results.grid[,7], 
                               lon_slope = results.grid[,8], 
                               lon_pval = results.grid[,9],
                               distance_ratio = results.grid[,10],
                               dratio_rand_mean = results.grid[,11],
                               dratio_rand_sd = results.grid[,12])

#assign shift directions
direction.lat <- c()
for(i in 1:35){
  slope <- results.grid.df$lat_slope[i]
  pval <- results.grid.df$lat_pval[i]
  if (slope > 0 & pval < 0.05) {
    direction.lat <- c(direction.lat, "north")
  } else if (slope < 0 & pval < 0.05) {
    direction.lat <- c(direction.lat, "south")
  } else 
    direction.lat <- c(direction.lat, "")
}
results.grid.df <- cbind(results.grid.df, direction.lat)

direction.lon <- c()
for(i in 1:35){
  slope <- results.grid.df$lon_slope[i]
  pval <- results.grid.df$lon_pval[i]
  if (slope > 0 & pval < 0.05) {
    direction.lon <- c(direction.lon, "east")
  } else if (slope < 0 & pval < 0.05) {
    direction.lon <- c(direction.lon, "west")
  } else 
    direction.lon <- c(direction.lon, "")
}
results.grid.df <- cbind(results.grid.df,direction.lon)
results.grid.df$shiftdir <- paste(results.grid.df$direction.lat,results.grid.df$direction.lon, sep = "")

popchange <- c()
for(i in 1:35){
  r <- results.grid.df$r[i]
  if (r > 0) {
    popchange <- c(popchange, "increasing")
  } else popchange <- c(popchange, "decreasing")
}
results.grid.df <- cbind(results.grid.df, popchange)
setwd("C:/Users/gdicecco/Desktop/git/bbs-centroid/bayesian-output/")
write.csv(results.grid.df, "centroids-by-grids-results.csv", row.names=F)
final.grids.df <- data.frame(species = huang_species$species, 
                             aou= huang_species$aou, 
                             shiftdistance = results.grid.df$shiftdist, 
                             velocity = results.grid.df$velocity, 
                             direction = results.grid.df$shiftdir, 
                             abundance = results.grid.df$popchange,
                             distanceratio = results.grid.df$distance_ratio)
write.csv(final.grids.df, "results-grids-summarized.csv", row.names=F)

####### Centroids by strata ########

#Calculate BCR/state strata centroids
regioncodes <- read.table("C:/Users/gdicecco/Desktop/git/bbs-centroid/regioncodes.txt", sep = "\t")
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
routes.short.centers <- left_join(routes.short, polys.merged.land, by = "statebcr")
counts.strata.centers <- right_join(routes.short.centers, counts.strata)

subs.routes4 <- counts.strata.centers %>%
  group_by(aou, statebcr) %>%
  distinct(stateroute) %>%
  summarize(total = n()) %>%
  filter(total > 4)

bcrroutes <- routes.short.centers %>%
  group_by(statebcr, year) %>%
  select(stateroute) %>%
  summarize(total.routes = n())

counts.subs.strata <- inner_join(counts.strata.centers, subs.routes4) %>%
  left_join(bcrroutes)

#mean centroid for each species in five year time windows per bcr/state strata
centroids.bcr <- counts.subs.strata %>%
  group_by(aou, t, statebcr) %>%
  summarize(centroid_lat = sum(y*abundind, na.rm = TRUE)/sum(abundind, na.rm = TRUE), centroid_lon = sum(x*abundind, na.rm = TRUE)/sum(abundind, na.rm = TRUE),
            mean_total_abund = mean(abundind))

#mean centroid weighted from bcr centroids
centroids.bcr2 <- centroids.bcr %>%
  group_by(aou, t) %>%
  summarize(centroid_lat = sum(centroid_lat*mean_total_abund, na.rm = TRUE)/sum(mean_total_abund, na.rm = TRUE), centroid_lon = sum(centroid_lon*mean_total_abund, na.rm = TRUE)/sum(mean_total_abund, na.rm = TRUE),
            mean_total_abund = mean(mean_total_abund))

colnames(centroids.bcr2)[2] <- "time.window"

setwd("C:/Users/gdicecco/Desktop/git/bbs-centroid/bayesian-output/")
write.csv(centroids.bcr2, "centroids-by-strata.csv", row.names=F)

#Plot centroid movement
#map(database="world",xlim = longs, ylim = lats) #plot just on states
#map(database = "state", add = TRUE)
plot(bcrshp[bcrshp$WATER == 3,], ylim = lats, xlim = longs, border = "gray73", col = "gray95") #plot centroids on BCR map
mtext("Centroids by strata",3,cex=2,line=.5)

#shifted distance, velocity, bearing of shift, population change, shift direction regression
results.bcr <- matrix(nrow = length(unique(centroids.bcr2$aou)), ncol = 12)
for(i in 1:length(unique(centroids.bcr2$aou))) {
  species <- unique(centroids.bcr2$aou)[i]
  results.bcr[i,1] <- species
  df <- centroids.bcr2 %>%
    filter(aou == species)
  
  results.bcr[i,2] <- distGeo(df[1,4:3], df[48,4:3])/1000
  results.bcr[i,3] <- results.bcr[i,2]/(2016-1969)
  results.bcr[i,4] <- bearing(df[1,4:3], df[48,4:3])
  results.bcr[i,5] <- log(mean(df$mean_total_abund[47:48])/mean(df$mean_total_abund[1:2]))
  distratios <- distance.ratio(centroids.bcr2, species, F)
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
  points(df$centroid_lon, df$centroid_lat, type = 'l', col = "gray45")
  points(df$centroid_lon[48], df$centroid_lat[48], pch = 17, col = 'blue')
  text(df$centroid_lon[48], df$centroid_lat[48]+1, species, cex = .75)
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

setwd("C:/Users/gdicecco/Desktop/git/bbs-centroid/bayesian-output/")
write.csv(results.bcr.df, "centroids-by-strata-results.csv", row.names=F)

final.bcr.df <- data.frame(aou= unique(centroids.bcr2$aou), 
                           shiftdistance = results.bcr.df$shiftdist, 
                           velocity = results.bcr.df$velocity, 
                           direction = results.bcr.df$shiftdir, 
                           abundance = results.bcr.df$popchange,
                           distanceratio = results.bcr.df$distance_ratio)
write.csv(final.bcr.df, "results-strata-summarized.csv", row.names=F)

####### Centroids by strata & weighted abundances #########

#A - ratio of stratum area over the area of all strata where species is present
#z - proportion of routes in specific strata on which spp was spotted
#Abundance index = A*z*total mean abund for stratum
#function abund.index: requires spp_abund_means_bcr, table with species abundances per route, year, statebcr
#requires polys.merged.land, df of statebcr strata containing Shape_Area for each strata
#requires routes.short.centers, df containing statebcrs and routes assigned to each strata
#returns abundance index weighted by scaling factors A and Z
abund.index <- function(year, spid, statebcr) {
  test <- as.character(statebcr)
  ai <- spp_abund_means_bcr %>%
    filter(time.window == year) %>%
    filter(aou == spid) %>%
    filter(statebcr == test)
  total.area <- sum(polys.merged.land$Shape_Area[polys.merged.land$statebcr %in% spp_abund_means_bcr$statebcr[spp_abund_means_bcr$aou == spid]])
  a <- polys.merged.land$Shape_Area[polys.merged.land$statebcr == strata]/total.area
  z <- length(unique(ai$stateroute))/length(unique(routes.short.centers$stateroute[routes.short.centers$statebcr == strata]))
  
  index <- a*z*mean(ai$avg_abund, na.rm = TRUE)
  return(index)
}

colnames(spp_abund_means_bcr)[9] <- "time.window"
abundance.indices <- data.frame(statebcr = 0, x = 0, y = 0, aou = 0, time.window = 0, abund.index = 0)

#Takes a really long time ~16 hours
init.time = Sys.time()
for(i in 1:length(unique(centroids.bcr2$aou))) {
  species <- unique(centroids.bcr2$aou)[i]
  bcrs <- spp_abund_means_bcr %>%
    filter(aou == species) 
  bcr.list <- unique(bcrs$statebcr)
  for(b in bcr.list) {
    for(year in seq(1, 48)) {
      if(b %in% bcrs$statebcr[bcrs$time.window == year]){
        tseries <- bcrs %>%
          filter(time.window == year, statebcr == b) 
        index <- abund.index(year, species, b)
        data <- data.frame(statebcr = b, x = tseries$x[tseries$statebcr == b], y = tseries$y[tseries$statebcr == b],
                           aou = species, time.window = year, abund.index = index)
        abundance.indices <- bind_rows(abundance.indices, unique(data))
      } else {
        data <- data.frame(statebcr = b, x = NA, y = NA,
                           bcr = NA,
                           aou = species, time.window = year, abund.index = NA)
        abundance.indices <- bind_rows(abundance.indices, unique(data))
      }
    }
  }
}
abundance.indices <- abundance.indices[-1,]
setwd("C:/Users/gdicecco/Desktop/git/bbs-centroid/bayesian-output/")
write.csv(abundance.indices, "weighted-abundance-indices-centroids.csv", row.names = F)

#mean centroid weighted from strata specific abund indices
centroids.weighted <- abundance.indices %>%
  group_by(aou, time.window) %>%
  summarize(centroid_lat = sum(y*abund.index, na.rm = TRUE)/sum(abund.index, na.rm = TRUE), 
            centroid_lon = sum(x*abund.index, na.rm = TRUE)/sum(abund.index, na.rm = TRUE),
            mean_total_ai = mean(abund.index, na.rm = TRUE))
write.csv(centroids.weighted, "centroids-weighted-abund.csv", row.names = F)

plot(bcrshp[bcrshp$WATER == 3,], ylim = c(26,60), xlim = c(-140,-60), border = "gray73", col = "gray95") #plot centroids on BCR map
mtext("Centroids by strata with weighted abundance index",3,cex=2,line=.5)

#shifted distance, velocity, bearing of shift, population change, shift direction regression
results.weighted <- matrix(nrow = length(unique(centroids.bcr2$aou)), ncol = 12)
for(i in 1:length(unique(centroids.bcr2$aou))) {
  species <- unique(centroids.bcr2$aou)[i]
  results.weighted[i,1] <- species
  df <- centroids.weighted %>%
    filter(aou == species)
  
  results.weighted[i,2] <- distGeo(df[1,4:3], df[48,4:3])/1000
  results.weighted[i,3] <- results.weighted[i,2]/(2016-1969)
  results.weighted[i,4] <- bearing(df[1,4:3], df[48,4:3])
  results.weighted[i,5] <- log(mean(df$mean_total_ai[47:48])/mean(df$mean_total_ai[1:2]))
  distratios <- distance.ratio(centroids.weighted, species, F)
  results.weighted[i,10] <- distratios$obs
  results.weighted[i,11] <- mean(distratios$rand)
  results.weighted[i,12] <- sd(distratios$rand)
  
  mod.test <- lm(df$centroid_lat ~ df$time.window)
  sum <- summary(mod.test)
  results.weighted[i,7] <- sum$coefficients[2,4]
  results.weighted[i,6] <- sum$coefficients[2,1]
  
  mod.test.2 <- lm(df$centroid_lon ~ df$time.window)
  sum.2 <- summary(mod.test.2)
  results.weighted[i,9] <- sum.2$coefficients[2,4]
  results.weighted[i,8] <- sum.2$coefficients[2,1]
  
  points(df$centroid_lon[1], df$centroid_lat[1], pch = 16, col = 'red')
  points(df$centroid_lon, df$centroid_lat, type = 'l', col = "gray45")
  points(df$centroid_lon[48], df$centroid_lat[48], pch = 17, col = 'blue')
  text(df$centroid_lon[48], df$centroid_lat[48]+1, species, cex = .75)
  
  
}
results.weighted.df <- data.frame(aou = results.weighted[,1], 
                                  shiftdist = results.weighted[,2], 
                                  velocity = results.weighted[,3], 
                                  bearing = results.weighted[,4], 
                                  r = results.weighted[,5],
                                  lat_slope = results.weighted[,6], 
                                  lat_pval = results.weighted[,7], 
                                  lon_slope = results.weighted[,8], 
                                  lon_pval = results.weighted[,9],
                                  distance_ratio = results.weighted[,10],
                                  dratio_rand_mean = results.weighted[,11],
                                  dratio_rand_sd = results.weighted[,12])

#assign shift directions
direction.lat <- c()
for(i in 1:34){
  slope <- results.weighted.df$lat_slope[i]
  pval <- results.weighted.df$lat_pval[i]
  if (slope > 0 & pval < 0.05) {
    direction.lat <- c(direction.lat, "north")
  } else if (slope < 0 & pval < 0.05) {
    direction.lat <- c(direction.lat, "south")
  } else 
    direction.lat <- c(direction.lat, "")
}
results.weighted.df <- cbind(results.weighted.df, direction.lat)

direction.lon <- c()
for(i in 1:34){
  slope <- results.weighted.df$lon_slope[i]
  pval <- results.weighted.df$lon_pval[i]
  if (slope > 0 & pval < 0.05) {
    direction.lon <- c(direction.lon, "east")
  } else if (slope < 0 & pval < 0.05) {
    direction.lon <- c(direction.lon, "west")
  } else 
    direction.lon <- c(direction.lon, "")
}
results.weighted.df <- cbind(results.weighted.df,direction.lon)
results.weighted.df$shiftdir <- paste(results.weighted.df$direction.lat,results.weighted.df$direction.lon, sep = "")

popchange <- c()
for(i in 1:34){
  r <- results.weighted.df$r[i]
  if (r > 0) {
    popchange <- c(popchange, "increasing")
  } else popchange <- c(popchange, "decreasing")
}
results.weighted.df <- cbind(results.weighted.df, popchange)

setwd("C:/Users/gdicecco/Desktop/git/bbs-centroid/bayesian-output/")
write.csv(results.weighted.df, "centroids-weighted-results.csv", row.names=F)

final.weighted.df <- data.frame( 
                                aou= results.weighted.df$aou, 
                                shiftdistance = results.weighted.df$shiftdist, 
                                velocity = results.weighted.df$velocity, 
                                direction = results.weighted.df$shiftdir, 
                                abundance = results.weighted.df$popchange,
                                distanceratio = results.weighted.df$distance_ratio)
write.csv(final.weighted.df, "results-strata-weighted-summarized.csv", row.names=F)

##### Summary plots #####

setwd("C:/Users/gdicecco/Desktop/git/bbs-centroid/bayesian-output/")
final.grids.df <- read.csv("results-grids-summarized.csv", stringsAsFactors = F)
final.bcr.df <- read.csv("results-strata-summarized.csv", stringsAsFactors = F)
final.weighted.df <- read.csv("results-strata-weighted-summarized.csv", stringsAsFactors = F)

#centroids by grid
final.grids.df$analysis <- rep(x = "By grid", times = nrow(final.grids.df))
#centroids by bcr
final.bcr.df$analysis <- rep(x = "By strata", times = nrow(final.bcr.df))
#centroids weighted
final.weighted.df$analysis <- rep(x = "By strata with weighted abund.", times = nrow(final.weighted.df))
#huang data
huang_species$distanceratio <- rep(NA)
huang_species$analysis <- rep(x = "Huang et al.", times = nrow(huang_species))

compiled.df <- rbind(final.grids.df[, -1], final.bcr.df, final.weighted.df, huang_species[, -1])

blank <- theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#compare shift distance
##boxplot
ggplot(compiled.df, aes(analysis, shiftdistance)) + theme_bw() + geom_boxplot(fill = "gray") + labs(x = "Analysis", y = "Shift Distance") + blank


#Polar plot
compiled.df$direction <- factor(compiled.df$direction, levels = c("north", "northeast", "east", "southeast", "south", "southwest", "west", "northwest"))

plot <- ggplot(compiled.df, aes(x = direction, fill = analysis))+stat_count()+theme_bw()+labs(x = "Direction of centroid shift", y = "Number of Species") +
  facet_wrap(~analysis, ncol = 4)+theme(axis.text.x = element_text(face = "bold", size = 7.5), legend.position = "none") + scale_fill_brewer(palette = "Set1")

plot + coord_polar(start = -pi/8, direction = 1)

ggplot(compiled.df, aes(x = aou, y = shiftdistance, color = analysis)) + theme_bw() + geom_point() + scale_color_brewer(palette = "Set1") + labs(x = "AOU", y = "Shift Distance", color = "Analysis")
