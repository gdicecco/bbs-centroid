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

#Remove routes that were not sampled at least once each time window between 1969 and 2016
##This only results in 29 usable routes

route.results <- matrix(nrow = 5065, ncol = 1)
stateroutes <- unique(RT1.routes$stateroute)

for(i in 1:5065) {
  route <- stateroutes[i]

  route2 <- subset(RT1.routes, stateroute == route, select = c("year"))

    results.route2 <- c()
  for(i in 1:length(route2$year)) {
    year <- route2$year[i]
    if (year >=1969 & year <=1973) {
      results.route2 <- c(results.route2, 1)
    } else if (year >=1974 & year <= 1978) {
      results.route2 <- c(results.route2, 2)
    } else if (year >= 1979 & year <= 1983) {
      results.route2 <- c(results.route2, 3)
    } else if (year >=1984 & year <= 1988) {
      results.route2 <- c(results.route2, 4)
    } else if (year >=1989 & year <= 1993) {
      results.route2 <- c(results.route2, 5)
    } else if (year >=1994 & year <= 1998) {
      results.route2 <- c(results.route2, 6)
    } else if (year >=1999 & year <= 2003) {
      results.route2 <- c(results.route2, 7)
    } else if (year >=2004 & year <= 2008) {
      results.route2 <- c(results.route2, 8)
    } else if (year >=2009 & year <= 2013) {
      results.route2 <- c(results.route2, 9)
    } else if (year >= 2014 & year <= 2016) {
      results.route2 <- c(results.route2, 10)
    } else 
    results.route2 <- c(results.route2, 11)
  }
  
    windows <- data.frame(window = unique(results.route2))
    onetoten <- data.frame(window = seq(1,10, by = 1))
    merge.windows <- merge(onetoten, windows, by = "window")

if (length(merge.windows$window) == 10 ) {
  route.results[i,1] <- route
}
}

stateroutes.subs <- na.omit(route.results)

#Average routes for each spp over five year time windows, calculate centroid, record shifted distance, velocity, shift bearing, population change
#install.packages("dplyr")
#install.packages("DBI")
#install.packages("lazyeval")
#install.packages("geosphere")
library("dplyr")
library("geosphere")

years <- seq(1966, 2016, by = 5)
final <- matrix(nrow = 35, ncol = 5)

setwd("C:/Users/gdicecco/Documents/bbs-centroid/results")

for(i in 1:35) {
  species <- huang_species[i,1]
  final[i,1] <- species

  spp1 <- subset(counts, aou == species, select = c("stateroute", "year", "aou", "speciestotal"))
  spp1.routes <- merge(spp1, RT1.routes[,c("stateroute", "latitude","longitude")], by = "stateroute", all.x = TRUE)

  results <- matrix(nrow = 10, ncol = 4) 
for (j in 1:10) {
yr <- years[j]  
results[j,1] <- yr
bobwhite <- distinct(filter(spp1.routes, year >= yr & year <= yr + 4, aou == species))
bobwhite.mean <- bobwhite %>%
  group_by(stateroute, latitude, longitude) %>%
  summarize(mean_abund = mean(speciestotal, na.rm = TRUE))
results[j, 4] <- mean(bobwhite.mean$mean_abund)
results[j,2] <- sum(bobwhite.mean$latitude*bobwhite.mean$mean_abund, na.rm=TRUE)/sum(bobwhite.mean$mean_abund, na.rm=TRUE)
results[j,3] <- sum(bobwhite.mean$longitude*bobwhite.mean$mean_abund, na.rm=TRUE)/sum(bobwhite.mean$mean_abund, na.rm = TRUE)
}
bobwhite.centroid <- data.frame(year = results[,1], centroid.lat = results[,2], centroid.lon = results[,3], mean_abund = results[,4])
write.csv(bobwhite.centroid, paste("AOU_", species, "_centroid_five_year_windows.csv", sep = ""))

shift <- distGeo(bobwhite.centroid[1,3:2], bobwhite.centroid[10,3:2])/1000
bring <- bearing(bobwhite.centroid[1,3:2], bobwhite.centroid[10,3:2])
velocity <- shift/((bobwhite.centroid[10,1]+4)-bobwhite.centroid[1,1])
delta_abund <- bobwhite.centroid[10,4] - bobwhite.centroid[1,4]

final[i, 2] <- shift
final[i, 3] <- velocity
final[i, 4] <- bring
final[i, 5] <- delta_abund
}

final.df <- data.frame(aou = final[,1], shifted_dist = final[,2], velocity = final[,3], bearing = final[,4], population_change = final[,5])
write.csv(final.df, 'centroid_shifts_all_spp.csv')

#determine direction and abundance change per Huang et al methods

files <- list.files("C:/Users/gdicecco/Documents/bbs-centroid/results")
results.2 <- matrix(nrow = 35, ncol = 5)
for(i in 1:35) {
  test<- read.csv(files[i], header = TRUE)
  
  mod.test <- lm(test$centroid.lat ~ test$year)
  sum <- summary(mod.test)
  results.2[i,2] <- sum$coefficients[2,4]
  results.2[i,1] <- sum$coefficients[2,1]
  
  mod.test.2 <- lm(test$centroid.lon ~ test$year)
  sum.2 <- summary(mod.test.2)
  results.2[i,4] <- sum.2$coefficients[2,4]
  results.2[i,3] <- sum.2$coefficients[2,1]

results.2[i,5] <- log(mean(test$mean_abund[9:10])/mean(test$mean_abund[1:2]))
}

results.df <- data.frame(species = huang_species$Species, aou = huang_species$ID, shift_distance = final.df$shifted_dist, slope_lat = results.2[,1], pvalue_lat = results.2[,2], slope_lon = results.2[,3], pvalue_lon = results.2[,4], r = results.2[,5])
write.csv(results.df, "centroid_abundance_shifts_all_spp.csv")
