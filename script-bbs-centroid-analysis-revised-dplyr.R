#install.packages("geosphere")
#install.packages("dplyr")
library("geosphere")
library("dplyr")

#mean species abundance for each route over time.windows
spp_abund_means <- counts.merged %>%
  group_by(stateroute, latitude, longitude, aou, time.window) %>%
  summarize(avg_abund = mean(speciestotal, na.rm = TRUE))

#mean centroid for each species in five year time windows
centroids.grids <- means %>%
  group_by(aou, time.window) %>%
  summarize(centroid_lat = sum(latitude*avg_abund, na.rm = TRUE)/sum(avg_abund, na.rm = TRUE), centroid_lon = sum(longitude*avg_abund, na.rm = TRUE)/sum(avg_abund, na.rm = TRUE),
            mean_total_abund = mean(avg_abund))
#setwd("C:/Users/gdicecco/Documents/bbs-centroid/results2/")
#write.csv(centroids, "centroids_five_year_windows_all_spp.csv")

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
  text(df$centroid_lon[10], df$centroid_lat[10]+.6, species, cex = .75)
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
direction <- c()
for(i in 1:35){
  slope <- results.df$lat_slope[i]
  pval <- results.df$lat_pval[i]
  if (slope > 0 & pval < 0.05) {
    direction <- c(direction, "north")
  } else if (slope < 0 & pval < 0.05) {
    direction <- c(direction, "south")
  } else 
    direction <- c(direction, "")
}
results.df <- cbind(results.df, direction)

direction2 <- c()
for(i in 1:35){
  slope <- results.df$lon_slope[i]
  pval <- results.df$lon_pval[i]
  if (slope > 0 & pval < 0.05) {
    direction2 <- c(direction2, "east")
  } else if (slope < 0 & pval < 0.05) {
    direction2 <- c(direction2, "west")
  } else 
    direction2 <- c(direction2, "")
}
results.df <- cbind(results.df,direction2)
results.df$shiftdir <- paste(results.df$direction,results.df$direction2, sep = "")

popchange <- c()
for(i in 1:35){
  r <- results.df$r[i]
  if (r > 0) {
    popchange <- c(popchange, "increasing")
  } else popchange <- c(popchange, "decreasing")
}
results.df <- cbind(results.df, popchange)

#setwd("C:/Users/gdicecco/Documents/bbs-centroid/results2/")
#write.csv(results.df, "spatial-temporal-coverage-results.csv")
