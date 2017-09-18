#Plot centroid on map of US

longs = c(-125,-60)
lats = c(25,50)

library(maps)
setwd("C:/Users/gdicecco/Documents/bbs-centroid/results4/")
centroids <- read.csv("centroids_five_year_windows_all_spp.csv")

#Blackcapped chickadee- 7350
blackcap <- centroids %>%
  filter(aou == 7350)
map(database="state",xlim = longs, ylim = lats)
points(blackcap$centroid_lon[1],blackcap$centroid_lat[1],cex=6,col=rgb(.5,.5,.5,.5),pch=16)
mtext("Blackcapped Chickadee 1969",3,cex=2,line=.5)

map(database="state",xlim = longs, ylim = lats)
points(blackcap$centroid_lon[10],blackcap$centroid_lat[10],cex=6,col=rgb(.5,.5,.5,.5),pch=16)
mtext("Blackcapped Chickadee 2014",3,cex=2,line=.5)

#Carolina chickadee- 7360
carolina <- centroids %>%
  filter(aou == 7360)
map(database="state",xlim = longs, ylim = lats)
points(carolina$centroid_lon[1],carolina$centroid_lat[1],cex=6,col=rgb(.5,.5,.5,.5),pch=16)
mtext("Carolina Chickadee 1969",3,cex=2,line=.5)

map(database="state",xlim = longs, ylim = lats)
points(carolina$centroid_lon[10],carolina$centroid_lat[10],cex=6,col=rgb(.5,.5,.5,.5),pch=16)
mtext("Carolina Chickadee 2014",3,cex=2,line=.5)
