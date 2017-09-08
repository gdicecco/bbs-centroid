#Plot centroid on map of US

longs = c(-125,-60)
lats = c(25,50)

library(maps)

map(database="state",xlim = longs, ylim = lats)
points(-87.86856473,35.57657859,cex=6,col=rgb(.5,.5,.5,.5),pch=16)
mtext("1966-1970",3,cex=2,line=.5)

map(database="state",xlim = longs, ylim = lats)
points(-87.2798307,33.30592728,cex=6,col=rgb(.5,.5,.5,.5),pch=16)
mtext("2011-2015",3,cex=2,line=.5)