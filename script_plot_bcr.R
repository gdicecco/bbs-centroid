library(rgdal)
bcrshp <- readOGR("C:/Users/gdicecco/Documents/bcr_terrestrial_shape/BCR_Terrestrial_master.shp")
longs = c(-125,-60)
lats = c(26,50)

plot(bcrshp, ylim = lats, xlim = longs)