# 
#install.packages("maps")
library(maps)

#install.packages("mapdata")
library(mapdata)

#install.packages("plotrix")
library(plotrix)
#
map('worldHires') ######The command gives you a plot of the World centered around the GMT timezone
# map('world2Hires') ####### The command gives you a plot of a map centered around the International Date line

# map('worldHires','Switzerland', resolution =1) ### We can select an individual country to plot
# map('worldHires','IR', resolution =2) ### We can select an individual country to plot
# map('worldHires', resolution =2) 
# map(database = "world", resolution =2)
# map(database = "world", regions= "Germany", resolution =2)
map(database = "world", xlim = c(-6,18), ylim = c(42,54), resolution =2, col= 9)
map.axes() #Add axes
map_info <- read.table('./pop_gen/popmap_2', header = TRUE)

?map

head(map_info) ### shows the first few rows of your file

col <- map_info$cluster
rbPal <- colorRampPalette(c('red3','green', 'blue3')) # 'blue3','green','red3'
color_ny <- rbPal(20)[as.numeric(cut(col,breaks = 20))]


# points(map_info$lon, map_info$lat, col='blue', pch=2) #### Shows the coordinates on the map. CAUTION: First longitude, second latitude! 
# points(map_info$lon, map_info$lat,col=rainbow(3), pch=2) #### Shows the coordinates on the map. CAUTION: First longitude, second latitude! 
points(map_info$lon, map_info$lat, col=color_ny, pch=1) #### Shows the coordinates on the map. CAUTION: First longitude, second latitude! 

#To add admixture plots based on the Structure results - here I used K = 3.
for (x in 1:nrow(map_info)) {floating.pie(map_info$lon[x],map_info$lat[x], 
                                          c(map_info$K1[x],map_info$K2[x],map_info$K3[x]),radius = 0.7,
                                          col=c("blue", "red", "green"))}
#
help("map")
help("points")
