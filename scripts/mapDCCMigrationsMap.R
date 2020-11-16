#Plots migration rates from BEAST on a world map, plots direction so assumes asymmetric phylogeography
#Argument 2 is the DCC to be plotted and determines the points plotted on the map, can be DCC1, DCC2 or DCC3
#Used to plot Figure 2C

args <- commandArgs(TRUE)

library(ggplot2)
library(maps)
library(grid)

#Import the supported migrations
migrationRates <- read.table(args[1],sep="\t",header=TRUE)

#Assign to the DCC being plotted
dcc <- args[2]

outFile <- args[3]

continentLocations <- data.frame(Continent = c("Asia","Europe","NorthAmerica","Oceania","SouthAmerica"), latitude = c(31.2304,51.5074,35.9132,-33.8688,-22.9068), longitude = c(121.4737,-0.1278,-79.0558,151.2093,-43.1729))

migrationRates$StartLatitude <- 0
migrationRates$StartLongitude <- 0
migrationRates$EndLatitude <- 0
migrationRates$EndLongitude <- 0

for (continent in 1:length(migrationRates[,1])) {
  migrationRates$StartLatitude[continent] <- continentLocations[which(continentLocations[,1] == toString(migrationRates[continent,1])),2]
  migrationRates$StartLongitude[continent] <- continentLocations[which(continentLocations[,1] == toString(migrationRates[continent,1])),3]
  migrationRates$EndLatitude[continent] <- continentLocations[which(continentLocations[,1] == toString(migrationRates[continent,2])),2]
  migrationRates$EndLongitude[continent] <- continentLocations[which(continentLocations[,1] == toString(migrationRates[continent,2])),3]}

#Plot different points based on the DCC
if (dcc == "DCC1") {
  uniqueContinents <- c("Asia", "Europe", "NorthAmerica", "Oceania")} else {
  if (dcc == "DCC2") {
    uniqueContinents <- c("Europe", "NorthAmerica", "Oceania")} else {
    uniqueContinents <- c("Asia","Europe","NorthAmerica","Oceania","SouthAmerica")}}
continentPoints <- data.frame(Continents = uniqueContinents,Latitude = 0,Longitude = 0)
for (continent in 1:length(uniqueContinents)) {
  continentPoints[continent,2] <- continentLocations[which(continentLocations[,1] == toString(uniqueContinents[continent])),2]
  continentPoints[continent,3] <- continentLocations[which(continentLocations[,1] == toString(uniqueContinents[continent])),3]}

worldMap <- map_data('world')
pdf(outFile)
baseWorld <- ggplot() + coord_fixed() + xlab('') + ylab('') + geom_polygon(data=worldMap,aes(x=long,y=lat,group=group),colour='tan',fill='tan') + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_rect(fill='white',colour='white'),axis.line=element_line(colour='white'),legend.position='none',axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())
mapData <- baseWorld + geom_curve(data = migrationRates,aes(x=StartLongitude,y=StartLatitude,xend=EndLongitude,yend=EndLatitude),arrow=arrow(type="closed"),colour=migrationRates$Colour) + geom_point(data = continentPoints,aes(x=Longitude,y=Latitude),size=2)
print(mapData)
dev.off()