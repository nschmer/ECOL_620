########################################################
########################################################
#Fletcher and Fortin 2019
#Chapter 8: Space Use and Resource Selection
########################################################
########################################################

#load packages

library(raster)           #for raster covariate data; version 2.6-7 used
library(reshape2)         #for re-formatting data; version 1.4.3 used
library(rgdal)            #for reading different types of GIS files; version 1.3-4 used
library(adehabitatLT)     #for trajectory information; version 0.3.23 used
library(adehabitatHR)     #Home range estimation; version 0.4.15 used
library(adehabitatHS)     #for selection ratios; version 0.3.13 used
library(survival)         #for conditional logit model; version 2.42-3 used
library(ggplot2)

#set working directory where data were downloaded
setwd("/Users/natalieschmer/Documents/GitHub/CSU/ECOL_620/scripts/lab7")

###################################################
#8.3.2 Prepping the data
###################################################
#label projection for later use
crs.land <- "+proj=aea +lat_0=24 +lon_0=-84 +lat_1=24 +lat_2=31.5 +x_0=400000 +y_0=0 +ellps=GRS80 +datum=WGS84 +units=m +no_defs"

#landcover source: fwc/fnai
land <- raster("/Users/natalieschmer/Documents/GitHub/CSU/ECOL_620/lab_7/data_for_lab7/panther_landcover")
plot(land)

#check projection
projection(land)

#Add panther data
panthers <- readOGR("/Users/natalieschmer/Documents/GitHub/CSU/ECOL_620/lab_7/data_for_lab7/panthers.shp")
panthers=spTransform(panthers, crs.land)
plot(panthers)
#the x and y variables are likely in a different CRS, so I have removed them to avoid confusion
panthers=panthers[,3:5]

#check projection
projection(panthers)

#inspect
summary(panthers)
unique(panthers$CatID) #the unique cat IDs, Julian date for BBMM model 
head(panthers)

#plot
plot(land)
points(panthers, col=panthers$CatID)
unique(land)
#load reclassification table for reclassifying map
classification <- read.table("/Users/natalieschmer/Documents/GitHub/CSU/ECOL_620/lab_7/data_for_lab7/landcover reclass.txt", header=TRUE)

#inspect
head(classification)
classification$Description    #original classification
classification$Description2   #re-class

#format for reclassify function;
class <- as.matrix(classification[,c(1,3)])
land_sub <- reclassify(land,rcl=class)

#plot
plot(land_sub)

#create forested wetlands layer
wetforest <- land_sub
values(wetforest) <- 0
wetforest[land_sub==9 | land_sub==11] <- 1

#create forested uplands layer
dryforest <- land_sub
values(dryforest) <- 0
dryforest[land_sub==10 | land_sub==12] <- 1

#5 km moving window to get neighborhood proportion, is animal selecting pixel or neighborhood around pixel? good for niche modeling
fw <- focalWeight(land_sub, 5000, 'circle')
dry.focal <- focal(dryforest, w=fw, fun="sum", na.rm=T)
wet.focal <- focal(wetforest, w=fw, fun="sum", na.rm=T)

#merge into a single raster stack
layers <- stack(land_sub, wet.focal, dry.focal)
names(layers) <- c("landcover", "wetforest", "dryforest")

#plot
plot(layers)

###################################################
#8.3.3 Home range analysis
###################################################

#------------------#
#mcp home range
#------------------#
# taking the 95 and 50% 
mcp95 <- mcp(panthers[,"CatID"], percent = 95)
mcp50 <- mcp(panthers[,"CatID"], percent = 50)

#inspect
class(mcp95)
head(mcp95@polygons)

#plot
plot(land_sub)
plot(panthers, add=TRUE, col=panthers$CatID)
plot(mcp95, add=TRUE)
plot(mcp50, add=TRUE, col="orange")

#------------------------------------#
#fixed bivariate kernel home range
#------------------------------------#

#kernel types: 
kernel.href.bivar <- kernelUD(panthers[,"CatID"], h="href", kern="bivnorm")
kernel.href.epa <- kernelUD(panthers[,"CatID"], h="href", kern="epa")

#plot
image(kernel.href.bivar)
image(kernel.href.epa)

#alternative plot for first cat
plot(kernel.href.bivar[[1]])
plot(kernel.href.epa[[1]])

#UD data
kernel.href.bivar[[1]]@data

#h value for bandwidth
kernel.href.bivar[[2]]@h
kernel.href.bivar[[2]]@h$h

#least-squares cross validation for h, this will yield a warning of convergence, can't converge on solution so can't use 
kernel.lscv.bivar <- kernelUD(panthers[,"CatID"], h="LSCV", kern="bivnorm")

#manually adjust h, the bandwith parameter, "hard coding"
kernel.bivar.h1000 <- kernelUD(panthers[,"CatID"], h=1000, kern="bivnorm")

#plot first cat
plot(kernel.bivar.h1000[[1]])
plot(kernel.href.bivar[[1]])

#contour maps of activity for fourth cat
plot(kernel.href.bivar[[4]])
contour.data <- as.image.SpatialGridDataFrame(kernel.href.bivar[[4]])
contour(contour.data, nlevels=5, add=TRUE)

#-------------------------------#
#local convex hull home range
#-------------------------------#

#subset
panther147 <- panthers[panthers$CatID==147, ]
panther100 <- panthers[panthers$CatID==100, ]

#plot
dev.off()
plot(panther147)
coordinates(panther147)

#initialize, focusing on 147, get start point to initialize values
k.int <- round(nrow(coordinates(panther147))^0.5,0)
a.int <- round(max(dist(coordinates(panther147))),0)
# searching across neighbors 
k.search <- seq(k.int, 10*k.int, by=5)
a.search <- seq(a.int, 2*a.int, by=3000)

#Parameter search for locoh-a: creating home ranges and going through the parameters 
LoCoH.a.range <- LoCoH.a.area(SpatialPoints(coordinates(panther147)), unout="km2", arange=a.search)

#Parameter search for locoh-k
LoCoH.k.range <- LoCoH.k.area(SpatialPoints(coordinates(panther147)), unout="km2", krange=k.search)

#plot, area relative to params, need a k  
plot(LoCoH.a.range)
plot(LoCoH.k.range) 

#inspect
a.search[5]
k.search[14] # k = 14 for all cats to converge 

#re-fit model
LoCoH.k.61 <- LoCoH.k(SpatialPoints(coordinates(panther147)), k=k.search[14])

#plot
plot(LoCoH.k.61) #stack of all the polygons formed in that process 

#re-fit model
LoCoH.a.100062 <- LoCoH.a(SpatialPoints(coordinates(panther147)), a=a.search[5])
class(LoCoH.a.100062)
#plot
plot(LoCoH.a.100062)

#-----------------------------#
#brownian bridge home range
#-----------------------------#
# need a time series to make this work bc random walk through time
#Re-format Juldate information:
#function for taking characters of a string from rightmost value
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#re-format
panthers$Juldate <- as.character(panthers$Juldate)
panther.date <- as.numeric(substrRight(panthers$Juldate, 3))
panthers$Date <-as.Date(panther.date, origin=as.Date("2006-01-01"))

#convert to POSIXct object
panthers$Date <- as.POSIXct(panthers$Date,"%Y-%m-%d", tz = "EST")

#convert to data frame
panther.df <- as.data.frame(panthers)

#make trajectory object
panther.ltraj <- as.ltraj(xy=coordinates(panthers), date=panthers$Date, id=panthers$CatID, typeII=T)

#inspect
head(panther.ltraj)
head(panther.ltraj[[1]], 2)
summary(panther.ltraj)

#plot, how much room is there to move across landscape 
plot(panther.ltraj)

#telemetry error (average)
sigma2 <- 450

#estimate sigma1, individual- specific metric of pace of movement, will vary
sigma1 <- liker(panther.ltraj, sig2 = sigma2, rangesig1 = c(2, 100))

#inspect
sigma1

#brownian bridge for Cat 147
bb.147 <- kernelbb(panther.ltraj[6], sig1 = 7.2, sig2 = sigma2, grid = 200)

#all panthers
sig1 <- c(sigma1[[1]]$sig1, sigma1[[2]]$sig1, sigma1[[3]]$sig1, sigma1[[4]]$sig1, sigma1[[5]]$sig1, sigma1[[6]]$sig1)
bb.panther <- kernelbb(panther.ltraj, sig1 = sig1, sig2 = sigma2, grid = 200)

#plot
plot(panther.ltraj[6])
plot(bb.147)

#----------------------------#
#contrast estimates
#----------------------------#

#home range area estimates
kernel.95 <- getverticeshr(kernel.href.bivar, percent=95)
bb.95 <- getverticeshr(bb.panther, percent=95)

#contrast area
mcp95$area
kernel.95$area
bb.95$area

#plot
par(mfrow=c(1,2))  #sets up the graph window to store two graphs
plot(land_sub)
plot(kernel.95, add=TRUE, col=kernel.95$id)

plot(land_sub)
plot(mcp95, add=TRUE, col=kernel.95$id)
dev.off()

#write to shapefile if needed
#writePolyShape(mcp95, "homerange")
# Here and below is the resource selection stuff 
########################################################
#8.3.4 Resource Selection
########################################################

###########################################
#8.3.4.1 Point selection functions
###########################################

#use data
use <- extract(layers, panthers)
use <- data.frame(use)

#inspect
head(use)
str(use)

#add CatID
use$CatID <- as.factor(panthers$CatID)

#reformat
useCatID <- dcast(use, CatID~landcover, length, value.var="CatID")

#inspect
useCatID

#add land-cover names
newclass.names <- unique(classification[,3:4])
names(useCatID) <- c("CatID", as.character(newclass.names[1:13,2]))

#inspect
useCatID

#---------------------------------------------------#
#design II availability: population availability
#---------------------------------------------------#

#get availability points
set.seed(8)
rand.II <- sampleRandom(layers, size=1000)
rand.II <- data.frame(rand.II)

#inspect
head(rand.II)
str(rand.II)

rand.II.land <- as.factor(rand.II$landcover)

#get counts of each landcover type
avail.II <- tapply(rand.II.land, rand.II.land, length)

#inspect
avail.II

#add land-cover names
names(avail.II) <- as.character(newclass.names[1:14,2])

#inspect
avail.II

#remove exotics, which were not observed in use sample
avail.II <- avail.II[c(-14)]

#--------------------------------------------------------------------------#
#design III availability: within home-range availability for each individual
#--------------------------------------------------------------------------#

cat.unique <- unique(panthers$CatID)
samples <- 200
rand.III <- matrix(nrow=0, ncol=4)

#loop for all individuals
for(i in 1:length(cat.unique)){

  id.i <- cat.unique[i]
  cat.i <- panthers[panthers$CatID==id.i,]
  mcp.i <- mcp(SpatialPoints(coordinates(cat.i)), percent = 99)
  rand.i <- spsample(mcp.i, type="random", n=samples)
  rand.i.sample <- extract(layers, rand.i)

  #make a matrix of CatID and rand samples
  cat.i <- as.numeric(rep(cat.unique[i], length(rand.i)))
  rand.cat.i <- cbind(cat.i, rand.i.sample)
  rand.III <- rbind(rand.III, rand.cat.i)
}

#inspect
head(rand.III)
class(rand.III)
str(rand.III)

#reshape data
rand.III <- data.frame(rand.III)
rand.III$cat.i <- as.factor(rand.III$cat.i)
avail.III <- dcast(rand.III, cat.i~landcover, length, value.var="cat.i")

names(avail.III)[2:14] <- as.character(newclass.names[1:13,2])
#inspect
avail.III

#---------------------------------------#
#selection ratios
#---------------------------------------#

#Design II:
sel.ratioII <- widesII(u=useCatID[,c(2:ncol(useCatID))], a=as.vector(avail.II), avknown=FALSE, alpha = 0.05)

#inspect
sel.ratioII
sel.ratioII$wi
sel.ratioII$se.wi

#plot
plot(sel.ratioII, errbar = c("CI"))

#Design III:
sel.ratioIII <- widesIII(u=useCatID[,c(2:ncol(useCatID))], a=avail.III[,2:14], avknown=FALSE, alpha = 0.05)

#inspect
sel.ratioIII
sel.ratioIII$wi
sel.ratioIII$se.wi
sel.ratioIII$ICwiupper
sel.ratioIII$ICwilower

#plot
plot(sel.ratioIII, errbar = c("CI"))

