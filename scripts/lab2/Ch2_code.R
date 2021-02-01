########################################################
########################################################
#Fletcher and Fortin 2019
#Chapter 2: Scale
########################################################
########################################################

#load packages
library(raster)      #for raster data; version 2.6-7 used
library(rgdal)       #for raster data, projections; version 1.3-4 used
library(rgeos)       #for buffer analysis; version 0.3-28 used

#set working directory where you downloaded the data
#setwd(choose.dir())
setwd("/Users/natalieschmer/Desktop/GitHub/ECOL_620/data/data_for_lab2")

#############################
#2.3.3 a simple example
#############################

# SETTING UP A RASTER W POISSON DIST
set.seed(16)
#xmx maximum x coordinate (right border). Modify the number to 10 and see what happens 

# 6x6 raster 
#rpois is to create values to fill in the raster (positive integer values with a mean of 3 aka lambda)
toy <- raster(ncol=6, nrow=6, xmn=0, xmx=6, ymn=0, ymx=6)
toy[] <- rpois(ncell(toy), lambda=3)

#plot
plot(toy, axes=F, box=F)
text(toy, digits=2)
res(toy) #resolution of the raster, i.e. grain 1 by 1

#check cell labeling/order
ncell(toy) #number of cells
toy2 <- toy
toy2[] <- 1:ncell(toy)

#plot
plot(toy2)
text(toy2, digits=2)

#aggregate = increase the grain size, by mean or modal
toy_mean <- aggregate(toy, fact=2, fun=mean) #mean value, descrease by factor of 2, 4 pixels --> 1 pixel
toy_maj <- aggregate(toy, fact=2, fun=modal) #majority rule

#plot mean rule
plot(toy_mean)
text(toy_mean,digits=1)

#plot majority rule (whatever is more, randomly picks if there is a tie)
plot(toy_maj)
text(toy_maj)

#contrast means/variances with cell stats for raster global metrics 
cellStats(toy, mean) #or: mean(toy[]) not quite 3 bc small sampling 
cellStats(toy, var) #variance #small sampling  so not quite 3 

#Note the reduction in variance, relative to the variance of toy
cellStats(toy_mean, mean)
cellStats(toy_mean, var) #bc coarser grain size so waaay smaller sampling

cellStats(toy_maj, mean) 
cellStats(toy_maj, var)

#decrease the grain = finer grain, what was 1 pixel is now 4 but still holds same value 
toy_dis2 <- disaggregate(toy, fact=2) #simple method, increase by a factor of 2 without any interpolation  
toy_dis2_bilinear <- disaggregate(toy, fact=2, method='bilinear') # increase by a factor of 2 using bilinear interpolation

#plot
plot(toy_dis2, axes=F, box=F)
plot(rasterToPolygons(toy_dis2), add=TRUE, border='gray50', lwd=1)
text(toy_dis2, cex=0.9)

#plot
plot(toy_dis2_bilinear, axes=F, box=F)
plot(rasterToPolygons(toy_dis2_bilinear), add=TRUE, border='gray50', lwd=1)
text(toy_dis2_bilinear, digits=1, cex=0.6)

#decrease the extent: narrowing boundary
e <- extent(2, 4, 2, 4)#first create new, smaller extent
toy_crop <- crop(toy, e) #clip boundary to new extent, 6x6 -> 2x2 

#plot
plot(toy, zlim=c(0,7))
rect(2, 2, 4, 4, border = "black", lwd = 2)
plot(toy_crop, zlim=c(0,7))

#increase the extent: widen boundary 
e <- extent(0, 7, 0, 7)#first create new, bigger extent
toy_big <- extend(toy,e) #adds a new column 

#plot
plot(toy)
plot(toy_big)

remove(toy,toy2, e, toy_big, toy_crop, toy_dis2, toy_dis2_bilinear, toy_maj, toy_mean)

####################################
#2.3.4.1 multi-scale analysis
####################################

#------------------#
#nlcd
#------------------#

nlcd<-raster("nlcd2011SE")

#inspect
proj4string(nlcd)   #from sp
projection(nlcd)    #alternative function from raster package
crs(nlcd)           #alternative function from raster package (replaces 'projection')

#set projection: just need to define it 
nlcd_proj <- projection(nlcd)

#inspect raster properties: resolution, number of cells, extent 
res(nlcd)
ncell(nlcd)
extent(nlcd)

#check raster values
levels(nlcd) # currently numeric 
nlcd <- as.factor(nlcd) #convert to factors, this may take a little while... (~1-2 mintues)
levels(nlcd) #cover types 
plot(nlcd) #a part of florida? 

#-------------------------------#
#site locations: shp file
#-------------------------------#

#site and reptile data
sites <- readOGR("reptiledata")

#inspect
class(sites)
proj4string(sites)
proj4string(sites) <- nlcd_proj #set projection of the sites to be the same as raster
summary(sites)
head(sites, 2)

#plot with custom color scheme: setting colors for diff colors 
my_col <- c("black","blue","darkorange","red","darkred","grey30","grey50", "lightgreen",
            "green", "darkgreen", "yellow", "goldenrod", "purple", "orchid","lightblue", "lightcyan")

#plot
plot(nlcd, col=my_col, axes=F, box=F) #plotting nlcd and putting our legend colors on 
plot(sites, add=T)

#subset points to remove corn land use
sites <- subset(sites, management!="Corn")
nrow(sites)

#crop raster to 10 km from sampling points: determine min/max coordinates for new extent
x.min <- min(sites$coords_x1) - 10000
x.max <- max(sites$coords_x1) + 10000
y.min <- min(sites$coords_x2) - 10000
y.max <- max(sites$coords_x2) + 10000

extent.new <- extent(x.min, x.max, y.min, y.max)
nlcd <- crop(nlcd, extent.new) #this may take ~20 seconds

#create a binary forest layer using nlcd as template
forest <- nlcd #copy the raster 
#sets all the values of forest to zero
values(forest) <- 0 #set to zero

#reclassify:
#with raster algebra; this is slow, ~2-3 mintues
forest[nlcd==41 | nlcd==42 | nlcd==43] <- 1  #locations with evergreen + mixed forest + deciduous forest, make as "1"

#reclassify with reclassify function is faster
levels(nlcd)[[1]]
reclass <- c(rep(0,7), rep(1,3), rep(0,6)) #first 7 values are 0, next 3 (forest types) are 1, last 6 are 0 
nlcd.levels <- levels(nlcd)[[1]]

#create reclassify matrix: first col: orginal; second: change to
reclass.mat <- cbind(levels(nlcd)[[1]], reclass) #cbinds by the ID? 
reclass.mat

#reclassify
forest <- reclassify(nlcd, reclass.mat) #formal class raster and the ID/binary 

#plot
plot(forest)
plot(sites, pch=21, col="white", add=T)

#define the buffer width (i.e., radius)
buf1km <- 1000 #in m 
buf5km <- 5000 #in m 

#buffer first site
buffer.site1.1km <- buffer(sites[1,], width=buf1km)
plot(buffer.site1.1km)

#buffer using rgeos, which is more flexible
buffer.site1.1km <- gBuffer(sites[1,], width=buf1km, quadsegs=10)
buffer.site1.5km <- gBuffer(sites[1,], width=buf5km, quadsegs=10)

#zoom in on plot for 5 km buffer at site 1
#can provide object to zoom on or click twice on layer

#plot
zoom(nlcd, buffer.site1.5km, col=my_col, box=F)
plot(buffer.site1.1km, border="red", lwd = 3, add=T)
plot(buffer.site1.5km, border="red", lwd = 3, add=T)
points(sites[1,], pch=19, cex=2)
plot(sites[1,], col="grey20", bg="black", pch=22, cex=1, add=T)
dev.off() #end zooming

#view just forest within buffer
zoom(forest, buffer.site1.1km, box=F)
plot(buffer.site1.1km, border="red", lwd = 3,add=T)
dev.off() #end zooming

#calculate forest area within buffer
buffer.forest1.1km <- crop(forest, buffer.site1.1km) #clip to the square extent of the buffer, speeds up the masking  
buffer.forest1.1km <- mask(buffer.forest1.1km, buffer.site1.1km) #the actual circle 

#plot forest within buffer
plot(buffer.forest1.1km)

#calculate percent forest cover
grainarea <- res(forest)[[1]]^2/10000#in ha, resetting resolution  
bufferarea <- (3.14159*buf1km^2)/10000#pi*r^2 #setting the buffer area
forestcover1km <- cellStats(buffer.forest1.1km, 'sum')*grainarea
percentforest1km <- forestcover1km/bufferarea*100
percentforest1km

#-----------------------------------------#
#Function that puts all the steps together
#requires:
#  points: one set of x,y coordinates
#  size: the buffer size (radius), in m
#  landcover: a binary land-cover map
#  grain: the resolution of the map
#-----------------------------------------#

BufferCover <- function(coords, size, landcover, grain){

  bufferarea.i <- pi*size^2/10000                             #size must be in m
  coords.i <- SpatialPoints(cbind(coords[i, 1],coords[i, 2])) #create spatial points from coordinates
  buffer.i <- gBuffer(coords.i, width=size)                   #buffer from rgeos
  crop.i <- crop(landcover, buffer.i)                         #crop with raster function
  crop.NA <- setValues(crop.i, NA)                            #empty raster for the rasterization
  buffer.r <- rasterize(buffer.i, crop.NA)                    #rasterize buffer
  land.buffer <- mask(x=crop.i, mask=buffer.r)                #mask by putting NA outside the boundary
  coveramount<-cellStats(land.buffer, 'sum')*grain            #calculate area
  percentcover<-100*(coveramount/bufferarea.i)                #convert to %

  return(percentcover)
}

#create empty vector for storing output
f1km <- rep(NA, length = nrow(sites))
f2km <- rep(NA, length = nrow(sites))

#with for loop (all five buffers: 910s; <=3km: 228s)
for(i in 1:nrow(sites)) {
  f1km[i] <- BufferCover(coords=sites,size=1000,landcover=forest,grain=grainarea)
  f2km[i] <- BufferCover(coords=sites,size=2000,landcover=forest,grain=grainarea)
  print(i)
}

#make a data frame
forest.scale <- data.frame(site=sites$site,
                         x=sites$coords_x1, y=sites$coords_x2,
                         f1km=f1km, f2km=f2km)

#plot
plot(f1km, f2km)

#correlation matrix
cor(forest.scale[,4:5])

####################################
#2.3.4.2 Scale of effect
####################################

#----------------------------------------#
#2.3.4.2 Buffer analysis
#----------------------------------------#

#herp data
flsk <- read.csv("reptiledata/reptiles_flsk.csv", header=T)
flsk <- merge(flsk, forest.scale, by="site", all=F)

#glms at 2 scales; see text for more scales considered
pres.1km <- glm(pres ~ f1km, family = "binomial", data = flsk)
pres.2km <- glm(pres ~ f2km, family = "binomial", data = flsk)

#summary information
summary(pres.1km)
summary(pres.2km)

#likelihoods
logLik(pres.1km)
logLik(pres.2km)

#accessing coefficients
pres.1km.ci <- confint(pres.1km)
pres.2km.ci <- confint(pres.2km)
pres.1km.ci
pres.2km.ci


